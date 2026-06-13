use serde::{Serialize, Deserialize};
use chrono::{DateTime, Utc};

use crate::api::watchers::model::ProductionWatcher;

use super::retention::{LifecycleTransition, RetentionClass, RetentionPolicy, StorageTier};
use super::schema::RegisterFileSchema;

/*
========================
File system and storage
========================
*/

#[derive(Debug, Clone, Serialize, Deserialize, clap::ValueEnum)]
pub enum FileType {
    ReadPaired,
    ReadSingle,
    QualityOutput,
    PanviralOutput,
    PathogenOutput,
    Other
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FileTag {
    #[serde(rename = "DNA")]
    Dna,
    #[serde(rename = "RNA")]
    Rna,
    #[serde(rename = "POS")]
    Pos,
    #[serde(rename = "NEG")]
    Neg,
    #[serde(rename = "NTC")]
    Ntc,
    #[serde(rename = "ENV")]
    Env,
    #[serde(rename = "HOST")]
    Host
}

pub type SeaweedFileId = String;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeaweedFile {
    pub id: String,
    pub date: String,
    pub name: String,
    pub hash: String,
    pub size: u64,
    pub fid: SeaweedFileId,
    pub tags: Vec<FileTag>,
    pub run_id: Option<String>,
    pub sample_id: Option<String>,
    pub ftype: Option<FileType>,
    pub watcher: Option<ProductionWatcher>,
    /// Filer object path, when the file was stored via the path-addressed filer
    /// rather than as a fid-addressed weed object. Preferred over `fid` for
    /// retrieval when present. Defaulted for documents registered before FS-2.
    #[serde(default)]
    pub path: Option<String>,
    /// Physical storage tier the file currently occupies.
    #[serde(default)]
    pub tier: StorageTier,
    /// Retention category assigned at registration.
    #[serde(default)]
    pub retention: RetentionClass,
    /// Absolute expiry computed from the retention policy in force at
    /// registration. `None` means "retain indefinitely".
    #[serde(default)]
    pub retain_until: Option<DateTime<Utc>>,
    /// When set, the file is exempt from expiry regardless of `retain_until`.
    #[serde(default)]
    pub legal_hold: bool,
    /// Observed replica count, populated by topology/health checks (FS-6).
    #[serde(default)]
    pub replicas: Option<u32>,
    /// Whether the object's data has been tier-moved to remote archival storage
    /// (S3 Glacier, Model B) and therefore requires a restore before it can be
    /// retrieved. Defaulted for documents registered before FS-4.
    #[serde(default)]
    pub archived: bool,
    /// When the result for this file's case was reported out. Anchors the
    /// retention clock and triggers the move to cold storage. Defaulted for
    /// documents registered before FS-7.
    #[serde(default)]
    pub reported_at: Option<DateTime<Utc>>,
}
impl SeaweedFile {
    pub fn from_schema(register_file_schema: &RegisterFileSchema) -> Self {
        Self {
            id: register_file_schema.id.clone(),
            run_id: register_file_schema.run_id.clone(),
            sample_id: register_file_schema.sample_id.clone(),
            date: register_file_schema.date.clone(),
            name: register_file_schema.name.clone(),
            hash: register_file_schema.hash.clone(),
            ftype: register_file_schema.ftype.clone(),
            fid: register_file_schema.fid.clone(),
            size: register_file_schema.size.clone(),
            watcher: register_file_schema.watcher.clone(),
            tags: Vec::new(),
            path: register_file_schema.path.clone(),
            tier: register_file_schema.tier,
            retention: register_file_schema.retention,
            retain_until: register_file_schema.retain_until,
            legal_hold: register_file_schema.legal_hold,
            replicas: register_file_schema.replicas,
            archived: register_file_schema.archived,
            reported_at: register_file_schema.reported_at,
        }
    }
    pub fn size_mb(&self) -> f64 {
        bytes_to_mb(self.size)
    }

    /// Identifier preferred for retrieval: the filer `path` when present and
    /// non-empty, otherwise the weed `fid`.
    ///
    /// This lets retrieval transparently handle both path-addressed (filer) and
    /// fid-addressed (weed) objects without the caller needing to know which
    /// backend stored the file.
    pub fn effective_identifier(&self) -> &str {
        match &self.path {
            Some(path) if !path.is_empty() => path.as_str(),
            _ => self.fid.as_str(),
        }
    }

    /// Whether the object must be restored from archival storage before it can
    /// be retrieved.
    ///
    /// True when the data has been tier-moved to remote Glacier storage
    /// (`archived`). Note that [`StorageTier::Cold`](crate::api::files::retention::StorageTier::Cold)
    /// alone does **not** imply this: in Model A the cold tier is local HDD and
    /// is directly readable, whereas in Model B it is S3 Glacier.
    pub fn requires_restore(&self) -> bool {
        self.archived
    }

    /// Apply the report-out lifecycle transition to this record in place.
    ///
    /// Records `reported_at`, re-anchors `retain_until` to that moment using the
    /// supplied [`RetentionPolicy`], and moves the record to its post-report
    /// tier: **warm** when `warm_available` (and a positive warm dwell is
    /// configured) so the data stays directly re-inspectable, otherwise
    /// **cold**. The returned [`LifecycleTransition`] is what a worker/server
    /// persists and acts on — including [`LifecycleTransition::cold_move_at`],
    /// the later warm→cold (S3) move. Legal hold is untouched and continues to
    /// override expiry.
    pub fn report_out(
        &mut self,
        reported_at: DateTime<Utc>,
        policy: &RetentionPolicy,
        warm_available: bool,
    ) -> LifecycleTransition {
        let transition = policy.report_out_transition(self.retention, reported_at, warm_available);
        self.reported_at = Some(transition.reported_at);
        self.retain_until = transition.retain_until;
        self.tier = transition.target_tier;
        transition
    }

    /// Whether the file is eligible for expiry at instant `now`.
    ///
    /// Always `false` while under [`SeaweedFile::legal_hold`]; otherwise `true`
    /// once `now` has reached `retain_until`. A file with no `retain_until`
    /// (indefinite retention) never expires.
    pub fn is_expired(&self, now: DateTime<Utc>) -> bool {
        if self.legal_hold {
            return false;
        }
        match self.retain_until {
            Some(until) => now >= until,
            None => false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SeaweedReads {
    pub reads_1: SeaweedFile,
    pub reads_2: Option<SeaweedFile>
}


/// Converts a size in bytes to megabytes (MB).
///
/// # Arguments
///
/// * `bytes` - The size in bytes.
///
/// # Returns
///
/// The size in megabytes (MB) as a floating point number.
fn bytes_to_mb(bytes: u64) -> f64 {
    const BYTES_PER_MB: f64 = 1024.0 * 1024.0; // 1 MB = 1024 * 1024 bytes
    bytes as f64 / BYTES_PER_MB
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;

    fn sample_file() -> SeaweedFile {
        SeaweedFile {
            id: "id".into(),
            date: "2025-01-01".into(),
            name: "reads_R1.fastq.gz".into(),
            hash: "abc".into(),
            size: 1024,
            fid: "3,01abcd".into(),
            tags: Vec::new(),
            run_id: None,
            sample_id: None,
            ftype: None,
            watcher: None,
            path: None,
            tier: StorageTier::default(),
            retention: RetentionClass::default(),
            retain_until: None,
            legal_hold: false,
            replicas: None,
            archived: false,
            reported_at: None,
        }
    }

    #[test]
    fn report_out_without_warm_moves_to_cold() {
        use super::super::retention::RetentionPolicy;
        let policy = RetentionPolicy { diagnostic_days: 365 * 4, intermediate_days: 365, transient_days: 30, warm_days: 90 };
        let reported_at = Utc.with_ymd_and_hms(2026, 6, 13, 0, 0, 0).unwrap();

        let mut f = sample_file(); // tier defaults to Hot, retention Diagnostic
        let transition = f.report_out(reported_at, &policy, false);

        assert_eq!(f.reported_at, Some(reported_at));
        assert_eq!(f.tier, StorageTier::Cold);
        assert_eq!(f.retain_until, Some(reported_at + chrono::Duration::days(365 * 4)));
        assert_eq!(transition.target_tier, StorageTier::Cold);
        // legal hold still governs expiry independently
        assert!(!f.is_expired(reported_at));
    }

    #[test]
    fn report_out_with_warm_moves_to_warm() {
        use super::super::retention::RetentionPolicy;
        let policy = RetentionPolicy { diagnostic_days: 365 * 4, intermediate_days: 365, transient_days: 30, warm_days: 90 };
        let reported_at = Utc.with_ymd_and_hms(2026, 6, 13, 0, 0, 0).unwrap();

        let mut f = sample_file();
        let transition = f.report_out(reported_at, &policy, true);

        assert_eq!(f.tier, StorageTier::Warm);
        assert_eq!(transition.cold_move_at, Some(reported_at + chrono::Duration::days(90)));
        // retention unchanged by the warm placement
        assert_eq!(f.retain_until, Some(reported_at + chrono::Duration::days(365 * 4)));
    }

    #[test]
    fn requires_restore_tracks_archived_flag() {
        let mut f = sample_file();
        assert!(!f.requires_restore());
        f.tier = StorageTier::Cold; // cold alone does not imply restore (Model A = local HDD)
        assert!(!f.requires_restore());
        f.archived = true; // tier-moved to Glacier (Model B)
        assert!(f.requires_restore());
    }

    #[test]
    fn effective_identifier_prefers_non_empty_path() {
        let mut f = sample_file();
        assert_eq!(f.effective_identifier(), "3,01abcd");
        f.path = Some(String::new());
        assert_eq!(f.effective_identifier(), "3,01abcd"); // empty path ignored
        f.path = Some("/run/sample/reads_R1.fastq.gz".into());
        assert_eq!(f.effective_identifier(), "/run/sample/reads_R1.fastq.gz");
    }

    #[test]
    fn legal_hold_blocks_expiry() {
        let now = Utc.with_ymd_and_hms(2030, 1, 1, 0, 0, 0).unwrap();
        let mut f = sample_file();
        f.retain_until = Some(Utc.with_ymd_and_hms(2025, 1, 1, 0, 0, 0).unwrap());
        assert!(f.is_expired(now)); // past retain_until, no hold
        f.legal_hold = true;
        assert!(!f.is_expired(now)); // hold overrides expiry
    }

    #[test]
    fn no_retain_until_never_expires() {
        let now = Utc.with_ymd_and_hms(2030, 1, 1, 0, 0, 0).unwrap();
        let f = sample_file();
        assert!(!f.is_expired(now));
    }
}
