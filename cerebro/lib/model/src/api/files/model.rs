use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};

use crate::api::watchers::model::ProductionWatcher;

use super::retention::{
    ExpiryState, LifecycleTransition, RestoreState, RetentionClass, RetentionPolicy, StorageTier,
};
use super::schema::RegisterFileSchema;

/*
========================
File system and storage
========================
*/

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
pub enum FileType {
    ReadPaired,
    ReadSingle,
    QualityOutput,
    PanviralOutput,
    PathogenOutput,
    /// A consensus genome assembled by the pipeline (diagnostic artefact).
    Consensus,
    /// The Cerebro model / released result document (primary diagnostic artefact).
    CerebroModel,
    /// A run provenance manifest (pipeline version, params, reference DBs, input
    /// hashes); a diagnostic chain-of-custody artefact.
    RunManifest,
    Other,
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
    Host,
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
    /// retrieval when present.
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
    /// Observed replica count, populated by topology/health checks.
    #[serde(default)]
    pub replicas: Option<u32>,
    /// Whether the object's data has been tier-moved to remote archival storage
    /// (S3 Glacier, Model B) and therefore requires a restore before it can be
    /// retrieved.
    #[serde(default)]
    pub archived: bool,
    /// Remote object key where the archived bytes live, set when the file is
    /// archived to the cold object store and cleared when it is restored to
    /// local. The dedicated relocate endpoint is the only writer, so an
    /// archived file always carries the key needed to fetch it back.
    #[serde(default)]
    pub archive_key: Option<String>,
    /// When the result for this file's case was reported out. Anchors the
    /// retention clock and triggers the move to cold storage.
    #[serde(default)]
    pub reported_at: Option<DateTime<Utc>>,
    /// Archival restore state for cold (Glacier) objects (Model B): tracks an
    /// in-flight or completed restore so retrieval can be gated correctly.
    /// Defaulted ([`RestoreState::NotArchived`]) for documents registered before
    /// Set at upload time only for archival objects, otherwise via the
    /// lifecycle update endpoint.
    #[serde(default)]
    pub restore_state: RestoreState,
    /// Soft-delete expiry state. [`ExpiryState::Active`] until a retention sweep
    /// quarantines the file. Defaulted for documents registered before this field was added.
    #[serde(default)]
    pub expiry_state: ExpiryState,
    /// When the file was quarantined by expiry (`None` while active). Anchors the
    /// quarantine grace window before an eligible hard purge.
    #[serde(default)]
    pub quarantined_at: Option<DateTime<Utc>>,
    /// In-flight tier move target claimed by a mover. `Some` while a
    /// physical move is in progress; cleared when the move commits, so the DB
    /// never reports a final `tier` whose bytes aren't verified at destination.
    #[serde(default)]
    pub pending_tier: Option<StorageTier>,
    /// When the file's `tier` was last committed by a verified move.
    #[serde(default)]
    pub tier_moved_at: Option<DateTime<Utc>>,
    /// When the current `pending_tier` claim was made by a mover.
    ///
    /// Stamped at claim time and cleared when the move commits or rolls back.
    /// Because [`tier_moved_at`](Self::tier_moved_at) is only written on a
    /// *successful* commit, a worker that dies mid-move would otherwise leave a
    /// `pending_tier` with no timestamp and nothing to age it out. `pending_since`
    /// gives the tier-move scan an unambiguous claim age so it can re-drive a
    /// stale claim (the mover is idempotent and CAS-guarded, so re-driving safely
    /// completes or rolls back).
    #[serde(default)]
    pub pending_since: Option<DateTime<Utc>>,
    /// When this file's stored bytes were last successfully verified against the
    /// catalogue BLAKE3 hash.
    ///
    /// `None` means never verified. The verify scan orders by this field ascending
    /// (nulls first) so a bounded per-run budget sweeps the *oldest-verified* files
    /// and achieves full estate coverage over time, rather than re-checking the
    /// same head of the list every run. Only a successful verify stamps it; a
    /// failure leaves it unchanged so the file stays at the front of the queue.
    #[serde(default)]
    pub verified_at: Option<DateTime<Utc>>,
    /// When an archival restore was requested.
    #[serde(default)]
    pub restore_requested_at: Option<DateTime<Utc>>,
    /// When a restore completed and the object became retrievable.
    #[serde(default)]
    pub restore_available_at: Option<DateTime<Utc>>,
    /// When the restored object's availability window elapses.
    #[serde(default)]
    pub restore_expires_at: Option<DateTime<Utc>>,
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
            archive_key: None,
            reported_at: register_file_schema.reported_at,
            restore_state: RestoreState::default(),
            expiry_state: ExpiryState::default(),
            quarantined_at: None,
            pending_tier: None,
            tier_moved_at: None,
            pending_since: None,
            verified_at: None,
            restore_requested_at: None,
            restore_available_at: None,
            restore_expires_at: None,
        }
    }

    /// Whether this file is eligible to be quarantined by an expiry sweep at
    /// `now`: it is past its `retain_until`, not under legal hold, and still
    /// active. Files with no `retain_until` (indefinite retention) are never
    /// eligible.
    pub fn is_expiry_eligible(&self, now: DateTime<Utc>) -> bool {
        matches!(self.expiry_state, ExpiryState::Active)
            && !self.legal_hold
            && matches!(self.retain_until, Some(until) if until <= now)
    }

    /// Whether this quarantined file is eligible for a hard purge at `now`: it is
    /// quarantined, not under legal hold, and its quarantine grace window
    /// (`grace_days`) has elapsed since `quarantined_at`.
    pub fn is_purge_eligible(&self, now: DateTime<Utc>, grace_days: i64) -> bool {
        matches!(self.expiry_state, ExpiryState::Quarantined)
            && !self.legal_hold
            && matches!(self.quarantined_at, Some(at) if at + chrono::Duration::days(grace_days) <= now)
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
    pub reads_2: Option<SeaweedFile>,
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
            archive_key: None,
            reported_at: None,
            restore_state: super::super::retention::RestoreState::NotArchived,
            expiry_state: super::super::retention::ExpiryState::Active,
            quarantined_at: None,
            pending_tier: None,
            tier_moved_at: None,
            pending_since: None,
            verified_at: None,
            restore_requested_at: None,
            restore_available_at: None,
            restore_expires_at: None,
        }
    }

    #[test]
    fn expiry_eligible_only_when_lapsed_active_and_unheld() {
        use super::super::retention::ExpiryState;
        let now = Utc.with_ymd_and_hms(2026, 6, 14, 0, 0, 0).unwrap();
        let past = now - chrono::Duration::days(1);
        let future = now + chrono::Duration::days(1);

        let mut f = sample_file();
        // No retain_until => indefinite => never eligible.
        assert!(!f.is_expiry_eligible(now));

        f.retain_until = Some(future);
        assert!(!f.is_expiry_eligible(now)); // not yet lapsed

        f.retain_until = Some(past);
        assert!(f.is_expiry_eligible(now)); // lapsed, active, unheld

        f.legal_hold = true;
        assert!(!f.is_expiry_eligible(now)); // hold overrides
        f.legal_hold = false;

        f.expiry_state = ExpiryState::Quarantined;
        assert!(!f.is_expiry_eligible(now)); // already quarantined
    }

    #[test]
    fn purge_eligible_after_grace_window() {
        use super::super::retention::ExpiryState;
        let now = Utc.with_ymd_and_hms(2026, 6, 14, 0, 0, 0).unwrap();
        let mut f = sample_file();

        f.expiry_state = ExpiryState::Quarantined;
        f.quarantined_at = Some(now - chrono::Duration::days(10));
        assert!(!f.is_purge_eligible(now, 30)); // within grace
        assert!(f.is_purge_eligible(now, 7)); // grace elapsed

        f.legal_hold = true;
        assert!(!f.is_purge_eligible(now, 7)); // hold overrides
    }

    #[test]
    fn report_out_without_warm_moves_to_cold() {
        use super::super::retention::RetentionPolicy;
        let policy = RetentionPolicy {
            diagnostic_days: 365 * 4,
            intermediate_days: 365,
            transient_days: 30,
            warm_days: 90,
        };
        let reported_at = Utc.with_ymd_and_hms(2026, 6, 13, 0, 0, 0).unwrap();

        let mut f = sample_file(); // tier defaults to Hot, retention Diagnostic
        let transition = f.report_out(reported_at, &policy, false);

        assert_eq!(f.reported_at, Some(reported_at));
        assert_eq!(f.tier, StorageTier::Cold);
        assert_eq!(
            f.retain_until,
            Some(reported_at + chrono::Duration::days(365 * 4))
        );
        assert_eq!(transition.target_tier, StorageTier::Cold);
        // legal hold still governs expiry independently
        assert!(!f.is_expired(reported_at));
    }

    #[test]
    fn report_out_with_warm_moves_to_warm() {
        use super::super::retention::RetentionPolicy;
        let policy = RetentionPolicy {
            diagnostic_days: 365 * 4,
            intermediate_days: 365,
            transient_days: 30,
            warm_days: 90,
        };
        let reported_at = Utc.with_ymd_and_hms(2026, 6, 13, 0, 0, 0).unwrap();

        let mut f = sample_file();
        let transition = f.report_out(reported_at, &policy, true);

        assert_eq!(f.tier, StorageTier::Warm);
        assert_eq!(
            transition.cold_move_at,
            Some(reported_at + chrono::Duration::days(90))
        );
        // retention unchanged by the warm placement
        assert_eq!(
            f.retain_until,
            Some(reported_at + chrono::Duration::days(365 * 4))
        );
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
