//! Storage tiering and retention policy types for Cerebro FS.
//!
//! These types describe *where* a file currently lives ([`StorageTier`]) and
//! *how long* it must be kept ([`RetentionClass`] resolved through a
//! [`RetentionPolicy`]). They are the data-model foundation that the deployment
//! models (ingle-server hot/cold, HPC + S3 cold) and the lifecycle
//! workers (Stage 3) act on, and that the audit/chain-of-custody layer
//! (Stage 2 · S2-6) records against.
//!
//! # Retention is configurable, not hardcoded
//!
//! NATA/NPAAC retention requirements vary by record type and jurisdiction, and
//! are a laboratory policy decision rather than something this crate can assert.
//! The durations in [`RetentionPolicy::default`] are therefore **placeholders**
//! intended only to give a safe, keep-rather-than-delete ordering out of the
//! box — they are not legal advice. Configure them to your laboratory's
//! accreditation scope.
//!
//! Two further design points support compliance:
//!
//! * A file stores **both** its [`RetentionClass`] and the absolute
//!   `retain_until` timestamp computed when it was registered. The class lets a
//!   schedule be re-evaluated later; the absolute stamp records the obligation
//!   that was in force at the time of registration.
//! * `legal_hold` overrides expiry entirely: a file under hold is never
//!   eligible for deletion regardless of `retain_until`.

use chrono::{DateTime, Duration, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Physical storage tier a file currently occupies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
pub enum StorageTier {
    /// Low-latency online storage (e.g. SSD); newly ingested and active data.
    Hot,
    /// Capacity online storage (e.g. HDD); ageing but still directly readable.
    Warm,
    /// Archival storage (e.g. S3 Glacier); retrieval may require a restore step.
    Cold,
}

impl Default for StorageTier {
    fn default() -> Self {
        StorageTier::Hot
    }
}

impl StorageTier {
    /// Ordinal rank from online to archival (Hot=0, Warm=1, Cold=2), used to
    /// classify a tier move's direction.
    pub fn rank(&self) -> u8 {
        match self {
            StorageTier::Hot => 0,
            StorageTier::Warm => 1,
            StorageTier::Cold => 2,
        }
    }

    /// True when moving to `target` is a demotion (ageing toward archival), e.g.
    /// Hot→Warm or Warm→Cold.
    pub fn is_demotion_to(&self, target: &StorageTier) -> bool {
        target.rank() > self.rank()
    }

    /// True when moving to `target` is a promotion (restoring toward online), e.g.
    /// Cold→Warm or Warm→Hot.
    pub fn is_promotion_to(&self, target: &StorageTier) -> bool {
        target.rank() < self.rank()
    }
}

impl std::fmt::Display for StorageTier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StorageTier::Hot => write!(f, "hot"),
            StorageTier::Warm => write!(f, "warm"),
            StorageTier::Cold => write!(f, "cold"),
        }
    }
}

/// Retention category assigned to a file at registration.
///
/// The category is durable; the *duration* it implies is resolved through a
/// [`RetentionPolicy`], so retained data can be re-evaluated if a laboratory's
/// retention schedule changes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
pub enum RetentionClass {
    /// Primary diagnostic evidence (raw reads, released results). Longest kept.
    Diagnostic,
    /// Working/intermediate outputs useful for review but not primary evidence.
    Intermediate,
    /// Scratch data with no retention obligation; safe to expire quickly.
    Transient,
}

impl Default for RetentionClass {
    /// Defaults to [`RetentionClass::Diagnostic`] — the safe, keep-by-default
    /// choice for a diagnostics platform.
    fn default() -> Self {
        RetentionClass::Diagnostic
    }
}

impl std::fmt::Display for RetentionClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RetentionClass::Diagnostic => write!(f, "diagnostic"),
            RetentionClass::Intermediate => write!(f, "intermediate"),
            RetentionClass::Transient => write!(f, "transient"),
        }
    }
}

/// Laboratory-configurable mapping from [`RetentionClass`] to a retention period
/// in days.
///
/// A non-positive number of days is treated as "no expiry" (the file is kept
/// indefinitely until an explicit administrative action removes it).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RetentionPolicy {
    /// Retention period (days) for [`RetentionClass::Diagnostic`].
    pub diagnostic_days: i64,
    /// Retention period (days) for [`RetentionClass::Intermediate`].
    pub intermediate_days: i64,
    /// Retention period (days) for [`RetentionClass::Transient`].
    pub transient_days: i64,
    /// How long reported-out data dwells on the **warm** tier (directly-readable
    /// HDD, for re-inspection) before ageing to the cold (S3) tier, in days.
    /// Only meaningful when the deployment has a warm tier (three-tier Model B);
    /// non-positive means "no warm dwell" (move straight to cold).
    #[serde(default = "default_warm_days")]
    pub warm_days: i64,
}

/// Default warm-tier dwell (re-inspection window) in days — placeholder, confirm
/// against your accreditation scope. Used by serde for configs written before the
/// `warm_days` field existed.
fn default_warm_days() -> i64 {
    365
}

/// Default quarantine grace window (days) between expiry-quarantine and an
/// eligible hard purge. Configurable via `CEREBRO_QUARANTINE_GRACE_DAYS`.
pub const DEFAULT_QUARANTINE_GRACE_DAYS: i64 = 30;

/// Quarantine grace window in days, from `CEREBRO_QUARANTINE_GRACE_DAYS`,
/// falling back to [`DEFAULT_QUARANTINE_GRACE_DAYS`].
///
/// This is the minimum dwell a soft-deleted (quarantined) file must spend before
/// it becomes eligible for a hard purge — a safety margin against accidental or
/// premature expiry.
pub fn quarantine_grace_days_from_env() -> i64 {
    let file = retention_config_file();
    resolve_config(&file, "CEREBRO_QUARANTINE_GRACE_DAYS")
        .and_then(|v| v.trim().parse::<i64>().ok())
        .unwrap_or(DEFAULT_QUARANTINE_GRACE_DAYS)
}

/// Soft-delete lifecycle state used by non-destructive expiry.
///
/// Expiry never deletes: a file past its retention (and not under legal hold)
/// is moved to [`ExpiryState::Quarantined`] and removed from active use. A
/// separate, explicitly-gated purge hard-deletes only quarantined files that
/// have cleared the quarantine grace window. New files are [`ExpiryState::Active`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ExpiryState {
    /// Live: within retention, or retained by a legal hold.
    Active,
    /// Past retention and soft-deleted; pending a gated hard purge.
    Quarantined,
}

impl Default for ExpiryState {
    fn default() -> Self {
        ExpiryState::Active
    }
}

impl std::fmt::Display for ExpiryState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExpiryState::Active => write!(f, "active"),
            ExpiryState::Quarantined => write!(f, "quarantined"),
        }
    }
}

impl Default for RetentionPolicy {
    /// Placeholder durations — **not legal advice**. Confirm and override these
    /// against your laboratory's NATA/NPAAC accreditation scope before relying
    /// on them. They are deliberately conservative (keep rather than delete).
    fn default() -> Self {
        Self {
            diagnostic_days: 365 * 4, // 4 years (NPAAC-aligned default; configurable)
            intermediate_days: 90,    // ~3 months (configurable)
            transient_days: 30,       // ~1 month (configurable)
            warm_days: 365,           // ~1 year warm re-inspection window (configurable)
        }
    }
}

/// Parse a day count from an optional environment value, falling back to
/// `default` when absent or unparseable.
fn parse_days(value: Option<String>, default: i64) -> i64 {
    value
        .and_then(|v| v.trim().parse::<i64>().ok())
        .unwrap_or(default)
}

/// Load `KEY=VALUE` pairs from the optional retention config file named by
/// `CEREBRO_RETENTION_CONFIG_FILE` (e.g. a Docker secret mounted at
/// `/run/secrets/retention_config`). Blank lines and `#` comments are ignored;
/// an unset or unreadable path yields an empty map.
fn retention_config_file() -> HashMap<String, String> {
    let mut map = HashMap::new();
    if let Ok(path) = std::env::var("CEREBRO_RETENTION_CONFIG_FILE") {
        if let Ok(contents) = std::fs::read_to_string(&path) {
            for line in contents.lines() {
                let line = line.trim();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                if let Some((key, value)) = line.split_once('=') {
                    map.insert(key.trim().to_string(), value.trim().to_string());
                }
            }
        }
    }
    map
}

/// Resolve a config value with precedence: retention config file, then the
/// process environment, then `None`.
fn resolve_config(file: &HashMap<String, String>, key: &str) -> Option<String> {
    file.get(key).cloned().or_else(|| std::env::var(key).ok())
}

/// Warm-tier availability (`CEREBRO_FS_WARM_AVAILABLE`), config-file-aware.
///
/// True when the deployment has a warm tier, so report-out targets warm before
/// cold. Resolved from the retention config file or the environment.
pub fn warm_available_from_env() -> bool {
    let file = retention_config_file();
    matches!(
        resolve_config(&file, "CEREBRO_FS_WARM_AVAILABLE").as_deref(),
        Some("true") | Some("1")
    )
}

impl RetentionPolicy {
    /// Build a policy from environment variables, falling back to the
    /// [`Default`] durations for any that are unset or unparseable.
    ///
    /// Reads `CEREBRO_RETENTION_DIAGNOSTIC_DAYS`,
    /// `CEREBRO_RETENTION_INTERMEDIATE_DAYS`, `CEREBRO_RETENTION_TRANSIENT_DAYS`
    /// and `CEREBRO_RETENTION_WARM_DAYS`. This is how a deployment supplies its
    /// accreditation-scoped periods instead of the placeholder defaults; set
    /// these to your confirmed NATA/NPAAC values.
    pub fn from_env() -> Self {
        let d = Self::default();
        let file = retention_config_file();
        Self {
            diagnostic_days: parse_days(
                resolve_config(&file, "CEREBRO_RETENTION_DIAGNOSTIC_DAYS"),
                d.diagnostic_days,
            ),
            intermediate_days: parse_days(
                resolve_config(&file, "CEREBRO_RETENTION_INTERMEDIATE_DAYS"),
                d.intermediate_days,
            ),
            transient_days: parse_days(
                resolve_config(&file, "CEREBRO_RETENTION_TRANSIENT_DAYS"),
                d.transient_days,
            ),
            warm_days: parse_days(
                resolve_config(&file, "CEREBRO_RETENTION_WARM_DAYS"),
                d.warm_days,
            ),
        }
    }

    /// Retention period in days configured for `class`.
    pub fn days_for(&self, class: RetentionClass) -> i64 {
        match class {
            RetentionClass::Diagnostic => self.diagnostic_days,
            RetentionClass::Intermediate => self.intermediate_days,
            RetentionClass::Transient => self.transient_days,
        }
    }

    /// Compute the expiry timestamp for a file registered at `from` under
    /// `class`.
    ///
    /// Returns `None` when the configured duration is non-positive, which is the
    /// canonical representation of "retain indefinitely".
    ///
    /// # Examples
    ///
    /// ```
    /// use chrono::{TimeZone, Utc};
    /// use cerebro_model::api::files::retention::{RetentionClass, RetentionPolicy};
    ///
    /// let policy = RetentionPolicy { diagnostic_days: 10, intermediate_days: 0, transient_days: 1, warm_days: 0 };
    /// let from = Utc.with_ymd_and_hms(2025, 1, 1, 0, 0, 0).unwrap();
    ///
    /// assert!(policy.retain_until(RetentionClass::Diagnostic, from).is_some());
    /// assert!(policy.retain_until(RetentionClass::Intermediate, from).is_none()); // 0 days => no expiry
    /// ```
    pub fn retain_until(
        &self,
        class: RetentionClass,
        from: DateTime<Utc>,
    ) -> Option<DateTime<Utc>> {
        let days = self.days_for(class);
        if days <= 0 {
            None
        } else {
            Some(from + Duration::days(days))
        }
    }

    /// Compute the lifecycle transition for a file reported out at `reported_at`
    /// under retention `class`.
    ///
    /// Anchors the retention clock to the report-out moment (the clinically and
    /// legally meaningful event) rather than to upload time, and decides the tier
    /// the data should move to:
    ///
    /// * **`warm_available` and `warm_days > 0`** (three-tier Model B): move to
    ///   the warm tier (directly-readable HDD) for re-inspection, and schedule a
    ///   later move to cold (S3) at `reported_at + warm_days`
    ///   ([`LifecycleTransition::cold_move_at`]).
    /// * **otherwise** (Model A, or no warm dwell configured): move straight to
    ///   the cold tier; there is no scheduled cold move.
    ///
    /// Whether "cold" is local HDD (Model A) or S3 Glacier (Model B, where the
    /// file also becomes `archived`) is resolved at execution time by the
    /// deployment-aware lifecycle worker (Stage 3). Retention (`retain_until`) is
    /// independent of placement and is unchanged by the later warm→cold move.
    pub fn report_out_transition(
        &self,
        class: RetentionClass,
        reported_at: DateTime<Utc>,
        warm_available: bool,
    ) -> LifecycleTransition {
        let (target_tier, cold_move_at) = if warm_available && self.warm_days > 0 {
            (
                StorageTier::Warm,
                Some(reported_at + Duration::days(self.warm_days)),
            )
        } else {
            (StorageTier::Cold, None)
        };
        LifecycleTransition {
            reported_at,
            retain_until: self.retain_until(class, reported_at),
            target_tier,
            cold_move_at,
        }
    }
}

/// The storage-lifecycle changes implied by reporting a result out.
///
/// When a diagnostic result is reported out, its data is no longer active and
/// the retention clock is anchored to that moment. This captures both decisions
/// in one value: the new expiry and the tier the data should move to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LifecycleTransition {
    /// Timestamp the result was reported out (the retention anchor).
    pub reported_at: DateTime<Utc>,
    /// Recomputed expiry: `reported_at` plus the retention period for the class
    /// (`None` for indefinite retention).
    pub retain_until: Option<DateTime<Utc>>,
    /// Tier the data should move to immediately on report-out: warm when the
    /// deployment has a warm tier (and a positive dwell), otherwise cold.
    pub target_tier: StorageTier,
    /// When the data should subsequently age from warm to the cold (S3) tier
    /// (`reported_at + warm_days`). `None` when it moves straight to cold.
    pub cold_move_at: Option<DateTime<Utc>>,
}

/// Lifecycle state of an archival (Glacier) restore for a cold-tiered object.
///
/// In Model B (HPC distributed), the cold tier is S3 Glacier: a file whose data
/// has been tier-moved there is `archived` and cannot be read directly. A
/// restore must first be initiated and completed, after which the object is
/// temporarily retrievable. This enum models that contract; see
/// [`SeaweedFile::requires_restore`](crate::api::files::model::SeaweedFile::requires_restore).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
pub enum RestoreState {
    /// No active restore. The object is directly retrievable (hot/warm) or an
    /// idle archived object with no restore in progress (see `archived`).
    NotArchived,
    /// A restore has been requested from archival storage, awaiting the provider.
    Requested,
    /// The archival provider (e.g. S3 Glacier) is restoring the object.
    InProgress,
    /// A restore completed; the object is temporarily retrievable until the
    /// availability window elapses.
    Restored,
    /// The most recent restore attempt failed.
    Failed,
    /// A completed restore's availability window has elapsed; the object is
    /// archived-only again.
    Expired,
}

impl Default for RestoreState {
    fn default() -> Self {
        RestoreState::NotArchived
    }
}

impl RestoreState {
    /// Whether a direct transition from `self` to `next` is valid (S2-11).
    ///
    /// The restore worker drives: `NotArchived → Requested → InProgress →
    /// Restored`, with `→ Failed` from the in-flight states (retryable back to
    /// `Requested`), `Restored → Expired` when the window lapses, and a return to
    /// `NotArchived` from terminal states. Self-transitions are handled as
    /// idempotent no-ops by the endpoint, not here.
    pub fn can_transition_to(&self, next: &RestoreState) -> bool {
        use RestoreState::*;
        matches!(
            (self, next),
            (NotArchived, Requested)
                | (Requested, InProgress)
                | (Requested, Restored)
                | (Requested, Failed)
                | (InProgress, Restored)
                | (InProgress, Failed)
                | (Restored, Expired)
                | (Restored, NotArchived)
                | (Failed, Requested)
                | (Failed, NotArchived)
                | (Expired, Requested)
                | (Expired, NotArchived)
        )
    }
}

/// Default restore availability window (days) — how long a restored archival
/// object stays retrievable. Configurable via `CEREBRO_RESTORE_AVAILABLE_DAYS`.
pub const DEFAULT_RESTORE_AVAILABLE_DAYS: i64 = 7;

/// Restore availability window in days, from `CEREBRO_RESTORE_AVAILABLE_DAYS`,
/// falling back to [`DEFAULT_RESTORE_AVAILABLE_DAYS`].
pub fn restore_available_days_from_env() -> i64 {
    std::env::var("CEREBRO_RESTORE_AVAILABLE_DAYS")
        .ok()
        .and_then(|v| v.trim().parse::<i64>().ok())
        .unwrap_or(DEFAULT_RESTORE_AVAILABLE_DAYS)
}

impl std::fmt::Display for RestoreState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RestoreState::NotArchived => write!(f, "not-archived"),
            RestoreState::Requested => write!(f, "requested"),
            RestoreState::InProgress => write!(f, "in-progress"),
            RestoreState::Restored => write!(f, "restored"),
            RestoreState::Failed => write!(f, "failed"),
            RestoreState::Expired => write!(f, "expired"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;

    #[test]
    fn defaults_are_keep_by_default() {
        assert_eq!(StorageTier::default(), StorageTier::Hot);
        assert_eq!(RetentionClass::default(), RetentionClass::Diagnostic);
    }

    #[test]
    fn tier_transition_direction() {
        assert!(StorageTier::Hot.is_demotion_to(&StorageTier::Warm));
        assert!(StorageTier::Warm.is_demotion_to(&StorageTier::Cold));
        assert!(StorageTier::Cold.is_promotion_to(&StorageTier::Hot));
        assert!(!StorageTier::Hot.is_demotion_to(&StorageTier::Hot));
        assert!(!StorageTier::Hot.is_promotion_to(&StorageTier::Cold));
    }

    #[test]
    fn restore_state_machine_transitions() {
        use RestoreState::*;
        assert!(NotArchived.can_transition_to(&Requested));
        assert!(Requested.can_transition_to(&InProgress));
        assert!(InProgress.can_transition_to(&Restored));
        assert!(Restored.can_transition_to(&Expired));
        assert!(Failed.can_transition_to(&Requested));
        assert!(Expired.can_transition_to(&Requested));
        // invalid jumps
        assert!(!NotArchived.can_transition_to(&Restored));
        assert!(!Restored.can_transition_to(&InProgress));
        assert!(!Requested.can_transition_to(&Requested)); // self handled as no-op, not a transition
    }

    #[test]
    fn display_matches_cli_values() {
        assert_eq!(StorageTier::Warm.to_string(), "warm");
        assert_eq!(RetentionClass::Transient.to_string(), "transient");
    }

    #[test]
    fn days_for_maps_each_class() {
        let p = RetentionPolicy {
            diagnostic_days: 100,
            intermediate_days: 50,
            transient_days: 5,
            warm_days: 90,
        };
        assert_eq!(p.days_for(RetentionClass::Diagnostic), 100);
        assert_eq!(p.days_for(RetentionClass::Intermediate), 50);
        assert_eq!(p.days_for(RetentionClass::Transient), 5);
    }

    #[test]
    fn retain_until_adds_configured_days() {
        let p = RetentionPolicy {
            diagnostic_days: 10,
            intermediate_days: 1,
            transient_days: 1,
            warm_days: 90,
        };
        let from = Utc.with_ymd_and_hms(2025, 1, 1, 0, 0, 0).unwrap();
        let until = p.retain_until(RetentionClass::Diagnostic, from).unwrap();
        assert_eq!(until, Utc.with_ymd_and_hms(2025, 1, 11, 0, 0, 0).unwrap());
    }

    #[test]
    fn non_positive_duration_means_no_expiry() {
        let p = RetentionPolicy {
            diagnostic_days: 0,
            intermediate_days: -1,
            transient_days: 1,
            warm_days: 0,
        };
        let from = Utc.with_ymd_and_hms(2025, 1, 1, 0, 0, 0).unwrap();
        assert!(p.retain_until(RetentionClass::Diagnostic, from).is_none());
        assert!(p.retain_until(RetentionClass::Intermediate, from).is_none());
        assert!(p.retain_until(RetentionClass::Transient, from).is_some());
    }

    #[test]
    fn report_out_without_warm_targets_cold() {
        // 4-year diagnostic retention, anchored at the report-out date.
        let policy = RetentionPolicy {
            diagnostic_days: 365 * 4,
            intermediate_days: 365,
            transient_days: 30,
            warm_days: 90,
        };
        let reported_at = Utc.with_ymd_and_hms(2026, 6, 13, 0, 0, 0).unwrap();
        let transition =
            policy.report_out_transition(RetentionClass::Diagnostic, reported_at, false);

        assert_eq!(transition.reported_at, reported_at);
        assert_eq!(transition.target_tier, StorageTier::Cold);
        assert_eq!(transition.cold_move_at, None);
        assert_eq!(
            transition.retain_until,
            Some(reported_at + chrono::Duration::days(365 * 4))
        );
    }

    #[test]
    fn report_out_with_warm_targets_warm_then_schedules_cold() {
        let policy = RetentionPolicy {
            diagnostic_days: 365 * 4,
            intermediate_days: 365,
            transient_days: 30,
            warm_days: 90,
        };
        let reported_at = Utc.with_ymd_and_hms(2026, 6, 13, 0, 0, 0).unwrap();
        let transition =
            policy.report_out_transition(RetentionClass::Diagnostic, reported_at, true);

        // Re-inspection tier first, with the cold (S3) move scheduled after the dwell.
        assert_eq!(transition.target_tier, StorageTier::Warm);
        assert_eq!(
            transition.cold_move_at,
            Some(reported_at + chrono::Duration::days(90))
        );
        // Retention is independent of placement and still anchored at report-out.
        assert_eq!(
            transition.retain_until,
            Some(reported_at + chrono::Duration::days(365 * 4))
        );
    }

    #[test]
    fn report_out_with_warm_but_zero_dwell_targets_cold() {
        let policy = RetentionPolicy {
            diagnostic_days: 365 * 4,
            intermediate_days: 365,
            transient_days: 30,
            warm_days: 0,
        };
        let reported_at = Utc.with_ymd_and_hms(2026, 6, 13, 0, 0, 0).unwrap();
        let transition =
            policy.report_out_transition(RetentionClass::Diagnostic, reported_at, true);
        assert_eq!(transition.target_tier, StorageTier::Cold);
        assert_eq!(transition.cold_move_at, None);
    }
}
