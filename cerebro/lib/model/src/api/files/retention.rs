//! Storage tiering and retention policy types for Cerebro FS.
//!
//! These types describe *where* a file currently lives ([`StorageTier`]) and
//! *how long* it must be kept ([`RetentionClass`] resolved through a
//! [`RetentionPolicy`]). They are the data-model foundation that the deployment
//! models (FS-3 single-server hot/cold, FS-4 HPC + S3 cold) and the lifecycle
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
//! accreditation scope (FS-5 surfaces them as deployment configuration).
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
}

impl Default for RetentionPolicy {
    /// Placeholder durations — **not legal advice**. Confirm and override these
    /// against your laboratory's NATA/NPAAC accreditation scope before relying
    /// on them. They are deliberately conservative (keep rather than delete).
    fn default() -> Self {
        Self {
            diagnostic_days: 365 * 7, // ~7 years, placeholder only
            intermediate_days: 365,   // ~1 year, placeholder only
            transient_days: 30,       // ~1 month, placeholder only
        }
    }
}

impl RetentionPolicy {
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
    /// let policy = RetentionPolicy { diagnostic_days: 10, intermediate_days: 0, transient_days: 1 };
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
    /// The object is directly retrievable; no restore is needed.
    NotRequired,
    /// A restore has been requested and is in progress.
    Pending,
    /// A restore has completed and the object is temporarily retrievable.
    Available,
    /// The most recent restore attempt failed.
    Failed,
}

impl Default for RestoreState {
    fn default() -> Self {
        RestoreState::NotRequired
    }
}

impl std::fmt::Display for RestoreState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RestoreState::NotRequired => write!(f, "not-required"),
            RestoreState::Pending => write!(f, "pending"),
            RestoreState::Available => write!(f, "available"),
            RestoreState::Failed => write!(f, "failed"),
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
    fn display_matches_cli_values() {
        assert_eq!(StorageTier::Warm.to_string(), "warm");
        assert_eq!(RetentionClass::Transient.to_string(), "transient");
    }

    #[test]
    fn days_for_maps_each_class() {
        let p = RetentionPolicy { diagnostic_days: 100, intermediate_days: 50, transient_days: 5 };
        assert_eq!(p.days_for(RetentionClass::Diagnostic), 100);
        assert_eq!(p.days_for(RetentionClass::Intermediate), 50);
        assert_eq!(p.days_for(RetentionClass::Transient), 5);
    }

    #[test]
    fn retain_until_adds_configured_days() {
        let p = RetentionPolicy { diagnostic_days: 10, intermediate_days: 1, transient_days: 1 };
        let from = Utc.with_ymd_and_hms(2025, 1, 1, 0, 0, 0).unwrap();
        let until = p.retain_until(RetentionClass::Diagnostic, from).unwrap();
        assert_eq!(until, Utc.with_ymd_and_hms(2025, 1, 11, 0, 0, 0).unwrap());
    }

    #[test]
    fn non_positive_duration_means_no_expiry() {
        let p = RetentionPolicy { diagnostic_days: 0, intermediate_days: -1, transient_days: 1 };
        let from = Utc.with_ymd_and_hms(2025, 1, 1, 0, 0, 0).unwrap();
        assert!(p.retain_until(RetentionClass::Diagnostic, from).is_none());
        assert!(p.retain_until(RetentionClass::Intermediate, from).is_none());
        assert!(p.retain_until(RetentionClass::Transient, from).is_some());
    }
}
