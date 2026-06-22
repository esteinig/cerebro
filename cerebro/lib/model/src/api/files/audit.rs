//! Tamper-evident audit trail / chain-of-custody (S2-6).
//!
//! Every lifecycle-affecting action on an artefact (upload, report-out, tier
//! move, restore, expiry, legal-hold change, tag change, delete) appends an
//! [`AuditEvent`] to a per-team, append-only chain. Each event is **sealed** with
//! `BLAKE3` over its canonical body — which includes the previous event's hash —
//! so the events form a hash chain: altering, reordering or removing any past
//! event breaks every subsequent hash and is detectable by [`verify_chain`].
//!
//! This realises the project decision that *Cerebro* owns artefact lifecycle and
//! that the lifecycle is auditable. The seal provides integrity and ordering
//! (a verifiable chain), not non-repudiation; cryptographic signing of events
//! would be an additive follow-on (as for the run manifest).

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use uuid::Uuid;

/// Genesis predecessor hash for the first event in a chain.
pub const AUDIT_GENESIS_HASH: &str =
    "0000000000000000000000000000000000000000000000000000000000000000";

/// Failure appending an event to the audit chain (S2-8).
///
/// Surfaced to the triggering operation when auditing is fail-closed (the
/// default); swallowed and logged when `CEREBRO_AUDIT_FAIL_OPEN` is set.
#[derive(Debug, thiserror::Error)]
pub enum AuditError {
    #[error("failed to read audit chain tip: {0}")]
    TipLookup(String),
    #[error("failed to append audit event: {0}")]
    Insert(String),
    #[error("audit append exhausted retries under sequence contention")]
    RetriesExhausted,
}

/// A lifecycle-affecting action recorded in the audit trail.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum AuditAction {
    Upload,
    ReportOut,
    TierMove,
    Restore,
    Expiry,
    LegalHoldChange,
    TagChange,
    Delete,
}

/// The authenticated principal that triggered an action (or the system).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct AuditActor {
    pub id: String,
    pub email: String,
}
impl AuditActor {
    pub fn system() -> Self {
        Self {
            id: String::from("system"),
            email: String::from("system"),
        }
    }
}

/// A single, sealed audit-trail entry.
///
/// All fields except `hash` are covered by the seal (`hash` is the seal itself).
/// `prev_hash` links to the previous event, forming the chain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AuditEvent {
    pub id: String,
    /// Monotonic per-team sequence number (chain position).
    pub sequence: u64,
    pub timestamp: DateTime<Utc>,
    pub action: AuditAction,
    pub file_id: Option<String>,
    pub run_id: Option<String>,
    pub sample_id: Option<String>,
    pub actor: AuditActor,
    /// Human-readable detail (e.g. `"legal_hold: false -> true"`).
    pub detail: String,
    /// Hash of the previous event (or [`AUDIT_GENESIS_HASH`] for the first).
    pub prev_hash: String,
    /// BLAKE3 seal over the canonical body. Empty until [`AuditEvent::seal`].
    pub hash: String,
}

impl AuditEvent {
    /// Build an unsealed event. Call [`AuditEvent::seal`] before persisting.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        sequence: u64,
        timestamp: DateTime<Utc>,
        action: AuditAction,
        file_id: Option<String>,
        run_id: Option<String>,
        sample_id: Option<String>,
        actor: AuditActor,
        detail: impl Into<String>,
        prev_hash: String,
    ) -> Self {
        Self {
            id: Uuid::new_v4().to_string(),
            sequence,
            timestamp,
            action,
            file_id,
            run_id,
            sample_id,
            actor,
            detail: detail.into(),
            prev_hash,
            hash: String::new(),
        }
    }

    /// Canonical bytes the seal covers: every field except `hash`.
    ///
    /// Deterministic (struct field order is stable; `prev_hash` is included, so
    /// the seal also covers the chain link).
    pub fn canonical_body(&self) -> Vec<u8> {
        let mut body = self.clone();
        body.hash = String::new();
        serde_json::to_vec(&body).unwrap_or_default()
    }

    /// The BLAKE3 seal this event should carry given its body and `prev_hash`.
    pub fn compute_hash(&self) -> String {
        let mut hasher = blake3::Hasher::new();
        hasher.update(&self.canonical_body());
        hasher.finalize().to_hex().to_string()
    }

    /// Seal the event by setting `hash` to the computed value.
    pub fn seal(&mut self) {
        self.hash = self.compute_hash();
    }

    /// Whether this event's own seal is intact (body not altered).
    pub fn is_sealed_valid(&self) -> bool {
        !self.hash.is_empty() && self.hash == self.compute_hash()
    }
}

/// Verify an ordered slice of events for internal chain integrity.
///
/// Checks, for the whole slice: every event's own seal is valid, the sequence
/// numbers strictly increment by one, and each event's `prev_hash` equals the
/// previous event's `hash`. When the slice starts at sequence 0 the first event's
/// `prev_hash` must be [`AUDIT_GENESIS_HASH`]. An empty slice is trivially valid.
pub fn verify_chain(events: &[AuditEvent]) -> bool {
    let mut prev: Option<&AuditEvent> = None;
    for event in events {
        if !event.is_sealed_valid() {
            return false;
        }
        match prev {
            None => {
                if event.sequence == 0 && event.prev_hash != AUDIT_GENESIS_HASH {
                    return false;
                }
            }
            Some(previous) => {
                if event.sequence != previous.sequence + 1 {
                    return false;
                }
                if event.prev_hash != previous.hash {
                    return false;
                }
            }
        }
        prev = Some(event);
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ev(seq: u64, prev_hash: String, detail: &str) -> AuditEvent {
        let mut e = AuditEvent::new(
            seq,
            Utc::now(),
            AuditAction::Upload,
            Some(format!("file{seq}")),
            Some("RUN".into()),
            None,
            AuditActor::system(),
            detail,
            prev_hash,
        );
        e.seal();
        e
    }

    fn chain(n: u64) -> Vec<AuditEvent> {
        let mut events = Vec::new();
        let mut prev = AUDIT_GENESIS_HASH.to_string();
        for seq in 0..n {
            let e = ev(seq, prev.clone(), "x");
            prev = e.hash.clone();
            events.push(e);
        }
        events
    }

    #[test]
    fn seal_then_verify() {
        let e = ev(0, AUDIT_GENESIS_HASH.into(), "first");
        assert!(e.is_sealed_valid());
    }

    #[test]
    fn valid_chain_verifies() {
        assert!(verify_chain(&chain(5)));
        assert!(verify_chain(&[]));
    }

    #[test]
    fn tampering_a_body_breaks_the_chain() {
        let mut events = chain(4);
        events[2].detail = "tampered".into(); // body changed, hash not recomputed
        assert!(!verify_chain(&events));
    }

    #[test]
    fn reordering_breaks_the_chain() {
        let mut events = chain(4);
        events.swap(1, 2);
        assert!(!verify_chain(&events));
    }

    #[test]
    fn broken_link_is_detected() {
        let mut events = chain(3);
        events[1].prev_hash = "deadbeef".into();
        events[1].seal(); // reseal so its own hash is valid, but link is wrong
        assert!(!verify_chain(&events));
    }
}
