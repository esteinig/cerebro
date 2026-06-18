# cerebro-fs — production readiness

**Current status: a strong, safe-by-design foundation that is pre-production for
authoritative clinical use.**

cerebro-fs is well past prototype — its architecture, its gating of destructive
actions, and its hash-anchored integrity model are production-grade. This document is
an honest, deliberately critical account of what still stands between the system and
holding the *only* copy of clinical data, written as a tracked roadmap. For what the
system does, see the [overview](overview.md) and [architecture](architecture.md); for
what is currently tested, see the [validation guide](validation.md).

## Assessment

cerebro-fs is suitable **today** in a **non-authoritative** role — running alongside a
managed, replicated store that remains the source of truth while cerebro-fs proves
itself. It is **not yet** recommended as the authoritative durability layer for
patient data until the go-live blockers below are closed. The three decisive gaps are
single-host durability for hot objects, the absence of end-to-end and recovery
testing, and the catalogue's high-availability / recovery-point story. A
hash-verified-reclaim fix and encryption at rest are the sharp secondary edges.

## How to read this

Items have stable IDs (`PRD-n`) and are grouped by priority tier. Severity reflects
worst-case data or compliance impact, not likelihood. Check an item off when its
recommended action has shipped **and** been validated in a running environment.

---

## Go-live blockers — required before authoritative clinical use

### PRD-1 — Multi-host replication + a durability invariant — *severity: critical*

- [ ] Multi-host replication for primary object data
- [ ] Durability-invariant checker (alert/refuse on unprotected objects and shared failure domains)

**Risk.** Non-archived primary objects are replicated only *within a single host's*
volumes; a host or array loss destroys every replica of hot data. Separately, the
requirement that the backup and cold stores be independently durable is an operator
assumption, not an enforced invariant — a misconfiguration (e.g. backup store on the
same disk as SeaweedFS) can silently create unrecoverable data.

**Action.** Add multi-host replication. Add a checker that alerts or refuses when any
object is neither replicated to the configured floor nor archived, and that flags a
backup or cold store sharing a failure domain with the primary.

**References.** [`fs-replication.md`](fs-replication.md) (multi-host upgrade),
[`disaster-recovery.md`](disaster-recovery.md) (durability matrix); code:
`cerebro/src/stack/deploy.rs` (replication configuration), `fs/`.

### PRD-2 — End-to-end and recovery testing — *severity: critical*

- [ ] Integration tests over the real I/O paths
- [ ] Fault-injection tests (lost volume, corrupt object, dropped catalogue)
- [ ] DR drills executed against a real stack, with timings recorded

**Risk.** Coverage today is unit tests over pure logic; the paths where production
actually fails are untested end-to-end — SeaweedFS I/O, `mongodump`/`mongorestore`,
worker↔API auth, the Faktory job lifecycle, and especially the real restore (cold
fetch → re-upload → relocate). A latent restore bug (a stale filer path after write)
was caught by inspection during hardening, not by a test, which is precisely the class
of defect that may still hide in the untested wiring. The disaster-recovery "drills"
are documented procedures, not executed proof.

**Action.** Build an integration + fault-injection suite, ideally in CI: kill a
volume, corrupt an object, drop the catalogue, and assert the documented recovery
restores byte-identical, hash-verified data. Turn the runbook drills into executed,
timed exercises.

**References.** [`validation.md`](validation.md),
[`disaster-recovery.md`](disaster-recovery.md) (drills); code: `worker/runners/`,
`worker/{backup,archive,reconcile}.rs`.

### PRD-3 — Catalogue resilience (HA + recovery point) — *severity: critical*

- [ ] MongoDB as a replica set
- [ ] Point-in-time recovery (oplog) for the catalogue
- [ ] Defined, alerting response to a failed audit-chain verification

**Risk.** The catalogue (MongoDB) is the single source of truth for identity,
location, hashes, retention, and the audit chain. Backups are a daily full
`mongodump`, so the catalogue recovery point is ~24 hours — up to a day of
provenance, retention, legal-hold, and audit mutations can be lost on a restore.
Whether MongoDB runs as a replica set is currently out of scope, so a single database
failure may be a full outage and inter-backup corruption is data loss. The audit-chain
check at backup time records a result but its failure response is undefined.

**Action.** Run MongoDB as a replica set and add oplog-based point-in-time recovery to
push the catalogue recovery point from ~24 h toward minutes. Define and alert on a
failed audit-chain verification (a broken append-only chain is a compliance event).

**References.** [`catalogue-backup.md`](catalogue-backup.md),
[`architecture.md`](architecture.md) (durability matrix); code: `worker/backup.rs`,
`worker/runners/catalogue_backup.rs`.

---

## Before authoritative use

### PRD-4 — Hash-verified reclaim — *severity: high*

- [ ] Verify the cold copy's integrity (not just presence) before deleting the last local copy

**Risk.** Local-copy reclaim deletes the *last local copy* of an archived object once
the cold copy is confirmed present and the grace window has passed. Presence is not
integrity: a truncated, zero-byte, or partially written cold object can pass an
existence check, after which the only good copy has been deleted on a false positive.

**Action.** Before reclaiming a local copy, require that the cold copy hash-verifies
(or was successfully verified recently), not merely that it exists.

**References.** [`archival.md`](archival.md) (local-copy reclamation),
[`verify-repair.md`](verify-repair.md); code: `worker/runners/archive_reclaim.rs`
(`is_reclaim_candidate`), `fs/client.rs` (`store_object_present`).

### PRD-5 — Observability and alerting — *severity: high*

- [ ] SLIs/SLOs for integrity, durability, and job health
- [ ] Shipped alert rules and a dashboard
- [ ] Reconcile and integrity findings routed to paging, not just to files

**Risk.** The safety model leans on operators reading reports and approving cleanups,
but reconcile findings land as files under a `reconcile/` prefix — files are not
alerts. If divergence is not actively surfaced, the operator gate degrades into either
unbounded orphan accumulation or rubber-stamping. Safe-automation jobs delete data
after grace and run unattended, so their health must be a guaranteed signal, not an
assumed practice.

**Action.** Define SLIs/SLOs; ship alert rules and a dashboard; route
reconcile/integrity findings to paging. Make safe-automation job health first-class.

**References.** [`maintenance.md`](maintenance.md) (where findings go),
[`administration.md`](administration.md) (monitoring); code: `worker/telemetry.rs`,
`worker/runners/reconcile.rs`, `server/api/telemetry.rs`.

### PRD-6 — Security hardening — *severity: high*

- [ ] Encryption at rest for object, cold, and backup stores
- [ ] Four-eyes (or strong audit) for mass-deletion operations
- [ ] Remove the invalid-TLS escape hatch from production paths

**Risk.** Encryption at rest is operator advice, not a built-in — below par for
PHI-adjacent data. The most destructive operation (`reconcile_reclaim` with an
operator-supplied key list) sits behind a single admin token with no approval step and
no confirmed infrastructure-level audit of who triggered a mass deletion. The
`CEREBRO_DANGER_ACCEPT_INVALID_TLS_CERTIFICATE` flag is a foot-gun that can leak into
production.

**Action.** Add encryption at rest across stores; require a second approval (or at
minimum a tamper-evident audit) for mass deletion; ensure the invalid-TLS flag cannot
be set on a production deployment.

**References.** [`administration.md`](administration.md) (secrets, configuration),
[`catalogue-backup.md`](catalogue-backup.md) (backup user); code:
`cerebro/src/stack/deploy.rs`, `server/api/files/handler.rs`, `server/api/jobs/`.

---

## Fast-followers — can run in parallel

### PRD-7 — Idempotency under redelivery — *severity: medium-high*

- [ ] Prove destructive jobs are safe under at-least-once redelivery
- [ ] Define read-path behaviour during an in-progress restore

**Risk.** Faktory delivers at least once and the poll-by-re-enqueue pattern leans into
that, so every destructive job must be safe to run twice, including after a partial
failure. A redelivered reclaim after a partial archive is where a double-execution
data-loss bug would hide. Separately, the behaviour of a read against a file that is
mid-restore is unspecified.

**Action.** Establish and test idempotency for tier-move, archive, and reclaim under
redelivery; specify and test the read path during restore.

**References.** [`maintenance.md`](maintenance.md) (per-object jobs); code:
`worker/context.rs` (re-enqueue), `worker/runners/tier_move.rs`, `worker/archive.rs`.

### PRD-8 — Scale validation — *severity: medium*

- [ ] Measure verify cadence vs. estate size against an integrity SLA
- [ ] Measure reconcile coverage vs. cadence × budget
- [ ] Measure catalogue backup time/size growth; consider incremental backup

**Risk.** Metagenomics implies multi-GB files and a growing object count. A rotating
verify only protects you if it completes within the integrity SLA; at scale it may
fall behind. Weekly budgeted reconcile covers the estate only if cadence × budget ≥
size, or orphans persist for many cycles. Full daily catalogue dumps grow unbounded.

**Action.** Measure these at representative scale; move to incremental catalogue backup
if full-dump cost climbs; tune cadences/budgets to meet stated SLAs.

**References.** [`validation.md`](validation.md), [`maintenance.md`](maintenance.md);
code: `worker/runners/verify.rs`, `worker/runners/reconcile.rs`, `worker/backup.rs`.

---

## Already solid — the baseline this builds on

- Destructive actions are gated on grace windows, presence checks, or explicit
  operator confirmation; detection is report-first and never mutates.
- Integrity is hash-anchored (BLAKE3) at every boundary; a corrupt copy cannot
  silently overwrite a good record.
- Backups carry a verifiable, checksummed, audit-chain-aware manifest.
- The cold archive store is independent of the primary object store by design.
- The core decision logic is cleanly factored and unit-tested.

## To confirm in code

Two items above are drawn from the system's design and are worth a short confirmation
before acting on them:

- **PRD-4** — confirm whether `is_reclaim_candidate` gates on cold-copy *presence*
  versus a hash or a recent successful verify. Code: `worker/runners/archive_reclaim.rs`.
- **PRD-6** — confirm whether admin and destructive actions (`reconcile_reclaim`,
  relocate, enqueue) are recorded in an infrastructure-level audit. Code:
  `server/api/files/handler.rs`, `server/api/jobs/`.
