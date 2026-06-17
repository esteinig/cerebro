# Verify-repair & restore-from-cold (Stage 4, S4-5)

S4-5 closes two loops left open earlier: it makes **restore from the cold store
real** (the direction deferred from S4-4), and it **consumes the S4-3 reconcile
reports** to confirm and escalate data loss.

## Restore from cold (closes the S4-4 loop)

The archival restore state machine (S3-3b: `NotArchived → Requested → InProgress →
Restored`) tracked progress but never actually brought bytes back — the
`InProgress → Restored` step just flipped state. S4-5 makes it real: when a restore
is ready, `restore_drive` re-materialises the object before declaring it
retrievable.

1. Fetch the bytes from the cold store at the file's `archive_key`
   (`ObjectStore::get`).
2. Verify them against the catalogue BLAKE3 (`hash_bytes`) — a corrupt cold copy
   can never silently replace a good catalogue entry.
3. Re-upload to SeaweedFS for a new fid (`write_object`, via the normal `weed
   upload` path, so the restored object lands replicated like any other).
4. Repoint the catalogue through the relocate endpoint (D3): `archived = false`,
   the new `fid`, `archive_key` cleared, landed on the Warm tier — all under a
   compare-and-set on the file's current (Cold) tier.

If any step fails, the file is marked `Failed` rather than declared restored, so a
non-retrievable object is never reported as available. With no cold store
configured the step is a no-op and the existing simulation seam
(`CEREBRO_RESTORE_SIMULATE_SECONDS`) drives dev/test as before.

`restore_object` (engine) and `write_object` (FS client) are the new primitives;
`RestoreDrive` gains the cold-store settings and the `rematerialize` step.

## Verify-repair (consumes reconcile reports)

`reconcile_scan` (S4-3) detects dangling references and writes a report,
report-first. `verify_repair` is the weekly pass that acts on it. For each dangling
reference it re-checks live state, which collapses to three outcomes:

- **Recovered** — the object reappeared since the scan (a transient miss, or
  SeaweedFS volume replication healed it). No action.
- **Archived since** — the file was tier-moved to cold after the scan, so its
  object is now *correctly* absent locally. No action.
- **Confirmed loss** — still non-archived and still missing. Escalated loudly
  (error-level, with the fid and tier) for operator action.

It mutates nothing. Its value is turning a passive weekly report into a confirmed,
transient-filtered data-loss signal — a clinical estate wants genuine loss
escalated, not buried in a report. Staggered well after `reconcile_scan` so it
reads a fresh report.

## Why confirmed loss has no automatic repair

A reconcile dangling reference is a *non-archived* file whose object is missing.
There is no in-system recovery source for it: the S4-1 replica is volume-level
within the same SeaweedFS fid (a 404 means replication did not save it), and the
S4-2 backups cover the *catalogue*, not objects. So confirmed loss is escalated,
not auto-repaired. The two genuine repair paths are elsewhere and remain active:

- **Lost redundancy** (a copy survives, replica count low) → SeaweedFS volume
  replication re-replicates automatically. Restoring redundancy is safe/auto (D6).
- **Archived objects** (the cold copy is authoritative) → the `restore_drive` loop
  above re-materialises on demand.

A natural follow-on is to feed **verify-scan** hash-mismatch findings (S3-3a) for
*archived* files into the same restore path — a corrupt local copy of an archived
file is repairable from cold — closing verify→repair for the integrity case too.

## What's tested here vs in your environment

`cargo test` covers the engine's pure pieces (`archive_key_for`; the restore hash
gate is exercised through `restore_object`'s verification branch). The byte
movement (`write_object` via `weed upload`, `read_object`, `ObjectStore` get/put),
the relocate repoint, the restore state-machine transitions, and the report
consumption are validated against a running stack. Correct-by-inspection — run
`cargo check` / `cargo test` after applying.
