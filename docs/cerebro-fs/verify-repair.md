# Verify-repair & restore-from-cold

This closes two integrity loops: it makes **restore from the cold store
real**, and it **consumes the reconcile
reports** to confirm and escalate data loss.

## Restore from cold

The archival restore state machine (`NotArchived → Requested → InProgress →
Restored`) tracked progress but never actually brought bytes back — the
`InProgress → Restored` step just flipped state. The restore pass makes it real: when a restore
is ready, `restore_drive` re-materialises the object before declaring it
retrievable.

1. Fetch the bytes from the cold store at the file's `archive_key`
   (`ObjectStore::get`).
2. Verify them against the catalogue BLAKE3 (`hash_bytes`) — a corrupt cold copy
   can never silently replace a good catalogue entry.
3. Write the bytes back to the file's **effective location**: a path-addressed
   filer object is overwritten in place (its path stays valid), otherwise a fresh
   replicated weed object is written and its fid returned.
4. Repoint the catalogue through the relocate endpoint: `archived = false`,
   the new `fid` (only when one was written), landed on the Warm tier — all under a
   compare-and-set on the file's current (Cold) tier. The `archive_key` is
   **retained** so the restored file keeps a durable cold backup.

If any step fails, the file is marked `Failed` rather than declared restored, so a
non-retrievable object is never reported as available. With no cold store
configured the step is a no-op and the existing simulation seam
(`CEREBRO_RESTORE_SIMULATE_SECONDS`) drives dev/test as before.

`restore_object` (engine), `write_object`, and `write_object_at_path` (FS client)
are the primitives; `RestoreDrive` gains the cold-store settings and the
`rematerialize` step.

## Verify-repair (consumes reconcile reports)

`reconcile_scan` detects dangling references and writes a report,
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
There is no in-system recovery source for it: the replica is volume-level
within the same SeaweedFS fid (a 404 means replication did not save it), and the
catalogue backups cover the *catalogue*, not objects. So confirmed loss is escalated,
not auto-repaired. The two genuine repair paths are elsewhere and remain active:

- **Lost redundancy** (a copy survives, replica count low) → SeaweedFS volume
  replication re-replicates automatically. Restoring redundancy is safe/auto.
- **Archived objects** (the cold copy is authoritative) → the `restore_drive` loop
  above re-materialises on demand.

## Integrity repair from cold

The dangling-reference case above has no in-system recovery source. The **integrity**
case does: a non-archived file that fails its BLAKE3 verify but has a retained cold
backup (`archive_key`). This closes verify→repair for it. When `verify_file`
finds a mismatch, before alerting it tries `repair_from_cold`:

- re-pull the bytes from cold and re-verify them against the catalogue hash (a
  corrupt cold copy can never overwrite the live entry);
- write them back to the file's effective location (in-place filer overwrite, or a
  fresh weed fid);
- repoint the catalogue, keeping the file live, on its current tier, with the cold
  backup retained; and
- if a fresh fid replaced the corrupt object, delete the now-orphaned old bytes.

Only a file with no cold backup, or a repair that errors, is escalated (metric +
error log) for operator action. The retained `archive_key` (above) is what makes a
restored file repairable on a later mismatch — the cold copy stays a durable backup
pointer for the file's whole life, not just until first restore.

The gate is presence + hash: repair proceeds only if the cold bytes match the
catalogue hash, so a cold copy that is itself corrupt fails closed (escalation),
never a bad overwrite.

## Running verify and repair on demand

Both passes are scheduled, but you can force them (admin token — see
[administration → operator CLI setup](administration.md#operator-cli-setup-authentication)):

```bash
# Re-hash one object against its catalogue hash; repairs from cold on mismatch.
id=$(cerebro-client jobs launch-job --kind verify_file \
  --args '{"file_id":"<FILE_ID>"}' --queue maintenance)
cerebro-client jobs job-status --id "$id"

# Re-check the latest reconcile report and escalate confirmed losses.
cerebro-client jobs launch-job --kind verify_repair --args '{}' --queue maintenance
```

Read the outcome in the worker logs and metrics: a success logs `integrity verified`
(or `integrity failure REPAIRED from cold backup`) and stamps `verified_at`; an
unrepairable mismatch logs `INTEGRITY FAILURE …` and increments the verify-failure
series of `cerebro_file_lifecycle_ops_total`. See
[disaster recovery DR-4](disaster-recovery.md#dr-4--object-integrity-failure-bit-rot)
for the decision tree on each outcome.

## What's tested here vs in your environment

`cargo test` covers the engine's pure pieces (`archive_key_for`; the restore hash
gate is exercised through `restore_object`'s verification branch). The byte
movement (`write_object` via `weed upload`, `read_object`, `ObjectStore` get/put),
the relocate repoint, the restore state-machine transitions, and the report
consumption are validated against a running stack. Correct-by-inspection — run
`cargo check` / `cargo test` after applying.
