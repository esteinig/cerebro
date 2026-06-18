# Maintenance guide

The cerebro-fs guarantees — integrity, durability, consistency, retention — are
enforced by background jobs the worker runs from a queue. This guide catalogues those
jobs, their cadence, and how to run them on demand. For the recovery procedures that
use them, see the [disaster-recovery runbook](disaster-recovery.md).

## How maintenance runs

The worker consumes jobs from a `maintenance` queue. Scheduled jobs are **seeded**
with stable identities so re-seeding only ever adds what is missing, and run on
staggered offsets so they do not all fire at once. Each scheduled **scan** is a
producer: it pages the catalogue and enqueues per-object work jobs, which lets a large
estate be swept incrementally within a per-run budget rather than all at once.

Jobs fall into three classes by how much they are trusted to act on their own:

- **Detection (read-only)** — find problems, record them, never mutate.
- **Safe automation** — act automatically, but only behind integrity and safety
  gates (hash checks, grace windows, confirmed-replacement checks).
- **Operator-gated** — never automatic; require an explicit confirmation and an
  explicit target list.

## Scheduled jobs

| Job | Cadence | Class | What it does |
|---|---|---|---|
| `retention_sweep` | daily | safe automation | Expires data whose retention has elapsed, honouring legal holds. |
| `tier_move_scan` | daily | safe automation | Finds files due for a tier change and enqueues `tier_move` for each. |
| `verify_scan` | daily | detection | Finds files due for re-verification and enqueues `verify_file`, rotating across the estate. |
| `purge_reclaim` | daily | safe automation | Reclaims space from purged/expired entries. |
| `restore_scan` | hourly | safe automation | Advances in-progress restores and enqueues `restore_drive`. |
| `catalogue_backup` | daily | safe automation | Dumps the catalogue to the backup store with a verifiable manifest. |
| `reconcile_scan` | weekly | detection | Compares the catalogue against the store; reports dangling references and orphan objects. |
| `verify_repair` | weekly | safe automation | Re-checks the reconcile findings and escalates confirmed loss; the integrity case is repaired by `verify_file` from a cold backup. |
| `archive_reclaim` | weekly | safe automation | Deletes the redundant local copy of an archived file past its grace window, once the cold copy is confirmed present. |

### Per-object work jobs

Enqueued by the scans above (or on demand):

- **`tier_move {file_id}`** — moves one file's bytes between tiers; archives to the
  cold store when moving to cold (if a cold store is configured).
- **`verify_file {file_id}`** — re-hashes one object against its catalogue hash; on a
  mismatch, repairs from the cold backup when one exists, otherwise alerts.
- **`restore_drive {file_id}`** — re-materialises an archived file from the cold store,
  hash-verified, and repoints the catalogue.

### Operator-gated

- **`reconcile_reclaim {confirm, keys}`** — deletes untracked store objects. Never
  seeded and never automatic: it requires `confirm: true` and an explicit list of
  keys, which it routes to the filer or volume delete as appropriate.

## Running a job on demand

Maintenance jobs can be triggered through the admin job endpoint (admin
authentication required). Enqueue a job:

```bash
curl -sS -X POST https://<host>/jobs/enqueue \
  -H "Authorization: Bearer $ADMIN_TOKEN" -H 'Content-Type: application/json' \
  -d '{"kind":"<job-kind>","args":{ ... },"queue":"maintenance"}'
```

Poll completion at `GET /jobs/enqueue/{id}`. Common uses:

- Force a backup now: `{"kind":"catalogue_backup","args":{}}`.
- Re-verify (and repair) a specific file:
  `{"kind":"verify_file","args":{"file_id":"<id>"}}`.
- Get a current consistency picture: `{"kind":"reconcile_scan","args":{}}`.
- Reclaim confirmed orphans (gated):
  `{"kind":"reconcile_reclaim","args":{"confirm":true,"keys":["<key>"]}}`.

## Where findings go

- **Reconcile reports** are written to the backup store under a `reconcile/` prefix;
  read the latest before acting on divergence.
- **Backups** land under the backup prefix as `{prefix}/{backup_id}/` with the
  archive and a `manifest.json` (hash + audit-chain-verified).
- **Integrity failures and confirmed losses** are surfaced on the worker metrics and
  in error-level logs for alerting.

## During an incident

The [disaster-recovery runbook](disaster-recovery.md#3-first-response--the-first-15-minutes)
opens by pausing the destructive maintenance jobs (`reconcile_reclaim`,
`archive_reclaim`, `tier_move`) so recovery never races deletion — the fastest safe
way is to stop the worker. Re-enable it once the scenario's verification passes.
