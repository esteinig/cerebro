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

Maintenance jobs are enqueued through the admin job API. The operator path is the
`cerebro-client` CLI; the raw endpoint is shown alongside for scripting. Both require
an **admin** token — set one up first (see
[administration → operator CLI setup](administration.md#operator-cli-setup-authentication)):

```bash
export CEREBRO_API_URL="https://api.<your-domain>"
export CEREBRO_API_TOKEN="$(cerebro-client login --email admin@cerebro --password '****')"
```

**Enqueue (CLI).** `launch-job` prints the job id (a UUID) on success:

```bash
cerebro-client jobs launch-job \
  --kind <job-kind> --args '<json-args>' --queue maintenance
# → 3f9a1c54-…   (job id; poll it below)
```

**Enqueue (raw API).** The CLI wraps this exact call; the server runs the job on
`cerebro-worker` and returns the id in `data`:

```bash
curl -sS -X POST "$CEREBRO_API_URL/jobs/enqueue" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" -H 'Content-Type: application/json' \
  -d '{"kind":"<job-kind>","args":{ },"queue":"maintenance"}'
# → {"status":"success","message":"Job enqueued","data":"3f9a1c54-…"}
```

**Poll completion.** Use the returned id (UUID or Faktory JID):

```bash
cerebro-client jobs job-status   --id 3f9a1c54-…   # completed=true status=succeeded result=…
cerebro-client jobs job-complete --id 3f9a1c54-…   # yes | no
# raw equivalent: GET $CEREBRO_API_URL/jobs/enqueue/3f9a1c54-…
```

### Common on-demand jobs

| Intent | Command |
|---|---|
| Force a catalogue backup now | `cerebro-client jobs launch-job --kind catalogue_backup --args '{}' --queue maintenance` |
| Re-verify (and repair) one file | `cerebro-client jobs launch-job --kind verify_file --args '{"file_id":"<id>"}' --queue maintenance` |
| Current consistency picture | `cerebro-client jobs launch-job --kind reconcile_scan --args '{}' --queue maintenance` |
| Reclaim confirmed orphans (**gated**) | `cerebro-client jobs launch-job --kind reconcile_reclaim --args '{"confirm":true,"keys":["<key>"],"max_delete":100}' --queue maintenance` |
| Tier-move / archive one file | `cerebro-client jobs launch-job --kind tier_move --args '{"file_id":"<id>","target":"Cold"}' --queue maintenance` |
| Reclaim local copies of archived files | `cerebro-client jobs launch-job --kind archive_reclaim --args '{}' --queue maintenance` |

Job-kind argument shapes are documented with each capability: `catalogue_backup`
([backup](catalogue-backup.md)), `verify_file`/`verify_repair`
([verify & repair](verify-repair.md)), `reconcile_scan`/`reconcile_reclaim`
([reconcile](consistency-reconcile.md)), `tier_move`/`archive_reclaim`
([archival](archival.md)), restore ([disaster recovery](disaster-recovery.md#restore-an-archived-file)).

### Scheduling and re-seeding

The seeded schedule (above) is created automatically. To add or override a schedule
on demand — e.g. a more frequent backup, or to re-seed after a fresh database — use
`schedule-job` (one-shot when `--interval-seconds` is omitted, recurring otherwise)
and inspect the schedule with `schedule-summary`:

```bash
# A recurring backup every 12h, first run at an explicit UTC time.
cerebro-client jobs schedule-job \
  --kind catalogue_backup --args '{}' --queue maintenance \
  --run-at 2026-06-24T15:00:00Z --interval-seconds 43200

cerebro-client jobs schedule-summary    # next/previous 10 scheduled runs
```

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
way is to stop the worker:

```bash
docker compose stop cerebro-worker     # halts job execution; API + data keep serving
# … recover …
docker compose up -d cerebro-worker    # resume once the scenario's Verify step passes
```

Producers (the seeded scans) keep queuing jobs harmlessly while the worker is stopped;
nothing executes until you restart it.
