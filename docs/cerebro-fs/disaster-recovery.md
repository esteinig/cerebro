# Cerebro disaster-recovery runbook

Operator procedures for detecting, containing, and recovering from data-loss and
integrity incidents on a Cerebro deployment. This is the **top-level** recovery
document; it orchestrates the machinery described in detail in:

- `fs-replication.md` — SeaweedFS single-host replication
- `catalogue-backup.md` — scheduled MongoDB backups to the cold store
- `consistency-reconcile.md` — dangling/orphan detection + gated reclaim
- `archival.md` — real cold-tier archival + local-copy reclaim
- `verify-repair.md` — restore-from-cold + integrity repair

> **Audience.** A Cerebro operator with shell access to the deployment host, the
> `docker compose` project, an admin API token, and read access to the backup and
> cold-store locations. Recovery touches clinical data — work deliberately and
> record what you do.

> **Golden rule.** Recovery is not complete when a command exits `0`. It is complete
> when the **Verify** step of the scenario passes. Every scenario below ends in one.

---

## 1. Severity triage — find your scenario

| Symptom | Likely class | Go to |
|---|---|---|
| A volume/disk failed; API still serves; some reads slow | One object copy lost | [DR-1](#dr-1--one-object-copy-lost-volume--disk-failure) |
| Catalogue queries fail / return wrong data; Mongo down or corrupt | Catalogue loss | [DR-2](#dr-2--catalogue-mongodb-loss-or-corruption) |
| Uploads/downloads by path fail; filer errors; listings empty | Filer metadata loss | [DR-3](#dr-3--filer-metadata-loss) |
| `verify` reports `INTEGRITY FAILURE`; a file's bytes don't match its hash | Object corruption | [DR-4](#dr-4--object-integrity-failure-bit-rot) |
| Reconcile reports dangling refs or orphans | Catalogue↔store divergence | [DR-5](#dr-5--cataloguestore-divergence) |
| Cold store unreachable/lost; restores fail; archived files unavailable | Cold-store loss | [DR-6](#dr-6--cold-object-store-loss) |
| Host is gone; rebuilding from scratch | Full rebuild | [DR-7](#dr-7--full-host-rebuild) |
| A reclaim/relocate ran wrongly and deleted/repointed data | Bad destructive job | [DR-8](#dr-8--bad-destructive-job) |

If you are not sure which class you are in, do [§3 First response](#3-first-response--the-first-15-minutes) first; it is safe in all cases.

---

## 2. What can be recovered, and from where

Every class of state has exactly one authoritative recovery source. Know yours
*before* an incident.

| Data class | Lives in | Recovery source | RPO (data age on recovery) | If the source is also gone |
|---|---|---|---|---|
| **Catalogue** (file records, lifecycle, audit chain) | MongoDB `cerebro-database` | Latest verified `mongodump` in the cold store | The backup interval | Catalogue is unrecoverable — see §3 prevention |
| **Filer metadata** (path → fid map) | MongoDB (filer store) | Same whole-instance dump, **if** the filer DB is co-located in `cerebro-database`; else SeaweedFS filer's own store | The backup interval | Rebuild path map by re-indexing volumes (vendor procedure) |
| **Live object bytes** (non-archived) | SeaweedFS volumes | The **replica** copy | None (current) | If archived: the cold store. If not: **lost** |
| **Archived object bytes** | The **cold store** (NFS/disk/S3) | The cold store *is* the backup | None (current) | Unrecoverable — the cold store must be independently durable |
| **Audit chain integrity** | MongoDB | The dump; its manifest records `audit_chain_verified` | The backup interval | A restore cannot launder a broken chain (manifest gate) |

Two consequences worth internalising:

- **A non-archived, single-replica object that loses both copies has no recovery
  source.** Replication and archival cadence are what keep objects out
  of this category. Prevention lives in those two knobs.
- **The cold store is a real backup target and must be durable on its own** (RAID +
  offsite/object-lock, or S3 with versioning). Cerebro does not replicate the cold
  store for you.

Set your targets in [§9 Objectives worksheet](#9-recovery-objectives-worksheet).

---

## 3. First response — the first 15 minutes

Do this on **any** suspected data-loss or integrity incident, before recovery.

1. **Declare and timestamp.** Note the time, the trigger, and who is responding. Open
   an incident log; you will record every command you run.

2. **Stop destructive automation.** Recovery must not race deletion. Pause the
   maintenance jobs that can delete or repoint data:

   - `reconcile_reclaim` — deletes store objects (operator-gated; never auto, but
     make sure no one runs it now).
   - `archive_reclaim` — deletes local copies of archived files (seeded weekly).
   - `tier_move` / `tier_move_scan` — moves bytes between tiers and archives.

   The fastest safe pause is to **stop the worker**, which halts all job execution
   while leaving the API and data serving:

   ```bash
   docker compose stop cerebro-worker
   ```

   Producers (the seeded scans) will queue jobs harmlessly; nothing executes until
   you restart the worker. (Re-enable in the scenario's final step.)

3. **Snapshot current state — do not overwrite anything.** Before any restore:
   - Capture a fresh catalogue dump *if Mongo is healthy* (so you can roll back a bad
     restore): see [§8 Trigger a catalogue backup](#trigger-a-catalogue-backup).
   - Do **not** `mongorestore --drop` onto a live, healthy instance.
   - Do **not** delete the failed volume/disk until recovery is verified.

4. **Classify** using the [§1 triage table](#1-severity-triage--find-your-scenario)
   and proceed to the scenario. When in doubt, run a **reconcile scan** (read-only)
   to get a current divergence picture: [§8 Run reconcile](#run-a-reconcile-scan-read-only).

---

## 4. System & state inventory

`docker compose` services (names as deployed):

| Service | Role | Holds durable state? |
|---|---|---|
| `cerebro-api` | REST API (actix) | No |
| `cerebro-app` | Web UI | No |
| `cerebro-database` | MongoDB: catalogue, audit chain, (often) filer metadata | **Yes** |
| `cerebro-fs-master` | SeaweedFS master (topology) | Cluster metadata |
| `cerebro-fs-primary` / `cerebro-fs-replica` | SeaweedFS volumes (replicated) | **Yes — object bytes** |
| `cerebro-fs-secondary` | Additional volume capacity | **Yes — object bytes** |
| `cerebro-fs-filer` | SeaweedFS filer (path addressing) | Metadata (in Mongo) |
| `faktory` | Job queue | Transient queue |
| `cerebro-worker` | Maintenance jobs (tiering, verify, reconcile, backup, archive, repair) | No |

External state (not in compose):

- **Backup store** — `CEREBRO_BACKUP_STORE_PATH`, prefix `CEREBRO_BACKUP_PREFIX`
  (default `catalogue`). Holds `{prefix}/{backup_id}/catalogue.archive.gz` +
  `manifest.json`, and reconcile reports under `reconcile/`.
- **Cold store** — `CEREBRO_ARCHIVE_STORE_PATH` (or `CEREBRO_ARCHIVE_S3_*`), prefix
  `CEREBRO_ARCHIVE_PREFIX` (default `archive`). Holds archived object bytes.

These two may be the same physical target with different prefixes, or separate. Both
must be **independently durable** from the SeaweedFS data.

---

## 5. Scenarios

Each scenario: **Detect → Assess → Recover → Verify → Abort.**

### DR-1 — one object copy lost (volume / disk failure)

**Detect.** A volume or disk failed. `cerebro-fs-master` reports a volume down;
some reads are slow or briefly fail; the API still serves.

**Assess.** With replication (`replication: "001"` or higher), every object has
a second copy. The cluster heals automatically; you confirm and, if needed, nudge it.

```bash
# Inspect volume health and replica placement
docker compose exec cerebro-fs-master weed shell <<'EOF'
volume.list
EOF
```

**Recover.**
1. Replace the failed disk / restart the volume service:
   ```bash
   docker compose up -d cerebro-fs-primary cerebro-fs-secondary
   ```
2. Force re-replication to restore the configured copy count if it has not already:
   ```bash
   docker compose exec cerebro-fs-master weed shell <<'EOF'
   volume.fix.replication
   EOF
   ```

**Verify.**
- `volume.list` shows every volume back at its target replica count.
- Restart the worker and run a **reconcile scan**
  ([§8](#run-a-reconcile-scan-read-only)); the report shows **0 dangling**
  references for the affected volumes. A reappeared object is the signal the copy
  healed.

**Abort.** If re-replication cannot find a surviving copy for some fids, those are in
DR-4/DR-5 territory (confirmed loss) — escalate via reconcile and recover archived
ones from cold; non-archived ones with no copy are unrecoverable.

---

### DR-2 — catalogue (MongoDB) loss or corruption

**Detect.** Catalogue queries fail or return inconsistent data; `cerebro-database`
is down, or a bad migration/operation corrupted records.

**Assess.** The object bytes are independent of the catalogue — this is a metadata
recovery, not a data-loss event, *provided* a recent backup exists. Identify the
newest **verified** backup before restoring.

```bash
# List backups (filesystem cold store shown; adapt for S3)
ls -1 "$CEREBRO_BACKUP_STORE_PATH/catalogue"/*/manifest.json
```

For the chosen `backup_id`, **verify the archive before trusting it**
([§8 Inspect & verify a backup](#inspect--verify-a-backup)): the archive's BLAKE3
must equal `manifest.archive_blake3`, and `manifest.audit_chain_verified` should be
`true`. A backup that fails either check is not a valid recovery point — step back to
the previous one.

**Recover.**
1. If Mongo is healthy but data is bad, **snapshot first** (§3.3) so you can roll
   back.
2. Stop services that write the catalogue:
   ```bash
   docker compose stop cerebro-worker cerebro-api
   ```
3. Restore the verified archive (mirror of the backup's `mongodump --archive --gzip`):
   ```bash
   mongorestore --gzip --archive=/path/to/catalogue.archive.gz \
     --uri="mongodb://<admin-user>:<pw>@localhost:27017/?authSource=admin" --drop
   ```
   `--drop` replaces the collections being restored. Only use it on an instance you
   intend to overwrite (a fresh/empty target, or after the snapshot in step 1).
4. Restart:
   ```bash
   docker compose up -d cerebro-api cerebro-worker
   ```

**Verify.**
- Catalogue queries return expected records; file counts match expectations.
- Re-run the audit-chain verification (the same check the backup manifest recorded);
  it must report the chain intact.
- Run a **reconcile scan**: dangling/orphan counts should reflect only changes since
  the backup point (objects uploaded/archived after the backup), not wholesale
  divergence. Large divergence means you restored the wrong point or `--drop` hit a
  live instance — abort.

**Abort.** If the restore makes things worse, restore the §3.3 snapshot. Never delete
the snapshot until verification passes.

---

### DR-3 — filer metadata loss

**Detect.** Path-addressed I/O fails (`download`/`upload` by path), filer listings
are empty, but volumes are healthy and the catalogue is intact.

**Assess.** The filer's path→fid map lives in MongoDB. If the filer DB is co-located
in `cerebro-database` (the common Cerebro deployment), the whole-instance backup
already contains it and **DR-2's restore recovers the filer too**. Confirm where the
filer store points before acting.

**Recover.**
- **Co-located filer DB:** follow [DR-2](#dr-2--catalogue-mongodb-loss-or-corruption);
  the filer metadata returns with the catalogue. Restart `cerebro-fs-filer`
  afterwards:
  ```bash
  docker compose restart cerebro-fs-filer
  ```
- **Separate filer store:** restore that store from its own backup. If no backup
  exists, the path map can be rebuilt by re-indexing volume contents (SeaweedFS
  vendor procedure) — but the object bytes themselves are safe in the volumes.

**Verify.**
- A known object is retrievable by its filer **path** again.
- `effective_identifier`-based retrieval works for path-addressed files (verify a
  sample by downloading and hashing against the catalogue hash).
- Reconcile orphan detection (filer mode) runs and reports a sane object count —
  a count near zero would indicate the filer is still not listing.

**Abort.** If path retrieval still fails after restore, fall back to fid-addressed
retrieval (objects remain reachable by fid via the volumes) while you escalate the
filer rebuild.

---

### DR-4 — object integrity failure (bit rot)

**Detect.** `verify_file`/`verify_scan` logs `INTEGRITY FAILURE on verify: stored
bytes do not match catalogue hash`, or the lifecycle verify metric fires an alert.

**Assess.** Cerebro repairs this automatically **when the file has a cold backup**
(`archive_key`). The retained-backup behaviour means any file that was archived
at least once keeps a durable cold copy for its whole life. Check the log line for
whether repair already happened:

- `integrity failure REPAIRED from cold backup` → already fixed; go to Verify.
- `INTEGRITY FAILURE ... (no cold backup to repair from)` → no cold source; manual.
- `repair from cold failed: ...` → repair attempted but the cold copy was bad/absent
  (it failed the hash gate, which is correct — a corrupt cold copy is never written
  over a live entry).

**Recover.**
- **Has a cold backup:** force a re-verify to trigger repair, or wait for the weekly
  scan. Trigger now:
  ```bash
  # enqueue a targeted verify — repairs from cold on mismatch (CLI; poll with job-status)
  cerebro-client jobs launch-job --kind verify_file \
    --args '{"file_id":"<FILE_ID>"}' --queue maintenance
  # raw API equivalent
  curl -sS -X POST "$CEREBRO_API_URL/jobs/enqueue" \
    -H "Authorization: Bearer $CEREBRO_API_TOKEN" -H 'Content-Type: application/json' \
    -d '{"kind":"verify_file","args":{"file_id":"<FILE_ID>"},"queue":"maintenance"}'
  ```
- **No cold backup, but the object is replicated:** the healthy replica copy is the
  source — SeaweedFS serves the good copy; re-replication overwrites the bad one
  (DR-1). If only one corrupt copy exists, the bytes are lost; restore from any
  external upstream (the original sequencer output) if retained, then re-register.
- **No cold backup and not recoverable:** this is confirmed loss; record it, alert
  the data owner, and re-ingest from the upstream source if available.

**Verify.**
- Re-run `verify_file` for the file; it must log `integrity verified` (or the repair
  line) and stamp `verified_at`.
- Download the object and confirm its BLAKE3 equals the catalogue hash.

**Abort.** If repair-from-cold keeps failing the hash gate, the cold copy is itself
corrupt — treat as DR-6 for that object and do not force any overwrite.

---

### DR-5 — catalogue↔store divergence

**Detect.** A reconcile report lists **dangling** references (catalogue entries whose
object is missing) and/or **orphans** (stored objects with no catalogue entry).

**Assess.** Reconcile is report-first and never deletes on its own. Read the latest
report from the backup store under `reconcile/`. `verify_repair` (weekly) re-checks
dangling refs and collapses them to *recovered* / *archived-since* / *confirmed
loss* — start from its escalations, not the raw scan, because many dangling refs are
transient.

**Recover.**
- **Dangling, confirmed loss, archived file:** trigger a restore (the cold copy is
  authoritative) — see [§8 Restore an archived file](#restore-an-archived-file).
- **Dangling, confirmed loss, non-archived file:** no in-system source; recover from
  replica (DR-1) if a copy survived, else re-ingest from upstream, else record as
  lost.
- **Orphans (stored object, no catalogue entry):** confirm they are truly orphaned
  (a just-uploaded object can look orphaned for a moment; the grace window filters
  most). When confirmed, reclaim them — **operator-gated, explicit keys only**:
  [§8 Reclaim orphans](#reclaim-orphans-gated).

**Verify.**
- Re-run the reconcile scan; the addressed dangling/orphan entries are gone and no
  new ones appeared.
- For each restored file, confirm retrievability + hash (DR-4 Verify).

**Abort.** If an "orphan" turns out to be referenced (e.g. a shared filer holding
non-Cerebro data, or a catalogue that was truncated by budget), **do not reclaim**.
Orphan detection only runs on a complete catalogue for this reason; if you raised the
budget, re-confirm completeness first.

---

### DR-6 — cold object store loss

**Detect.** The cold store is unreachable or lost; restores fail to fetch bytes;
`archive_reclaim` logs `archived file missing its cold copy`; new archival errors.

**Assess.** Impact depends on what was *only* in cold:
- Files that are **archived and whose local copy was reclaimed** have the cold
  store as their only copy — those are at risk.
- Files still holding a local copy are safe; the cold loss only removes their backup.
- The **catalogue backups** also live here — losing the cold store can also mean
  losing your catalogue recovery point (DR-2). Treat both.

**Recover.**
1. **Stop reclaim immediately** (§3.2) so `archive_reclaim` cannot delete any more
   local copies while their cold backups are gone.
2. Restore or replace the cold store from its **own** durable backup (this is why the
   cold store must be independently durable). Re-point `CEREBRO_ARCHIVE_STORE_PATH` /
   `CEREBRO_ARCHIVE_S3_*` to the recovered location.
3. For files whose local copy still exists but whose cold backup was lost,
   **re-archive** to rebuild the cold copy: let `tier_move` re-run, or re-trigger
   archival for the affected files.
4. For files that were cold-only and the cold store is unrecoverable: confirmed loss
   — record and re-ingest from upstream if available.

**Verify.**
- A sample archived file restores end-to-end and matches its catalogue hash.
- `archive_reclaim` (dry-run by leaving the worker stopped, then a single bounded
  run) reports `missing_cold = 0` for files past grace.
- Catalogue backups are present and verify (DR-2 Assess).

**Abort.** Keep the worker's reclaim paused until a sample restore verifies; do not
resume `archive_reclaim` until the cold store is confirmed healthy.

---

### DR-7 — full host rebuild

**Detect.** The deployment host is lost; you are rebuilding from infrastructure +
backups.

**Assess.** You will rebuild the stack, restore the catalogue, re-attach or rebuild
the object volumes, and reconcile. Order matters: catalogue and filer metadata first,
then verify objects against them.

**Recover.**
1. Re-provision the host and bring up infrastructure only:
   ```bash
   docker compose up -d cerebro-database faktory cerebro-fs-master
   ```
2. **Re-attach the SeaweedFS volume data** if the disks survived (mount them back to
   `cerebro-fs-primary`/`-replica`/`-secondary`), then bring those up. If the disks
   are gone, objects come from the cold store (archived) or are lost (non-archived) —
   plan re-ingestion for the latter.
3. **Restore the catalogue** from the newest verified backup — [DR-2 Recover](#dr-2--catalogue-mongodb-loss-or-corruption). This also restores the
   co-located filer metadata.
4. Bring up the rest:
   ```bash
   docker compose up -d cerebro-fs-filer cerebro-api cerebro-app cerebro-worker reverse-proxy
   ```
5. Re-seed schedules if needed (idempotent; only adds missing) and let the maintenance
   jobs resume.

**Verify.**
- API serves; catalogue counts match the backup point.
- Run a **reconcile scan**: dangling refs identify objects the catalogue expects but
  the rebuilt store lacks — these are your re-ingestion / cold-restore worklist.
- Restore one archived file and one live file end-to-end and hash-check both.
- Confirm the audit chain verifies.

**Abort.** If the catalogue restore point predates significant ingestion, you will see
large dangling counts for recent files — recover those from cold (if archived) or
upstream, and document the RPO gap.

---

### DR-8 — bad destructive job

**Detect.** A reclaim or relocate ran incorrectly — local copies deleted, fids
repointed, or objects removed — e.g. a `reconcile_reclaim` with wrong keys, or an
`archive_reclaim`/`tier_move` misconfiguration.

**Assess.** Identify the blast radius from the job logs and the audit trail (every
relocate/tier-move/delete writes an audit event). The gates limit damage: reclaim is
key-explicit, archive_reclaim requires a confirmed-present cold copy, and repair/
restore re-verify hashes. Most "bad job" outcomes are recoverable from cold or replica.

**Recover.**
1. **Stop the worker** (§3.2) to halt further damage.
2. For **deleted local copies of archived files** (`archive_reclaim`): no recovery
   needed — the cold copy is intact by the gate's guarantee; the file restores on
   demand. Verify the cold copy exists for the affected keys.
3. For **wrongly reclaimed orphans** that were *not* truly orphaned: restore from
   replica (if a copy survived) or from cold (if archived); otherwise re-ingest.
4. For **wrong repoints** (`relocate`): the audit trail records the previous state;
   re-relocate to the correct fid/tier. If the repoint pointed at corrupt bytes,
   `verify` + repair-from-cold (DR-4) corrects it.

**Verify.**
- The affected files retrieve and hash-check.
- A reconcile scan shows no residual dangling/orphan entries for the affected set.
- Review the audit trail to confirm the corrective actions are recorded.

**Abort.** If recovery sources are exhausted for some files, record confirmed loss and
move to re-ingestion. Then address the root cause (job args, schedule, config) before
resuming the worker.

---

## 6. Routine recovery drills

> An untested backup is not a backup. These rehearsals prove the recovery paths work
> *before* you need them, and produce the evidence a clinical audit will ask for.

| Drill | Cadence | What it proves |
|---|---|---|
| **Catalogue restore-to-scratch** | Monthly | The newest backup restores into a throwaway Mongo and the audit chain verifies. |
| **Archived-file restore** | Monthly | A random archived file restores from cold and matches its catalogue hash. |
| **Integrity repair** | Quarterly | A deliberately corrupted *test* object is auto-repaired from cold (DR-4). |
| **Replica heal** | Quarterly | Taking one volume offline triggers re-replication with no dangling refs (DR-1). |
| **Full game-day** | Annually | A from-backup rebuild (DR-7) on a staging host, timed against the RTO target. |

**Drill procedure (restore-to-scratch example):**
1. Pick the newest backup; verify its BLAKE3 + `audit_chain_verified`
   ([§8](#inspect--verify-a-backup)).
2. `mongorestore --gzip --archive=<copy> --uri=<scratch-instance>` (never the live
   instance).
3. Run the audit-chain verification against the scratch instance.
4. Record the result in the log below. Destroy the scratch instance.

**Drill sign-off log:**

| Date | Drill | Backup/object id | Result | RPO/RTO observed | Operator |
|---|---|---|---|---|---|
|  |  |  |  |  |  |

---

## 7. Prevention checklist

The cheapest recovery is the one you never need. Confirm periodically:

- [ ] SeaweedFS replication is `≥ 001` so every live object has a second copy
      (`fs-replication.md`).
- [ ] Catalogue backups run on schedule and the **latest** verifies
      (`audit_chain_verified = true`, BLAKE3 matches).
- [ ] The backup store and cold store are **independently durable** (separate disk +
      offsite/object-lock) from the SeaweedFS data.
- [ ] Archival cadence keeps important data with a cold copy (so DR-4/DR-6 have a
      source).
- [ ] `reconcile_scan` runs and recent reports show low, explained divergence.
- [ ] The local-copy grace (`CEREBRO_ARCHIVE_LOCAL_GRACE_DAYS`) is long enough that
      a verify pass runs before a local copy is reclaimed.

---

## 8. Operations reference

> All on-demand jobs go through the admin-only `POST /jobs/enqueue` endpoint and are
> executed by `cerebro-worker`. Poll completion at `GET /jobs/enqueue/{id}`. The
> worker must be running for jobs to execute (re-start it after a §3 pause).
>
> The operator path is the `cerebro-client` CLI (it wraps these endpoints); the raw
> `curl` form is kept for scripting. Set an admin token up first:
>
> ```bash
> export CEREBRO_API_URL="https://api.<your-domain>"
> export CEREBRO_API_TOKEN="$(cerebro-client login --email admin@cerebro --password '****')"
> ```
> Enqueue with `cerebro-client jobs launch-job --kind <k> --args '<json>' --queue maintenance`
> and poll with `cerebro-client jobs job-status --id <id>`. See the
> [maintenance guide](maintenance.md#running-a-job-on-demand) and the
> [operations walkthrough](walkthrough.md) for the full tour.

### Trigger a catalogue backup
```bash
# CLI (operator path) — prints the job id; poll it with job-status.
cerebro-client jobs launch-job --kind catalogue_backup --args '{}' --queue maintenance

# raw API equivalent
curl -sS -X POST "$CEREBRO_API_URL/jobs/enqueue" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" -H 'Content-Type: application/json' \
  -d '{"kind":"catalogue_backup","args":{},"queue":"maintenance"}'
```

### Inspect & verify a backup
A backup is two keys under `{prefix}/{backup_id}/`: `catalogue.archive.gz` and
`manifest.json`. Before restoring, **verify** the checksum and the chain flag:

```bash
STORE="$CEREBRO_BACKUP_STORE_PATH"; PREFIX="${CEREBRO_BACKUP_PREFIX:-catalogue}"
ls -1 "$STORE/$PREFIX"/*/manifest.json                     # list backups
BID=$(ls -1 "$STORE/$PREFIX" | sort | tail -n1)            # newest id
jq '{archive_blake3, audit_chain_verified, created_at}' "$STORE/$PREFIX/$BID/manifest.json"
test "$(b3sum --no-names "$STORE/$PREFIX/$BID/catalogue.archive.gz")" \
   = "$(jq -r .archive_blake3 "$STORE/$PREFIX/$BID/manifest.json")" \
   && echo "checksum OK" || echo "CHECKSUM MISMATCH — do not restore this backup"
```

Require `audit_chain_verified == true`. If it is `false`/`null`, the chain was not
intact at backup time — prefer an earlier good backup.

### Verify the live audit chain
After any catalogue restore (and on a schedule), confirm the live chain verifies. The
audit endpoint is team-scoped and returns the chain status in `data.verified`:

```bash
curl -sS "$CEREBRO_API_URL/audit?team=<TEAM>" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" | jq '.data.verified'   # must be true
```

If the live chain and a backup's `audit_chain_verified` disagree, stop and investigate
— this is the guard against a restore silently laundering a broken chain.

### Find a file id / list files
The verify and restore actions take a file id (and team). List the catalogue to find
them (team-scoped; `limit=0` returns all):

```bash
curl -sS "$CEREBRO_API_URL/files?team=<TEAM>&page=0&limit=0" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" \
  | jq -r '.data[] | [.id, .name, .tier, .archived, .restore_state] | @tsv'
```

### Restore the catalogue
The operator path is `restore-catalogue.sh`, which **verifies the archive's BLAKE3
against the manifest and prints the audit-chain state before touching MongoDB** — use
it rather than a bare `mongorestore`:

```bash
# Validate the latest backup is loadable, without writing anything.
restore-catalogue.sh --uri "mongodb://admin:***@cerebro-database:27017/?authSource=admin" \
  --backup latest --dry-run
# Restore a specific backup (asks for typed confirmation; --drop for a clean restore).
restore-catalogue.sh --uri "mongodb://admin:***@cerebro-database:27017/?authSource=admin" \
  --backup 20260617T031500Z --drop
```

The underlying command it runs is `mongorestore --gzip --archive=<file> --uri=<uri>`
(add `--drop` only against a target you mean to overwrite). Requires `jq`, `b3sum`, and
`mongorestore` on the host. After loading, re-verify the live audit chain (below).

### Restore an archived file
Archived files are brought back by the restore state machine. **Request** the restore
by transitioning the file to `Requested`; the hourly `restore_scan` advances it and
`restore_drive` re-materialises the bytes from cold, hash-verifies them, writes them
back to the file's effective location, and keeps the cold backup. The restore endpoint
is team-scoped (`?team=`), so pass the owning team:

```bash
# Request the restore (FILE_ID + TEAM from the catalogue listing / reconcile report).
curl -sS -X POST "$CEREBRO_API_URL/files/<FILE_ID>/restore?team=<TEAM>" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" -H 'Content-Type: application/json' \
  -d '{"target":"Requested"}'

# Optional: push the pass now instead of waiting for the hourly scan.
cerebro-client jobs launch-job --kind restore_scan --args '{}' --queue maintenance
```

Then confirm: list the file and check `restore_state` is `Restored` (the serialized
value, with `archived` back to `false`), download it, and verify its BLAKE3 equals the
catalogue hash (see [DR-4 Verify](#dr-4--object-integrity-failure-bit-rot)).

### Run a reconcile scan (read-only)
```bash
# CLI (operator path). Raise the budget to cover a large catalogue + enable orphan detection.
cerebro-client jobs launch-job --kind reconcile_scan \
  --args '{"budget":100000,"grace_days":1}' --queue maintenance

# raw API equivalent
curl -sS -X POST "$CEREBRO_API_URL/jobs/enqueue" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" -H 'Content-Type: application/json' \
  -d '{"kind":"reconcile_scan","args":{},"queue":"maintenance"}'
```
The report is written to the backup store under `reconcile/`. Orphan detection runs
only in filer mode and only when the catalogue is complete within budget. Reading the
report is shown in [consistency-reconcile.md](consistency-reconcile.md#reconcile_scan-scheduled-safe).

### Reclaim orphans (gated)
**Destructive, operator-gated.** Requires `confirm: true` and the *explicit* keys to
delete — reclaim never enumerates-and-deletes on its own:
```bash
# CLI (operator path).
cerebro-client jobs launch-job --kind reconcile_reclaim \
  --args '{"confirm":true,"keys":["<key1>","<key2>"],"max_delete":100}' --queue maintenance

# raw API equivalent
curl -sS -X POST "$CEREBRO_API_URL/jobs/enqueue" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" -H 'Content-Type: application/json' \
  -d '{"kind":"reconcile_reclaim","args":{"confirm":true,"keys":["<key1>","<key2>"]},"queue":"maintenance"}'
```
Path keys route to the filer delete, fids to the volume delete, automatically.

### Maintenance schedule (seeded)
Producers run weekly/daily on the `maintenance` queue (offsets stagger them):
`retention_sweep`, `tier_move_scan`, `verify_scan`, `restore_scan`, `purge_reclaim`,
`catalogue_backup`, `reconcile_scan`, `verify_repair`, `archive_reclaim`. The
operator-gated `reconcile_reclaim` is **not** seeded — it is run by hand only.

### Environment pointers
Backup: `CEREBRO_BACKUP_MONGO_URI`, `CEREBRO_BACKUP_STORE_PATH`,
`CEREBRO_BACKUP_PREFIX`, `CEREBRO_BACKUP_KEEP_LAST`, `CEREBRO_BACKUP_MIN_AGE_DAYS`.
Archive: `CEREBRO_ARCHIVE_STORE_PATH` / `CEREBRO_ARCHIVE_S3_*`,
`CEREBRO_ARCHIVE_PREFIX`, `CEREBRO_ARCHIVE_LOCAL_GRACE_DAYS`. Full tables live in
`catalogue-backup.md` and `archival.md`.

---

## 9. Recovery objectives worksheet

Fill these in for your lab; they are policy, not defaults Cerebro can choose.

| Data class | Target RPO | Target RTO | Backup/replica cadence set to meet it | Owner |
|---|---|---|---|---|
| Catalogue | | | `catalogue_backup` interval = | |
| Live objects | | | replication = | |
| Archived objects | | | archival cadence + cold durability = | |

Escalation contacts:

| Role | Name | Reach |
|---|---|---|
| Primary on-call |  |  |
| Data owner / lab lead |  |  |
| Infrastructure / storage |  |  |

---

## 10. Review

Review this runbook after any incident, after any topology change, and at least
twice a year. Re-run the §6 drills on cadence and keep the sign-off log current — the
log is the evidence that recovery works.
