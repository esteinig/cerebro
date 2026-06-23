# Operations walkthrough

A hands-on, task-by-task how-to for running a cerebro-fs deployment. Where the other
documents explain *what* each capability is, this one walks through *how* to drive it
from the operator's seat — the actual commands, their arguments, and how you confirm
each one worked. It doubles as a runnable example and as the operational guide for the
routine administrative and maintenance tasks.

Work through it top to bottom the first time (stand up → connect → operate); after
that, jump to the task you need. Every task follows the same shape: **what/why**, the
**command(s)**, and **how you know it worked**.

> **Audience & prerequisites.** A Cerebro operator with shell access to the deployment
> host, the `docker compose` project, and an **admin** account on the API. Two CLIs do
> almost everything here: `cerebro` (deploy/stack) and `cerebro-client` (API
> operations). SeaweedFS volume operations use `weed shell` inside the master
> container; the catalogue restore uses the shipped `restore-catalogue.sh`.

## Conventions used below

The `cerebro-client` examples read the API URL and token from the environment, set
once per session. The `docker compose …` examples assume you are in the stack's
`docker/` directory (or that you pass `-f <outdir>/docker/docker-compose.yml`). Service
names are those in the [recovery runbook §4](disaster-recovery.md#4-system--state-inventory).

```bash
export CEREBRO_API_URL="https://api.<your-domain>"     # default: http://localhost:8080
cd /opt/cerebro/stack/docker                            # where deploy rendered the compose file
```

---

## Part 1 — Stand up and connect

### 1.1 Render and bring up the stack

**What/why.** `cerebro stack deploy` renders the compose file, env files, and secrets
into an output directory; `docker compose` runs them. Both steps are idempotent.

```bash
cerebro stack deploy \
  --name prod --outdir /opt/cerebro/stack \
  --config localhost-fs \
  --fs-model single-server-replicated \
  --fs-hot /mnt/ssd0/hot         --fs-cold /mnt/hdd0/cold \
  --fs-replica-hot /mnt/ssd1/hot --fs-replica-cold /mnt/hdd1/cold

docker compose -f /opt/cerebro/stack/docker/docker-compose.yml up -d
```

**Verify.** Every service is up:

```bash
docker compose ps        # all services Up; none restarting
```

The deployment models and disk-layout rules (and what `001` protects against) are in
[storage & replication](fs-replication.md); the full configuration reference is in the
[administration guide](administration.md).

### 1.2 Authenticate the CLI

**What/why.** Maintenance/job commands need an admin token. `login` prints the token to
stdout; capture it into the variable the CLI reads.

```bash
export CEREBRO_API_TOKEN="$(cerebro-client login --email admin@cerebro --password '****')"
```

**Verify.** The token is accepted:

```bash
cerebro-client ping-server      # logs: Cerebro API status: ok
```

Alternatives: `--token-file token.json` (then point the CLI at it with `-f`, or set
`CEREBRO_API_TOKEN_FILE`); add `--bot` for a service `Role::Bot` token.

### 1.3 First health check

**What/why.** Before relying on automation, confirm the API, the worker, and the object
store are healthy.

```bash
curl -sS http://<worker-host>:9464/health                       # → ok
docker compose exec cerebro-fs-master weed shell <<'EOF'
volume.list
EOF
```

**Verify.** `/health` returns `ok`, and `volume.list` shows every volume at its target
replica count (two copies under the replicated model).

---

## Part 2 — The everyday pattern

Almost every maintenance action is "enqueue a job, then poll it". Learn this once and
the rest of the tasks are variations on it.

### 2.1 Enqueue a job and wait for it

**What/why.** `launch-job` enqueues a worker job on the `maintenance` queue and prints
the job id; `job-status` reports completion and the result/error.

```bash
id=$(cerebro-client jobs launch-job --kind <job-kind> --args '<json-args>' --queue maintenance)
cerebro-client jobs job-status   --id "$id"     # completed=true status=succeeded result=…
cerebro-client jobs job-complete --id "$id"     # yes | no   (handy in scripts)
```

The same call as a raw API request (the CLI wraps this) is in the
[maintenance guide](maintenance.md#running-a-job-on-demand). The worker must be running
for jobs to execute.

### 2.2 Inspect the schedule

**What/why.** See what the seeded producers will run next, and when they last ran.

```bash
cerebro-client jobs schedule-summary        # next / previous 10 scheduled runs
```

To add or override a schedule (e.g. a more frequent backup, or to re-seed after a fresh
database), use `schedule-job` — see
[maintenance → scheduling](maintenance.md#scheduling-and-re-seeding).

### 2.3 Find a file in the catalogue

**What/why.** Targeted verify, archive, and restore actions take a file id (and the
owning team). List the catalogue to find them (team-scoped; `limit=0` returns all,
`limit>0` paginates with `page`).

```bash
curl -sS "$CEREBRO_API_URL/files?team=<TEAM>&page=0&limit=0" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" \
  | jq -r '.data[] | [.id, .name, .tier, .archived, .restore_state] | @tsv'
```

`cerebro-fs list -r <run-id>` lists a run's files from the same catalogue if you
prefer the FS client; the API form above is easiest when you want ids to feed into the
job commands below.

---

## Part 3 — Maintenance & administration tasks

### 3.1 Force a catalogue backup (and verify it)

**What/why.** Besides the daily schedule, take a backup before risky changes — e.g.
before a restore, so you can roll back.

```bash
id=$(cerebro-client jobs launch-job --kind catalogue_backup --args '{}' --queue maintenance)
cerebro-client jobs job-status --id "$id"
```

**Verify.** A new backup directory appears and its checksum + chain flag are good:

```bash
STORE="$CEREBRO_BACKUP_STORE_PATH"; PREFIX="${CEREBRO_BACKUP_PREFIX:-catalogue}"
BID=$(ls -1 "$STORE/$PREFIX" | sort | tail -n1)
jq '{archive_blake3, audit_chain_verified, created_at}' "$STORE/$PREFIX/$BID/manifest.json"
test "$(b3sum --no-names "$STORE/$PREFIX/$BID/catalogue.archive.gz")" \
   = "$(jq -r .archive_blake3 "$STORE/$PREFIX/$BID/manifest.json")" && echo "checksum OK"
```

`audit_chain_verified` must be `true` and the checksum must match. Detail:
[catalogue backup](catalogue-backup.md).

### 3.2 Verify and repair a file's integrity

**What/why.** Re-hash one object against its catalogue hash. On a mismatch the worker
repairs from the cold backup automatically when the file has an `archive_key`.

```bash
id=$(cerebro-client jobs launch-job --kind verify_file --args '{"file_id":"<FILE_ID>"}' --queue maintenance)
cerebro-client jobs job-status --id "$id"
```

**Verify.** The worker logs `integrity verified` (or `integrity failure REPAIRED from
cold backup`) and stamps `verified_at`; an unrepairable mismatch logs `INTEGRITY
FAILURE …` and increments the verify-failure series of
`cerebro_file_lifecycle_ops_total`. The decision tree per outcome is
[DR-4](disaster-recovery.md#dr-4--object-integrity-failure-bit-rot); mechanism detail is
in [verify & repair](verify-repair.md).

### 3.3 Run a consistency reconcile and read the report

**What/why.** Detect catalogue↔store divergence (dangling references, orphan objects).
Report-first: the scan never deletes. Raise the budget to cover a large catalogue and
enable orphan detection.

```bash
id=$(cerebro-client jobs launch-job --kind reconcile_scan \
  --args '{"budget":100000,"grace_days":1}' --queue maintenance)
cerebro-client jobs job-status --id "$id"
```

**Verify / read the result.** The report lands under the backup store's `reconcile/`
prefix:

```bash
RID=$(ls -1 "$CEREBRO_BACKUP_STORE_PATH/reconcile" | sort | tail -n1)
jq '{store_enumerated, dangling: (.dangling|length), orphans: (.orphans|length)}' \
   "$CEREBRO_BACKUP_STORE_PATH/reconcile/$RID"
```

If `store_enumerated` is `false`, orphans were not computed (volume mode, or the
catalogue exceeded the budget — raise it and re-run). Full semantics:
[consistency reconcile](consistency-reconcile.md).

### 3.4 Reclaim orphans (destructive, operator-gated)

**What/why.** Delete stored objects that have no catalogue entry. This is the one
destructive maintenance action that is **never** automatic: it requires `confirm:true`
and the *explicit* keys you confirmed from the reconcile report.

```bash
cerebro-client jobs launch-job --kind reconcile_reclaim \
  --args '{"confirm":true,"keys":["3,01abcd…","5,07ef12…"],"max_delete":100}' \
  --queue maintenance
```

**Verify.** Re-run a reconcile scan (3.3) and confirm the reclaimed keys are gone and no
new orphans appeared. If a key turns out to be referenced, **do not** reclaim it — see
the abort guidance in [DR-5](disaster-recovery.md#dr-5--cataloguestore-divergence).

### 3.5 Archive a file and reclaim its local copy

**What/why.** Move a file's bytes to the cold store (archival), then later reclaim the
now-redundant local copy once it is safe.

```bash
# Archive now: move to Cold (copies to cold store + repoints the catalogue).
cerebro-client jobs launch-job --kind tier_move \
  --args '{"file_id":"<FILE_ID>","target":"Cold"}' --queue maintenance

# Reclaim local copies of archived files past the grace window (bounded; cold copy must exist).
cerebro-client jobs launch-job --kind archive_reclaim \
  --args '{"budget":1000,"max_delete":100}' --queue maintenance
```

**Verify.** After archival the file lists as `archived=true` with an `archive_key`
(3.3 listing). `archive_reclaim` only deletes a local copy where the cold copy is
confirmed present; it logs each reclaim. Detail: [archival](archival.md).

### 3.6 Restore an archived file

**What/why.** Bring an archived object back. You **request** the restore; the hourly
`restore_scan` advances it and `restore_drive` re-materialises the bytes from cold,
hash-verified, keeping the cold backup. The restore endpoint is team-scoped.

```bash
curl -sS -X POST "$CEREBRO_API_URL/files/<FILE_ID>/restore?team=<TEAM>" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" -H 'Content-Type: application/json' \
  -d '{"target":"Requested"}'

# Optional: push the restore pass now instead of waiting for the hourly scan.
cerebro-client jobs launch-job --kind restore_scan --args '{}' --queue maintenance
```

**Verify.** List the file (3.3) and confirm `restore_state` is `Restored` (the
serialized value the listing returns, alongside `archived` flipping back to `false`);
download it and confirm its BLAKE3 equals the catalogue hash
([DR-4 Verify](disaster-recovery.md#dr-4--object-integrity-failure-bit-rot)).

### 3.7 Rotate a secret

**What/why.** Replace a rendered credential and restart its consumer. Example: the
read-only catalogue-backup MongoDB URI.

```bash
printf 'mongodb://cerebro_backup:<new-password>@cerebro-database:27017/?authSource=admin' \
  > /opt/cerebro/stack/mongodb/backup_mongo_uri.secret
docker compose restart cerebro-worker
```

**Verify.** Trigger a catalogue backup (3.1) and confirm it succeeds with the new
credential. The same overwrite-and-restart pattern applies to any rendered secret
([administration → secrets](administration.md#secrets)).

### 3.8 Check and fix replication health

**What/why.** A silently under-replicated volume is a latent single copy. `volume.list`
is the authoritative view; `volume.fix.replication` re-replicates once a failed
disk/server has been replaced.

```bash
docker compose exec cerebro-fs-master weed shell <<'EOF'
volume.list
EOF
# After replacing failed hardware, nudge re-replication:
docker compose exec cerebro-fs-master weed shell <<'EOF'
volume.fix.replication
EOF
```

**Verify.** `volume.list` shows every volume back at its target replica count; a
follow-up reconcile scan (3.3) reports no dangling references for the affected volumes.
Monitoring guidance: [storage & replication](fs-replication.md#monitoring-under-replication).

---

## Part 4 — Handling an incident

When something has failed, **do not** improvise. The recovery runbook is the
authority; this is the first-15-minutes shape so you reach for the right page fast.

1. **Pause destructive automation** so recovery never races deletion — stop the worker:

   ```bash
   docker compose stop cerebro-worker
   ```

2. **Snapshot before you change anything.** If MongoDB is healthy, take a catalogue
   backup (3.1) so you can roll back a bad restore. Do **not** `mongorestore --drop`
   onto a live, healthy instance, and do **not** delete a failed disk until recovery
   verifies.

3. **Triage and follow the scenario.** Match the symptom in the
   [severity-triage table](disaster-recovery.md#1-severity-triage--find-your-scenario)
   and work that scenario to its **Verify** step. When unsure, run a read-only reconcile
   scan (3.3) for a current divergence picture — it is safe in every case.

4. **Resume the worker** only once the scenario's verification passes:

   ```bash
   docker compose up -d cerebro-worker
   ```

The catalogue restore itself is the `restore-catalogue.sh` operator path (it verifies
the archive's checksum and prints the audit-chain state before touching MongoDB):

```bash
restore-catalogue.sh --uri "mongodb://admin:***@cerebro-database:27017/?authSource=admin" \
  --backup latest --dry-run                      # validate, write nothing
restore-catalogue.sh --uri "mongodb://admin:***@cerebro-database:27017/?authSource=admin" \
  --backup 20260617T031500Z --drop               # restore (typed confirmation)
```

After loading, re-verify the live audit chain:

```bash
curl -sS "$CEREBRO_API_URL/audit?team=<TEAM>" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" | jq '.data.verified'   # must be true
```

---

## Part 5 — Drills (prove recovery before you need it)

An untested backup is not a backup. Run the [recovery
drills](disaster-recovery.md#6-routine-recovery-drills) on cadence and keep the sign-off
log current. The quickest, the monthly **catalogue restore-to-scratch**, never touches
the live instance:

```bash
# Verify the newest backup, then load it into a THROWAWAY Mongo (never the live URI).
restore-catalogue.sh --uri "mongodb://admin:***@scratch-host:27017/?authSource=admin" \
  --backup latest --dry-run
# … then a real restore into the scratch instance, verify the chain, and destroy it.
```

---

## Command quick-reference

| Task | Command |
|---|---|
| Deploy / render the stack | `cerebro stack deploy --name … --outdir … --config localhost-fs --fs-model …` |
| Bring the stack up / down | `docker compose up -d` · `docker compose stop cerebro-worker` |
| Authenticate | `export CEREBRO_API_TOKEN="$(cerebro-client login --email … --password …)"` |
| Health | `curl …:9464/health` · `docker compose exec cerebro-fs-master weed shell <<<'volume.list'` |
| Enqueue a job / poll | `cerebro-client jobs launch-job --kind … --args '…' --queue maintenance` · `… job-status --id …` |
| Schedule / inspect schedule | `cerebro-client jobs schedule-job …` · `cerebro-client jobs schedule-summary` |
| List the catalogue | `curl "$CEREBRO_API_URL/files?team=…&page=0&limit=0"` |
| Catalogue backup | `cerebro-client jobs launch-job --kind catalogue_backup --args '{}' --queue maintenance` |
| Verify / repair a file | `cerebro-client jobs launch-job --kind verify_file --args '{"file_id":"…"}' --queue maintenance` |
| Reconcile (scan / reclaim) | `… --kind reconcile_scan …` · `… --kind reconcile_reclaim --args '{"confirm":true,"keys":[…]}'` |
| Archive / reclaim local copy | `… --kind tier_move --args '{"file_id":"…","target":"Cold"}'` · `… --kind archive_reclaim …` |
| Restore an archived file | `curl -X POST "$CEREBRO_API_URL/files/<id>/restore?team=…" -d '{"target":"Requested"}'` |
| Restore the catalogue | `restore-catalogue.sh --uri … --backup latest --dry-run` |
| Verify the audit chain | `curl "$CEREBRO_API_URL/audit?team=…" \| jq '.data.verified'` |

## Where to next

- The complete per-capability semantics: [maintenance](maintenance.md),
  [catalogue backup](catalogue-backup.md), [consistency reconcile](consistency-reconcile.md),
  [archival](archival.md), [verify & repair](verify-repair.md),
  [storage & replication](fs-replication.md).
- The full recovery procedures and drills: [disaster-recovery runbook](disaster-recovery.md).
- Deployment and configuration: [administration guide](administration.md).
- What is proven vs. what your environment must confirm: [validation guide](validation.md).
