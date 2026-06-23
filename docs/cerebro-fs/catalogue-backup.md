# Catalogue & audit backup / restore

The MongoDB control plane — every team's file catalogue, lifecycle state,
provenance, users, schedules, and the tamper-evident audit chain — is
irreplaceable. Lose it and the objects in `cerebro-fs` are orphaned and the audit
history is gone. cerebro-fs adds a scheduled, checksum-verified, audit-chain-aware
backup of that control plane, plus a scripted restore.

## What runs

A `catalogue_backup` worker job (Faktory `kind`) runs daily, seeded by the
scheduler (`seed:catalogue_backup`, queue `maintenance`). Each run:

1. **Verifies the audit chain first.** It asks the API for the audit trail and
   records whether the chain verified intact, in the backup manifest. This is
   best-effort — if the check can't run it records "unknown" rather than failing
   the backup — but it means a restore can detect a backup taken over an already
   broken chain.
2. **Dumps** the control plane with `mongodump --archive --gzip`. By default it
   dumps the whole instance (the catalogue spans every team's database), against a
   **read-only** MongoDB user.
3. **Checksums** the archive with BLAKE3 (the same hash `cerebro-fs` uses).
4. **Uploads** the archive and a `manifest.json` to the backup object store.
5. **Prunes** old backups: keep the last *N* (default 14), and only delete
   older-than-*N* backups once they pass an age floor (default 7 days), so a burst
   of manual backups can't evict the daily history.

The job is heavy I/O and lives in the worker, never in an API request handler.

## Triggering and verifying a backup on demand

Besides the daily schedule, an operator can force a backup (admin token required —
see [administration → operator CLI setup](administration.md#operator-cli-setup-authentication)):

```bash
# Enqueue a backup and wait for it.
id=$(cerebro-client jobs launch-job --kind catalogue_backup --args '{}' --queue maintenance)
cerebro-client jobs job-status --id "$id"     # completed=true status=succeeded …
```

Then confirm the new backup is well-formed — this is the same check a restore makes,
so run it after any backup you intend to rely on:

```bash
STORE="$CEREBRO_BACKUP_STORE_PATH"; PREFIX="${CEREBRO_BACKUP_PREFIX:-catalogue}"
BID=$(ls -1 "$STORE/$PREFIX" | sort | tail -n1)            # newest backup id
jq '{archive_blake3, audit_chain_verified, created_at}' "$STORE/$PREFIX/$BID/manifest.json"
# The archive's BLAKE3 must equal manifest.archive_blake3, and audit_chain_verified must be true:
test "$(b3sum --no-names "$STORE/$PREFIX/$BID/catalogue.archive.gz")" \
   = "$(jq -r .archive_blake3 "$STORE/$PREFIX/$BID/manifest.json")" && echo "checksum OK"
```

A backup whose checksum does not match, or whose `audit_chain_verified` is not `true`,
is not a valid recovery point — prefer the previous good backup.

## Object store backends

The backup engine is written against a small `ObjectStore` trait, so the target
is swappable:

- **Filesystem (default).** Objects are files under `CEREBRO_BACKUP_STORE_PATH`
  (mounted at `/backups`). For an *independent failure domain* this must be a disk
  physically separate from the SeaweedFS data, or — for host-failure durability —
  an NFS / remote mount. The default compose named volume `cerebro_backups` lives
  on the same Docker host, so override the `/backups` mount for production.
- **S3-compatible.** AWS S3 / MinIO / SeaweedFS-S3 land behind the same
  trait, with no change to the backup or retention logic.

## Configuration

Backup is **opt-in**: the worker only runs it when both a backup Mongo URI and a
store path are set.

| env | meaning | default |
|-----|---------|---------|
| `CEREBRO_BACKUP_MONGO_URI` / `_FILE` | read-only MongoDB URI for `mongodump` | — (required) |
| `CEREBRO_BACKUP_STORE_PATH` | filesystem store root | — (required) |
| `CEREBRO_BACKUP_MONGO_DB` | scope to one DB; empty = whole instance | empty |
| `CEREBRO_BACKUP_PREFIX` | key prefix under the store | `catalogue` |
| `CEREBRO_BACKUP_KEEP_LAST` | most-recent backups always kept | `14` |
| `CEREBRO_BACKUP_MIN_AGE_DAYS` | age floor before older backups are pruned | `7` |

The compose template wires these on the worker and mounts the store at `/backups`.

### The backup user + secret (auto-provisioned)

`mongodump` needs MongoDB credentials, and the catalogue backup must never run with
write access. This is **auto-provisioned** by `cerebro stack deploy` — no
manual step:

- `mongo-init` creates a `cerebro_backup` user (configurable via
  `mongodb.backup_username`) with MongoDB's built-in **`backup`** role on `admin`,
  which grants exactly the read access `mongodump` needs to dump the whole instance
  and nothing more.
- The password defaults to a generated 32-character alphanumeric string (URI-safe,
  so no percent-encoding is needed) unless you set `mongodb.backup_password`.
- deploy renders the connection URI to the Docker secret the compose `secrets:`
  block references:
  ```
  {{ outdir }}/mongodb/backup_mongo_uri.secret
  # contents, one line:
  mongodb://cerebro_backup:<generated>@cerebro-database:27017/?authSource=admin
  ```

`mongo-init` only runs on first database initialisation, so this is the credential
the backup authenticates with from then on. To use an external or differently-scoped
backup user, overwrite `backup_mongo_uri.secret` after deploy.

The worker image installs `mongodb-database-tools` (for `mongodump`) in
`Dockerfile.worker`; the pinned version/arch there may need adjusting for your
platform.

## Restoring

`templates/stack/scripts/restore-catalogue.sh` performs an operator-gated restore
from the filesystem store. It is intentionally **not** an automatic job — a
restore overwrites live data.

```bash
# Validate a backup is loadable, without writing (checksum + mongorestore --dryRun):
restore-catalogue.sh --uri "mongodb://admin:***@cerebro-database:27017/?authSource=admin" \
  --backup latest --dry-run

# Restore the latest backup (asks for confirmation; --drop for a clean restore):
restore-catalogue.sh --uri "mongodb://admin:***@cerebro-database:27017/?authSource=admin" \
  --backup 20260617T031500Z --drop
```

The script:

1. Resolves `latest` (ids are sortable timestamps) and locates the archive +
   manifest.
2. **Verifies the archive against the manifest BLAKE3 checksum before touching
   MongoDB** — a corrupt or tampered archive aborts the restore.
3. Prints the manifest's `audit_chain_verified` and warns if the chain was not
   confirmed intact when the backup was taken.
4. For a real restore, requires typed confirmation, then runs `mongorestore
   --gzip --archive` (optionally `--drop`).
5. Reminds you to **re-verify the audit chain after loading** (chain-verify-after):
   fetch the audit trail and confirm `data.verified` is `true`. If the live chain
   and the manifest disagree, stop and investigate — this is the guard against a
   restore silently laundering a broken chain.

Restore requires `jq`, `b3sum`, and `mongorestore` on the host running it.

## What's testable here vs in your environment

Unit tests cover the parts that don't need a live database: object-store
round-trip, the manifest JSON shape, `mongodump` argument construction, backup-ref
parsing, and the retention selection math. The actual `mongodump`/`mongorestore`,
the chain verify API call, and the image's `mongodb-database-tools` install are
validated against a running stack.

## Limits / follow-ons

- **Point-in-time recovery.** Scheduled `mongodump` gives daily granularity.
  Continuous PITR needs a MongoDB **replica set** with oplog capture — documented
  as an HA upgrade, not enabled here.
- **Encryption at rest.** Backups are checksummed but not encrypted by the engine;
  encrypt the underlying store (or the S3 bucket) for confidentiality.
- **Backup-user auto-provisioning** from the deploy CLI is **delivered**: the
  configure/render step generates the `cerebro_backup` user and writes the backup
  URI secret, so backups are active by default. Override by replacing the rendered
  secret.
