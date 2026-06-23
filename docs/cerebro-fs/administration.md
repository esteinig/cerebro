# Administration guide

How to deploy, configure, and administer a cerebro-fs deployment. For the system
model see the [architecture](architecture.md); for day-to-day maintenance see the
[maintenance guide](maintenance.md); for recovery see the
[disaster-recovery runbook](disaster-recovery.md).

## Deployment

cerebro-fs runs as a `docker compose` project. The deploy CLI renders the compose
file, service environment, and secrets from a single `stack.toml` configuration, then
brings the stack up. Re-running the configure/render step is idempotent — it only
generates secrets that are missing — so a deployment can be reconfigured in place.

The services and their roles are listed in the [architecture](architecture.md#components).
The state-bearing services — MongoDB, the SeaweedFS volumes, and the external backup
and cold stores — are the ones whose durability matters; the rest are stateless and
can be recreated.

### Render and bring up the stack

The deploy CLI (`cerebro stack deploy`) renders the compose file, env files, and
secrets into `--outdir`; you then bring the project up with `docker compose`. Both
steps are safe to re-run (render only fills in what is missing).

```bash
# 1. Render the stack into an output directory (single-host replicated FS model).
#    Replication 001 needs two volume servers on physically separate disks.
cerebro stack deploy \
  --name prod --outdir /opt/cerebro/stack \
  --config localhost-fs \
  --fs-model single-server-replicated \
  --fs-hot /mnt/ssd0/hot         --fs-cold /mnt/hdd0/cold \
  --fs-replica-hot /mnt/ssd1/hot --fs-replica-cold /mnt/hdd1/cold

# 2. Bring the compose project up (the rendered file lives under <outdir>/docker/).
docker compose -f /opt/cerebro/stack/docker/docker-compose.yml up -d

# 3. Confirm every service is healthy.
docker compose -f /opt/cerebro/stack/docker/docker-compose.yml ps
```

`--config` selects a base template (`localhost`, `localhost-insecure`, `localhost-fs`,
`web`); `-f <stack.toml>` loads a prepared configuration file instead. The FS model and
disk paths are detailed in [storage & replication](fs-replication.md). For a single-copy
(non-durable) layout use `--fs-model single-server`; the durability gate refuses a
replication code the rendered topology cannot satisfy.

> The compose commands below assume you are in the stack's `docker/` directory (or
> pass `-f <outdir>/docker/docker-compose.yml` each time). All `docker compose …`
> examples in this documentation set use the service names from
> [§4 of the recovery runbook](disaster-recovery.md#4-system--state-inventory).

### Operator CLI setup (authentication)

Day-to-day operations go through the `cerebro-client` CLI against the API. Point it at
the API and obtain a token once; the maintenance/job commands require an **admin**
token (the admin account is created at deploy time).

```bash
export CEREBRO_API_URL="https://api.<your-domain>"     # default: http://localhost:8080
# Log in; prints the access token to stdout, which we capture into the env var the
# CLI reads. (Use a real admin email/password from your deployment.)
export CEREBRO_API_TOKEN="$(cerebro-client login --email admin@cerebro --password '****')"
cerebro-client ping-server                              # confirms the token works
```

Alternatives: write the token to a file with
`cerebro-client --token-file token.json login --email admin@cerebro` (connection flags
precede the subcommand) and point later calls at it via `-f/--token-file` (or
`CEREBRO_API_TOKEN_FILE`); add `--bot` to the `login` subcommand to obtain a `Role::Bot`
token for service automation. A full task-by-task tour is in the
[operations walkthrough](walkthrough.md).

## Secrets

The deploy step generates and renders the credentials the stack needs, so an operator
does not hand-manage them:

- **Service and worker bot credentials** — the internal accounts the API and worker
  use to authenticate to each other.
- **Backup database URI** — a dedicated, least-privilege MongoDB account for the
  catalogue backup, written as a rendered secret. Catalogue backups are therefore
  active by default; override the credential by replacing the rendered secret.

Secrets are written into the deployment's secret tree and referenced by the services;
they are not committed to configuration. Rotate a secret by replacing its rendered
file and restarting the consuming service:

```bash
# Example: rotate the catalogue-backup MongoDB URI (read-only backup user).
# Replace the rendered secret (one line, no trailing newline) …
printf 'mongodb://cerebro_backup:<new-password>@cerebro-database:27017/?authSource=admin' \
  > /opt/cerebro/stack/mongodb/backup_mongo_uri.secret
# … then restart the service that consumes it.
docker compose restart cerebro-worker
```

The same pattern applies to any rendered secret: overwrite the file, restart the
consumer. See [catalogue backup](catalogue-backup.md#the-backup-user--secret-auto-provisioned)
for the backup user specifically.

## Configuration reference

Configuration is environment-driven. The variables below are grouped by subsystem;
the per-topic docs ([catalogue backup](catalogue-backup.md),
[archival](archival.md)) carry the exhaustive defaults for their areas.

### Storage client (cerebro-fs)

| Variable | Meaning | Default |
|---|---|---|
| `CEREBRO_FS_URL` | SeaweedFS master URL | `http://localhost` |
| `CEREBRO_FS_PORT` | SeaweedFS master port | `9333` |
| `CEREBRO_FS_FILER_URL` | SeaweedFS filer base URL | filer default |
| `CEREBRO_FS_ACCESS` | Object I/O mode: `filer` (path-addressed) or `weed` (fid-addressed) | per deployment |
| `CEREBRO_DANGER_ACCEPT_INVALID_TLS_CERTIFICATE` | Accept invalid TLS (dev only) | `false` |

### Worker

| Variable | Meaning |
|---|---|
| `CEREBRO_API_URL` | API base URL the worker calls |
| `CEREBRO_API_TOKEN` / `CEREBRO_API_TOKEN_FILE` | API token (inline or file) |
| `CEREBRO_API_BOT_EMAIL` / `CEREBRO_API_BOT_PASSWORD` | Bot login credentials |
| `CEREBRO_TEAM` | Team scope for the worker, when applicable |
| `CEREBRO_WORKER_QUEUES` | Faktory queues the worker consumes |
| `CEREBRO_WORKER_VERIFY_ON_MOVE` | Re-verify an object's hash after a tier move |
| `CEREBRO_WORKER_METRICS_ADDR` | Address to expose worker metrics on |
| `CEREBRO_WORKER_RELOGIN_SECS` | Re-authentication interval |
| `CEREBRO_RESTORE_SIMULATE_SECONDS` | Dev/test seam: simulate restore latency |

### Catalogue backup

| Variable | Meaning | Default |
|---|---|---|
| `CEREBRO_BACKUP_MONGO_URI` | Read-only backup DB URI (rendered as a secret) | required to enable |
| `CEREBRO_BACKUP_STORE_PATH` | Backup store root | required to enable |
| `CEREBRO_BACKUP_MONGO_DB` | DB to dump; empty = whole instance | empty (whole instance) |
| `CEREBRO_BACKUP_PREFIX` | Key prefix for backups | `catalogue` |
| `CEREBRO_BACKUP_KEEP_LAST` | Backups to retain | see [catalogue-backup.md](catalogue-backup.md) |
| `CEREBRO_BACKUP_MIN_AGE_DAYS` | Age floor before a backup may be pruned | see [catalogue-backup.md](catalogue-backup.md) |

### Cold archive

| Variable | Meaning | Default |
|---|---|---|
| `CEREBRO_ARCHIVE_STORE_PATH` | Filesystem/NFS cold store root | selects filesystem backend |
| `CEREBRO_ARCHIVE_S3_ENDPOINT` / `_BUCKET` / `_REGION` / `_ACCESS_KEY` / `_SECRET_KEY` | S3-compatible cold store (needs the `s3` build feature) | selects S3 backend |
| `CEREBRO_ARCHIVE_PREFIX` | Key prefix for archived objects | `archive` |
| `CEREBRO_ARCHIVE_LOCAL_GRACE_DAYS` | Days after archival before the local copy is reclaimed | `7` |

If neither an archive store path nor an S3 endpoint is set, tiering is logical only
(no real cold copies), and the archival and reclaim jobs are no-ops.

## Durability prerequisites

Three settings determine what the system can recover; confirm them per deployment:

1. **Replication factor** — set so every live object has more than one copy (see
   [storage & replication](fs-replication.md)). This is what protects non-archived
   data from a single disk/volume loss.
2. **Backup store durability** — `CEREBRO_BACKUP_STORE_PATH` must be on a target
   durable independently of the SeaweedFS data (separate disk, offsite, or object
   lock).
3. **Cold store durability** — the cold archive store is the *only* copy of a
   reclaimed archived object, so it too must be independently durable.

The [validation guide](validation.md) and the
[disaster-recovery prevention checklist](disaster-recovery.md#7-prevention-checklist)
cover how to confirm these.

## Routine administration

- **Monitoring.** The worker exposes a health probe and Prometheus metrics on
  `CEREBRO_WORKER_METRICS_ADDR` (default `0.0.0.0:9464`). Watch the integrity-verify
  failure signal and per-job success/failure rates, and alert on a volume that drops
  below its replica count or a backup that fails to verify.

  ```bash
  # Liveness, and the two metric families that matter most.
  curl -sS http://<worker-host>:9464/health           # → ok
  curl -sS http://<worker-host>:9464/metrics | grep -E \
    'cerebro_file_lifecycle_ops_total|cerebro_worker_jobs_total'
  # Integrity-verify failures are the op=verify, outcome=failure series; job health is
  # cerebro_worker_jobs_total{outcome="failure"}. SeaweedFS volume health comes from the
  # master metrics (cerebro-fs-master:9324) — see fs-replication.md#monitoring-under-replication.
  ```

- **On-demand operations.** Maintenance jobs are triggered with the CLI
  (`cerebro-client jobs launch-job …`, admin token required) and polled with
  `cerebro-client jobs job-status`. The full catalogue of jobs, arguments, and the
  equivalent raw API calls is in the [maintenance guide](maintenance.md); a worked
  end-to-end tour is in the [operations walkthrough](walkthrough.md).

  ```bash
  # Force a catalogue backup now and wait for it to finish.
  id=$(cerebro-client jobs launch-job --kind catalogue_backup --args '{}' --queue maintenance)
  cerebro-client jobs job-status --id "$id"
  ```

- **Recovery.** All recovery procedures, and the drills that prove them, are in the
  [disaster-recovery runbook](disaster-recovery.md).
