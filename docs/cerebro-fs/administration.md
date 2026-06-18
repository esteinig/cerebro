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
file and restarting the consuming service.

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

- **Monitoring.** The worker exposes metrics (`CEREBRO_WORKER_METRICS_ADDR`); watch
  the integrity-verify failure signal and job success/failure rates, and alert on a
  volume that drops below its replica count or a backup that fails to verify.
- **On-demand operations.** Maintenance jobs can be triggered on demand through the
  admin job endpoint — see the [maintenance guide](maintenance.md).
- **Recovery.** All recovery procedures, and the drills that prove them, are in the
  [disaster-recovery runbook](disaster-recovery.md).
