# Cerebro production docs

Operator and developer documentation for cerebro-fs — the storage, lifecycle, and
durability layer of Cerebro, the clinical metagenomic diagnostics platform. This is
the entry point; it orients you and points at the rest of the documentation set.

## Where to start

- **New to cerebro-fs** → the [overview](overview.md), then the
  [architecture](architecture.md).
- **Deploying / configuring** → the [administration guide](administration.md).
- **Running maintenance** → the [maintenance guide](maintenance.md).
- **Operating / recovering a deployment** → the
  [disaster-recovery runbook](disaster-recovery.md) (triage → scenario → commands →
  verify).
- **Validating a deployment** → the [validation guide](validation.md).
- **Planning the road to production** → the
  [production-readiness assessment](production-readiness.md).

## The stack

A `docker compose` deployment (service names in `disaster-recovery.md §4`):

- **API / app** (`cerebro-api`, `cerebro-app`) — REST + web UI (stateless).
- **MongoDB** (`cerebro-database`) — the catalogue, audit chain, and (typically) the
  filer metadata.
- **SeaweedFS** (`cerebro-fs-master`, `-primary`, `-replica`, `-secondary`, `-filer`)
  — object storage, replicated across volumes; the filer gives objects stable paths.
- **Faktory** + **worker** (`faktory`, `cerebro-worker`) — the maintenance jobs
  (tiering, verification, reconcile, backup, archival, restore, repair).

State that must be durable: MongoDB, the SeaweedFS volumes, and the external
**backup** and **cold** stores. The last two must be durable independently of the
SeaweedFS data.

## The file/data lifecycle

The spine that ties the docs together — a file moves through:

1. **Upload + register.** Bytes land in SeaweedFS (filer path or weed fid); a
   catalogue record is created with a BLAKE3 hash.
2. **Lifecycle + retention.** The record carries tier, retention, legal hold, and an
   append-only audit trail.
3. **Tiering.** Scheduled jobs move bytes between hot/warm/cold tiers.
4. **Archival.** Cold-tier files are copied to the **cold store** and repointed
   (`archived`, `archive_key`); after a grace window the redundant local copy is
   reclaimed — [`archival.md`](archival.md).
5. **Integrity + repair.** A rolling verify re-hashes objects; a mismatch on a file
   with a cold backup is **repaired from cold** — [`verify-repair.md`](verify-repair.md).
6. **Consistency.** A reconcile pass finds catalogue↔store divergence (dangling
   references, orphan objects) for gated cleanup —
   [`consistency-reconcile.md`](consistency-reconcile.md).
7. **Backup.** The whole catalogue is dumped to the backup store on a schedule, with
   a verifiable manifest — [`catalogue-backup.md`](catalogue-backup.md).
8. **Recovery.** When something fails, the [`disaster-recovery.md`](disaster-recovery.md)
   runbook ties replication, backups, reconcile, archival, restore, and repair into
   recovery drills.

## Documentation map

**Overview**

| Doc | Covers |
|---|---|
| [`overview.md`](overview.md) | What cerebro-fs is, the problem it solves, capabilities, guarantees |
| [`architecture.md`](architecture.md) | Components, storage and catalogue model, data lifecycle, durability model |

**Operations**

| Doc | Covers |
|---|---|
| [`administration.md`](administration.md) | Deploying, configuring (configuration reference), secrets, durability prerequisites |
| [`maintenance.md`](maintenance.md) | The scheduled-job catalogue and on-demand operations |
| [`fs-replication.md`](fs-replication.md) | SeaweedFS replication; multi-host durability |
| [`catalogue-backup.md`](catalogue-backup.md) | Scheduled catalogue backups; manifest; backup user |
| [`consistency-reconcile.md`](consistency-reconcile.md) | Dangling/orphan detection; store enumeration; gated reclaim |
| [`archival.md`](archival.md) | Cold-tier archival; local-copy reclaim |
| [`verify-repair.md`](verify-repair.md) | Restore-from-cold; integrity repair |
| [`disaster-recovery.md`](disaster-recovery.md) | Operator recovery runbook + drills |

**Reference**

| Doc | Covers |
|---|---|
| [`validation.md`](validation.md) | Build, test, smoke checks, the test boundary |
| [`production-readiness.md`](production-readiness.md) | Honest readiness assessment, tracked gaps, and the hardening roadmap |

## Conventions

- Each component doc ends with a **"tested here vs your environment"** note — the
  honest boundary between what unit tests cover and what only a running stack can
  validate. [`validation.md`](validation.md) consolidates this.
- Destructive operations are **operator-gated** or **grace + presence-gated**; the
  runbook's first step pauses them during an incident.
