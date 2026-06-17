# Cerebro production docs

Operator and developer documentation for the `feat/production` hardening of Cerebro,
the clinical metagenomic diagnostics platform. This is the entry point; it orients
you and points at the per-topic docs.

## Where to start

- **Operating / recovering a deployment** → [`disaster-recovery.md`](disaster-recovery.md)
  is the top-level runbook (triage → scenario → commands → verify).
- **Understanding a subsystem** → the component docs below.
- **Validating after applying patches or before go-live** →
  [`validation.md`](validation.md) (build, test, smoke checks, drills).

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
   append-only audit trail (Stage 2).
3. **Tiering.** Scheduled jobs move bytes between hot/warm/cold tiers (Stage 3).
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

## Component docs

| Doc | Covers | Stage |
|---|---|---|
| [`fs-replication.md`](fs-replication.md) | SeaweedFS single-host replication; multi-host durability | S4-1 |
| [`catalogue-backup.md`](catalogue-backup.md) | Scheduled `mongodump` to the cold store; manifest; backup user | S4-2, H1 |
| [`consistency-reconcile.md`](consistency-reconcile.md) | Dangling/orphan detection; store enumeration; gated reclaim | S4-3, H2 |
| [`archival.md`](archival.md) | Real cold-tier archival; local-copy reclaim | S4-4, H3 |
| [`verify-repair.md`](verify-repair.md) | Restore-from-cold; integrity repair | S4-5, H4 |
| [`disaster-recovery.md`](disaster-recovery.md) | Operator recovery runbook + drills | S4-7 |
| [`validation.md`](validation.md) | Build, test, smoke checks, the test boundary | S5 |

## Stage history

- **Stage 1** — `cerebro-fs`: SeaweedFS object storage integration.
- **Stage 2** — lifecycle, retention, provenance, append-only audit, telemetry.
- **Stage 3** — Faktory workers: tiering, verification, the restore state machine,
  retention sweep.
- **Stage 4** — backup & recovery: replication (S4-1), catalogue backup (S4-2),
  consistency reconcile (S4-3), cold archival (S4-4), verify-repair / restore (S4-5),
  hardening (S4-6: H1 backup-user provisioning, H2 orphan enumeration, H3 local-copy
  reclaim, H4 repair-from-cold), and the DR runbook (S4-7).
- **Stage 5** — documentation & testing hardening: unit tests over the pure logic,
  drift fixes, this index, and the validation guide.

## Conventions

- Each component doc ends with a **"tested here vs your environment"** note — the
  honest boundary between what unit tests cover and what only a running stack can
  validate. [`validation.md`](validation.md) consolidates this.
- Destructive operations are **operator-gated** or **grace + presence-gated**; the
  runbook's first step pauses them during an incident.
