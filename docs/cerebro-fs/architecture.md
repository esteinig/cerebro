# cerebro-fs — architecture

A higher-level technical view of how cerebro-fs is put together: its components, how
data flows through them, the storage and catalogue model, and the maintenance
subsystem that enforces the system's guarantees. For the conceptual summary see the
[overview](overview.md); for operating it see the [administration](administration.md)
and [maintenance](maintenance.md) guides.

## Components

```
            ┌──────────────────────────────────────────────┐
            │                Web app / REST API            │
            │     (clinical users, programmatic clients)   │
            └───────────────┬──────────────────────────────┘
                            │
        ┌───────────────────┼─────────────────────┬──────────────────┐
        ▼                   ▼                     ▼                  ▼
 ┌──────────────┐   ┌──────────────────┐   ┌────────────┐   ┌────────────────┐
 │   MongoDB    │   │    SeaweedFS     │   │  Faktory   │   │  (clients call │
 │  catalogue   │   │  master /        │   │  job queue │   │   the API)     │
 │  audit chain │   │  volumes /       │   └─────┬──────┘   └────────────────┘
 │  filer meta  │   │  filer           │         │
 └──────────────┘   └──────────────────┘   ┌─────▼──────┐
                            ▲               │   Worker   │
                            │               │ maintenance│
                            └───────────────┤   jobs     │
                                            └─────┬──────┘
                                ┌─────────────────┴──────────────┐
                                ▼                                ▼
                        ┌───────────────┐                ┌───────────────┐
                        │  Backup store │                │  Cold store   │
                        │  (catalogue   │                │  (archived    │
                        │   dumps)      │                │   objects)    │
                        └───────────────┘                └───────────────┘
```

- **Web app / REST API** — the stateless front door. Clients upload and retrieve
  files and manage their lifecycle; administrators trigger maintenance and recovery
  operations.
- **MongoDB** — the durable system of record: the file **catalogue**, the append-only
  **audit chain**, and (typically) the SeaweedFS **filer metadata**.
- **SeaweedFS** — the object store: a **master** tracks topology, **volume** servers
  hold replicated object bytes, and the **filer** gives objects stable paths backed by
  MongoDB.
- **Faktory + worker** — the job queue and the process that runs scheduled and
  on-demand maintenance: tiering, verification, reconciliation, backup, archival,
  reclamation, restore, and repair.
- **Backup store** and **cold store** — two independently durable targets, one for
  catalogue backups and one for archived object bytes. They may be the same physical
  system under different prefixes, or separate, but must be durable independently of
  the SeaweedFS data.

## Storage model

Every object is stored once in SeaweedFS and described once in the catalogue.

- **Addressing.** An object is reached by either a **filer path** (a stable,
  human-meaningful path such as `/run/sample/reads_R1.fastq.gz`) or a **volume
  identifier** (an opaque fid). The catalogue records whichever applies; retrieval
  prefers the path when present and falls back to the fid.
- **Replication.** Volumes are replicated so that every live object has more than one
  copy. The replication factor is a configured property of the cluster; losing one
  copy heals from another.
- **Integrity.** A BLAKE3 content hash is recorded for every object at ingestion and
  is the reference for all later verification and repair.
- **Tiers.** Each file carries a logical storage tier (hot / warm / cold). Tiering is
  a catalogue property: a non-archived file is directly retrievable regardless of its
  tier. The cold tier becomes *physical* only when a file is **archived** — its bytes
  copied to the cold store and the catalogue marked with an archive key.

## Catalogue model

The catalogue record for a file binds together its identity and its policy:

- **Identity** — id, name, content hash, location (path and/or fid), and provenance
  (run, sample, file type, tags).
- **Lifecycle** — storage tier, archival state (whether archived and where its cold
  copy lives), restore state, and timestamps for the last tier move and the last
  integrity verification.
- **Policy** — retention class, retain-until date, and legal hold.
- **Accountability** — an append-only audit trail recording every consequential
  action (upload, tier move, archival/restore, retention, legal-hold and tag changes,
  deletion). The chain is verifiable and its integrity is checked at backup time.

## Data lifecycle

```
 upload → register → [ tier moves ] → archive → reclaim local copy
   │         │            │              │             │
   └─ hash   └─ catalogue └─ hot/warm/   └─ cold store  └─ space reclaimed,
      recorded   entry       cold label     copy + key     cold copy authoritative

 throughout:  verify (re-hash) ──► repair from cold on mismatch
              reconcile (catalogue ↔ store) ──► gated cleanup
              retention sweep ──► expiry under policy + legal hold
              catalogue backup ──► verifiable dump to the backup store
```

1. **Ingestion.** A file is uploaded, its content hashed, and a catalogue record
   created.
2. **Tiering.** Scheduled jobs move bytes between tiers according to policy.
3. **Archival.** Cold files are copied to the cold store and repointed; after a grace
   window their redundant local copy is reclaimed, leaving the cold copy
   authoritative and reclaiming primary-store space.
4. **Verification and repair.** A rolling pass re-hashes stored objects; a mismatch on
   a file that has a cold backup is repaired by re-fetching the verified-good bytes.
5. **Reconciliation.** A periodic pass compares the catalogue against the store and
   reports divergence — catalogue entries whose object is missing, and stored objects
   with no catalogue entry — for operator-gated cleanup.
6. **Retention.** A sweep expires data whose retention has elapsed, honouring legal
   holds.
7. **Backup.** The whole catalogue is dumped to the backup store on a schedule, with a
   manifest that records the dump's hash and confirms the audit chain was intact.

## Maintenance subsystem

The guarantees above are enforced by jobs the worker runs from the queue. They are
**seeded** to run on a staggered schedule and can also be triggered **on demand** by
an administrator. They divide into:

- **Detection (read-only):** reconcile scan, verify scan — find problems, record them,
  never mutate.
- **Safe automation:** tiering, archival, local-copy reclaim (grace + presence
  gated), integrity repair (hash gated), catalogue backup, retention sweep.
- **Operator-gated:** consistency reclaim (explicit confirmation and an explicit list)
  — the deletion of untracked store objects is never automatic.

The [maintenance guide](maintenance.md) catalogues each job, its cadence, and how to
run it; the [disaster-recovery runbook](disaster-recovery.md) covers the recovery
procedures that use them.

## Durability and recovery model

Each class of state has one authoritative recovery source:

| State | Protected by | Recovered from |
|---|---|---|
| Live object bytes | Volume replication | A surviving replica |
| Archived object bytes | The cold store (independently durable) | The cold store |
| Catalogue + audit chain | Scheduled verifiable backups | The latest verified backup |
| Filer metadata | Co-located in the catalogue backup (typical) | The catalogue restore |

The recovery procedures, their verification steps, and the rehearsal drills that prove
them are in the [disaster-recovery runbook](disaster-recovery.md).
