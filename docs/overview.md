# cerebro-fs — overview

cerebro-fs is the storage, lifecycle, and durability layer of the Cerebro clinical
metagenomic diagnostics platform. It manages the sequencing data and derived files a
diagnostic laboratory produces — from ingestion through tiering, archival, and
eventual retention — while guaranteeing that every object stays **intact**,
**durable**, **recoverable**, and **auditable** for as long as policy requires.

## The problem it solves

Clinical metagenomics generates large sequencing read sets and analysis outputs that
must be kept under strict clinical-data obligations: nothing may be silently lost or
corrupted, every change must be accountable, retention and legal holds must be
honoured, and the laboratory must be able to recover from hardware, data, and
operator failures with confidence. Doing this by hand across a file store and a
database does not scale and is error-prone. cerebro-fs makes these obligations a
property of the system rather than a manual discipline.

## What it provides

- **Object storage** for files of any size, replicated across volumes, addressed by a
  stable path or identifier, with a content hash recorded for every object.
- **A catalogue** that records each file's identity, location, content hash, storage
  tier, retention, legal hold, and an append-only audit trail of everything that
  happened to it.
- **Lifecycle management** — scheduled movement of data between hot, warm, and cold
  tiers, archival of cold data to an independent object store, and reclamation of
  redundant local copies once it is safe.
- **Integrity assurance** — a rolling re-verification of stored bytes against their
  recorded hash, with automatic repair from a cold backup when one exists.
- **Consistency reconciliation** — periodic detection of divergence between the
  catalogue and the store (missing objects, untracked objects), surfaced for review
  and gated cleanup.
- **Backup** of the whole catalogue to an independent store on a schedule, with a
  verifiable manifest and an audit-chain integrity check.
- **Disaster recovery** — a documented runbook and rehearsal drills that tie
  replication, backups, reconciliation, archival, restore, and repair into proven
  recovery procedures.

## Design principles

- **Correctness over convenience.** Integrity is checked with content hashes at every
  boundary; a corrupt copy can never silently overwrite a good record.
- **Destructive actions are gated.** Operations that delete or repoint data require an
  explicit operator confirmation, or proceed only after a grace window *and* a
  positive confirmation that a replacement copy exists.
- **Report first, act second.** Detection passes record findings without mutating;
  repair and reclamation act on confirmed, transient-filtered signals.
- **Independent durability.** The backup store and the cold archive store are
  expected to be durable on their own, separate from the primary object store, so a
  single failure domain cannot take out both the data and its backup.
- **Operator control and auditability.** Every consequential action is recorded in an
  append-only audit trail, and operators retain explicit control over destructive
  recovery steps.

## Guarantees, and their limits

cerebro-fs is honest about what it can and cannot recover:

- A live object is protected by **replication**; losing one copy heals from another.
- A cold-archived object is protected by the **cold store**, which must be
  independently durable.
- The **catalogue** is protected by scheduled, verifiable backups; recovery is to the
  most recent verified backup point.

The one combination with no recovery source is an object that was never replicated
*and* never archived and then lost — prevention there lives in the replication factor
and the archival cadence, which the administration guide explains how to set.

## Audience and entry points

- **Operators** running or recovering a deployment start with the
  [disaster-recovery runbook](disaster-recovery.md) and the
  [administration guide](administration.md).
- **Engineers** understanding the system start here, then the
  [architecture](architecture.md) and the per-topic operational docs.
- **Anyone validating a deployment** uses the [validation guide](validation.md).

See the [documentation index](README.md) for the full map.
