# Real archival to the cold object store (Stage 4, S4-4)

Before S4-4 the Cold tier was a catalogue *label*: a file marked `Cold` still had
its only copy in SeaweedFS. S4-4 makes the cold tier real — a file moved to Cold
has its bytes copied to a cold object store, the catalogue repointed, and
`archived` flipped — so hot storage can be reclaimed and the cold copy lives in an
independent failure domain.

## The cold store is an `ObjectStore` (D1)

Archival is backend-agnostic: it writes through the same `ObjectStore` trait used
for catalogue backups (S4-2).

- **Filesystem backend (default, compile-safe).** A mounted cold disk or NFS
  export — slower, cheaper storage in a separate failure domain. Set
  `CEREBRO_ARCHIVE_STORE_PATH`. No extra dependencies.
- **S3 backend (feature-gated).** AWS S3 / MinIO / SeaweedFS-S3 / Glacier. Build
  the worker with `--features s3` and set `CEREBRO_ARCHIVE_S3_ENDPOINT`,
  `CEREBRO_ARCHIVE_S3_BUCKET`, `CEREBRO_ARCHIVE_S3_REGION`, and the access/secret
  keys (`CEREBRO_ARCHIVE_S3_ACCESS_KEY`/`_SECRET_KEY`, `_FILE` variants honoured).

Both satisfy the same trait, so the archival engine and the tier-move integration
are identical regardless of where cold bytes land.

> The S3 backend uses the `rust-s3` blocking API. Its method names/signatures
> shift across minor versions, and it is excluded from the default build, so it is
> the one module that cannot be compile-checked here. Verify it against the pinned
> `rust-s3` version when you enable `--features s3`.

## Configuration

`ArchiveSettings::from_env` is opt-in. With nothing set, tiering keeps its prior
label-only behaviour. An S3 endpoint selects the S3 backend; otherwise a store
path selects the filesystem backend. `CEREBRO_ARCHIVE_PREFIX` (default `archive`)
namespaces keys.

## The archive flow

When `tier_move` commits a file to Cold and a cold store is configured:

1. Read the object's bytes from SeaweedFS (`read_object`, a direct `GET`).
2. Put them to the cold store at a deterministic key,
   `<prefix>/<team>/<file_id>` (`archive_key_for`).
3. Repoint the catalogue through the dedicated relocate endpoint (below): set
   `archived = true` and record `archive_key`.

Archival runs **after** the tier commit and is **non-fatal**: if it fails, the
file is Cold-committed and still directly retrievable, and the next move retries.
The repoint is the moment the file becomes officially archived, so a crash between
copy and repoint just leaves a harmless extra copy in the cold store (the key is
deterministic, so the retry overwrites it).

## The dedicated relocate endpoint (D3)

`POST /files/{id}/relocate` (`FileRelocateSchema`) repoints a file between local
and remote storage in one compare-and-set, kept deliberately separate from the
lifecycle update so the fid repoint — which atomically changes `tier`, `archived`,
`fid`, and `archive_key` — is never tangled with routine edits. It applies only if
the file's current `tier` equals `expected_tier` (a `409` means another mover got
there first), and records a `TierMove` audit event with an `archive`/`restore`
detail. The client wrapper is `CerebroClient::relocate_file`.

The schema carries both directions: archival sends `archived = true` with an
`archive_key`; restore sends `archived = false`, the freshly re-uploaded local
`fid`, and `clear_archive_key`.

## Restore (paired with S4-5)

The relocate endpoint and the catalogue already model the restore direction, and
the restore *state machine* (S3-3b: `NotArchived → Requested → InProgress →
Restored`) is in place. Wiring the real fetch — `ObjectStore::get(archive_key)` →
re-upload to SeaweedFS for a new fid → relocate back to a local tier — is done with
S4-5 (verify-repair), which re-materialises bytes from a replica or the archive
using the same primitives. Until then, archived objects are retrieved through the
existing restore flow / simulation seam.

## Local-copy reclamation (D4, live in H3)

S4-4 archives and repoints but keeps the local copy; the `archive_reclaim` runner
(H3) deletes it once it is safe. A weekly, bounded pass reclaims the local copy of
an archived file when **all** hold: the file is `archived` with an `archive_key`;
its archival (`tier_moved_at`) is older than the grace window
(`CEREBRO_ARCHIVE_LOCAL_GRACE_DAYS`, default 7); the **cold copy is confirmed
present** (`ObjectStore::exists` — the safety gate, so a local copy is never
deleted without its replacement); and a local copy still exists.

It deletes by filer path when present (cleaning filer metadata and data together),
else by fid, and is active only when a cold store is configured. After deletion the
file stays `archived`, so reconcile skips it and retrieval goes through the restore
path; the stale fid is harmless and is replaced on the next restore.

The gate confirms cold-copy *presence*, not integrity. The archival wrote verified
bytes and restore re-verifies the BLAKE3 on the way back, so the residual risk is a
cold object that exists but is corrupt; operators wanting maximum safety can run a
verify pass before enabling reclaim, or raise the grace window.

## What's tested here vs in your environment

Unit tests cover the pure pieces: `archive_key_for` (namespacing, sanitisation,
determinism). The relocate endpoint mirrors the audited, CAS-guarded lifecycle
handler and is correct-by-inspection. The byte movement (`read_object`, the
`ObjectStore` put), the S3 backend, and the end-to-end archive-on-tier-move are
validated against a running stack — and the S3 backend additionally needs
`--features s3` plus a `rust-s3` version check.
