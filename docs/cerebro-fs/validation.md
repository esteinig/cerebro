# Validation guide

How to validate the `feat/production` changes after applying the stage patches:
build, run the unit tests, understand what those tests do and do **not** cover, and
run the environment-side smoke checks and recovery drills that close the gap.

> Why this matters: the code in these stages was written to be **correct by
> inspection** and delivered as patches; the author could not compile or run it. The
> unit tests below are your first gate, and the environment smoke checks + DR drills
> are the second. Nothing is "validated" until both have run on your stack.

## 1. Build

From the workspace root:

```bash
cargo check --workspace
# the S3 cold-store backend is feature-gated and pulls rust-s3 (pin per platform):
cargo check -p cerebro-worker --features s3
```

`mongodb-database-tools` (for `mongodump`/`mongorestore`) is installed in the worker
image, not required to build.

## 2. Run the unit tests

```bash
cargo test --workspace
# or per crate:
cargo test -p cerebro-model
cargo test -p cerebro-fs
cargo test -p cerebro-worker
```

The unit tests are hermetic: no MongoDB, SeaweedFS, Faktory, or network. Filesystem
tests use a unique `temp_dir()` subtree and clean up after themselves.

## 3. What the unit tests cover (the inventory)

Roughly 128 test functions across the workspace. The tests deliberately target the
**pure** logic — parsing, predicates, key derivation, retention selection — where a
small deterministic test is decisive. The table below maps the production-relevant
modules.

| Module | Unit-tested logic | Tests |
|---|---|---|
| `model/api/files/model.rs` | `effective_identifier` (path vs fid), legal-hold/expiry rules, fixtures | 8 |
| `fs/src/filer.rs` | object-URL joining; **filer-listing parser** (dir bit, mtime parse, pagination, defaults) | 7 |
| `fs/src/client.rs` | `is_filer_path`, remote-path building, segment sanitisation | 6 |
| `worker/backup.rs` | `path_for` traversal-safety, FS store round-trip, **`exists`** (FS + trait default), backup-ref parsing, retention selection, manifest round-trip | 7 |
| `worker/reconcile.rs` | `find_dangling`, `find_orphans`, `within_grace` (the reconcile engine) | 3 |
| `worker/archive.rs` | `archive_key_for` (namespacing, sanitisation, determinism) | 2 |
| `worker/runners/archive_reclaim.rs` | **`is_reclaim_candidate`** gate (archived + key + grace) | 3 |
| `worker/runners/verify.rs` | `retrievable`, **`is_repairable`**, `is_verify_due` rotation, arg defaults | 4 |

## 4. The boundary: unit-tested vs environment-validated

Unit tests **cannot** cover the parts that talk to infrastructure. Those are
correct-by-inspection here and must be validated on a running stack:

- **SeaweedFS I/O** — weed/filer upload, download, list, delete, `exists`, and volume
  operations. The listing *parser* is unit-tested; the *HTTP walk* is not.
- **MongoDB** — `mongodump`/`mongorestore` execution, the catalogue/audit queries,
  CAS-guarded lifecycle and relocate updates. The `mongodump` *arg vector* is
  unit-tested; the *process* is not.
- **The S3 cold-store backend** — feature-gated (`--features s3`); needs a real
  S3/MinIO endpoint and a `rust-s3` version check.
- **Actix handlers** — the relocate endpoint and `POST /jobs/enqueue` are
  correct-by-inspection over the audited lifecycle pattern; exercise them against a
  running API.
- **End-to-end job flows** — `tier_move`, `restore_drive`, `archive_reclaim`,
  `verify_repair`, `catalogue_backup` wire the pure engines to the clients; the
  engines are tested, the wiring is validated by the smoke checks below.

Each component doc restates its own slice of this boundary in its "tested here vs
your environment" section.

## 5. Environment smoke checks

After deploying the built images, confirm the happy paths before relying on the
maintenance automation:

1. **Health.** All compose services up; the API answers; `weed shell` `volume.list`
   shows volumes at their target replica count.
2. **Round-trip + integrity.** Upload a file, confirm it registers with a hash, then
   enqueue `verify_file` for it (`POST /jobs/enqueue`) and confirm it logs
   `integrity verified`.
3. **Backup + verify.** Enqueue `catalogue_backup`; confirm a new
   `{prefix}/{backup_id}/` appears with `catalogue.archive.gz` + `manifest.json`,
   that the archive's BLAKE3 matches `manifest.archive_blake3`, and that
   `audit_chain_verified` is `true`.
4. **Reconcile.** Enqueue `reconcile_scan`; confirm a report lands under `reconcile/`
   with sane catalogue/store counts and low, explained divergence.
5. **Archival round-trip** (if a cold store is configured). Let a file archive, then
   restore it and confirm the restored bytes match the catalogue hash. This exercises
   archival, the retained backup pointer, and restore together.

The commands for each step (admin token set per
[administration → operator CLI setup](administration.md#operator-cli-setup-authentication)):

```bash
# 1. Health.
docker compose ps
curl -sS http://<worker-host>:9464/health
docker compose exec cerebro-fs-master weed shell <<'EOF'
volume.list
EOF

# 2. Round-trip + integrity (FILE_ID from the catalogue listing).
id=$(cerebro-client jobs launch-job --kind verify_file --args '{"file_id":"<FILE_ID>"}' --queue maintenance)
cerebro-client jobs job-status --id "$id"      # then check logs for "integrity verified"

# 3. Backup + verify.
id=$(cerebro-client jobs launch-job --kind catalogue_backup --args '{}' --queue maintenance)
cerebro-client jobs job-status --id "$id"      # verify the manifest as in catalogue-backup.md

# 4. Reconcile.
cerebro-client jobs launch-job --kind reconcile_scan --args '{"budget":100000}' --queue maintenance

# 5. Archival round-trip — request a restore, then hash-check (see disaster-recovery.md §8).
curl -sS -X POST "$CEREBRO_API_URL/files/<FILE_ID>/restore?team=<TEAM>" \
  -H "Authorization: Bearer $CEREBRO_API_TOKEN" -H 'Content-Type: application/json' \
  -d '{"target":"Requested"}'
```

## 6. Recovery drills

The recovery paths are validated by rehearsal, not by unit tests. Run the drills in
[`disaster-recovery.md §6`](disaster-recovery.md#6-routine-recovery-drills) on
cadence — catalogue restore-to-scratch, archived-file restore, integrity repair,
replica heal, and the annual game-day — and keep the sign-off log current. That log
is the evidence that recovery works.

## 7. Definition of done (per deployment)

- [ ] `cargo check --workspace` (and `--features s3` if used) is clean.
- [ ] `cargo test --workspace` passes.
- [ ] §5 smoke checks pass on the deployed stack.
- [ ] At least the catalogue restore-to-scratch and archived-file restore drills
      (§6) have a recent sign-off.
