# Consistency reconcile

The catalogue (MongoDB) and the object store (SeaweedFS) can drift apart:

- **Dangling reference** — a catalogue entry whose backing object is missing.
  This is latent data loss: the system believes a file exists and will hand out
  its metadata, but the bytes are gone.
- **Orphan object** — a stored object with no catalogue entry. Wasted capacity,
  and deleting one is irreversible.

cerebro-fs detects both and reports them. It is **report-first**: detection is automatic
and safe; the one destructive action (deleting an orphan) is operator-gated.
Repair of dangling references — re-replicating from a replica or backup — is handled by the verify-repair pass,
which consumes these reports.

## `reconcile_scan` (scheduled, safe)

A weekly worker job (`seed:reconcile_scan`, queue `maintenance`):

1. Enumerates the catalogue via the API, paginated, up to a budget.
2. For each **non-archived** entry, cheaply probes object presence with a
   `HEAD {fs}/{fid}` (`object_exists`). Archived entries are skipped — their object
   lives in remote cold (S3 Glacier), legitimately not in the local store, and is
   retrieved through the restore path.
3. Computes dangling references: entries whose object is definitely absent, older
   than a grace window (default 1 day, to avoid racing an in-flight write), and not
   archived. A probe that *errors* is treated as unknown, never as missing — a
   transient fault can't manufacture a false dangling finding.
4. Emits a report: it logs each dangling reference, records metrics, and (when the
   backup object store is configured) writes a JSON report under `reconcile/` in
   that store.

It never mutates catalogue or object state. Tune per-run with job args
`{ "budget": N, "grace_days": D }`.

**Run it on demand and read the report** (admin token — see
[administration → operator CLI setup](administration.md#operator-cli-setup-authentication)):

```bash
# Enqueue a scan (raise the budget to cover a large catalogue and enable orphan detection).
id=$(cerebro-client jobs launch-job --kind reconcile_scan \
  --args '{"budget":100000,"grace_days":1}' --queue maintenance)
cerebro-client jobs job-status --id "$id"

# Read the newest report from the backup store's reconcile/ prefix.
RID=$(ls -1 "$CEREBRO_BACKUP_STORE_PATH/reconcile" | sort | tail -n1)
jq '{store_enumerated, dangling: (.dangling|length), orphans: (.orphans|length)}' \
   "$CEREBRO_BACKUP_STORE_PATH/reconcile/$RID"
jq -r '.orphans[]?.key' "$CEREBRO_BACKUP_STORE_PATH/reconcile/$RID"   # candidate keys for reclaim
```

If `store_enumerated` is `false`, orphans were not computed (volume mode, or the
catalogue exceeded the budget); raise the budget and re-run to enable it.

## `reconcile_reclaim` (operator-gated, destructive)

Deleting an orphan is irreversible, so reclaim is a separate job that is **not**
seeded. An operator enqueues it after reviewing a report:

```json
{ "confirm": true, "keys": ["3,01abcd…", "5,07ef12…"], "max_delete": 100 }
```

It refuses unless `confirm` is `true`, deletes only the explicit `keys` (capped by
`max_delete`), and logs each deletion. It acts on the operator's reviewed list, so
the grace decision is made at report time; it does not re-enumerate.

Run it with the **explicit keys** you confirmed from the report — never a wildcard:

```bash
cerebro-client jobs launch-job --kind reconcile_reclaim \
  --args '{"confirm":true,"keys":["3,01abcd…","5,07ef12…"],"max_delete":100}' \
  --queue maintenance
# Path keys route to the filer delete, fids to the volume delete, automatically.
```

## Orphan detection (live in filer mode)

Both directions are now live. Orphan detection enumerates the object store
independently and finds stored objects with no catalogue entry:

- In **filer mode** (the Cerebro deployment's mode) the scan walks the filer tree
  (`enumerate_objects`, a bounded depth-first listing) and compares the object
  **paths** against the catalogue `path` set.
- A safety rule prevents false positives: orphan detection runs only when the
  catalogue path set is **complete** — i.e. the catalogue enumeration finished
  within the budget. A truncated catalogue would make real objects look orphaned,
  so if the catalogue exceeds the budget the scan skips orphans and logs a note
  (raise the budget to enable it). A truncated *store* walk is safe — it only
  misses some orphans, never invents them, since every key found is checked
  against the full catalogue.
- Found orphans are reported (logged + in the persisted report), never deleted by
  the scan. Reclaim stays operator-gated: `reconcile_reclaim` with `confirm=true`
  and an explicit key list, now routing each key to the filer delete (paths) or the
  volume delete (fids) automatically.

In **weed/volume mode** the store cannot be cheaply enumerated, so the scan leaves
`store_enumerated: false` and reports dangling references only; `reconcile_reclaim`
still works on an operator-supplied key list. The `find_orphans` engine is unchanged
and unit-tested; the scan simply feeds it a real store listing.

## What's tested here vs in your environment

Unit tests cover the engine: dangling selection (present/absent/archived/grace
cases), orphan selection (store−catalogue with grace), the grace predicate, and
the catalogue-date parsing. The live catalogue enumeration, the `object_exists`
HEAD probe against SeaweedFS, report persistence, and the gated `delete_file`
reclaim are validated against a running stack.

## Design notes

- **Why HEAD, not hash.** Existence is a cheap `HEAD`; hashing (the verify scan) is a separate, heavier integrity check. Reconcile answers "does the object
  exist", verify answers "are its bytes intact".
- **Grace window.** Registration happens after a successful upload, so the
  catalogue→object race is small, but the grace window still suppresses the brief
  window around very recent writes on both sides.
- **Report-first split.** Detection is safe and automatic; repair and
  destructive reclaim are deliberately separate and gated, so a scan can run
  unattended without ever risking data.
