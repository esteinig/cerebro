# Consistency reconcile (Stage 4, S4-3)

The catalogue (MongoDB) and the object store (SeaweedFS) can drift apart:

- **Dangling reference** — a catalogue entry whose backing object is missing.
  This is latent data loss: the system believes a file exists and will hand out
  its metadata, but the bytes are gone.
- **Orphan object** — a stored object with no catalogue entry. Wasted capacity,
  and deleting one is irreversible.

S4-3 detects both and reports them. It is **report-first**: detection is automatic
and safe; the one destructive action (deleting an orphan) is operator-gated (D6).
Repair of dangling references — re-replicating from a replica or backup — is S4-5,
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

## `reconcile_reclaim` (operator-gated, destructive)

Deleting an orphan is irreversible, so reclaim is a separate job that is **not**
seeded. An operator enqueues it after reviewing a report:

```json
{ "confirm": true, "keys": ["3,01abcd…", "5,07ef12…"], "max_delete": 100 }
```

It refuses unless `confirm` is `true`, deletes only the explicit `keys` (capped by
`max_delete`), and logs each deletion. It acts on the operator's reviewed list, so
the grace decision is made at report time; it does not re-enumerate.

## Orphan detection — current status

Dangling detection is fully live. Orphan *detection* additionally needs an
**independent enumeration of the object store** (what objects actually exist),
which the SeaweedFS client does not yet expose cheaply — enumerating fid space
means walking the master's volume list and each volume's needles (or, in filer
mode, walking the filer tree). The reconcile engine already computes orphans
(`find_orphans`) from any provided store listing and is unit-tested; until the
enumeration lands, `reconcile_scan` sets `store_enumerated: false` and reports no
orphans, while `reconcile_reclaim` still works on an operator-supplied key list
(e.g. from a manual `weed shell volume.list` audit). Wiring the store enumeration
into the scan is the S4-3 follow-on.

## What's tested here vs in your environment

Unit tests cover the engine: dangling selection (present/absent/archived/grace
cases), orphan selection (store−catalogue with grace), the grace predicate, and
the catalogue-date parsing. The live catalogue enumeration, the `object_exists`
HEAD probe against SeaweedFS, report persistence, and the gated `delete_file`
reclaim are validated against a running stack.

## Design notes

- **Why HEAD, not hash.** Existence is a cheap `HEAD`; hashing (the verify scan,
  S3-3a) is a separate, heavier integrity check. Reconcile answers "does the object
  exist", verify answers "are its bytes intact".
- **Grace window.** Registration happens after a successful upload, so the
  catalogue→object race is small, but the grace window still suppresses the brief
  window around very recent writes on both sides.
- **Report-first split.** Detection is safe and automatic; repair (S4-5) and
  destructive reclaim are deliberately separate and gated, so a scan can run
  unattended without ever risking data.
