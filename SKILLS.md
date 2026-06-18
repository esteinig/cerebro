# SKILLS.md — Cerebro feature development

A living playbook for building features in the Cerebro platform. It captures how we
work (a staged, patch-driven, evidence-honest method), what Cerebro is
(architecture), how Rust is written here, and the posture that matters most: clinical
production and regulatory rigour at high throughput. Read it before starting a
feature; update it when the platform or the method moves on.

- **Audience.** Engineers and AI assistants doing Cerebro feature work.
- **Scope.** Development method, architecture orientation, code style, and the
  clinical/regulatory/throughput posture. Not a substitute for the per-feature
  operational docs under `docs/`.
- **How to use.** Start a feature with *The staged method*, hold every change to *Code
  style* and *Clinical posture*, and close with a production-readiness assessment.
- **Default lens.** High-throughput clinical/diagnostic operation. The
  research / small-lab profile is a deliberate, narrow relaxation — never of data
  safety, only of operational burden.

---

## 1. Cerebro at a glance

Cerebro is a clinical metagenomic diagnostics platform: it turns sequencing data into
pathogen-detection results suitable for clinical decision-making, and manages the
data, provenance, and lifecycle around that with the durability and auditability
clinical work demands.

### Architecture (kept generic on purpose)

The platform is a small set of cooperating services. Treat the *roles* below as
stable; treat specific technologies as current choices that live behind abstractions
and may be swapped.

- **REST API server** — the stateless front door for clients and the web app; owns
  request handling, authorization, and orchestration of work.
- **Catalogue & audit database** (MongoDB) — the durable system of record: domain
  records, the append-only audit chain, and service metadata. This is the single
  source of truth; protect it accordingly.
- **Object store** (SeaweedFS: master, volumes, filer) — replicated storage for files
  of any size, addressed by a stable path or an opaque identifier, with a content hash
  recorded for every object.
- **Job queue + workers** (Faktory + a Rust worker) — asynchronous, scheduled and
  on-demand background work (lifecycle, integrity, backup, recovery). Producers scan
  and enqueue; consumers do bounded, idempotent units of work.
- **Web application** (SvelteKit) — the operator/clinician UI over the API.
- **Pipeline orchestration** (Nextflow, plus the in-repo tower/orchestration crate) —
  runs the analysis pipelines that produce results and files, and emits run
  provenance.
- **Deployment tooling** (a CLI that renders a `docker compose` stack and its secrets
  from a single configuration) — one idempotent path from config to a running stack.

### Code organisation

A Rust workspace plus pipelines and the web app:

- `lib/` — shared domain crates (models, API schemas, domain types). The model crate
  is the contract every service depends on; changes here ripple, so they are
  backwards-compatible by default.
- `stack/` — the services (`server`, `worker`, `fs`, `client`, `tower`, `app`, …),
  each a focused crate.
- `cerebro/` — the top-level CLI and deploy/stack orchestration.
- `docs/` — operational and reference documentation, written as a system (see §6).

### Data flow (the mental model)

Sequencing run → pipeline analysis → results + files → API registers them in the
catalogue and stores bytes in the object store → lifecycle, integrity, backup, and
recovery run continuously over them → reporting and retrieval. Provenance and an audit
record follow the data the whole way.

---

## 2. The staged method

This is the working style. It exists to make large, high-stakes changes land safely,
reviewably, and with honest evidence — and to keep development scaffolding out of the
shipped product.

### Stages and sub-stages

- A **stage** is a coherent unit of work with one goal and explicit exit criteria
  (e.g. "backup & recovery", "documentation consolidation"). Capture it in a
  `STAGE-NN-PLAN.md`: goal, decisions, sub-stages, ordering, exit, out-of-scope.
- A **sub-stage** is the smallest independently useful change — one focused commit.
  Each is built on the real branch and delivered as one git patch.

### Locked decisions

Before building, write down the handful of design decisions the stage turns on as
**numbered, locked choices**, each with its rationale and the alternative considered.
This creates a decision record and stops mid-build re-litigation. If reality forces a
change, change the decision explicitly — don't drift.

### One patch per sub-stage, verified to apply

- Build on the actual feature branch; produce one commit per sub-stage with a clear
  message (ID + concise subject + a *what and why* body that references the relevant
  docs/code).
- **Verify mechanically.** Every patch must apply cleanly to its parent
  (`git apply --check` then `git am`) — ideally in a throwaway worktree so the working
  branch never detaches. For a batch, replay the whole sequence onto the base and
  confirm the result.
- For text-only edits (comments, docs), run a balance check (braces, quotes) on every
  edited file; a comment edit must never change code.

### Risk-first, behaviour-preserving discipline

- Order sub-stages so the riskiest or most foundational work comes first.
- Keep behaviour-changing work separate from pure refactors and documentation. When a
  change is meant to preserve behaviour, say so and keep it true.

### Evidence honesty (the boundary)

Never claim more validation than was performed. State plainly what is **unit-tested**
versus what only a **running environment** can validate. If something is written
correct-by-inspection because it can't be compiled/run here, say so. If a claim is
inferred from memory of the code rather than confirmed, flag it for confirmation
rather than asserting it. This honesty is load-bearing in a clinical context.

### Tests: hermetic, house-pattern, no surprise dependencies

- Unit-test the **pure decision logic**. Factor decisions into pure functions
  (predicates, parsers, selectors) that run without I/O — this is the single biggest
  testability lever (see §3).
- Tests are deterministic and self-cleaning (temp dirs + unique ids + teardown; pure
  (de)serialisation; in-test fakes for I/O). Add no new dependencies for tests.
- Integration tests are a deliberate choice, not a default. If they're out of scope
  for a stage, record that — and record it as a production-readiness gap, because the
  I/O and recovery paths are where production breaks.

### Close every feature with two artifacts

- A **production-readiness assessment** (`docs/production-readiness.md` style): an
  honest, critical list of remaining gaps as tracked items with priority tiers, each
  with risk, recommended action, and references. Be a hard marker.
- A short, concise **PR summary** referencing the relevant docs and code.

### Keep scaffolding out of the product

Stages, sub-stage IDs, and decision numbers belong in plans and commit history —
**not** in shipped docstrings, comments, or operational docs. The product reads as a
*system*, not a development timeline. (Removing that scaffolding before release is
itself a worthwhile stage.)

---

## 3. Rust code style

Conventions here favour correctness, testability, and clarity over cleverness.

### Formatting & linting

- `rustfmt` and `clippy` clean; warnings are not left to rot.
- Prefer small, focused modules. Public surfaces are intentional and documented.

### Documentation

- `//!` module docs state the module's job and its invariants; `///` item docs explain
  **what and why**, plus gotchas and edge cases — not the obvious. Document the
  reasoning a future reader will lack.
- Doc examples that can be doctests are welcome where they clarify usage.
- **Doctests compile and run as standalone programs — give them their imports.**
  Bring in everything the example uses with visible `use` lines, or with hidden
  `# `-prefixed setup lines (e.g. `# use crate_name::module::Thing;`) that are compiled
  but not rendered; don't rely on the surrounding module's scope. A doctest missing its
  imports fails to build and so is silently not tested.

### Types & errors

- Push domain meaning into the **shared model crate**; avoid stringly-typed data; use
  newtypes for identifiers and enums for states.
- Typed, enumerated errors (a `thiserror`-style error enum per crate); return
  `Result`; attach context. **No `unwrap`/`expect` on production paths** — only in
  tests or behind an invariant that is documented and locally provable.
- Make illegal states unrepresentable where practical; validate at the boundary.

### Abstraction for I/O and testability

- Put external systems behind **traits** (e.g. an object-store trait, a client
  abstraction) so backends are swappable and code is testable with fakes. Feature-gate
  optional backends.
- **Extract pure logic** out of I/O. A function that decides *whether* to reclaim,
  repair, or move — given plain inputs — is unit-testable and reviewable; the I/O
  wrapper around it stays thin. This is the canonical pattern in this codebase.

### Concurrency & jobs

- Async for I/O and protocol; offload blocking/CPU work to a blocking pool — don't
  stall the async executor.
- Jobs must be **idempotent under at-least-once delivery** and safe to re-run after a
  partial failure. Long waits use poll-by-re-enqueue rather than parking a worker.

### Serialization & compatibility

- Domain records evolve **backwards-compatibly**: new fields use `#[serde(default)]`
  (or an explicit default) so existing documents still deserialize. Never silently
  break the on-disk/on-wire shape.

### Configuration, secrets, logging

- Configuration is explicit and environment-driven; no hidden globals. Document every
  variable in the consolidated configuration reference.
- Never put secrets in code or logs; least-privilege credentials, rendered as secrets,
  rotatable.
- Structured logging and metrics with a shared taxonomy; emit the signals operators
  need to see (see §4).

---

## 4. Clinical, regulatory & production posture (the default lens)

Assume every object may be patient-linked and every action may be audited. The
following are non-negotiable design properties, not features to bolt on later.

### Integrity

- Content-hash every object at ingestion; that hash is the reference for all later
  verification and repair.
- Verify before you trust. A corrupt or unverified copy must never overwrite or
  replace a good record. Re-verify on a rolling schedule.

### Provenance & audit

- Every consequential action (create, move, archive/restore, retention, legal-hold and
  tag changes, deletion, and **destructive operator actions**) is recorded in an
  **append-only, tamper-evident audit trail**.
- The chain's integrity is itself verifiable; define and alert on the response to a
  broken chain — a broken audit chain is a compliance event.

### Retention & legal hold

- Enforce retention and legal hold at **every** delete path, not only the retention
  sweep. Legal hold overrides expiry, always.

### Destructive-action safety

- Never delete the **last copy** of data on a weak signal. Gate destructive actions on
  a grace window **plus an integrity-verified replacement**, or on **explicit operator
  confirmation** with an explicit target list.
- Detection is **report-first**: scans find and record, they do not mutate. Consider
  four-eyes (a second approval) for mass deletion.

### Durability & recovery

- Maintain **independent failure domains**: the primary store, the backup store, and
  the cold archive must not share a disk or host. Don't assume it — check it.
- Hold a **durability invariant**: every object is either replicated to the configured
  floor or archived to an independently durable store. Surface violations.
- Define **RPO/RTO** explicitly and make recovery procedures **tested**, not merely
  documented. A runbook drill that has never been executed is a hypothesis.

### Regulatory mindset

Design so an auditor can follow the trail end-to-end. Even without claiming a specific
certification, build to support the expectations of clinical-laboratory regimes (for
example ISO 15189 and CAP/CLIA-style requirements): traceability of who/what/when,
change control through small reviewable patches, validation evidence, reproducibility,
and clear separation of validated production behaviour from anything experimental.

### Observability is part of safety

Safe-automation jobs delete data after grace and run unattended; their health must be
a **guaranteed signal**. Ship SLIs/SLOs, alert rules, and dashboards. Route integrity
and consistency findings to **paging**, not to a file nobody reads.

### High-throughput operations (the scale the design targets)

Clinical metagenomics means multi-gigabyte files, large object counts, and sustained
ingest. Build for it:

- **Stream, don't buffer.** Hash and verify by streaming; avoid temp-disk copies of
  large objects.
- **Incremental, budgeted, staggered batch work.** Producers page the catalogue and
  enqueue per-object jobs within a per-run budget, on staggered schedules, so a large
  estate is swept without thundering herds.
- **Bounded concurrency, backpressure, idempotent retries.** Workers do bounded units
  and survive redelivery.
- **Measure cadence against the estate.** A rolling verify only protects you if it
  meets its integrity SLA at your scale; budgeted reconcile only covers the estate if
  cadence × budget ≥ size. Measure these; tune them; alert when they slip.
- **Avoid unbounded operations.** Prefer incremental backups and pruning over
  ever-growing full dumps as the catalogue grows. Scale stateless services
  horizontally; partition work.

---

## 5. Research / small-lab profile (the deliberate minor case)

Some deployments are smaller research or low-volume settings that don't need the full
operational apparatus. Support them as a **lighter operational profile**, not a
different system:

- Relax **operational burden**: single-node components may be acceptable, drill
  cadence and on-call expectations can be lower, retention policies simpler, and the
  durability posture may be explicitly non-authoritative.
- **Never relax the data-safety primitives.** Integrity hashing, the audit trail, and
  legal-hold/retention enforcement stay on. The small-lab path lowers operational
  overhead, not data correctness or accountability.
- Make the profile a **configuration choice with a clear, documented trade-off**, so
  an operator knows exactly what they have and have not signed up for. When in doubt,
  default to the clinical posture.

---

## 6. Documentation & runbooks (the common docs tree)

Operational documentation is a **first-class deliverable**, not an afterthought. Every
feature ships its docs into the **shared `docs/` tree**, indexed and cross-linked
alongside every other feature's — one place an operator learns to understand, run,
verify, and recover the system. The runbooks proved among the most valuable artifacts
of this work; treat them as part of "done".

The documentation answers four questions: **understand it, operate it, verify it,
recover it.**

### Where it lives and how it's organised

- A single `docs/` tree, indexed by `docs/README.md` as a **categorised map**
  (Overview / Operations / Reference), with every doc cross-linked.
- Written as a **system, not a development history** — no stage scaffolding (see §2),
  professional and present-tense.
- Each doc states its **tested-vs-validated boundary** honestly, so a reader knows
  what is proven versus what their environment must confirm.
- Docs are kept current with the code: documentation drift is a bug, not cosmetic.

### The standard set (modelled on what worked here)

- **Overview layer** — a short, high-level pair: an *overview* (what the feature is,
  the problem it solves, its guarantees and their limits) and an *architecture* doc
  (components, the data and storage model, the lifecycle).
- **Operations layer** — the runbooks an operator lives in:
  - an *administration* guide (deploy, configure, secrets, a consolidated
    configuration reference, prerequisites);
  - a *maintenance* guide (the scheduled-job catalogue with cadences, and how to run
    operations on demand);
  - per-subsystem *topic runbooks* (one per major capability);
  - a *disaster-recovery runbook* (see the anatomy below).
- **Reference layer** — a *validation* guide (how to build, test, and smoke-check, and
  the test boundary) and a *production-readiness* assessment (tracked gaps + roadmap).

### Anatomy of an operational runbook

The recovery runbook pattern that worked, reusable for any operational procedure:

1. **First response / triage.** What to check first, and — critically — *pause
   destructive automation* so recovery never races deletion.
2. **Scenario index.** A table of failure scenarios with **stable IDs**, each mapping
   cause → recovery steps → verification. Stable IDs let people refer to "scenario
   N" for years.
3. **Concrete, copy-pasteable commands** — real commands and endpoints, not prose
   descriptions of them.
4. **A verification step for every procedure** — how you *know* it worked
   (hash-verified bytes, object counts reconcile, health green). A procedure without a
   verification step is incomplete.
5. **Drills.** Rehearsals that prove the procedure against a real stack, run
   periodically. An unexecuted runbook is a hypothesis, not a capability (see §4).
6. **A recovery matrix** — what is recoverable, from where, and the RPO/RTO for each
   class of state.

### Operate and verify (the heart of it)

- **Operate.** Document routine operation and every on-demand action with the actual
  command or endpoint, the arguments, and the expected result — including the
  operator-gated destructive actions and exactly how they are confirmed.
- **Verify.** Give the reader a way to confirm health at every level: a quick
  smoke-check after deploy or after applying changes; the integrity/consistency checks
  that prove data is intact; the metrics and signals to watch; and the validation
  guide for building and running the tests. "How do I know it's working?" always has a
  documented answer.

### Done when

A feature's documentation is complete when a competent operator who has never seen the
code can, from `docs/` alone: understand what it does, stand it up and configure it,
run its routine and on-demand operations, confirm it is healthy, and recover it from
each failure in the scenario index — and when those recovery procedures have been
drilled at least once.

---

## 7. Putting it together — feature workflow

1. **Frame.** State the feature, who relies on it, and its clinical and throughput
   stakes. Identify the data-safety and scale risks up front.
2. **Plan.** Write `STAGE-NN-PLAN.md`: goal, locked decisions (with alternatives),
   sub-stages in risk-first order, exit criteria, out-of-scope.
3. **Build.** Sub-stage by sub-stage: pure logic first and unit-tested; behaviour and
   refactors kept separate; one verified patch each; the tested-vs-validated boundary
   stated.
4. **Document as a system.** Provide the full documentation set into the shared
   `docs/` tree — overview, architecture, administration, maintenance, topic runbooks,
   and a recovery runbook with drills (see §6) — and keep stage scaffolding out of
   shipped artifacts.
5. **Assess.** Close with an honest, critical production-readiness assessment as
   tracked, prioritized gaps.
6. **Summarise.** Write a concise PR summary referencing the docs and code.

### Definition of done (per sub-stage)

- One focused commit, message with ID + *what/why*, patch verified to apply.
- Pure logic unit-tested; tests hermetic and dependency-free; behaviour preserved
  where claimed (balance-checked for text-only edits).
- Docs updated; no development-stage markers in shipped comments or product docs.
- Data-safety properties (§4) upheld; throughput implications considered; honest
  validation boundary stated.

---

## 8. Anti-patterns to avoid

- Claiming validation that wasn't done, or stating an inferred fact as confirmed.
- Auto-deleting the last copy of data on a presence/existence check rather than a
  verified-integrity check.
- Enforcing retention or legal hold in only one delete path.
- Burying integrity/consistency findings in a file instead of alerting on them.
- Letting destructive operations run without a gate, a grace window, or an audit
  record.
- Buffering whole multi-GB objects into memory or temp disk where streaming would do.
- Unbounded sweeps or full dumps that don't scale with the estate.
- Leaking development-stage scaffolding (S/H/Stage/decision numbers) into shipped
  docstrings, comments, or operational docs.
- Large, multi-concern commits that can't be reviewed or reverted cleanly.
- Breaking the serialized shape of a domain record without a compatibility path.

---

## Maintaining this document

This playbook is intentionally generic so it ages well. Update it when an architectural
role or technology changes, when a new cross-cutting safety requirement is adopted, or
when the method itself improves. Keep technology specifics behind their abstractions,
keep the clinical posture as the default, and keep the examples grounded in the actual
codebase rather than abstract ideals.
