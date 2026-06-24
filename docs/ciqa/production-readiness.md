# ciqa — production readiness

An honest, CIQA-scoped readiness summary: the gaps that stand between this documentation's
operational present tense and a clinically deployed, fully validated CIQA + META-GPT subsystem.
Tiers: **P0** blocks clinical production; **P1** required before broad rollout; **P2** important
hardening; **P3** nice-to-have. This page summarises and points at the full, cross-subsystem
register in the integration plan; for how each item is checked, see [validation](validation.md).

> **The dominant caveat.** Parts of the surface this set documents **exist today** (the CIQA CLI
> commands, the training subsystem, the GPU container); parts are **delivered by the integration**
> (the `local` prefetch source, QC datasets/baselines/regression, the `ciqa/*` API, the Nextflow
> GPU stage). The largest class of remaining work is **environment validation** — the GPU, stack,
> and download paths are drilled, not unit-tested.

## P0 — Blocks clinical production

- **Build & decision-logic evidence.** A clean `cargo build --workspace --locked` (CPU) and
  `--features local` (GPU), plus the pure-logic test suites, are the entry gate
  ([validation](validation.md)).
- **`local` vs `stack` prefetch parity.** The clinical input to the model must be identical
  whichever source built it; the differential check must pass on a representative panel, with the
  one documented difference accepted or closed ([prefetch & filtering](prefetch-filtering.md)).
- **Prevalence-counter equivalence.** The in-pipeline prevalence set must match the stack's on a
  known collection before it is trusted.
- **Evaluator number-preservation.** The shared evaluator must reproduce the existing training and
  CIQA numbers exactly (cross-check tests) — a regulatory metric must not move under a refactor.
- **GPU end-to-end.** The diagnose/regression stage must be drilled on a real GPU host
  ([nextflow-gpu](nextflow-gpu.md)); until then the headline capability is designed, not delivered.

## P1 — Required before broad rollout

- **Audit + provenance** for the consequential actions — baseline promotion, dataset delete,
  regression submit — written to the append-only chain, with promotion operator-gated (four-eyes
  recommended) ([datasets](datasets.md), [regression-testing](regression-testing.md)).
- **Integrity hashing** of derived artifacts (prefetch, diagnoses, manifest, baseline, report),
  verified before a baseline is trusted.
- **Model-output coupling.** A round-trip test of the `DiagnosticResult` shape against a pinned
  real model output, re-run whenever the model dependency is re-pinned.
- **Model download mechanism** defined per deployment, with hash-verified weights treated as a
  clinical-path dependency ([administration](administration.md)).
- **Integration tests** for the `ciqa/*` endpoints and the standalone `metagpt` workflow against an
  ephemeral stack — the I/O paths the unit tests do not cover.

## P2 — Important hardening

- **Library-grouping shape.** Confirm whether the pipeline emits one model per library or a
  combined model, and encode the DNA/RNA grouping explicitly
  ([prefetch & filtering](prefetch-filtering.md)).
- **No-regex guarantee.** Keep the legacy directory-name parser quarantined and deprecated; a CI
  grep-gate asserts the production and regression paths never call it.
- **Throughput** of `local` prefetch + prevalence over large runs — stream taxids rather than
  materialising all taxa; bound concurrency.
- **Error taxonomy** for the new failure modes (missing model dir, GPU init failure,
  baseline-not-found, manifest parse error) — typed errors, no panics on production paths.
- **Statistical-power reporting** in every regression report (panel size, discordant counts, chosen
  adjustment), with the absolute floor as the primary gate
  ([regression-testing](regression-testing.md)).
- **Idempotent resume** of the manifest/regression steps under Nextflow `-resume`.

## P3 — Nice-to-have

- A **web UI** for QC datasets / baselines / regression, mirroring the training UI.
- **Multi-node GPU** orchestration beyond a single host.
- The optional **contamination-history enrichment** in `local` mode, to close the one documented
  `local`/`stack` difference if a use case appears.

## Where the full register lives

This page is the CIQA-scoped view. The complete, cross-subsystem readiness register — with the
same tiers, the per-item risk and recommended action, and the roadmap that sequences clearing them
— is maintained alongside the integration plan (`05-PRODUCTION-READINESS.md`). A deployment is
ready for clinical use when the P0 items are cleared with evidence and the P1 items are scheduled;
the P2 hardening precedes broad rollout.
