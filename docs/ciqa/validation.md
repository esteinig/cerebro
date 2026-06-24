# ciqa — validation

How to validate a CIQA deployment, and — the important part — the honest boundary between what the
unit tests prove and what only a running stack or a real GPU can confirm. This page consolidates
the "tested here vs your environment" notes from the rest of the set. For the gaps and the
hardening roadmap, see [production-readiness](production-readiness.md).

## Build

CIQA is a workspace member built with the rest of Cerebro:

```bash
# default (CPU) build — no model, no GPU
cargo build --workspace --locked

# the local-GPU diagnosis path (CUDA; on a GPU host)
cargo build -p cerebro-ciqa --release --locked --features local
```

The default build is feature-free and CPU-only; the `local` feature pulls the local META-GPT
generator and is GPU-oriented. A clean `--locked` build from a fresh checkout is the first gate —
it confirms the workspace resolves with path dependencies and a single lockfile.

## Test — the pure core

The logic that matters for correctness is pure and hermetic; run it without any stack or GPU:

```bash
cargo test --workspace
# focused:
cargo test -p cerebro-pipeline taxa::          # filtering + prevalence
cargo test -p cerebro-model    prefetch_local  # self-contained prefetch assembly
cargo test -p cerebro-model    diagnostics     # the shared evaluator (sens/spec)
```

These cover:

- **Tiered filtering, aggregation, prevalence counting** — the prefetch's decisions.
- **Self-contained prefetch assembly** — that `local` builds the same `PrefetchData` shape.
- **The diagnostic evaluator** — decisions → sensitivity / specificity / PPV / NPV, cross-checked
  so the training and QC numbers are identical for the same inputs.
- **The regression comparison and McNemar** — the maths of the gate.
- **Manifest and `DiagnosticResult` (de)serialization** — round-tripped, the latter against a real
  model-output fixture.

If these pass, the *decision logic* is sound. They do **not** prove the system runs in your
environment — that is the next two sections.

## Smoke checks — your environment

After deploying, confirm the parts the unit tests cannot:

1. **Self-contained prefetch (no stack).** With the stack down:

   ```bash
   cerebro-ciqa prefetch --source local --models <run-dir> --outdir prefetch/
   ```

   → one `{sample}.prefetch.json` per sample.

2. **`local` vs `stack` parity.** Against a dev stack, build a reference sample's prefetch both
   ways and confirm they are equal except for the documented contamination-history enrichment (the
   differential check). This is the proof that the self-contained output matches production.

3. **Prevalence-counter equivalence.** Confirm the prevalence set from the in-pipeline function
   matches the stack's contamination endpoint on a known collection.

4. **GPU diagnosis.** On a GPU host, run the diagnose step (or the Nextflow `metagpt` entry); the
   benchmark sidecars must show **non-zero peak VRAM** — confirmation the model ran on the GPU.

5. **Regression round-trip.** Register a QC dataset, promote a baseline, submit a run manifest, and
   retrieve the report; confirm `regressed = false` against a matching run and that the baseline is
   unchanged. Confirm an audit entry for the promotion.

6. **Off-by-default.** Diff a normal pipeline run's outputs with the META-GPT stage disabled — they
   must be identical to before the feature.

## The tested-vs-environment boundary

| Concern | Proven by unit tests | Needs your environment |
|---|---|---|
| Tiered filtering / prevalence / prune | ✔ pure | parity with the stack on real samples |
| Self-contained prefetch assembly | ✔ pure | reading real output models; library grouping shape |
| Diagnostic evaluator (sens/spec) | ✔ pure, cross-checked | — |
| McNemar + adjustment, regression maths | ✔ pure | statistical power of your panel |
| `DiagnosticResult` / manifest shapes | ✔ round-trip | round-trip vs the real model output |
| Model load / GPU execution / VRAM | — | ✔ GPU host only |
| Model download | — | ✔ your download source + hash |
| `ciqa/*` persistence + audit | — | ✔ running stack |
| Nextflow GPU wiring (`--nv`/`--gpus`) | — | ✔ GPU host only |

## Honesty note

The pure-logic tests are the regulatory core and are designed to be runnable and trustworthy
without infrastructure. The infrastructure-dependent paths — the GPU, the stack, model download —
are **drilled, not unit-tested**, and the runbooks ([`regression-testing.md`](regression-testing.md),
[`nextflow-gpu.md`](nextflow-gpu.md)) are how a deployment confirms them. The remaining gaps,
prioritised, are in [production-readiness](production-readiness.md) and the integration plan's full
readiness register.

## Report rendering (Typst/WASM)

The regression report and the gated META-GPT appendix are validated by an integration test
(`cerebro-report` `tests/render.rs`, CI job `render`): it renders the templates for PASS and
REGRESSED fixtures to PDF and asserts the verdict/stats content, and checks that "Appendix C"
appears only when a `meta_gpt` block is present (the back-compat guarantee, rendered). The native
job is the gate; a non-blocking `wasm-pack` job smoke-builds the browser `cdylib`. The render
needs the Typst toolchain + fonts, so it runs in CI rather than locally. See the
[Typst/WASM render plan](../../11-STAGE-TYPST-WASM-RENDER-PLAN.md) for the full scope.
