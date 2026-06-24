# ciqa

Operator and developer documentation for **CIQA** — the continuous integration of quality
assurance for Cerebro, the clinical metagenomic diagnostics platform — and the **META-GPT**
generative diagnostic practitioner it drives. This is the entry point; it orients you and
points at the rest of the documentation set.

CIQA takes the taxonomic evidence a Cerebro run produces, filters it into a tiered
**prefetch**, asks META-GPT (a local or hosted language model) to render a **diagnosis** per
sample, and scores those diagnoses against a known truth — as **sensitivity / specificity /
PPV / NPV** — so that a model or configuration change can be **regression-tested** against a
frozen **baseline** before it reaches clinical use. It is the laboratory's change-control
instrument for the generative step.

## Where to start

- **New to CIQA / META-GPT** → the [overview](overview.md), then the
  [architecture](architecture.md).
- **Deploying / configuring** → the [administration guide](administration.md) (the full
  configuration reference, the GPU image build, and model download live here).
- **Running it day to day** → the [maintenance guide](maintenance.md).
- **Understanding the filtered input** → [prefetch & filtering](prefetch-filtering.md).
- **Running the model** → [META-GPT diagnosis](meta-gpt-diagnosis.md).
- **Working with the training & QC sets** → [datasets](datasets.md).
- **Regression-testing a model/config** → the
  [regression-testing runbook](regression-testing.md) (triage → scenario → commands → verify).
- **Running it on a local GPU in the pipeline** → the
  [Nextflow GPU runbook](nextflow-gpu.md).
- **Validating a deployment** → the [validation guide](validation.md).
- **Planning the road to production** → the
  [production-readiness assessment](production-readiness.md).

## The data-flow spine

The spine that ties the docs together — evidence for one run moves through:

1. **Output models.** A Cerebro run emits per-sample/library **Cerebro model** JSON carrying
   the taxonomic profile (`taxa`). These are the input to everything CIQA does —
   [`architecture.md`](architecture.md).
2. **Prefetch (tiered + prevalence filtering).** The taxa are filtered into three tiers
   (**primary / secondary / target**) and split from prevalence-contaminants, producing a
   `PrefetchData` per sample. This is produced **either from the live stack** or
   **self-contained from the run's output models** (a user setting) —
   [`prefetch-filtering.md`](prefetch-filtering.md).
3. **Diagnosis (META-GPT).** Each `PrefetchData` is handed to the generative practitioner,
   which renders a `DiagnosticResult` (an infectious / non-infectious / tumour call plus a
   pathogen candidate). On a local GPU this runs inside the **NVIDIA CUDA container** —
   [`meta-gpt-diagnosis.md`](meta-gpt-diagnosis.md).
4. **Evaluation (sens/spec).** Each `DiagnosticResult` is compared to the sample's reference
   truth to classify it **TP / TN / FP / FN / Excluded**, and the set is reduced to
   **sensitivity / specificity / PPV / NPV** — [`regression-testing.md`](regression-testing.md).
5. **Regression vs baseline.** The run's statistics are compared to a stored, versioned
   **baseline** (absolute thresholds + a paired **McNemar** test) to decide whether the
   change **regressed** — report-first, never auto-mutating —
   [`regression-testing.md`](regression-testing.md).
6. **Datasets.** The labelled prefetch items that feed interactive **training** and the
   reference-truthed **QC** sets share one shape (a record + `PrefetchData` in GridFS) and one
   evaluator — [`datasets.md`](datasets.md).

The one fixed point across all of this is the **stable seam**:

```
PrefetchData (in)  ──►  META-GPT (run_local)  ──►  DiagnosticResult (out)
```

Everything else is built around that seam without changing it.

## Documentation map

**Overview**

| Doc | Covers |
|---|---|
| [`overview.md`](overview.md) | What CIQA + META-GPT are, the problem they solve, capabilities, the stable seam, guarantees |
| [`architecture.md`](architecture.md) | Components, the data flow, the crate/dependency layering, where pure logic lives, the GPU container |

**Operations**

| Doc | Covers |
|---|---|
| [`administration.md`](administration.md) | Deploying, the full configuration reference, the GPU image build, model download, secrets |
| [`maintenance.md`](maintenance.md) | On-demand operations, baseline cadence, dataset upkeep |
| [`prefetch-filtering.md`](prefetch-filtering.md) | Tiered filtering, prevalence filtering, the `stack` vs `local` source setting |
| [`meta-gpt-diagnosis.md`](meta-gpt-diagnosis.md) | The generative practitioner, local-GPU execution, model download, the output contract |
| [`datasets.md`](datasets.md) | The training sets and the QC sets — shared shape, different lifecycle |
| [`regression-testing.md`](regression-testing.md) | The regression runbook: compute sens/spec, compare to baseline, McNemar + thresholds |
| [`nextflow-gpu.md`](nextflow-gpu.md) | Running prefetch → diagnose → regression on a local GPU in Nextflow |

**Reference**

| Doc | Covers |
|---|---|
| [`validation.md`](validation.md) | Build, test, smoke checks, the tested-vs-environment boundary |
| [`production-readiness.md`](production-readiness.md) | Honest readiness assessment, tracked gaps, and the hardening roadmap |

## Conventions

- Each operational doc ends with a **"tested here vs your environment"** note — the honest
  boundary between what unit tests cover and what only a running stack or a real GPU can
  validate. [`validation.md`](validation.md) consolidates this.
- **Authoritative and destructive operations are operator-gated and audited.** Promoting a
  baseline (the truth a future change is judged against) and deleting a QC dataset both
  require explicit operator action and are written to the append-only audit trail.
- **Report first, act second.** Regression *detects and records*; it never promotes a baseline
  or deletes data on a weak signal. Promotion is a separate, deliberate step.
- The **META-GPT operating model and data ingestion are stable** — `PrefetchData` in,
  `DiagnosticResult` out. Everything documented here produces or consumes those two artifacts.
