# ciqa — overview

CIQA is the quality-assurance layer for the **generative diagnostic** step of the Cerebro
clinical metagenomic diagnostics platform. Where the core pipeline turns sequencing reads into
a taxonomic profile, CIQA turns that profile into a **diagnosis** — by way of META-GPT, a
language model acting as a diagnostic practitioner — and then holds that diagnosis to account:
it measures how often the model is right, against known truth, and refuses to let a model or
configuration change ship until it has been **regression-tested** against a frozen baseline.

## The problem it solves

A generative model in a clinical path is only as trustworthy as the evidence that it has not
silently regressed. Swapping a model, changing a quantization, adjusting a filter threshold, or
editing a prompt can all move sensitivity and specificity — sometimes invisibly, sometimes only
on a particular sample type. Doing this safely requires three things that are tedious and
error-prone by hand: a **reproducible filtered input** to the model, a **truth set** to score
against, and a **disciplined comparison** to a known-good baseline with enough statistical care
to tell a real change from noise. CIQA makes those three a property of the system: the prefetch
is built the same way every time, the truth lives in versioned datasets, and the comparison is a
single command that produces an auditable report.

It also solves a deployment problem. Running the model historically required a live Cerebro
stack to produce the filtered input. CIQA can now produce that input **self-contained from a
run's output models**, so META-GPT can run **inside the pipeline on a local GPU** with no stack
attached — the operating profile a smaller laboratory or an offline validation host needs.

## What it provides

- **Tiered prefetch filtering.** Per sample, the run's taxa are filtered into **primary /
  secondary / target** tiers and separated from prevalence-driven contaminants, yielding a
  `PrefetchData` — the model's input. Produced from the live stack **or** self-contained from
  the run's output models, selected by configuration.
- **Prevalence filtering in the pipeline.** Cross-sample contaminant taxa are identified over
  the samples of a run (or fetched from the stack), and removed, under a configuration toggle.
- **Generative diagnosis.** META-GPT renders a `DiagnosticResult` per sample — an
  infectious / non-infectious / tumour call and a pathogen candidate — locally on a GPU or via
  a hosted endpoint, with the model **downloaded by configuration**.
- **One diagnostic evaluator.** A single, typed routine classifies each diagnosis as
  **TP / TN / FP / FN / Excluded** against reference truth and computes **sensitivity /
  specificity / PPV / NPV**. The same evaluator backs both the interactive **training** review
  and the automated **QC** scoring, so the numbers mean the same thing everywhere.
- **Training and QC datasets.** Labelled prefetch items for interactive model/prompt
  development (training) and reference-truthed, versioned sets for automated quality control
  (QC) share one storage shape (a record plus the `PrefetchData` in GridFS).
- **Regression testing against a baseline.** A run's statistics are compared to a stored,
  **versioned baseline** — an absolute sensitivity/specificity floor plus a paired **McNemar**
  test (with multiple-comparison adjustment when several configurations are compared) — to
  decide whether a change has regressed, emitting a report and a gate signal.
- **Pipeline integration.** The whole prefetch → diagnose → regression sequence runs as a
  config-gated stage of the Nextflow pipeline on a local GPU, inside the NVIDIA CUDA container.

## Design principles

- **The model's contract is stable.** META-GPT consumes a `PrefetchData` and produces a
  `DiagnosticResult`; CIQA never changes those shapes or the model's operating logic. Every
  feature is built around that seam.
- **One source of filtering and one source of evaluation.** The prefetch reuses the same
  filtering functions the stack uses, and scoring uses one evaluator shared with training, so a
  self-contained run and a stack-backed run agree by construction.
- **Report first, act second.** Regression detection records and alerts; it never promotes a
  baseline or deletes data on its own. Changing the truth a future run is judged against is a
  deliberate, gated, audited act.
- **Reproducibility over convenience.** Model/configuration provenance is carried in an
  explicit, typed **run manifest** rather than recovered from directory names, and every
  persisted artifact carries a content hash.
- **Honest about validation.** The system distinguishes what its unit tests prove (the pure
  filtering, evaluation, and comparison logic) from what only a running stack or a real GPU can
  confirm; see [`validation.md`](validation.md).

## The stable seam

```
                 (tiered + prevalence filtering)
 run output models ─────────────────────────────►  PrefetchData
                                                        │
                                                        ▼
                                          META-GPT  (run_local)   ← stable: not changed by CIQA
                                                        │
                                                        ▼
                                                  DiagnosticResult
                                                        │
                          (one evaluator vs reference truth)
                                                        ▼
                              sensitivity / specificity / PPV / NPV
                                                        │
                                     (vs versioned baseline + McNemar)
                                                        ▼
                                              regression report + gate
```

## Guarantees, and their limits

CIQA is honest about what its measurements do and do not mean:

- **Sensitivity / specificity are only as good as the truth set.** They are computed against
  the reference results and candidates recorded in the dataset; a gap or error in the truth set
  is a gap or error in the metric. Datasets are versioned so a baseline is always tied to the
  exact truth it was measured against.
- **A non-significant McNemar result is not proof of no change.** Small panels have low
  statistical power; the report always states the panel size and discordant counts, and the
  absolute threshold floor is the primary gate, with significance as a secondary signal.
- **`local` and `stack` prefetch agree by construction, but parity is verified, not assumed.**
  The one documented difference (an optional stack-only contamination-history enrichment) is
  called out where it applies; a differential check confirms equality on a reference sample.
- **A baseline is a human decision.** CIQA computes and reports; an operator promotes. The
  system will not silently adopt a new baseline, even a better-looking one.

## Audience and entry points

- **Assay/method developers and QA** regression-testing a model or configuration start with the
  [regression-testing runbook](regression-testing.md) and [datasets](datasets.md).
- **Operators** running the model in the pipeline start with the
  [Nextflow GPU runbook](nextflow-gpu.md) and the [administration guide](administration.md).
- **Engineers** understanding the system start here, then the
  [architecture](architecture.md) and the per-topic docs.
- **Anyone validating a deployment** uses the [validation guide](validation.md).

See the [documentation index](index.md) for the full map.
