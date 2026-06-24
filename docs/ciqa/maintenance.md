# ciqa — maintenance

The day-to-day and periodic operations that keep CIQA's quality-control function meaningful:
running the workflow on demand, keeping datasets and baselines current, and the cadence of
regression. For deploying and configuring, see [administration](administration.md); for the
step-by-step runbooks, see [regression-testing](regression-testing.md) and
[nextflow-gpu](nextflow-gpu.md).

## The command surface

CIQA is a single binary (`cerebro-ciqa`) whose subcommands map onto the data-flow spine. The
ones that matter operationally:

| Command | Does |
|---|---|
| `prefetch` | Build the tiered + prevalence-filtered `PrefetchData` per sample (source `stack` or `local`) |
| `diagnose` | Run META-GPT (the generative practitioner) over the prefetch on a local GPU; emit `DiagnosticResult` + run manifest |
| `review` | Score diagnostic outputs against the reference truth → sensitivity/specificity |
| `mcnemar` / `mcnemar-adjust` | Compare two reviews, or a series against a baseline, with McNemar (and multiple-comparison adjustment) |
| `upload` | Register/stage a CI/QA dataset to Cerebro Fs and index it for the pipeline |
| `plate` / `plate-table` | Create/inspect the reference plate (the truth layout) |
| `*-plot`, `summarize`, `compute-summary` | Plots and summaries for review and reporting |

The plotting, summary, and TUI commands are for human review and reporting; the operational
core is `prefetch → diagnose → review`/`mcnemar-adjust`, and in the pipeline the equivalent
Nextflow processes.

## On-demand operations

### Re-run a diagnosis over an existing run

To try a different model or configuration over a run you have already profiled — without
re-profiling — use the standalone path: build the prefetch from the run's output models
(`prefetch --source local`), then `diagnose` with the new model/config. This is the cheap inner
loop of model/prompt development and is also exposed as a standalone Nextflow entry
([`nextflow-gpu.md`](nextflow-gpu.md)).

### Score and compare

`review` turns a set of diagnostic outputs into sensitivity/specificity against the truth set;
`mcnemar-adjust` compares a series of configurations against a baseline with multiple-comparison
adjustment. For the production gate, prefer the regression workflow
([`regression-testing.md`](regression-testing.md)), which combines the absolute threshold floor
with the paired test and produces an auditable report.

### Idempotency and resume

The diagnose step skips samples whose result already exists unless forced, so an interrupted run
resumes cheaply; in Nextflow the processes are safe under `-resume`. Re-running a regression is
always safe — it reads outputs and a baseline and writes a report; it never mutates the baseline.

## Periodic upkeep

### Refresh QC datasets

A QC dataset is the truth a baseline is measured against. Keep it current as the assay's truth
evolves: add new reference-truthed samples, and **version** the dataset when its truth changes
(a baseline is tied to a specific dataset version). Adding samples or correcting truth without a
version bump would silently change what a baseline means; version instead. See
[`datasets.md`](datasets.md).

### Promote and rotate baselines

A baseline is the known-good statistics a future change is judged against. The cadence:

1. After a model/config you trust passes review on the current dataset version, **promote** its
   statistics as the baseline for that dataset version. Promotion is an explicit, audited,
   operator-gated action (four-eyes recommended for a clinical baseline).
2. Keep prior baselines — they are immutable and versioned, so historical comparisons remain
   reproducible.
3. When the dataset version changes, establish a new baseline against the new version before
   gating changes on it.

CIQA never promotes a baseline for you, even a better-looking one; promotion is always a human
decision ([`regression-testing.md`](regression-testing.md)).

### Re-pin the model

The generative model is pinned by configuration. When you intend to move to a new model or
quantization, that *is* a change to regression-test: run it against the current baseline,
inspect the report, and only then decide whether to adopt it and re-baseline.

## Health checks

- A regression run completes and writes a report; `regressed = false` and
  `passed_threshold = true` against a matching baseline indicates the gate is working.
- The diagnose step's NVML benchmark sidecars show non-zero peak VRAM — confirmation the model
  actually ran on the GPU.
- Consequential actions (baseline promotion, dataset delete) appear in the audit chain.

## Tested here vs your environment

The command behaviours described here are exercised by unit tests for the pure logic
(filtering, evaluation, comparison) and are otherwise confirmed on your stack and GPU host. The
cadence advice (versioning, promotion, rotation) is operational policy, not enforced by code
beyond the gating and audit described; adapt it to your laboratory's change-control regime. See
[`validation.md`](validation.md) for the boundary.
