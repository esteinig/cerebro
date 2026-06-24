# ciqa — META-GPT diagnosis

How CIQA runs **META-GPT**, the generative diagnostic practitioner, over the prefetch to produce
a per-sample diagnosis — on a local GPU, with the model downloaded by configuration — and the
**stable contract** that keeps the model's operating logic fixed while everything around it
evolves. For the input, see [prefetch & filtering](prefetch-filtering.md); for scoring the
output, see [regression-testing](regression-testing.md).

## What META-GPT does

META-GPT is a language model acting as a diagnostic practitioner. Given a sample's prefetch — its
tiered, filtered taxonomic evidence and context — it renders a **`DiagnosticResult`**: a
**diagnosis** and a **pathogen candidate**. The diagnosis is one of:

| `Diagnosis` | Mapped reference result |
|---|---|
| `Infectious`, `InfectiousReview` | `Positive` |
| `NonInfectious`, `NonInfectiousReview` | `Negative` |
| `Tumor` | (neither — excluded from infectious sensitivity) |

The `*Review` variants mark a call the model flags for human review. The pathogen candidate is the
predicted species, normalised for comparison against the truth set.

CIQA does not change how META-GPT reasons or what it emits; it supplies the input and consumes the
output.

## The stable seam

Everything CIQA builds is arranged around one fixed contract:

```
PrefetchData (in)  ──►  META-GPT (run_local)  ──►  DiagnosticResult (out)
```

- **Input:** a `PrefetchData` (the prefetch for one sample) plus the sample context, optional
  clinical notes, the assay context, the agent primer, and the post-filter configuration.
- **Output:** a `DiagnosticResult`, written as `{sample}.model.json`.

This seam does not change — not when the prefetch source changes, not when scoring or regression
changes, not when the run moves into Nextflow. New capabilities produce a `PrefetchData` or
consume a `DiagnosticResult`; they never reach inside the model step.

## Running on a local GPU

The local diagnosis (the `diagnose` command, built with the `local` feature) runs the model on
the host's GPUs:

1. It loads the local model from the **model directory** (the weights), once per GPU.
2. It splits the samples across the available GPUs and runs them in parallel, monitoring VRAM via
   NVML and writing benchmark sidecars (peak VRAM, runtime) — useful evidence that the model
   actually ran on the GPU.
3. For each sample it reads `{sample}.prefetch.json`, builds the diagnostic agent, applies the
   post-filter, and calls `run_local`, writing `{sample}.model.json` and the agent's state log.
4. Samples whose result already exists are skipped unless forced, so an interrupted run resumes
   cheaply.

```bash
# self-contained: prefetch from the run's models, then diagnose on the GPU
cerebro-ciqa prefetch --source local --models <run-dir> --outdir prefetch/
cerebro-ciqa diagnose \
  --plate plate.json \
  --prefetch prefetch/ \
  --outdir diagnose/ \
  --model qwen3-14b \
  --model-dir /models/qwen3-14b \
  --num-gpu 1
# → diagnose/{sample}.model.json + state logs + benchmark sidecars
```

Inside the pipeline this is the `metaGpt` GPU process, in the NVIDIA CUDA container; see
[`nextflow-gpu.md`](nextflow-gpu.md).

## Model download by configuration

The model directory either **already contains** the weights (mounted) or is **populated by a
configured download** before diagnosis. The model id, the directory, and the download toggle are
configuration ([administration](administration.md)). Treat downloaded weights as a
provenance-tracked, hash-verified dependency — they are part of the clinical path.

## The run manifest

Alongside the `{sample}.model.json` files, the diagnose step emits a typed **run manifest**
recording exactly what produced this set of diagnoses: the model id, quantization, parameters,
the clinical flag, the replicate, a hash of the configuration used, and each sample's prefetch and
result paths with content hashes. The manifest — not a parsed directory name — is what the
regression step consumes to attribute results to a model/configuration. It is the provenance
record for the generative step.

## Hosted vs local

The local-GPU path above is the self-contained default for in-pipeline and offline use. A hosted
model endpoint can drive the same seam (the agent's `run` path) where a deployment prefers it; the
input and output contracts are identical, so scoring and regression are unaffected by which the
diagnosis came from. The manifest records which model produced the result either way.

## Tested here vs your environment

The agent construction, post-filter, manifest assembly, and the mapping from `DiagnosticResult` to
a scoreable prediction are **unit-tested pure logic**. The model load, GPU execution, VRAM
behaviour, and download are **environment-validated only** — they require a real GPU host and are
confirmed by the GPU run in [`nextflow-gpu.md`](nextflow-gpu.md) and the checks in
[`validation.md`](validation.md). The on-disk `DiagnosticResult` shape is held stable and is
round-trip tested against a real model output fixture.

## Running the diagnosis (CLI)

The diagnose step reads the prefetch, runs the agent on a local GPU, and writes one
`{sample}.model.json` per sample. `--emit-manifest` additionally writes the
`MetaGptRunManifest` the regression step consumes. The seam (`run_local`) is unchanged.

```bash
cerebro-ciqa diagnose-local \
  --prefetch prefetch --plate plate.json --outdir diagnose \
  --model qwen3-8b-q8-0 --model-dir /weights --num-gpu 1 --post-filter \
  --emit-manifest diagnose/run.manifest.json \
  --run-id RUN-2026-06-001 --model-id qwen3-8b-q8-0 --replicate 1
```

Peak VRAM per GPU is recorded to `diagnose/state_logs/gpu*.bench.json`. For the full flag set,
the regression and report commands, and the equivalent endpoints, see
[Commands & API](commands.md). To run the whole prefetch → diagnose → regression path on a GPU
host under Nextflow, see [Nextflow GPU diagnosis & regression](nextflow-gpu.md).
