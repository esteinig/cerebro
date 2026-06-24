# Nextflow GPU diagnosis & regression (CIQA)

This page is the operator runbook for the local-GPU META-GPT diagnosis and the
regression-against-baseline capability wired into the pipeline (Stage 4 of the CIQA
integration). It assumes Stages 1–3 are in place (the binaries build in-workspace, the
self-contained prefetch exists, and the manifest/regression types are available).

> **Off by default.** Everything here is gated by
> `params.pathogenDetection.metaGpt.enabled` (default `false`). A normal pathogen run is
> byte-for-byte unchanged unless this is set to `true`.

## Data path

```
CreateCerebroModel (run Cerebro models, *.json)
        │
        ▼
MetaGptPrefetch        cerebro-ciqa prefetch --prefetch-source local --run-models … --plate-json …
        │  *.prefetch.json   (self-contained tiered + prevalence prefetch — no live stack)
        ▼
MetaGptDiagnose        cerebro-ciqa diagnose-local …  (GPU; the stable seam: PrefetchData → run_local → DiagnosticResult)
        │  *.model.json  +  run.manifest.json (MetaGptRunManifest)  +  state_logs/ (NVML VRAM)
        ▼
MetaGptRegression      cerebro-ciqa regression --manifest … --baseline … --dataset …  (optional, report-first)
        │  regression.report.json   (threshold floor + paired McNemar; exit code reflects the gate)
        ▼
(optional) UploadMetaGptResults
```

The diagnosis seam itself (`agent.run_local`) is unchanged; the only additions are the
self-contained prefetch source (Stage 2) and the run manifest the diagnose step emits (Stage 3).

## The image (S4-D1)

The GPU image is `templates/stack/docker/Dockerfile.gpt.dgx` (NVIDIA CUDA 12.2.2, builds
`--features local`, ships `cerebro-ciqa` + `cerebro-client`, creates the NVML symlink at
entrypoint). Build it from the repository base:

```bash
docker build -f templates/stack/docker/Dockerfile.gpt.dgx -t cerebro-metagpt:<tag> .
# apptainer host:
apptainer build cerebro-metagpt.sif docker-daemon://cerebro-metagpt:<tag>
```

It is referenced like any other tool via `params.resources.apptainer.metaGpt`. Model **weights
are not baked in** — `modelDir` is the contract: either mount it with the weights, or set
`params.pathogenDetection.metaGpt.download = true` with a `fetchCommand` that populates it (D11).

## GPU request (S4-D2)

The GPU is requested on the `metaGpt` process label only (`process { withLabel: metaGpt }`):
`containerOptions = params.pathogenDetection.metaGpt.containerOptions` (`--nv` for apptainer, or
`--gpus all` for docker) and `accelerator = …gpus`. The `dgx` profile sets `--nv` and keeps the
existing `--bind /raid:/raid`.

## Configuration

All under `params.pathogenDetection.metaGpt` (a peer of `taxonomicProfile`/`metagenomeAssembly`):

| Param | Default | Meaning |
|---|---|---|
| `enabled` | `false` | master gate (S4-D4) |
| `prefetchSource` | `"local"` | `local` (run models, no stack) or `stack` (Stage 2) |
| `model` | `"qwen3-8b-q8-0"` | generator model id |
| `modelDir` | `null` | weights dir (mounted or populated) |
| `download` / `fetchCommand` | `false` / `null` | optional weight-download hook |
| `gpus` / `containerOptions` | `1` / `"--nv"` | GPU request on the metaGpt label |
| `plate` | `null` | reference plate JSON (samples to diagnose) |
| `postFilter` / `clinicalNotes` | `true` / `false` | diagnosis options |
| `runId` / `replicate` | `null` | explicit manifest attribution (not parsed) |
| `modelsDir` | `null` | standalone entry: dir of an existing run's Cerebro models |
| `regression.enabled` | `false` | run the in-pipeline regression |
| `regression.baseline` / `regression.dataset` | `null` | baseline + QC truth for regression |
| `upload` | `false` | upload results + manifest to the stack |

## Standalone entry (S4-D5)

Re-run diagnosis/regression over a finished run without re-profiling:

```bash
nextflow run main.nf -entry metagpt -profile dgx \
  --pathogenDetection.metaGpt.enabled true \
  --pathogenDetection.metaGpt.modelsDir <run>/results/models \
  --pathogenDetection.metaGpt.plate plate.json \
  --pathogenDetection.metaGpt.modelDir <weights>
```

## Validation drill (must be run on a GPU host)

Nextflow/container/GPU behaviour is only meaningfully validated on a real NVIDIA-GPU host;
none of it is exercised by the unit tests. Run this drill once and record the outcome:

1. **Off-by-default proof.** A `pathogen` run with `metaGpt.enabled=false` produces the same
   published outputs as before (diff them).
2. **GPU run.** `-entry metagpt -profile dgx … --pathogenDetection.metaGpt.enabled true`
   produces `*.prefetch.json` (self-contained), `*.model.json`, `run.manifest.json`, and (if
   enabled) `regression.report.json`; the NVML sidecars (`state_logs/gpu*.bench.json`) show
   non-zero peak VRAM.
3. **Self-contained proof.** The same run with the stack **down** still completes through
   diagnosis (because `prefetchSource=local`).
4. **Regression gate.** A run against a matching baseline reports `regressed=false`; a degraded
   model reports `regressed=true` with an alert; the baseline is unchanged.

### Failure scenarios

| ID | Symptom | Recovery |
|---|---|---|
| GPU-OOM | CUDA out-of-memory during diagnosis | reduce `gpus`/batch, use a smaller quantization; `errorStrategy = "ignore"` keeps the run alive |
| MODEL-DL | weights missing in `modelDir` | mount the weights, or set `download=true` with a working `fetchCommand`; retry |
| PREFETCH-MISS | a sample has no prefetch | diagnose skips it (already handled); the manifest simply omits it |
| STACK-DOWN | `stack` prefetch source unreachable | switch `prefetchSource` to `local` |

An unexecuted runbook is a hypothesis — record the drill's result (and the host) in production
readiness once run.
