# Commands & API (CLI and curl)

Concrete, copy-pasteable commands and endpoints for the CIQA subsystem — the prefetch →
diagnose → regression path and the report generation, as actually implemented. Prose
descriptions of the data path live in [architecture](architecture.md); this page is the
command reference.

Conventions used below:

```bash
export CEREBRO_API_URL="https://cerebro.example.org/api"
export CEREBRO_API_TOKEN="…"          # bearer token (or use --token-file)
export TEAM="MyTeam" DB="MyDb" PROJECT="MyProject"
```

All authenticated API calls take `Authorization: Bearer $CEREBRO_API_TOKEN` and pass team /
database / project as **query parameters** (`?team=…&db=…&project=…`) — the same scheme the
`cerebro-client` uses internally.

---

## 1. Prefetch — build tiered, prevalence-filtered `PrefetchData`

The prefetch step produces one `{sample}.prefetch.json` per sample. It runs from a live stack
(`--prefetch-source stack`) or **self-contained** from a run's own Cerebro models
(`--prefetch-source local`, no API needed).

```bash
# Self-contained (local): taxa come from the run's output models — no stack required.
cerebro-ciqa prefetch \
  --prefetch-source local \
  --run-models results/models/*.json \
  --plate-json plate.json \
  --outdir prefetch

# With the contam_history rescue (Task B): rescue prevalence-flagged contaminants that are
# genuinely elevated in the sample, using a run-local regression. off | run | stack | run-plus-stack.
cerebro-ciqa prefetch \
  --prefetch-source local --run-models results/models/*.json \
  --plate-json plate.json --outdir prefetch \
  --contam-history run

# From the live stack instead (queries the API for taxa + prevalence):
cerebro-ciqa prefetch \
  --prefetch-source stack \
  --plate-json plate.json --outdir prefetch \
  --url $CEREBRO_API_URL --token $CEREBRO_API_TOKEN \
  --team $TEAM --db $DB --project $PROJECT
```

Useful flags: `--contamination <config.json>`, `--contamination-min-rpm <f64>`,
`--tiered-filter <config.json>`, `--disable-prevalence-control`, `--disable-filter`,
`--samples <id…>` (subset), `--force` (overwrite), `--threads <n>`.

**Equivalent stack endpoints (curl).** The `stack` source is built on these read endpoints:

```bash
# Filtered, tiered taxa for a sample set (the stack prefetch source):
curl -H "Authorization: Bearer $CEREBRO_API_TOKEN" \
  "$CEREBRO_API_URL/cerebro/taxa?team=$TEAM&db=$DB&project=$PROJECT"

# Cross-sample taxon history used by the contam_history rescue (per-taxon vs Homo sapiens):
curl -H "Authorization: Bearer $CEREBRO_API_TOKEN" \
  "$CEREBRO_API_URL/cerebro/taxa/history?taxon_label=s__Escherichia%20coli&host_label=s__Homo%20sapiens&team=$TEAM&db=$DB&project=$PROJECT"
```

---

## 2. Diagnose — local-GPU META-GPT (the stable seam)

Reads `{sample}.prefetch.json`, runs the agent on a local GPU, and writes
`{sample}.model.json` (a `DiagnosticResult`). With `--emit-manifest` it also writes a
`MetaGptRunManifest` for the regression step. The diagnosis seam itself is unchanged.

```bash
cerebro-ciqa diagnose-local \
  --prefetch prefetch \
  --plate plate.json \
  --outdir diagnose \
  --model qwen3-8b-q8-0 --model-dir /weights \
  --num-gpu 1 \
  --post-filter \
  --emit-manifest diagnose/run.manifest.json \
  --run-id RUN-2026-06-001 --model-id qwen3-8b-q8-0 --replicate 1
```

Key flags: `--clinical-notes`, `--sample-context`, `--temperature`, `--sample-len`,
`--num-gpu <n>` (batches samples across GPUs), `--force`. Requires the `local` feature build
(the DGX image ships it). VRAM is captured to `diagnose/state_logs/gpu*.bench.json`.

---

## 3. Regression — score a run against a baseline (report-first)

Evaluates a run's `DiagnosticResult`s (via the manifest) against a QC dataset's reference
truth, compares to a stored baseline, runs McNemar, and writes a report. The exit code
reflects the gate; the **report is the primary artifact**.

```bash
cerebro-ciqa regression \
  --manifest diagnose/run.manifest.json \
  --baseline baseline.json \
  --dataset qc_dataset.json \
  --report regression.report.json
# exit 0 = within tolerance; non-zero = regressed (the report carries the reasons either way)
```

---

## 4. Reports — render PDFs

`cerebro-report` renders the regression verdict and (optionally) the META-GPT assessment as
PDFs through the Typst pipeline.

```bash
# Regression report (Stage 08):
cerebro-report render \
  --report-type Regression \
  --config templates/regression_report.json \
  --format pdf --output regression.pdf

# Pathogen report with the generative Appendix C (a config carrying a "meta_gpt" block):
cerebro-report render \
  --report-type PathogenDetection \
  --config sample.report.json \
  --format pdf --output report.pdf
```

---

## 5. Upload — register results on the stack

```bash
# Cerebro models (and, by extension, the diagnosis outputs alongside them):
cerebro-client --token $CEREBRO_API_TOKEN --url $CEREBRO_API_URL \
  --team $TEAM --db $DB --project $PROJECT \
  upload-models --models results/models/*.json
```

Equivalent endpoint (model insert):

```bash
curl -X POST -H "Authorization: Bearer $CEREBRO_API_TOKEN" \
  -H "Content-Type: application/json" \
  --data @model.json \
  "$CEREBRO_API_URL/cerebro?team=$TEAM&db=$DB&project=$PROJECT"
```

---

## 6. Nextflow — run the whole thing on a GPU host

Gated by `params.pathogenDetection.metaGpt.enabled` (off by default). See
[Nextflow GPU diagnosis & regression](nextflow-gpu.md) for the full runbook.

```bash
# Standalone: re-run diagnosis/regression over an existing run's models, no re-profiling.
nextflow run main.nf -entry metagpt -profile dgx \
  --pathogenDetection.metaGpt.enabled true \
  --pathogenDetection.metaGpt.modelsDir <run>/results/models \
  --pathogenDetection.metaGpt.plate plate.json \
  --pathogenDetection.metaGpt.modelDir <weights>
```

---

## 7. How do I know it worked?

| After… | Check |
|---|---|
| prefetch | one `{sample}.prefetch.json` per sample; tiers + `*_contamination` populated; with `--contam-history`, previously-excluded elevated taxa reappear in signal |
| diagnose | one `{sample}.model.json` per sample; `run.manifest.json` lists them with `content_hash`; `state_logs/gpu*.bench.json` shows non-zero peak VRAM |
| regression | `regression.report.json` with `regressed` true/false + reasons; exit code matches |
| report | a PDF with the verdict box / the generative appendix |

> Validation boundary: the GPU/Nextflow/Typst steps are validated on a real GPU host (see
> [validation](validation.md)); the pure prefetch/regression logic is unit-tested.
