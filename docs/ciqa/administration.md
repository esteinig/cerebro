# ciqa — administration

How to deploy, configure, and feed CIQA + META-GPT: the configuration reference, the GPU image
build, model download, secrets, and the prefetch source setting. For what the system is, see the
[overview](overview.md) and [architecture](architecture.md); for running it day to day, see
[maintenance](maintenance.md) and the [Nextflow GPU runbook](nextflow-gpu.md).

## Prerequisites

- **For `local` prefetch + local diagnosis (self-contained):** a host with an NVIDIA GPU, the
  container runtime (apptainer with `--nv`, or Docker with the NVIDIA Container Toolkit), the
  CIQA GPU image (below), and the **model weights** mounted or downloadable. **No Cerebro stack
  is required.**
- **For `stack` prefetch and for persisting QC datasets / baselines / reports:** a reachable
  Cerebro stack (API URL + token, and a team / database / project), as for any Cerebro client.
- **For regression against a baseline:** a registered QC dataset with reference truth and a
  promoted baseline (see [datasets](datasets.md), [regression-testing](regression-testing.md)).

## Client connection (global options)

CIQA shares Cerebro's client options. The common ones:

| Option | Env | Meaning |
|---|---|---|
| `--api-url` | `CEREBRO_API_URL` | Cerebro API base URL (only needed for `stack` mode and persistence) |
| `--token` / `--token-file` | `CEREBRO_API_TOKEN` | API token (prefer the env/file form; do not put tokens on the command line) |
| `--team`, `--db`, `--project` | — | Team, database, and project scoping for requests that need it |
| `--fs-*` | — | SeaweedFS master address/port for dataset upload/staging |
| `--danger-accept-invalid-certs` | — | Skip TLS verification — **danger**, non-production only |

In a self-contained local run none of these are required; the connection options are only
consulted when the prefetch source is `stack` or when uploading/persisting to the stack.

## The prefetch source setting

The single switch that decides where the model's input comes from
([`prefetch-filtering.md`](prefetch-filtering.md)):

| `prefetch.source` | Builds `PrefetchData` from | Needs a live stack? | Default |
|---|---|---|---|
| `stack` | the API, after the run's models are uploaded | yes | deployed stack |
| `local` | the run's output models, filtered in-process | no | in-pipeline / offline |

Both sources emit the identical `PrefetchData`. Choose `local` to run the model inside the
pipeline or on an offline validation host; choose `stack` when a stack is already the source of
truth and you want the API's view (including the optional contamination-history enrichment).

## Configuration reference — the diagnostic config

The per-sample configuration that drives filtering and diagnosis is the **`MetaGpConfig`**.
Its fields:

| Field | Type | Meaning |
|---|---|---|
| `sample` | string | Biological sample identifier |
| `sample_type` | `SampleType` | One of `Csf`, `Eye`, `Tis`, `Ntc`, `Pos`, `Lod`, `Env` (CSF, eye, tissue, no-template control, positive control, limit-of-detection, environmental) |
| `test_result` | `TestResult?` | Reference result for the sample (`Positive` / `Negative`) when known — used for scoring |
| `candidates` | `[string]?` | Reference pathogen candidate(s) (normalised species, e.g. `s__…`) when known |
| `exclude_lod` | bool? | Exclude this sample from sensitivity at the limit of detection |
| `identifiers` | schema | The Cerebro identifiers resolving the sample's libraries |
| `filter_configs` | `TieredFilterConfig` | The three-tier filter (below) |
| `contamination` | `PrevalenceContaminationConfig` | Prevalence filtering (below) |

### Tiered filter — `TieredFilterConfig`

Three independent `TaxonFilterConfig`s applied to build the three prefetch tiers:

| Tier | Purpose (typical) |
|---|---|
| `primary` | The broad first-pass evidence set presented to the model |
| `secondary` | A second, usually stricter, evidence view |
| `target` | A focused set (e.g. specific targets / controls) |

Each tier is a full `TaxonFilterConfig` (the same type the core pipeline and the stack use to
filter taxa — abundance, evidence, rank, domain, and tag-based criteria). Tier configs can be
supplied inline or via a JSON file; presets exist for default, development, and "no filtering"
profiles. Because the same `apply_filters` is used here and in the stack, a tier behaves
identically whether the prefetch source is `local` or `stack`.

### Prevalence filtering — `PrevalenceContaminationConfig`

Identifies and removes cross-sample contaminant taxa:

| Field | Type | Meaning |
|---|---|---|
| `threshold` | f64 | Minimum cross-sample prevalence (fraction of samples a taxon appears in) to be called a contaminant |
| `min_rpm` | f64 | Minimum reads-per-million for a taxon to count toward prevalence |
| `sample_type` | string? | Restrict the prevalence computation to a sample type, when set |
| `outliers` | `PrevalenceOutliers` | Per-tier outlier overrides applied alongside the prevalence set |

In `local` mode the prevalence set is computed **over the run's samples**, honouring
`threshold` and `min_rpm`; in `stack` mode it is fetched from the API
(`cerebro/taxa/contamination`). The toggles `--disable-prevalence-control` and
`--disable-prevalence-outlier` turn the two parts off; `--contamination-min-rpm` sets `min_rpm`
at the command line. See [`prefetch-filtering.md`](prefetch-filtering.md) for the full behaviour.

### Post-filter — `PostFilterConfig`

Applied to the model's view of candidates after filtering, to keep the presented set clean:

| Field | Type | Meaning |
|---|---|---|
| `collapse_variants` | bool | Collapse strain/variant rows into their species |
| `best_species` | bool | Keep only the best-supported species |
| `best_species_min` | usize | Minimum support to qualify as "best species" |
| `best_species_domains` | `[string]` | Restrict "best species" to these domains |
| `best_species_base_weight` | f64? | Base weight in the best-species scoring |
| `exclude_phage` | bool | Drop phage from the presented set |
| `exclude_phage_list` | set | Explicit phage exclusions |

## The META-GPT GPU image

Build the production GPU image from the repository base (the Dockerfile expects the repo as
build context):

```bash
docker build -f templates/stack/docker/Dockerfile.gpt.dgx -t cerebro-metagpt:<tag> .
```

It is an NVIDIA CUDA 12.2 image that builds `cerebro-ciqa` with the `local` feature and ships
the `cerebro-ciqa` and `cerebro-client` binaries, with the NVML symlink handled at entrypoint.
If the run uses apptainer, convert it:

```bash
apptainer build cerebro-metagpt.sif docker-daemon://cerebro-metagpt:<tag>
```

Point Nextflow at the image via the META-GPT resource parameter (see
[`nextflow-gpu.md`](nextflow-gpu.md)). The `Dockerfile.gpt.dev` variant is a development shell,
not for production.

## Model download

META-GPT loads a local model from a directory at run time; **weights are never baked into the
image** (D11 in the integration plan). Two supported shapes, chosen by configuration:

- **Mounted weights.** Mount a directory of weights into the container and point the model
  directory at it. Nothing is downloaded.
- **Download by configuration.** When download is enabled, the weights are fetched into the
  model directory before diagnosis (e.g. by model id). Treat downloaded weights as a
  provenance-tracked artifact: record the source and verify a hash before use, as for any
  externally sourced clinical-path dependency.

The model id, the model directory, and the download toggle are configuration
([`nextflow-gpu.md`](nextflow-gpu.md) lists the Nextflow params; the CLI exposes the
equivalent `--model` / `--model-dir` options on the diagnose command).

## Secrets and configuration hygiene

- API tokens come from `CEREBRO_API_TOKEN` or a token file, never the command line or a
  committed file.
- Model-weight sources (and any credentials they need) are configuration, documented per
  deployment; downloaded weights are hash-verified.
- TLS verification is on in production; `--danger-accept-invalid-certs` is for local
  development only.
- All persisted CIQA artifacts (datasets, baselines, reports) inherit the catalogue's
  durability and audit; see [`cerebro-fs`](../cerebro-fs/administration.md) for the underlying
  storage configuration.

## Tested here vs your environment

The configuration *shapes* above are read from the code; the *operational* behaviour of the GPU
image, model download, and stack persistence is environment-specific and is confirmed by the
checks in [`validation.md`](validation.md) and the runbooks. In particular, the exact
reads-per-million predicate used by prevalence filtering, and the per-deployment model-download
mechanism, are confirmed against your stack and GPU host, not asserted here.

> **Commands & endpoints.** For copy-pasteable CLI and curl for every step, see [Commands & API](commands.md).
