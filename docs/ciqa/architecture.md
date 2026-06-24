# ciqa — architecture

A higher-level technical view of how CIQA and META-GPT are put together: the components, how
data flows from a run's output models to a regression verdict, the crate layering that
determines where logic lives, and the GPU container that runs the model. For the conceptual
summary see the [overview](overview.md); for operating it see the
[administration](administration.md) and [maintenance](maintenance.md) guides.

## Components

```
   ┌──────────────────────────────────────────────────────────────────────────┐
   │                          Cerebro run (Nextflow)                           │
   │   reads → quality control → taxonomic profile → ProcessOutput → models    │
   └───────────────────────────────────┬──────────────────────────────────────┘
                                        │  output models (Cerebro JSON, carry `taxa`)
                                        ▼
        ┌───────────────────────────────────────────────────────────────┐
        │                      CIQA  (cerebro-ciqa)                       │
        │                                                                 │
        │   prefetch ──► PrefetchData ──► diagnose ──► DiagnosticResult   │
        │   (tiered +                      (META-GPT,      + run manifest │
        │    prevalence)                    local GPU)                    │
        │        │                                            │           │
        │        │  source = local │ stack                    │           │
        │        ▼                                            ▼           │
        │   apply_filters / aggregate / prevalence       evaluate vs      │
        │   (reused from cerebro-pipeline)                reference truth  │
        │                                                     │           │
        └─────────────────────────────────────────────────────┼──────────┘
                       ▲ (source = stack)                       ▼
            ┌──────────┴───────────┐                  ┌──────────────────────┐
            │   Cerebro stack API   │                  │  regression vs        │
            │  (taxa, prevalence,   │                  │  versioned baseline   │
            │   ciqa/* datasets,    │◄─────────────────┤  (McNemar + floor)    │
            │   baselines, reports) │   persist/report │  → report + gate      │
            └──────────┬───────────┘                  └──────────────────────┘
                       ▼
            ┌──────────────────────┐
            │  MongoDB + GridFS    │  QC datasets, baselines, regression reports,
            │  (catalogue + audit) │  PrefetchData blobs; consequential actions audited
            └──────────────────────┘
```

- **Cerebro run (Nextflow)** — produces the per-sample/library **Cerebro model** JSON. Each
  model carries the taxonomic profile in its `taxa` field; these are CIQA's universal input.
- **CIQA (`cerebro-ciqa`)** — the command-line and in-pipeline tool. It builds the prefetch,
  drives META-GPT, evaluates the result, and runs regression. It is a workspace member that
  depends on the in-tree crates by path.
- **META-GPT** — the generative practitioner, an external crate consumed by CIQA behind the
  stable `run_local` seam. On a local GPU it runs inside the NVIDIA CUDA container; it is not
  modified by CIQA.
- **Cerebro stack API** — optional. In `stack` prefetch mode CIQA fetches filtered taxa and
  prevalence-contaminants from the API; in `local` mode it does not. The API also persists QC
  datasets, baselines, and regression reports via the `ciqa/*` endpoints.
- **MongoDB + GridFS** — the durable system of record for QC datasets, baselines, and reports;
  `PrefetchData` blobs live in GridFS exactly as the training prefetch records do. Consequential
  actions are written to the append-only audit chain.

## Data flow

One run, end to end:

1. **Output models → taxa.** CIQA reads the run's Cerebro model JSON. Each model's `taxa`
   (`Vec<Taxon>`) and its DNA/RNA tags are the raw material.
2. **Tiered filtering.** For each of the three tiers (**primary / secondary / target**), the
   taxa are aggregated across the sample's libraries and filtered by that tier's
   `TaxonFilterConfig`. The filtering is the **same function the stack uses** (`apply_filters`),
   so the two prefetch sources agree.
3. **Prevalence split.** Cross-sample contaminant taxids — computed over the run's samples (or
   fetched from the stack) — are separated from each tier into the prefetch's contamination
   vectors. The combined result is pruned of cross-tier and contaminant duplicates into a
   `PrefetchData`.
4. **Diagnosis.** Each `PrefetchData` is handed to META-GPT via `run_local`, which returns a
   `DiagnosticResult` (a `Diagnosis` — infectious / non-infectious / tumour, with review
   variants — plus a pathogen candidate). The diagnose step also emits a typed **run manifest**
   recording the model id, quantization, parameters, clinical flag, replicate, and per-sample
   output paths with content hashes.
5. **Evaluation.** Each `DiagnosticResult` is mapped to a predicted (result, candidate) and
   compared to the sample's reference truth by the **one shared evaluator**, classifying it
   **TP / TN / FP / FN / Excluded** and reducing the set to **sensitivity / specificity / PPV /
   NPV**.
6. **Regression.** The run's statistics are compared to a stored **baseline** for the same
   dataset version: an absolute floor from the baseline's thresholds, plus a paired **McNemar**
   test on per-sample correctness (with Bonferroni / Holm / Benjamini-Hochberg adjustment when
   several configurations are compared). The result is a **regression report** and an exit-code
   gate; nothing about the baseline is mutated.

## Crate / dependency layering — where logic lives

The placement of new logic follows the workspace's dependency direction, which is a clean
stack (lower layers know nothing of higher ones):

```
cerebro-pipeline   base: Taxon, TaxonFilterConfig, apply_filters, aggregate, prevalence
      ▲                  (pure filtering primitives; no internal cerebro deps)
cerebro-model      PrefetchData, MetaGpConfig, TieredFilterConfig, the Cerebro model w/ taxa,
      │            the shared diagnostic evaluator, training + QC dataset records
      ▼            (depends on cerebro-pipeline)
cerebro-client / cerebro-fs    REST + object-store clients
      ▼
cerebro-server     the API: cerebro/taxa, cerebro/taxa/contamination, ciqa/* endpoints
      ▼
cerebro-ciqa       the tool: prefetch, diagnose (META-GPT), evaluate, regression; the CLI
```

Consequences that matter:

- **Filtering primitives** (`apply_filters`, `aggregate`, and the prevalence counter) live in
  **`cerebro-pipeline`**, the base layer. The server and CIQA both call them, so there is one
  implementation of "what gets filtered".
- **The self-contained prefetch builder** lives in **`cerebro-model`**, because it needs to see
  both `PrefetchData` (model) and the pipeline filters. It takes taxa as plain data, so it is
  unit-testable with synthetic input and no I/O.
- **The diagnostic evaluator** lives in **`cerebro-model`** and is shared by the training review
  and the QC/regression path — one definition of TP/FP/TN/FN and of sensitivity/specificity.
- **The model and CIQA never depend on the server**; the prefetch and evaluation logic is pure
  and runs with or without a stack.

## The META-GPT container

Local-GPU diagnosis runs inside an NVIDIA CUDA container built from
`templates/stack/docker/Dockerfile.gpt.dgx`:

- A multi-stage build on `nvidia/cuda:12.2.2-cudnn8-devel-ubuntu20.04` → `…-runtime-…`.
- Builds `cerebro-ciqa` with the `local` feature (which wires the local META-GPT generator) and
  ships the `cerebro-ciqa` and `cerebro-client` binaries.
- Handles the NVML library symlink (`libnvidia-ml.so`) at entrypoint so GPU VRAM monitoring
  works inside the container.
- **Model weights are not baked in** — they are mounted or downloaded at run time, depending on
  configuration ([`administration.md`](administration.md), [`meta-gpt-diagnosis.md`](meta-gpt-diagnosis.md)).

In Nextflow the container is bound to a `metaGpt` process label that requests the GPU
(apptainer `--nv` / docker `--gpus`), exactly as other tools are bound to their containers; see
[`nextflow-gpu.md`](nextflow-gpu.md).

## Where state lives, and what must be durable

| State | Lives in | Durable how |
|---|---|---|
| Run output models (`taxa`) | The run's work/results directory | The pipeline's own outputs; transient input to CIQA |
| `PrefetchData`, `DiagnosticResult`, run manifest | The diagnose output directory; optionally uploaded | Files first; persisted to GridFS when uploaded |
| QC datasets, baselines, regression reports | MongoDB + GridFS via `ciqa/*` | The catalogue's durability (see `cerebro-fs` docs) |
| Audit of baseline promotion / dataset delete | MongoDB append-only audit chain | Same as the catalogue |
| Model weights | Mounted volume or download cache | External to CIQA; provenance tracked (see administration) |

For how the catalogue and GridFS themselves are kept durable and recoverable, see the
[`cerebro-fs`](../cerebro-fs/index.md) documentation; CIQA's records are ordinary catalogue
state subject to the same backup, verification, and recovery machinery.
