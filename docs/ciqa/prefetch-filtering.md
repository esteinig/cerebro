# ciqa — prefetch & filtering

How CIQA turns a run's taxonomic profile into the **prefetch** — the tiered, prevalence-filtered
evidence set that is META-GPT's input — and how the same prefetch is produced **either from the
live stack or self-contained from the run's output models**. For the configuration fields, see
[administration](administration.md); for what happens to the prefetch next, see
[META-GPT diagnosis](meta-gpt-diagnosis.md).

## What the prefetch is

A `PrefetchData` is the model's view of one sample's evidence, organised into three **tiers** and
separated from contaminants:

- **primary / secondary / target** — three filtered taxon sets, each built by applying that
  tier's `TaxonFilterConfig` to the sample's taxa.
- **primary / secondary / target contamination** — the taxa that were removed from each tier as
  prevalence-driven contaminants, kept alongside so the decision is transparent.
- the **`MetaGpConfig`** that produced it (sample, type, the tier configs, the prevalence config,
  and the reference truth when known).

After the tiers and contamination are assembled the prefetch is **pruned** of cross-tier and
contaminant duplicates, so the model sees each taxon once, in its most specific tier.

## How filtering works

Two filters compose to build the prefetch.

### Tiered filtering

For each tier the sample's taxa are first **aggregated** across the sample's libraries (DNA and
RNA), then **filtered** by that tier's `TaxonFilterConfig`. The filter is the **same function the
core pipeline and the stack use** — abundance, evidence, rank, domain, and tag-based criteria —
so a tier means exactly the same thing here as it does in the stack's filtered-taxa view. The
three tiers let the model be shown a broad set, a stricter set, and a focused/target set
simultaneously.

### Prevalence filtering

Some taxa appear across many samples because they are reagent or environmental contaminants
rather than true findings. Prevalence filtering identifies them by **how often a taxon appears
across the samples**, above a reads-per-million floor, and separates them into the contamination
vectors:

- a taxon whose cross-sample prevalence is at or above `threshold` (counting only taxa above
  `min_rpm`) is treated as a contaminant;
- per-tier **outlier** overrides can adjust this alongside the prevalence set.

This is the prevalence control the [administration](administration.md) page documents; it is
toggled by `--disable-prevalence-control` / `--disable-prevalence-outlier` and tuned by
`--contamination-min-rpm`.

## The two sources — `stack` and `local`

The prefetch is identical whichever way it is built; the difference is **where the taxa and the
prevalence set come from**.

### `stack` — from the live API

The classic path. CIQA resolves the sample's libraries through the API, POSTs each tier's filter
to the stack, and the **server** aggregates and filters the stored taxa and returns the result;
the prevalence-contaminant taxids come from the stack's collection-wide contamination endpoint.
Requires a reachable, populated stack (the run's models must already be uploaded).

```bash
cerebro-ciqa --api-url "$CEREBRO_API_URL" --token-file token \
  prefetch --prefetch-source stack --plate-json plate.json --outdir prefetch
```

### `local` — self-contained from the run's output models

The self-contained path. CIQA reads the run's **Cerebro model JSON** (the "output models",
each carrying `taxa`), aggregates and filters them **in-process** with the same functions the
server uses, and computes the prevalence set **over the run's samples**. **No stack is required.**

```bash
cerebro-ciqa prefetch --prefetch-source local --run-models <run>/results/models/*.json --plate-json plate.json --outdir prefetch
# prevalence computed over the run; toggles apply:
#   --disable-prevalence-control   --contamination-min-rpm <x>
```

This is what lets META-GPT run inside the pipeline on a local GPU with nothing else attached
([`nextflow-gpu.md`](nextflow-gpu.md)).

### Why the two agree — and the `contam_history` rescue

The `local` path reuses the **exact filtering and aggregation functions** the stack calls, and
the prevalence logic is one shared implementation, so the two produce a byte-identical
`PrefetchData` on the same input — verified by a differential check on a reference sample, not
merely asserted.

The stack path can also optionally **rescue** prevalence-flagged contaminants that are genuinely
elevated in *this* sample, using a per-taxon regression of taxon-RPM vs host (Homo sapiens) RPM
across samples. This `contam_history` rescue is now available **self-contained in `local` mode**
too (the regression history is built from the run's own models), reusing the same `RpmAnalyzer`
and the same rescue rule as the stack — so the decisions match given the same history. It is
selected with `--contam-history`:

```bash
# off (default) | run (self-contained) | stack | run-plus-stack
cerebro-ciqa prefetch \
  --prefetch-source local --run-models <run>/results/models/*.json \
  --plate-json plate.json --outdir prefetch \
  --contam-history run
```

`off` reproduces the plain `local` behaviour (no rescue). `run` is fully self-contained;
`stack` / `run-plus-stack` enrich with the collection history and need a configured API. See the
full command reference in [Commands & API](commands.md).

## Worked example

Building a prefetch self-contained for every sample in a finished run:

```bash
# 1. point at the run's output models; build the prefetch with default tiers and prevalence
cerebro-ciqa prefetch \
  --source local \
  --models results/run-2024-06/cerebro_models/ \
  --plate plate.json \
  --outdir results/run-2024-06/prefetch/

# 2. inspect: one {sample}.prefetch.json per sample, each with primary/secondary/target
ls results/run-2024-06/prefetch/
# sampleA.prefetch.json  sampleB.prefetch.json  ...
```

The same command with `--prefetch-source stack` (and API options) produces the same files from the stack. For the full flag set and curl equivalents, see [Commands & API](commands.md).

## Tested here vs your environment

The tiering, aggregation, prevalence counting, and prune are **pure logic with hermetic unit
tests** — they run and are checked without a stack. **Parity** between `local` and `stack` on a
real sample, and the exact reads-per-million predicate, are confirmed in your environment by the
differential check ([`validation.md`](validation.md)). Whether a run emits one model per library
or a combined model (which the grouping depends on) is confirmed against your pipeline's output.
