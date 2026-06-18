# Cerebro

**Metagenomic and metatranscriptomic host–pathogen diagnostics for clinical and
public-health production environments.**

Cerebro turns metagenomic sequencing data into auditable, clinically reportable
pathogen-detection results. It pairs a set of Nextflow analysis pipelines with an
application stack (API, database, object storage, workers, and a collaborative
reporting web app) so a diagnostic laboratory can run metagenomics as a service —
from sequencing run to clinical report — with the integrity, provenance, retention,
and recoverability clinical work demands.

> **Research today, production on this branch.** The pipelines and tools below are in
> active research and reference-laboratory use now. The `feat/production` branch
> extends them into a hardened, continuously operated production service: a durable
> storage and lifecycle subsystem (`cerebro-fs`), governed maintenance and recovery,
> and an automated, tiered-threshold candidate-review path for pathogen detection.
> Sections marked *(production, on this branch)* describe that extension.

---

<details>
<summary>🩸 <b>Metagenomic diagnostic core functions</b></summary>
<br>

**Core**

- Multi-classifier taxonomic profiling, metagenome assembly, and alignment in Nextflow
- Optimised pangenome host and background depletion with [`Scrubby`](https://github.com/esteinig/scrubby)
- Viral infection detection, pan-viral enrichment protocols, and syndrome-specific
  subtyping panels using [`Vircov`](https://github.com/esteinig/vircov)
- Differential host tumour-DNA diagnostics via segmental CNV detection and methylation
  classification (Sturgeon)

**Support**

- Species identification with GTDB reference assemblies and the GTDB-Tk toolkit
- MAG recovery from enriched culture and sample co-assembly; unclassified viral-bin and
  novel-virus prediction (geNomad, RdRP)
- Custom database and index construction, grafted taxonomies, and genome cleaning with
  `Cipher`

</details>

<details>
<summary>📰 <b>Collaborative clinical reporting (Bug Board)</b></summary>
<br>

- Collaborative, auditable pathogen determination from metagenomic sequencing results
- Multi-tenant Svelte application and API with secure local or web-server deployment
- Scalable stack deployment with configurable data-security and collaboration models
- Stack configuration and deployment integrated into the primary command-line interface
- Clinical reporting with [`Typst`](https://typst.app)-formatted templates linked to the
  underlying multi-classifier evidence
- Secure, `wasm`-enabled in-browser report generation for sensitive reports, with
  interactive data visualisations
- Auditable team comments and results discussion for expert-panel review ("Bug Board")

</details>

<details>
<summary>🏥 <b>Clinical and public-health production environments</b></summary>
<br>

- In-silico syndromic reference panels for signal- and read-level simulation with `Cipher`
- Quality assurance and regression evaluation of simulated and patient datasets
- Contamination assessment (background / sample-site / kitome) for clinical environments
- Distributed sequence and analysis storage, file-system lifecycle, and data-retention
  policies through [`SeaweedFS`](https://github.com/seaweedfs/seaweedfs) (the `cerebro-fs`
  subsystem)
- Standard operating procedures and runbooks for continuous operation as a clinical
  diagnostic service

</details>

---

## Subsystems at a glance

| Subsystem | Role | Status |
|---|---|---|
| **Pipelines** (Nextflow) | Analysis: QC, taxonomic profiling, assembly, subtyping, CNV | Research; production entry points on this branch |
| **API + database** | System of record: results, catalogue, append-only audit | Research; hardening on this branch |
| **Reporting app** (Svelte) | Collaborative determination and clinical reports (Bug Board) | Research; hardening on this branch |
| **cerebro-fs** | Storage, lifecycle, durability, backup, recovery | Production focus of this branch — see [documentation](https://esteinig.github.io/cerebro/) |
| **Workers** (Faktory) | Tiering, verification, reconcile, backup, archival, restore | Production focus of this branch |
| **CIQA + META-GPT** | Diagnostic classification, candidate prioritisation, QA/regression | Research; automation planned (see below) |

Full operations documentation for the production subsystems is published at
**[esteinig.github.io/cerebro](https://esteinig.github.io/cerebro/)** (Home →
Subsystems → cerebro-fs).

## Getting started

> [!NOTE]
> The full Docker application stack is **not** required for the core metagenomic
> pipelines and report generation. The Nextflow pipelines can be run on their own, and
> the Cerebro CLI used for processing pipeline outputs and generating clinical reports.

**Minimum requirements**

- Linux
- Nextflow ≥ v24
- Conda / Mamba / Docker (or Apptainer)

Computational requirements range from a single workstation for the application stack to
large server infrastructure for pipelines and the web application, depending on the
number of laboratories, sequencing throughput, storage, and the chosen infrastructure,
data-security, and collaboration model.

## Pipelines

All pipelines are entry points of the main workflow, selected with `-entry`:

```bash
nextflow run -r v1.0.0 https://github.com/esteinig/cerebro \
  -profile <resources>,mamba \
  -entry <pathogen|panviral|culture|quality|aneuploidy|production> \
  --outputDirectory out/ \
  --databaseDirectory db/
```

| Entry | Pipeline | Primary use |
|---|---|---|
| `pathogen` | Pathogen detection | Low-biomass sterile-site diagnostics (needle-in-a-haystack) |
| `panviral` | Pan-viral enrichment | Probe-capture viral panels, genome recovery, subtyping |
| `culture` | Culture & MAG assembly | Hybrid assembly, prokaryotic MAG typing, novel-virus detection |
| `quality` | Quality control | Read QC, host/background depletion, control accounting |
| `aneuploidy` | Host CNV / aneuploidy | Differential host tumour-DNA diagnostics |
| `production` | Integrated production | App-staged, continuous clinical operation |

### Pathogen detection — `-entry pathogen`

The primary diagnostic pipeline for high-human, low-microbial-biomass sterile-site
samples (validated for ocular fluid and cerebrospinal fluid), where distinguishing a
true pathogen from contamination is the core challenge.

It runs **quality control and host/background depletion** ([`Scrubby`](https://github.com/esteinig/scrubby)),
then **multi-classifier taxonomic profiling** — k-mer and alignment-based classifiers
(Kraken2 + Bracken, Metabuli, Sylph, KMCP, Ganon) and reference-alignment coverage
([`Vircov`](https://github.com/esteinig/vircov)) — and a parallel **metagenome assembly**
track (metaSPAdes / MEGAHIT) with binning (CONCOCT, MetaBAT2, SemiBin2) and contig-level
classification (BLAST / DIAMOND). All evidence is aggregated into a single Cerebro model
and uploaded to the API for collaborative review.

*(Production, on this branch)* A **prefetch** step then applies a **tiered-threshold
configuration** to the aggregated evidence, prioritising candidate taxa into evidence
tiers (e.g. strong / supporting / weak) across classifiers, coverage, and assembly. This
tiered-threshold prefetch is the designated **automation decision point** at which
candidate review is triggered — surfaced to clinicians on the Bug Board today, and the
insertion point for the automated review described under **CIQA & META-GPT** below.

```bash
nextflow run -r v1.0.0 https://github.com/esteinig/cerebro \
  -profile dgx,large,mamba,cns,cipher \
  -entry pathogen \
  --outputDirectory out/ \
  --databaseDirectory db/ \
  --fastqPaired 'fastq/*_{R1_001,R2_001}.fastq.gz'
```

This is **not** intended for high-biomass sample types (respiratory, environmental),
where an abundant background microbiome is the main challenge (haystack-full-of-needles).

### Pan-viral enrichment — `-entry panviral`

For probe-hybridisation capture panels (e.g. Agilent or Twist) targeting viral
pathogens. After QC and depletion, [`Vircov`](https://github.com/esteinig/vircov)
evaluates per-reference alignment coverage to call confident viral diagnoses from low
titres and high host background, recovers consensus genomes, and applies
species-to-lineage **subtyping schemes** for automated genotyping. Short-read and
Nanopore inputs are supported.

```bash
nextflow run -r v1.0.0 https://github.com/esteinig/cerebro \
  -profile large,mamba \
  -entry panviral \
  --outputDirectory out/ --databaseDirectory db/ \
  --fastqPaired 'fastq/*_{R1,R2}.fastq.gz'
```

### Culture & metagenome-assembled genomes (MAG) — `-entry culture`

A distinct assembly-centric pipeline for bacterial culture identification, enriched
culture, and sample co-assembly. It performs **hybrid (short- + long-read) assembly**,
recovers **metagenome-assembled genomes** through binning, and assigns **prokaryotic
species-level taxonomy with the GTDB-Tk toolkit** against GTDB reference assemblies. The
same assembly graph feeds **novel-virus and unclassified viral-bin detection** (geNomad,
RdRP), extending diagnosis beyond reference-contained taxa.

```bash
nextflow run -r v1.0.0 https://github.com/esteinig/cerebro \
  -profile large,mamba \
  -entry culture \
  --outputDirectory out/ --databaseDirectory db/ \
  --fastqPaired 'fastq/*_{R1,R2}.fastq.gz' \
  --fastqNanopore 'fastq/*.nanopore.fastq.gz'
```

### Quality control — `-entry quality`

Read quality control and adapter/quality trimming, host and background depletion
([`Scrubby`](https://github.com/esteinig/scrubby)), and accounting of spike-in and
ERCC/sequencing controls — the foundation every diagnostic entry shares, runnable on its
own for QC reporting. Illumina and Nanopore inputs are supported (`--nanopore`).

```bash
nextflow run -r v1.0.0 https://github.com/esteinig/cerebro \
  -profile medium,mamba \
  -entry quality \
  --outputDirectory out/ --databaseDirectory db/ \
  --fastqPaired 'fastq/*_{R1,R2}.fastq.gz'
```

### Host CNV / aneuploidy — `-entry aneuploidy`

Differential host tumour-DNA diagnostics from the host background of short-read
metagenomic data: large segmental copy-number-variation (CNV) detection across human
chromosomes with CNVkit, complemented by methylation-based classification (Sturgeon).
Reads are aligned with minimap2 against a host reference and compared to a normal
control.

```bash
nextflow run -r v1.0.0 https://github.com/esteinig/cerebro \
  -profile medium,mamba \
  -entry aneuploidy \
  --outputDirectory out/ \
  --referenceFasta host.fasta \
  --normalControlBam normal.bam \
  --fastqPaired 'fastq/*_{R1,R2}.fastq.gz'
```

### Integrated production — `-entry production`

The continuously operated mode used by the application stack. Sample files are staged by
the stack (via `cerebro-fs`), the pipeline branch is selected per sample
(pan-viral / pathogen / genome-assembly), and results are written back to the API and
catalogue with full provenance. This entry is normally driven by the stack rather than
invoked by hand; see the [documentation](https://esteinig.github.io/cerebro/) for
deployment and operation.

## Application stack

The stack is configured and deployed through the Cerebro CLI, which renders a
`docker compose` deployment and its secrets from a single configuration. Its services:

- **API + database** — the system of record for results, the file catalogue, and an
  append-only audit trail.
- **Reporting web app** — the multi-tenant Svelte application for collaborative
  determination and clinical reporting (Bug Board, Typst templates, in-browser `wasm`
  report generation).
- **cerebro-fs** — the storage, lifecycle, and durability subsystem: replicated object
  storage, tiering, archival, integrity verification, catalogue backup, and recovery.
  This is the production focus of this branch; see its
  [documentation](https://esteinig.github.io/cerebro/) (overview, architecture,
  administration, maintenance, disaster-recovery runbook, and production-readiness
  governance).
- **Workers** — scheduled and on-demand maintenance (tiering, verification, reconcile,
  backup, archival, restore) behind a detection / safe-automation / operator-gated
  taxonomy.

## CIQA & META-GPT — assisted diagnostic review and quality assurance

> *This section describes a designed and partially implemented capability that
> integrates the **CIQA** quality-assurance harness and the **META-GPT** diagnostic
> reasoning component with Cerebro's evidence model. It operates as decision support
> under expert oversight, not autonomous diagnosis, and is under active development.*

Cerebro aggregates every classifier, coverage, and assembly signal for a sample into a
single evidence model. **META-GPT** is the diagnostic-classification component that
reasons over that aggregated evidence the way an expert reviewer reads the Bug Board: it
consumes the **prefetch** evidence — the candidate taxa and their supporting signals,
organised by the **tiered-threshold configuration** — and produces a prioritised,
explained candidate determination.

- **Diagnostic classification and candidate prioritisation.** META-GPT ranks and
  short-lists candidate pathogens from the tiered prefetch evidence, prioritising review
  toward the strongest-supported taxa and explaining the basis for each call. Its output
  feeds candidate prioritisation in the clinical report and the Bug Board as decision
  support for the reviewing scientist.
- **Clinical-notes-anchored, local-only review.** Models are anchored with de-identified
  clinical notes to emulate the clinical review of the tiered-threshold evidence, and run
  **locally only** — no patient-derived data leaves the deployment. This keeps assisted
  review inside the laboratory's data-governance boundary.
- **Quality assurance and regression testing (CIQA).** Curated prefetch datasets carry
  **reference candidates** (the known-correct determination). CIQA replays them through
  the pipeline and the review model and scores the outcome against the reference —
  whether the determination matched the expected candidate — turning pipeline, database,
  and model changes into measurable **regression tests**. This reuses Cerebro's training
  subsystem, in which a review produces a `Decision` and a candidate selection that is
  compared to the reference candidates to record whether any expected candidate was
  matched.
- **Tiered-threshold prefetch as the automation point** *(to be implemented)*. The
  tiered-threshold prefetch step in the pathogen pipeline is the designated point at which
  this assisted review is invoked automatically: evidence crossing the configured tiers
  triggers a META-GPT review pass whose prioritised, explained candidates are presented
  for expert sign-off. Today this step surfaces the tiered evidence for manual review;
  the automated invocation is the planned production extension.

## Databases and taxonomy

Reference databases, indices, grafted taxonomies, and in-silico syndromic simulation
panels are constructed and maintained with `Cipher`, which also supports genome cleaning
and syndromic diversity injection for quality assurance.

## Status

Under active development for production release. The `feat/production` branch hardens the
storage, lifecycle, and operational subsystems for continuous clinical use; subsystem
documentation states, where relevant, the boundary between what is validated and what a
given deployment must confirm in its own environment.

This is a preliminary public release of code for the viral-enrichment branch of the
pipeline used in:

> Michael A Moso, George Taiaroa, Eike Steinig, *et al.* **Non-SARS-CoV-2 respiratory
> viral detection and whole genome sequencing from COVID-19 rapid antigen test devices:
> a laboratory evaluation study.** *The Lancet Microbe* (2024).
> [10.1016/S2666-5247(23)00375-0](https://doi.org/10.1016/S2666-5247(23)00375-0)

## Documentation

Operations and subsystem documentation is published at
**[esteinig.github.io/cerebro](https://esteinig.github.io/cerebro/)**. To build it
locally, see [`DOCS.md`](DOCS.md).

## License

See [`LICENSE`](LICENSE).
