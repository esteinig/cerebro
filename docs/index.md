# Cerebro

A metagenomic diagnostics pipeline and collaborative reporting stack for pathogen
detection, species identification, host-genome analysis, and quality assurance —
designed for deployment in clinical and public-health production environments.

Cerebro turns metagenomic sequencing data into auditable, clinically reportable
pathogen-detection results, and provides the application stack, storage, and
governance needed to run that process as a service in a diagnostic laboratory.

## What Cerebro does

**Metagenomic diagnostics.** Multi-classifier taxonomic profiling, metagenome assembly,
and alignment in Nextflow; optimised host and background depletion; viral enrichment
and syndrome-specific subtyping; and supporting workflows for species identification,
MAG recovery, and custom database construction.

**Collaborative clinical reporting.** Auditable, multi-tenant pathogen determination
from sequencing results, with a web application and API, configurable data-security and
collaboration models, templated clinical reports linked to the underlying evidence, and
team review of results ("Bug Board").

**Clinical & public-health production.** In-silico syndromic simulation panels for
quality assurance, contamination assessment, distributed sequence and analysis storage
with retention policies, and standard operating procedures for continuous operation as
a diagnostic service.

## Getting started

!!! note
    The full Docker application stack is **not** required for the core metagenomic
    pipelines and report generation. The Nextflow pipelines can be run on their own and
    the Cerebro CLI used for processing pipeline outputs and generating clinical
    reports.

**Minimum requirements**

- Linux
- Nextflow v24
- Conda / Mamba / Docker

Computational requirements range from a single workstation for the application stack to
large server infrastructure for pipelines and the web application, depending on the
number of laboratories, sequencing throughput, storage, and the chosen infrastructure,
data-security, and collaboration model.

**Pathogen detection — quick start**

Paired-end Illumina reads from metagenomic sequencing of sterile-site samples
(validated for ocular fluid and cerebrospinal fluid):

```bash
nextflow run -r v1.0.0 https://github.com/esteinig/cerebro \
  -profile dgx,large,mamba,cns,cipher \
  -entry pathogen \
  --outputDirectory outputTest/ \
  --databaseDirectory db/ \
  --fastqPaired 'fastq/*_{R1_001,R2_001}.fastq.gz'
```

## Subsystems

Cerebro's production capabilities are organised into subsystems, each with its own
operations documentation and governance (assurance and production-readiness).

**[cerebro-fs](cerebro-fs/index.md)** — the storage, lifecycle, and durability layer.
It manages sequencing data and derived files from ingestion through tiering, archival,
and retention, while guaranteeing that every object remains intact, durable,
recoverable, and auditable for as long as policy requires. See its
[overview](cerebro-fs/overview.md), [architecture](cerebro-fs/architecture.md),
operations runbooks, and [production-readiness governance](cerebro-fs/production-readiness.md).

## Status

!!! warning "Under active development for production release"
    Cerebro is under active development. Subsystem documentation states, where relevant,
    the boundary between what is validated and what a given deployment must confirm in
    its own environment.

This documentation accompanies the viral-enrichment branch of the pipeline used in
Moso et al., *Non-SARS-CoV-2 respiratory viral detection and whole genome sequencing
from COVID-19 rapid antigen test devices: a laboratory evaluation study*, **The Lancet
Microbe** (2024), [10.1016/S2666-5247(23)00375-0](https://doi.org/10.1016/S2666-5247(23)00375-0).
