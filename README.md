# Cerebro

Metagenomic diagnostics pipeline and collaborative reporting stack for pathogen detection, species identification, host genome analysis, quality assurance and deployment in clinical and public health production environments.

<details>
<summary>🩸 Metagenomic diagnostic core functions </summary>
<br>
 
Main

- Multi-classifier taxonomic profiling, metagenome assembly and alignment in Nextflow pipelines
- Optimized pangenome host depletion and background depletion with [`Scrubby`]() and `Metabuli`/`Strobealign`
- Viral infections, pan-viral enrichment protocols and syndrome-specific subtyping panels using [`Vircov`]()
- Differential host tumor DNA diagnostics using segmental CNV detection and methylation classifiers ([`Sturgeon`]())

Support
 
- Species identification pipelines with [`GTDB`]() for prokaryotic ONT/Illumina reference level genomes
- MAG recovery from enriched culture and sample co-assembly, unclassified viral bin prediction ([`geNomad`, `RdRP`]())
- Custom database and index construction, grafted taxonomies, genome cleaning and syndromic diversity injection with [`Cipher`]()

</details>

<details>
<summary>📰 Collaborative clinical reporting (Bug Board) </summary>
<br>
 
- [Collaborative and auditable pathogen determination]() from metagenome sequencing results
- Multi-tenant Svelte application and API with secure local or web-server deployment configs
- Scalable application stack deployment with different data security and collaboration models
- Stack configuration and deployment integrated into the primary command-line interface ([Cerebro CLI]()) 
- Clinical reporting with [`Typst`]() formatted templates linked into the database of evidence from multi-classifier/databases
- Secure [`wasm` enabled report generation]() in-browser for sensitive reports, interactive data visualizations
- Auditable team member comments and results discussion for expert panel reviews of data ([online "Bug Board"]())

</details>

<details>
<summary>🏥 Clinical and public health production environments </summary>
<br>
 
- Simulations using in silico syndromic reference panels for ONT/Illumina signal-level and read-level data with [`Cipher`]()
- Evaluation of simulation and patient datasets for continous integration of quality assurance with [`Cipher`]() and [`Cerebro`]()
- Background/sample site/kitome contamination issues in general clinical or public health environments via the Cerebro API
- Distributed sequence and analysis storage, file system and data retention policies, cloud storage etc. through [`SeaweedFS`]() integration
- [Standard operating procedures]() for continous operation of `Cerebro` as a service for clinical diagnostic reporting
- Experimental protocols for reference labs for optimisation of the [UMI-adapter DNA/RNA protocol]() for low abundance clinical sample types

</details>

## Getting started

Let's step through some common tasks and core functions of `Cerebro` and its data application and reporting stack. This section provides some examples of how to get started quickly with `Cerebro`. For more details and how to deploy and operate the full application in production please see the [documentation](). 

Minimum requirements:

* Linux OS
* Nextflow v2024.04
* Conda/Mamba/Docker

Computational resource requirements are variable and range from a standard laptop for the application stack to full nation-wide server infrastructure for pipelines and web-application (if you were so inclined). This is because the application stack for data and reporting can be deployed with various [infrastructure, data security and collaboration models]() in mind and depends on the number of laboratories, collaborators, sequencing throughput, data storage and many other considerations.

> [!NOTE]
You do not need the `Docker` stack for core metagenome diagnostic pipelines and report generation - you can run the [Nextflow pipelines]() separately and use the [`Cerebro CLI`](#command-line-client) for data manipulation, processing of pipeline outputs and clinical report generation.

## Nextflow pipeline 

### Quick start

```

```

## Cerebro CLI

### Quick start

```

```


## Clinical reporting

### Quick start

```

```


## Application stack

### Quick start

```

```


## Cerebro API and FS

### Quick start

```

```

## Status

Under active development for production release. Not recommended for deployment at this stage. 

This is a preliminary public release of code for the viral enrichment branch of the pipeline used in:

> Michael A Moso, George Taiaroa, Eike Steinig, Madiyar Zhanduisenov, Grace Butel-Simoes, Ivana Savic, Mona L Taouk, Socheata Chea, Jean Moselen, Jacinta O’Keefe, Jacqueline Prestedge, Georgina L Pollock, Mohammad Khan, Katherine Soloczynskyj, Janath Fernando, Genevieve E Martin, Leon Caly, Ian G Barr, Thomas Tran, Julian Druce, Chuan K Lim, Deborah A Williamson - **Non-SARS-CoV-2 respiratory viral detection and whole genome sequencing from COVID-19 rapid antigen test devices: a laboratory evaluation study** - Lancet Microbe (2024) -[10.1016/S2666-5247(23)00375-0](https://doi.org/10.1016/S2666-5247(23)00375-0)

