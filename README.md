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

Let's step through some common tasks and core functions of `Cerebro` and the application and reporting stack. This section provides some examples of how to get started quickly with `Cerebro`. For more details and how to deploy and operate the full application in production please see the [documentation](). 

Minimum requirements:

* Linux OS
* Nextflow v24
* Conda/Mamba/Docker

Computational resource requirements are variable and range from a standard laptop for the application stack to full nation-wide server infrastructure for pipelines and web-application (if you were so inclined). This is because the application stack for data and reporting can be deployed with various [infrastructure, data security and collaboration models]() in mind and depends on the number of laboratories, collaborators, sequencing throughput, data storage and many other considerations.

> [!NOTE]
You do not need the `Docker` stack for core metagenome diagnostic pipelines and report generation - you can run the [Nextflow pipelines]() separately and use the [`Cerebro CLI`](#command-line-client) for data manipulation, processing of pipeline outputs and clinical report generation.

## Nextflow pipeline 

### Quick start

Pathogen detection with PE Illumina reads from metagenomic sequencing of sterile-site samples (validated for ocular fluid and cerebrospinal fluid):

```bash
nextflow run -r v1.0.0 https://github.com/esteinig/cerebro \
  -profile dgx,large,mamba,cns,cipher \
  -entry pathogen \
  --outputDirectory outputTest/ \
  --databaseDirectory db/ \
  --fastqPaired 'fastq/*_{R1_001,R2_001}.fastq.gz'
```

This will run the default quality control, taxonomic profiling and metagenome assembly configuration for pathogen identification in low microbial biomass sample types (such as ocular fluids or CSF) where distinction from contamination and incidental background organisms is the main challenge for diagnostics (`needle-in-a-haystack`). 

> [!WARNING]
This pipeline is not suitable for high microbial biomass sample types (such as respiratory or environmental samples) where a diverse and abundant background microbiome is the main challenge for diagnostics (`haystack-full-of-needles`).

For the default configuration you will need the `Cipher` diagnostic database in the `--databaseDirectory`, which is an amalgamation of archaeal/bacterial (GTDB), eukaryotic (EuPath, Wormbase) and viral (ICTV) reference genome collections and taxonomies. Depending on sequencing protocol, you will also need a human reference index for alignment-based depletion and additional background sequences (like synthetic spike-ins or internal phage controls) for depletion in the quality control module.

Download the `Cipher v2` database:

```

```

If you are using the associated short-read sequencing protocol also download the following files:

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


## Cerebro API

### Quick start

```

```

## Cerebro FS

### Quick start

```

```

## Databases and taxonomy

### Quick start

```

```


## Production

With the stack running (see below) production mode will upload outputs directly to the team database and 

```
{sample_id}__{tag1}__{tag2}__{tag3}_R1_001.fastq.gz
```

The most important tags are (forward read examples, matches reverse read file names):

* `DNA` or `RNA` for example: `DW-63-V01__DNA_R1_001.fastq.gz` and a matching `DW-63-V01__RNA_R1_001.fastq.gz` read file
* `NTC` and `ENV` for negative template and environmental controls for example: `DW-63-V420__DNA__NTC_R1_001.fastq.gz` for the DNA negative template control library
* `POS` for a positive control mock sample for extraction and sequencing controls for example `DW-63-V07__DNA__NTC_R1_001.fastq.gz`  for the DNA positive control library

Any other tag outside of the above reserved ones can also be added but has not specific functions in the stack.

## Status

Under active development for production release. Not recommended for deployment at this stage. 

This is a preliminary public release of code for the viral enrichment branch of the pipeline used in:

> Michael A Moso, George Taiaroa, Eike Steinig, Madiyar Zhanduisenov, Grace Butel-Simoes, Ivana Savic, Mona L Taouk, Socheata Chea, Jean Moselen, Jacinta O’Keefe, Jacqueline Prestedge, Georgina L Pollock, Mohammad Khan, Katherine Soloczynskyj, Janath Fernando, Genevieve E Martin, Leon Caly, Ian G Barr, Thomas Tran, Julian Druce, Chuan K Lim, Deborah A Williamson - **Non-SARS-CoV-2 respiratory viral detection and whole genome sequencing from COVID-19 rapid antigen test devices: a laboratory evaluation study** - Lancet Microbe (2024) -[10.1016/S2666-5247(23)00375-0](https://doi.org/10.1016/S2666-5247(23)00375-0)

## Supplementary Pipelines

### Aneuploidy from host background

Dependencies:

* Conda/Mamba
* Nextflow v24

Large segmental copy number variation (CNV) detection across human chromosomes from host background of short-read metagenomic sequencing data with `CNVkit` ([Talevich et al. 2016](https://doi.org/10.1371/journal.pcbi.1004873)).

This is a supplementary pipeline not intended for production - to execute please clone the repository first:

```
git clone https://github.com/esteinig/cerebro
```

Execute with relevant parameters:

```
nextflow run ./cerebro/lib/standalone/aneuploidy/main.nf -profile mamba,medium \
  --pairedReads '*_{R1_001,R2_001}.fastq.gz' \
  --outdir test_aneuploidy \
  --referenceFasta resources/CHM13v2.fasta \
  --normalControlBam resources/HG007.5x.QC.bam \
  --resources.threads.minimap2 32
```

Other parameters can be found in the [`nextflow.config`](https://github.com/esteinig/cerebro/blob/main/lib/standalone/aneuploidy/nextflow.config). Resource dependencies are the CHM13v2 human reference genome and a sub-sampled (5x) reference alignment of [HG007 (ChineseTrio, mother)](https://github.com/genome-in-a-bottle/giab_data_indexes). We tested this default configuration on Detroit cell-lines which derive from a female pharyngeal cancer patient and show strong patterns of segmental aneuploidy across chromosomes when compared to a known healthy patient sample.


