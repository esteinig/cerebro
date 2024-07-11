# Cerebro

Metagenomic diagnostic pipelines and collaborative reporting for pathogen detection, species identification, and differential host genome analysis in clinical and public health production. 

🩸 **Metagenomic diagnostic core functions** 

- Multi-classifier taxonomic profiling, metagenome assembly and alignment in Nextflow pipelines
- Optimized pangenome host depletion and background depletion with [`Scrubby`]()
- Viral infections, pan-viral enrichment protocols and syndrome specific subtyping panels using [`Vircov`]()

- Differential host tumor DNA diagnostics using segmental CNV detection and methylation classifiers ([`Sturgeon`]())
- Species identification pipelines with [`GTDB`]() for prokaryotic ONT/Illumina reference level genomes
- MAG recovery from enriched culture and sample co-assembly, unclassified viral bin prediction ([`geNomad`, `RdRP`]())
- Custom database and index construction, grafted taxonomies, genome cleaning and syndromic diversity injection with [`Cipher`]()

 📰 **Collaborative clinical reporting (Bug Board)**

- [Collaborative and auditable pathogen determination]() from metagenome sequencing results
- Multi-tenant Svelte application and API with secure local or web-server deployment configs
- Scalable application stack deployment integrated into the primary command-line interface ([Cerebro CLI]()) 
- Clinical reporting with [`Typst`]() formatted templates linked into the database of evidence from multi-classifier/databases
- Secure [`wasm` enabled report generation]() in-browser for sensitive reports, interactive data visualizations
- Auditable team member comments and results discussion for expert panel reviews of data ([online "Bug Board"]())

🏥 **Clinical and public health production operations** 

- Simulations using in silico syndromic-specific reference panels for ONT/Illumina signal-level and read-level data with [`Cipher`]()
- Evaluation simulation and patient datasets for continous integration of quality assurance with [`Cerebro`]()
- Background/sample site/kitome contamination issues in general clinical or public health environments
- Distributed sequence and analysis storage, file system and data retention policies, cloud storage support etc. through [`SeaweedFS`]()
- Experimental protocols for reference labs for optimisation of the [UMI-adapter DNA/RNA protocol]() for low abundance clinical sampel types

## Getting started

Let's step through some common tasks and core functions of `Cerebro` and its Docker stack. 

### Nextflow pipelines 

❗ You do not need the `Docker` stack if you want to run only the Nextflow pipelines, or the command-line interface binary for data manipulation and core utility tasks around pipeline outputs. 

Nextflow pipelines are purposefully decoupled from the API/DB/App so that data can be analysed and summarized into the database models within a secure network (i.e. where primary samples are sequenced). These models can then be transfered outside the secure network e.g. by uploading to a web-deployment of the Docker stack, as it does not contain potentially identifiable information like incidental host sequences which may be abundant in some clinical sample types.


## Status

Under active development for production release. Not recommended for deployment at this stage. 

This is a preliminary public release of code for the viral enrichment branch of the pipeline used in:

> Michael A Moso, George Taiaroa, Eike Steinig, Madiyar Zhanduisenov, Grace Butel-Simoes, Ivana Savic, Mona L Taouk, Socheata Chea, Jean Moselen, Jacinta O’Keefe, Jacqueline Prestedge, Georgina L Pollock, Mohammad Khan, Katherine Soloczynskyj, Janath Fernando, Genevieve E Martin, Leon Caly, Ian G Barr, Thomas Tran, Julian Druce, Chuan K Lim, Deborah A Williamson - **Non-SARS-CoV-2 respiratory viral detection and whole genome sequencing from COVID-19 rapid antigen test devices: a laboratory evaluation study** - Lancet Microbe (2024) -[10.1016/S2666-5247(23)00375-0](https://doi.org/10.1016/S2666-5247(23)00375-0)

# Table of contents

1. Background
2. Install
3. CSF protocol 
4. Nextflow pipelines and core functions
5. The Rust command-line client and stack 
6. Docker stack deployments with web application and Rust API
7. Svelte web application for collaborative reporting

