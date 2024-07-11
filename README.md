# Cerebro

Metagenomic diagnostic pipelines and collaborative reporting for pathogen detection, species identification, and differential host genome analysis in clinical and public health production. 

## Overview

Cerebro is being developed primarily for neurological infections using a [UMI-barcoded deep Illumina DNA/RNA short-read sequencing protocol](). Its implementation as outlined below has focused on the particular challenges of low abundance infections in sterile fluid high host background sequence content (cerebrospinal fluid, ocular fluid, brain tissue) with a wide variety of applications that are may be of use for anyone operating clinical metagenomic sequencing,

🩸 **Metagenomic diagnostic core functions** 

- Viral infections in general and pan-viral enrichment protocols, with automatic consensus genome recovery and syndrome-representative subtyping using `Vircov`
- Optimized pangenome host depletion and background depletion with `Scrubby`, differential host diagnostic arm of the main Nextflow pipeline
- Segmental copy-number variation for tumor DNA detection with Illumina and methylation profile prediction with the `STURGEON` classfier with ONT 
- Multi-classifier taxonomic profiling, metagenome assembly and alignment approaches in the `Cerebro` Nextflow pipelines
- Species identification pipelines with `GTDB` for prokaryotic ONT/Illumina reference level genomes and MAG recovery from enriched culture
- Custom metagenomic database and index construction with `Cipher` with auditable database updates and taxonomy changes, genome cleaning and diversity injection

 📰 **Clinical collaborative reporting**

- Collaborative and auditable decision making on pathogen determination from metagenome samples with multi-tenant Svelte application
- Scalable application stack deployment integrated into the primary command-line interface (Cerebro CLI) with preconfigured deployment profiles for local
- Clinical reporting linked into the database of evidence from multi-classifier/databases and to team member comments / discussion (online "Bug Board")
- Low-level background/sample site/kitome contamination issues in general clinical or public health environments

🏥 **Clinical production operations** 

- Simulations using in silico syndromic-specific reference panels for ONT/Illumina signal-level and read-level data with  `Cipher`
- Evaluation simulation and patient datasets for continous integration of quality assurance with `Cerebro`
- Distributed primary sequence data and pipeline output storage, file system and data retention policies, cloud storage support etc. through `SeaweedFS`
- Experimental protocols for reference labs for optimisation of the UMI-adapter DNA/RNA protocol we use for CNS infections

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

