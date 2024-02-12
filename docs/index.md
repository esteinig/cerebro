
`Cerebro` is a metagenomics platform for clinical diagnostics and reporting in production environments. 

`Cerebro` uses an ensemble-like approach with multiple classification and profiling strategies, tools and databases targeting **pathogen detection**.

## Application

Our primary workflow is optimized for [pathogen detection in low microbial biomass samples]() using short-read deep sequencing assays including our open-source protocol.

`Cerebro` is **not** designed for complex community analysis or methods that require high microbial biomass samples, such as recovering metagenome-assembled genomes and strain typing.

!!! info
    
    `Cerebro` can be re-configured for other applications using `Nextflow` configuration profiles to recover [viral consensus genomes and genotypes](), perform [metagenome co-assembly]() to assess background contamination or run simple best-practice implementations like the pathogen detection [protocol for Kraken]().

## Components

We provide a collection of integrated modules that aim to ease taxonomic profiling, data management, interpretation and reporting in clinical production environments. 

* `Nextflow`: metagenomic- and transcriptomic pathogen detection pipeline 
* `Rust`: command-line application for pipeline integration and quality assurance
* `Rust`: database server and application programming interface for data management 
* `Svelte`: web application for collaborative exploration and interpretation of results 
* `Rust`: templating engine for data integration in auditable clinical reports

!!! info

    Pipeline and command-line utilities are independent of the full-stack application. Server and web-interface setups are recommended for ease of interpretation and reporting.

## Documentation

* [Setup](index.md)  
Installation and configuration of `Cerebro`

* [Pipeline](pipeline.md)  
`Cerebro` pipeline configuration and execution 

* [Deployment](deployment.md)  
Deploying the `Cerebro` application stack

* [Interface](interface.md)  
`Cerebro` interface configuration and navigation 

* [Reports](reports.md)  
Clinical report compilation with `Cerebro`

* [Security](security.md)  
Data privacy and security in `Cerebro`

* [Development](development.md)  
Setting up development environments for `Cerebro`

* [Dependencies](dependencies.md)  
Description and citations of tools used by `Cerebro`


### Tutorials

* [Taxonomic profiling for clinical applications](dependencies.md)  
How to setup, execute and interpret outputs from taxonomic profiling

* [Panviral enrichment and consensus assembly](dependencies.md)  
Viral detection and consensus genome recovery from panviral enrichment assays

* [Data interpretation and reporting for diagnostics](dependencies.md)  
How to interpret the data outputs and generate reports in the user interface

* [Local production deployment as diagnostic service](dependencies.md)  
Deploy one or multiple stacks in a production environment

* [Stack administration for multiple users and projects](dependencies.md)  
Learn how to manage teams, databases and data collections in production


### Schematics

Pipeline schematic:



Platform schematic:

## Considerations

### Data privacy

Metagenomic sequencing inevitably captures host nucleic acid. We have concerns about data privacy using well-advertised, philanthropically- or commerically-funded metagenomics platforms where sequence- and meta-data are uploaded to servers that are under third party control and not accountable to the public. We strongly care about the ethical and privacy implications of this kind of data. 

`Cerebro` aims to:

* Retain patient sequence data withing local networks where it was generated
* Decouple sequence analysis from diagnostic results accessible in the database and interface
* Produce clinical reports without sending sensitive information over the network
* Enable quality assurance and routine operation as clinical or public health service
* Provide all tools necessary for accreditation and validation forproduction environments

Our priority was developing pipelines and reporting tools that allow users to retain complete control over source code, deployment and operation on their own infrastructure. We hope that this conceptualisation is useful for clinical metagenomic services where there may be concerns about how sequence data is handled and evaluated by third parties.


### Data security

!!! info

    You can read more about data security and privacy considerations in the [security](#security.md) section.


`Cerebro` de-couples pipeline execution from the analytical stack, and does not store sequence or patient (meta-)data in the database. Once the pipeline outputs are generated, they are generally safe to be transmitted on a network. 

However, there are residual risks for accidental (or purposeful) data linkage, which are described in the [data security]() section of the deployment guide. Special consideration should be given to how sensitive patient information can be included in the [clinical report headers]().


## Support 

### Contributors


`Cerebro` is being developed as part of an open-source protocol for central-nervous system infections in the state of Victoria, Australia. Our work contributes to the [`META-GP`]() program funded by the Medical Research Future Fund (MRFF). 

`META-GP` central nervous system infection team at [The Peter Doherty Institute for Infection and Immunity]():

* Deborah Williamson
* Tim Stinear
* Prashanth Ramachandran
* Eike Steinig 
* Chhay Lim
* Marcelina Krysiak
* Janath Fernando
* George Taiaroa
* Kirti Deo

Translational diagnostics group at the [Victorian Infectious Diseases Reference Laboratory]():

* Chuan Lim
* Leon Caly
* Jean Moselen
* Ammar Aziz

Microbiology department at the [Royal Melbourne Hospital]():

* Katherine Bond
* Michael Moso

### Publications

Viral detection and genome recovery from rapid antigen devices using `Cerebro`:

> Moso, Taiaroa, Steinig et al. (2023) - Non-SARS-CoV-2 respiratory viral detection and whole genome sequencing from COVID-19 rapid antigen test devices: A laboratory evaluation study - The Lancet Microbe

### Citation

If you use `Cerebro` for your research, please cite:

> Steinig et al. (2024) - Clinical metagenomic diagnostics and reporting in production environments


If you use the panviral enrichment configuration, please also cite:

> Moso, Taiaroa, Steinig et al. (2023) - Non-SARS-CoV-2 respiratory viral detection and whole genome sequencing from COVID-19 rapid antigen test devices: A laboratory evaluation study - The Lancet Microbe

### Contact

If you have a feature request or find a bug, please feel free to [open an issue](https://github.com/esteinig/cerebro) in the repository.

If you would like to contribute to `Cerebro` development (you're awesome) please read the [Development](#development.md) section of the documentation and get into contact through the repository.

If you would like to get in contact to discuss extensions to `Cerebro` or local production setups for public health or clinical applications please contact [Prashanth Ramachandran]() or [Eike Steinig](https://github.com/esteinig).