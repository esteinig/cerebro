
<!-- <p align="left" style="margin-top: 2rem; margin-bottom: 2rem">
    <img src="assets/logo.png#only-light"  width="300" height="200">
    <img src="assets/logo.png#only-dark" width="300" height="200">
</p> -->

`Cerebro` is a metagenomics platform for clinical diagnostics and reporting in production environments. 

`Cerebro` is broadly applicable for short-read (for now) metagenomics and -transcriptomics. It is optimized for sample types with low microbial biomass (for now).

`Cerebro` has been primarily tested on central nervous system infections and cerebrospinal fluid samples in combination with our [clinical assay protocol](). Viral detection and consensus genome assembly modules have been validated for panviral enrichment using clinical samples extracted from [rapid antigen test devices]().


## Components

We provide a collection of integrated tools that aim to ease taxonomic classification, data management, interpretation and reporting in clinical production environments. 

* Metagenomic and -transcriptomic [pathogen detection pipeline]()
* Database server and application programming interface (API) 
* [Command-line application]() for pipeline integration, quality assurance and API
* Interface for (multi-user) exploration of taxonomic classifications and clinical reports


!!! info

    Pipeline and command-line utilities are independent of the platform. Server and interface setups are slightly more involved, but are recommended for ease of interpretation and reporting.


## Documentation

* [Setup](index.md)  
Installation and configuration of `Cerebro`

* [Pipeline](pipeline.md)  
Pipeline configuration and execution 

* [Server](server.md)  
Server configuration and API specification

* [Interface](interface.md)  
User interface configuration and navigation 

* [Reports](reports.md)  
Clinical report compilation with `Cerebro`

* [Security](security.md)  
Data privacy and security in `Cerebro`

* [Deployment](deployment.md)  
Deploying or hosting the `Cerebro` stack

* [Dependencies](dependencies.md)  
Description and citations of tools used in `Cerebro`

## Tutorials



## Data Privacy

Metagenomic sequencing inevitably captures host nucleic acid. We strongly care about the ethical and privacy implications of this kind of data. 

`Cerebro` aims to:

* Retain patient sequence data on infrastructure where it was generated
* Decouple sequence analysis from results accessible in the database and interface
* Produce clinical reports without sending sensitive information over the network
* Enable quality assurance for accredited production services

We have concerns about data privacy using well-advertised, philanthropically- or commerically-funded metagenomics platforms where patient sequence- and meta-data are uploaded to servers that are not under our control. 

Our priority was developing pipelines and reporting tools that allow users to retain complete control over source code, deployment and operation on their own infrastructure. We hope that this conceptualisation is useful for clinical metagenomic services where there may be concerns about how sequence data is handled and evaluated by third parties.


## Data Security

!!! info

    You can read more about data security and privacy considerations in the [security](#security.md) section.


`Cerebro` de-couples pipeline execution from the analytical stack, and does not store sequence or patient data in the database. Once the pipeline outputs are generated, they are generally safe to be transmitted on a network. However, there are residual risks for accidental (or purposeful) data linkage, which are described in the [data security]() section of the deployment guide. Special consideration should be given to how sensitive patient information can be included in the [clinical report headers]().

Note that in the current version access to the database and interface are implemented through scoped token authentication with logging and rate limiting on critical endpoints, verification and password-reset emails for user management and traffic encryption through a reverse-proxy (SSL/TLS). All user and database management runs through the administrator of the application who controls the container stack. 

One major reason for the current implementation is one of scale and purpose. We primarily deploy the database, server and interface on a secure local network to which a limited number of users from the production team have access. We also routinely spin up the stack on a laptop to quickly view and interpret workflow results after transfer of pipeline outputs through SSH over a VPN.

!!! warning

    While user authentication follows the `Authorization Code Flow` defined in [OAuth 2.0 RFC 6749, section 4.1](https://datatracker.ietf.org/doc/html/rfc6749#section-4.1) it is **not** a replacement for certified authentication/authorization flows. Do not deploy `Cerebro` on the open internet unless you are familiar with the implications for data security.


## Schematic

Pipeline schematic:



Platform schematic:


## Contributors


`Cerebro` is being developed as part of an open-source protocol for central-nervous system infections in the state of Victoria, Australia. Our work contributes to the `META-GP` program funded by the Medical Research Future Fund (MRFF) - a collaboration between Australia’s leading researchers in pathogen genomics and clinical infectious diseases. 

`META-GP` central nervous system infection team led by [Prof. Deborah Williamson]() at The Peter Doherty Institute for Infection and Immunity:

* Prashanth Ramachandran
* Marcelina Krysiak
* Janath Fernando
* Eike Steinig 
* Chhay Lim
* Kirti Deo

Victorian Infectious Diseases Reference Laboratory:

* Chuan Lim
* Leon Caly
* Jean Moselen

Royal Melbourne Hospital (RMH):

* Katherine Bond

## Citations

If you use `Cerebro` for your research please cite:

> Steinig (2023) - Clinical metagenomic diagnostics and reporting of central nervous system infections with Cerebro

## Contact

[Eike Steinig](https://github.com/esteinig)