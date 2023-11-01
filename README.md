# Cerebro

Central nervous system specific metagenomic diagnostics workflow for protocols at the Victorian Infectious Diseases Reference Laboratory (VIDRL) and the CNS Meta-GP group (see below) at the Peter Doherty Institute. 

Table of contents

1. [Features](#features)
2. [Basic usage](#basic-usage)
3. [Nextflow pipeline](#cerebro-pipeline)
4. [Cerebro UI and API stack](#cerebro-stack)
5. [Command-line interface](#cerebro-client)
6. [Deployment scenarios](#)
7. [Contributions]()

## Features

- `Nextflow` pipeline for meta-metagenomics taxonomic profiling targeting quantitative pathogen detection using Illumina/Nanopore data offering in depth quality control taxonomic classifier options. Includes custom tools for the classification pipeline focusing on human sequence removal (`scrubby`) and viral alignment filters and auto-mated reference selection for coverage remapping subworkflows (`vircov`)

- `Rust` library to run the database and user interace API server including authenticated routes and authorization middleware with `ActixWeb`; providing endpoints for workflow output processing (quality control, classification, workflow and sample configurations) and access to structured storage. Includes command-line client to run supplementary tools in the main pipeline support production execution of the pipeline (create sample sheets, run automated input framework for production, ...)

- `Docker` and `docker-compose` configs templated from the `cerebro` binary for deployment of a `Svelte UI` with `Skeleton UI` component library + `MongoDB` and `Redis` databases to manage, interpret and collaboratively decide on candidate pathogens and generate auditable clinical reports with the `Tectonic` LaTeX to PDF engine and report template.

- Provisioning of traceable and versioned metagenomic reference databases and quality or diagnostic genome feature masks; continuous integration testing of full diagnostic workflow against clinical samples with orthogonal diagnostic data or synthetic simulated datasets (`cipher`)

## Basic usage

Requirements:
- `Linux OS`

`Cerebro` consists of two parts: 

1. `Nextflow` pipeline for sequence data (Illumina PE) ingestion and processing, including integration testing and assay validation workflows.
2. `Rust` and `docker-compose` application stack (`MongoDB` + `Redis` + `Sveltekit` + `Rust` library that implements the `ActixWeb` API and pipeline integration via command-line interface).

The `cerebro` command-line interface is available as cross-compiled binary release (Linux, macOS, Windows not tested). It exposes structured subcommands for interaction with the application stack and pipeline (binary assumed to be on user `$PATH`)

```bash
# get cerebro cli binary to support cerebro workflows
curl https://github.com/esteinig/cerebro/releases/download/latest/cerebro-latest-Linux_x86_64.tar.xz -o - | tar -xzO > cerebro

# assume cerebro cli on user path
cerebro --help
```

## `cerebro` pipeline

Requirements for local execution: 

- `Conda` or `Mamba`
- `Nextflow >= v23.10.04`

```bash
# pull latest from github, show help menu, use mamba envs
nextflow run -r latest esteinig/cerebro -profile mamba --help

# provision with latest cipher kmer db build, may take some time
nextflow run -r latest esteinig/cerebro -profile mamba -entry cipher \
    --revision latest \
    --representation full \
    --outdir cipher_db/

# default qc and tax profile on input read dir (pe illumina)
nextflow run esteinig/cerebro -r latest -profile mamba \
    --fastq "fastq/*_{R1,R2}.fq.gz" \
    --databases cipher_db/

# kmer tax profile on input read dir (pe illumina)
nextflow run esteinig/cerebro -r latest -profile mamba,kmer \
    --fastq "fastq/*_{R1,R2}.fq.gz" \   
    --databases cipher_db/

# production: cerebro client to create input sample sheet
cerebro pipeline sample-sheet 
    --input fastq/ \
    --output sample_sheet.csv \
    --run-id production_test \
    --glob "*_{R1,R2}.fq.gz"

# production: tax profile on sample sheet input (pe illumina)
nextflow run esteinig/cerebro -r latest -profile mamba \
    --production true \
    --sample_sheet sample_sheet.csv \
    --databases cipher_db/

# production: db upload on successful completion 
nextflow run esteinig/cerebro -r latest -profile mamba \
    --production true \
    --sample_sheet sample_sheet.csv \    
    --databases cipher_db/ \
    --cerebro.api.enabled true \
    --cerebro.api.url $CEREBRO_API_URL \ 
    --cerebro.api.token $CEREBRO_API_TOKEN \
    --cerebro.api.upload true
```

#### Nextflow command-line configurations 

```bash
# runtime profiles
nextflow run -r latest esteinig/cerebro \
    # local conda,mamba env created for each process
    -profile conda,mamba
    # kmer, alignment, assembly classifier modules
    -profile kmer,alignment,assembly
    # process-configured resource profiles
    -profile small,medium,large,galactic
    # specific protocol configs
    -profile cns@v1                                # central nervous system meta-gp protocol
    -profile panviral@v1                           # panviral enrichment protocol rat sequencing
    -profile kraken@v1                             # kraken2 nature protocols pathogen detection 
    -profile aneuploidy@v1                         # wgs copy number variation for consented patients
    # taxonomy + database provisions
    -profile cipher-db@v1
    # workflow integration tests
    -profile cipher-cns@v1
    # workflow dev
    -profile io,dev

# nested module params in `nextflow.config`
nextflow run -r latest esteinig/cerebro -profile mamba \

    # input/output and workflow provision
    --db "cipher_db/" \
    # pe read input
    --fastq "fastq/*_{R1,R2}.fq.gz" \
    # production input
    --production \
    --sample_sheet "sample_sheet.csv" \
    # output directory
    --outdir "test_run" \

    # qc read processing module
    --qc.enabled \
    --qc.deduplication.enabled \
    --qc.deduplication.method "umi-naive" \
    --qc.reads.fastp.enabled \
    --qc.controls.ercc.enabled \
    --qc.controls.phage.enabled \
    --qc.host.depletion.enabled \
    --qc.background.mask.enabled \

    # taxon profiling with references and taxonomy files in --db
    --taxa.enabled \
    --taxa.kmer.enabled \
    --taxa.kmer.kraken2uniq.enabled \
    --taxa.kmer.kraken2bracken.enabled \
    --taxa.kmer.kmpc.enabled \
    --taxa.kmer.metabuli.enabled \
    --taxa.kmer.postalign.enabled \               # k-mer id reads align v db genome
    --taxa.alignment.enabled \
    --taxa.alignment.subset.enabled \             # mash screen tax pre-subset + db re-index 
    --taxa.alignment.minimap2.enabled true \      
    --taxa.alignment.strobealign.enabled false \  # one of
    --taxa.alignment.bowtie2.enabled false \      
    --taxa.assembly.enabled \                     # metaspades + align lca - ncbi nt/nr
    --taxa.assembly.blastn.enabled \
    --taxa.assembly.diamond.enabled \

    # api: requires --sample_sheet and --production
    --cerebro.api.enabled true \  
    --cerebro.api.url $CEREBRO_API_URL \ 
    --cerebro.api.token $CEREBRO_API_TOKEN \

    # api: live workflow status and sample tracking
    --cerebro.api.run.status.enabled \                # status report logging of pipeline updates
    --cerebro.api.run.status.slack.enabled \          # status report logging to slack channel
    --cerebro.api.run.report.enabled \                # create run report for sample tracking + qc
    --cerebro.api.run.report.slack.enabled \          # post run report to slack channel

    # api: upload to team collection
    --cerebro.api.upload.enabled \
    --cerebro.api.upload.team "VIDRL" \
    --cerebro.api.upload.database "PRODUCTION" \
    --cerebro.api.upload.collection "MGP-CNS-20231012"

    # resource labels for each process
    --resources.minimap2_align.cpus 32 \
    --resources.minimap2_align.memory "32 GB" \
    --resources.minimap2_align.conda "envs/minimap2.replacement.yml" \
    --resources.minimap2_align.container "biocontainers/minimap2:latest"
```

#### Example pipeline configurations

Example commands for executing the `Nextflow` pipeline and processing local output with the `cerebro` command-line interface:

1. Reference data and index provision for pipeline setup ([`Cipher`]())
2. [Kraken pathogen detection]() protocol and post-processing
3. Multi-classifer, multi-database taxonomic classification ([meta-metagenomics protocol v1]())
4. Taxonomic profiling for [central nervous system assay v1]() (META-GP)
5. [Viral genome recovery for panviral enrichment]() from rapid antigen tests (RAT)
6. `Cerebro` [pipeline integration testing]() with `Cipher` reference datasets and actions
7. Testing impact of host read depletion on classification performance with `Cipher` and `Scrubby`

## `Cerebro` stack

The `Cerebro` stack is a `docker-compose` file that runs:

* `Cerebro` API server and the `MongoDB` database holding the user models, stack admin models and team access structured databases and project collections of [workflow output models]()
* `Redis` databases for registration and session authentication of token authentication middleware in the API for access via `Sveltekit` user interface with `Skeleton UI` components and styling
* `Cerebro` command-line client that support data and user management functions, and allows for automation of routine sequencing and workflow runs with auditable diagnostic reporting in clinical production environments (e.g. `Tectonic` rendering engine for `LaTeX` templating of diagnostic reporting templates).

The stack slightly more involved in setting up and there are security caveats that the __administrator__ of the stack should consider when deploying the stack in their specific production environment or context. We have made scenario specific deployment templates available for admin consiguration as discussed in the next few sections

Requirements: 

- `Docker` and `docker-compose`

```bash
# deploy app on local dev machine with local security profile
# docker stack components are prefixed `stack_local`
cerebro stack deploy --config local --docker-prefix stack_local --outdir stack_local --http

# ... on first deployment admin users and passwords need to be setup ...

# reverse proxy network is included in the docker-compose config 
# to allow local http requests on `api.cerebro.localhost` and  
# to access he app on `app.cerebro.localhost` 

# start the stack, suffix with `-d` flag to run detached
cd cerebro_stack_local && docker-compose up

# go to app landing page in browser on http://app.cerebro.localhost/
```

#### Local user interface

`Cerebro` web-application stack deployment without sequence data processing or integration testing using `Nextflow` pipeline runs. Requires minimal resources but recommended around 32GB RAM, 8 CPU, 100GB - 1TB storage depending on scale, no GPU. Can be configured as web server (see documentation). Runs the `Docker` stack with `MongoDB` and `Redis` databases.

Can run on laptops with fast local upload of disk-stored models ( `.json`) of workflow outputs for the `cerebro api upload` command-line task and the local stack database to load these models from disk on demand.

```

```

#### Local development setup

`Cerebro` frontend with test outputs from workflows for `dev` and development system can run on a normal home system (16 CPU, 64GB RAM, 1 NVIDIA GPU) and is specifically configured (database size, test datasets) to run quick integration "smoke" tests with [`cipher`] for continuous integration of quality assurance during development.

```

```

#### Local production setup

Routine automated runs of the `Cerebro` profiling workflow are provisioned ideally
(NextSeq 2000 150 PE DNA + RNA + NTC libraries) around 256 CPU, 2 TB RAM, 25TB storage.

```

```

We run nanopore basecalling via `dorado` on 6 NVIDIA-A100 GPUs.

```

```

## `Cerebro` client

```bash
# get cerebro binary for support tasks
curl https://github.com/esteinig/cerebro/releases/download/latest/cerebro-latest-x86_64-unknown-linux-musl.tar.gz -o - | tar xf && mv cerebro-latest-x86_64-unknown-linux-musl cerebro

# assume `cerebro` on $PATH 
cerebro --help
```

The `cerebro-utils` client is a supplementary Python package responsible for most of the plotting and data-wrangling applications.

```bash
pip install cerebro-utils
```

The following sections outline how the client can be used in a manual workflow to deploy the stack, parse and process pipeline outputs, upload data via the `cerebro api` subcommand and work with the `cerebro api` and `cerebro-utils` clients to interrogate taxon profiles,

#### Docker stack deployment

```bash
# localhost http proxy stack config by user
cerebro stack deploy \
    --interactive \
    --config local \  # or stack_local.toml without --interactive
    --docker-prefix stack_local \
    --outdir stack_local \
    --http  # traefik rev proxy for http://{api,app}.cerebro.localhost/

# integrated traefik rev proxy network deployed
cd stack_local && docker-compose up -d

# go to app landing page in browser on http://app.cerebro.localhost/
```

#### Pipeline processing

```bash
# sample sheet for production run from fastq/
cerebro pipeline sample-sheet 
    --input fastq/ \
    --output sample_sheet.csv \
    --glob "*_{R1_001,R2_001}.fastq.gz" \
    --run-id MGP-20230420 \
    --run-date 2023-10-07 \
    --sample-group vidrl-prod \
    --sample-type csf \
    --symlinks

# nextflow run -r latest esteinig/cerebro -profile mamba,large \
#   --db cipher_db/ --production --sample_sheet sample_sheet.csv \
#   --outdir $WF_OUTDIR

# transform pipeline outputs to json models with pipeline taxonomy
cerebro pipeline process -i $WF_OUTDIR/results -o samples/ -t cipher_db

# quality control summary csv from sample models
cerebro pipeline quality -i samples/*.json -o qc.tsv -H
```

#### API interaction

When using the `cerebro api` command-line interfaces, two environment variables are considered for the API base URL (`CEREBRO_API_URL` ) and auhthentication token (`CEREBRO_API_TOKEN`). Access token expiration is set during stack deployment (default is two hours). 

When deploying locally the initial stack admin login email and password are those set in the `cerebro stack deploy` config specification fields (`cerebro_admin_email` and `cerebro_admin_password`) the main stack administration credentials.

```bash
# assume permission for data manipulation (Role::Data)
export CEREBRO_API_URL="http://api.cerebro.localhost"

# api ok
cerebro api status

export CEREBRO_API_TOKEN=$(cerebro api login -e $cerebro_admin_email -p $cerebro_admin_passowrd)

# 2023-10-29T00:20:46Z [INFO] - Login successful. Welcome back, $USER!

# api routes ok
cerebro api ping

# create a new team with default database
cerebro api team create \
    --team-name "vidrl" \
    --team-description "reference laboratory in victoria, australia"

# upload results from production workflow run
# to the default team data collection 

cerebro api upload \
    --input samples/*.json \
    --sample-sheet $WF_OUTDIR/sample_sheet.csv \
    --pipeline-config $WF_OUTDIR/config.json \ 
    --team-name "vidrl" \
    --project-name "default" \
    --output models/

# get modified default classifier evidence filter json
cerebro api utils taxa-filter-config \
    --alignment_min_ref_length 10000 \
    --kmer_min_reads 10 \
    --assembly_min_contig_length 150 \
    --output "test.filter.json"

# access all taxa per library with filter json
cerebro api taxa 
    --team-name "vidrl" \
    --project-name "default" \
    --output "taxa.filtered.csv" \
    --run-ids MGP-20230420 \
    --filter-config "test.filter.json"

# plot summary of taxa called in libraries across this run
cerebro-utils taxa plot-run-summary "taxa.filtered.csv" \
    --outdir taxa_run_filtered/ \ 
    --min_total_rpm 5 \                 # k-mer and alignment sum threshold
    --min_contigs 0 \                   # no lca diamond/blastn contig allowed
    --negative_control_tag "NTC" \      # {SAMPLE_UID}__{DNA|RNA}__{NTC|OTHER_TAG_LABEL}_{SUFFIX}
    --group_by "sample_group"           # plot panels by sample group from prod sample sheet
```


#### Clinical report templating