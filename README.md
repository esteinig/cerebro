# Cerebro

Central nervous system specific metagenomic diagnostics workflow for protocols at the Victorian Infectious Diseases Reference Laboratory (VIDRL) and the CNS Meta-GP group (see below) at the Peter Doherty Institute. 

## Features

- `Nextflow` pipeline for meta-metagenomics taxonomic profiling targeting quantitative pathogen detection using Illumina/Nanopore data offering in depth quality control taxonomic classifier options. Includes custom tools for the classification pipeline focusing on human sequence removal (`scrubby`) and viral alignment filters and auto-mated reference selection for coverage remapping subworkflows (`vircov`)

- `Rust` library to run the database and user interace API server including authenticated routes and authorization middleware with `ActixWeb`; providing endpoints for workflow output processing (quality control, classification, workflow and sample configurations) and access to structured storage. Includes command-line client to run supplementary tools in the main pipeline support production execution of the pipeline (create sample sheets, run automated input framework for production, ...)

- `Docker` and `docker-compose` configs templated from the `cerebro` binary for deployment of a `Svelte UI` with `Skeleton UI` component library + `MongoDB` and `Redis` databases to manage, interpret and collaboratively decide on candidate pathogens and generate auditable clinical reports with the `Tectonic` LaTeX to PDF engine and report template.

- Provisioning of traceable and versioned metagenomic reference databases and quality or diagnostic genome feature masks; continuous integration testing of full diagnostic workflow against clinical samples with orthogonal diagnostic data or synthetic simulated datasets (`cipher`)

## Usage

Requirements:
- `Linux OS`

`Cerebro` consists of two parts: 

1. `Nextflow` pipeline for sequence data (Illumina PE) ingestion and processing, including integration testing and assay validation workflows.
2. `Rust` and `docker-compose` application stack (`MongoDB` + `Redis` + `Sveltekit` + `Rust` processinglibrary that implements an `ActixWeb` API server).

The `cerebro` command-line interface is available as binary release and implements supportive functions and interactions with the application stack and pipeline. 

```bash
# get cerebro binary for support tasks
curl https://github.com/esteinig/cerebro/releases/download/latest/cerebro-latest-x86_64-unknown-linux-musl.tar.gz -o - | tar xf && mv cerebro-latest-x86_64-unknown-linux-musl cerebro

# assume `cerebro` on $PATH 
cerebro --help
```

### `Cerebro` pipeline

Requirements for local execution: 

- `Conda` or `Mamba`
- `Nextflow >= v23.10.04`

```bash
# pull latest from github, show help menu, use mamba local env
nextflow run -r latest esteinig/cerebro -profile mamba --help

# k-mer profiling example: run standard qc and kmer tax classifiers with mamba envs and cipher kmer db build directory using a large resource profile

# provision with latest cipher kmer db build, may take some time
nextflow run -r latest esteinig/cerebro -profile mamba \
    --cipher_download \
    --cipher_revision latest \
    --cipher_modules kmer \
    --outdir cipher_kmer_db/

# kmer tax profile on input read dir (pe illumina)
nextflow run -r latest esteinig/cerebro -profile mamba,large,kmer \
    --db cipher_kmer_db/ \
    --outdir kmer_test_run/ \
    --fastq "fastq/*_{R1,R2}.fq.gz"

# without kmer profile classifiers are manually enabled to demonstrate how nested params in `nextflow.config` are modified on the nextflow cli
nextflow run -r latest esteinig/cerebro -profile mamba,large \
    --db cipher_kmer_db/ \
    --outdir kmer_test_run/ \
    --fastq "fastq/*_{R1,R2}.fq.gz" \
    --taxa.alignment.enabled false \  # disable alignment module
    --taxa.assembly.enabled false \   # disable assembly module
    --taxa.kmer.enabled true \
    --taxa.kmer.kraken2uniq.enabled true \
    --taxa.kmer.kmpc.enabled true \
    --taxa.kmer.metabuli.enabled true 
```

Example commands for executing the `Nextflow` pipeline and processing local output with the `cerebro` command-line interface:

1. Reference data and index provision for pipeline setup ([`Cipher`]())
2. [Kraken pathogen detection]() protocol and post-processing
3. Multi-classifer, multi-database taxonomic classification ([meta-metagenomics protocol v1]())
4. Taxonomic profiling for [central nervous system assay v1]() (META-GP)
5. [Viral genome recovery for panviral enrichment]() from rapid antigen tests (RAT)
6. `Cerebro` [pipeline integration testing]() with `Cipher` reference datasets and actions
7. Testing impact of host read depletion on classification performance with `Cipher` and `Scrubby`

### `Cerebro` stack

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


### `Cerebro` client

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

# ... pipeline run using --production and --sample_sheet ...

# transform pipeline outputs to json models with pipeline taxonomy
cerebro pipeline process -i $WF_OUTDIR/results -o samples/ -t $WF_OUTDIR/cipher_db

# quality control summary csv from sample models
cerebro pipeline quality -i samples/*.json -o qc.tsv -H
```

#### API interaction

Global args on the `cerebro` client supercede the default environment variables `CEREBRO_API_URL` for the API base URL and `CEREBRO_API_TOKEN` for the access token we can get and set by logging in to an authenticated session. Access token expiration is set during stack deployment (default is two hours) and when deploying locall the account email and password are those set in the `cerebro stack deploy` config specification (`cerebro_admin_email` and `cerebro_admin_password`)

```bash
# assume account with permission for data manipulation (Role::Data)

export CEREBRO_API_URL="http://api.cerebro.localhost"

# api ok
cerebro api status

export CEREBRO_API_TOKEN=$(cerebro api login -e $CEREBRO_ACCOUNT_EMAIL -p $CEREBRO_ACCOUNT_PASSWORD)

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
    --group_by "sample_group"           # sample group col from prod sample sheet for plot panels
```

### Deployment scenarios

Potential ways of operating the `Cerebro` stack or pipeline on different resourced infrastructure.


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
