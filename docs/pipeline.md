
#### Cerebro pipeline command-line configs 

```bash
# main profiles
nextflow run -r latest esteinig/cerebro \
    # local conda,mamba env created for each process
    -profile conda,mamba
    # kmer, alignment, assembly classifier modules
    -profile kmer,alignment,assembly
    # process-configured resource profiles
    -profile small,medium,large,galactic
    # specific protocol configs
    -profile cns_assay,panviral,aneuploidy
    # data provisioning
    -profile cipher_db
    # workflow testing
    -profile cipher_test
    # workflow dev
    -profile io_mode,qc_mode,dev_mode

# nested module params in `nextflow.config`
nextflow run -r latest esteinig/cerebro -profile mamba \
    # all reference db + idx + tax
    --db "cipher_db/" \
    # paired read input
    --fastq "fastq/*_{R1,R2}.fq.gz" \
    # production input tracked files
    --production \
    --sample_sheet "sample_sheet.csv" \
    # output directory
    --outdir "module_test" \
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
    --taxa.assembly.diamond.enabled
    # post workflow processing and api interaction
    --cerebro.quality.enabled \                   # create run summary qc table
    --cerebro.sample.enabled \                    # create run aggregated cerebro sample json
    # require: --sample_sheet and --production
    --cerebro.api.enabled true \  
    --cerebro.api.url $CEREBRO_API_URL \ 
    --cerebro.api.token $CEREBRO_API_TOKEN \
    --cerebro.api.status.enabled \                # status report logging of pipeline updates
    --cerebro.api.status.slack.enabled \          # status report logging to configured slack channel
    --cerebro.api.report.enabled \                # create run report (sample, qc tracking)
    --cerebro.api.report.slack.enabled \          # post run report to configured slack channel
    --cerebro.api.upload.enabled \
    --cerebro.api.upload.team "VIDRL" \
    --cerebro.api.upload.database "PRODUCTION" \
    --cerebro.api.upload.collection "MGP-CNS-20231012"


# label-specific resource config with nested params
# labels are defined in `lib/configs/resources.config` and
# applied to process definitions in `lib/processes/*.nf`

# cpus, memory, time, conda, container (docker)
nextflow run -r latest esteinig/cerebro -profile mamba \
    --resources.kraken2uniq.cpus 64 \
    --resources.kraken2uniq.memory "256 GB" \
    --resources,minimap2_align.cpus 32 \
    --resources.minimap2_align.memory "32 GB" \
    --resources.minimap2_align.conda "envs/minimap2.dropin.yml"
```

## Cerebro pipeline application examples


### Reference data and index provision

Provision the `Cerebro` pipeline with the `Cipher` reference databases, database masks and taxonomy based on processed NCBI and GTDB genomes icludinng human pangenome host databases.

```bash
# run db provision task, may take some time
nextflow run -r latest esteinig/cerebro -profile mamba,cipher_db \ 
    --cipher_revision latest --outdir cipher_db
```

### Kraken pathogen detection protocol and post-processing

Run the `Kraken2Uniq` pathogen detection protocol as described in Nature Protolcs. Include automated remapping of classified reads at species level again a representative genome auto-selected from the reference index sequences.


```

```

### Multi-classifier, multi-database taxonomic classification


```

```


### Taxonomic profiling for central nervous system assay (Meta-GP)

Local execution with large resource profile with `mamba` dependencies on a specific central nervous system assay based configuration which included NEBNext UMI index deduplication, ERCC and T4 phage assay controls and taxonomic profiling withoutput generated for just 

```bash
nextflow run -r latest esteinig/cerebro -profile mamba,cns_assay \
    --outdir cns_test --db cipher_db/
```

Create sample sheet and run in production mode

```bash
nextflow run -r latest esteinig/cerebro -profile mamba,cns_assay \
    --outdir cns_test --db cipher_db/ --sample_sheet sample_sheet.csv
```

Local execution with


#### Viral genome recovery for pan-viral enrichment from clinical samples

Protocol for data derived from pan-viral enrichment as for our publicaiton on sequencing a range of respiratory viruses from rapid antigen tests (RAT). 

Tested on other clinical sample types treated with the enrichment protocol (brain tissue, ocular fluid, nasal swab) at VIDRL and Williamson Research Group. Read more in our recent publication o n RAT sequencing for virus detection and whole genome subtyping.


#### Experimental host genome analysis from background nucleic acid


## Deployment scenarios and resource requirements

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