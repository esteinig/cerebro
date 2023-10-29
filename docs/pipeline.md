## Examples

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