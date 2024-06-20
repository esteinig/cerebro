# Cerebro

Under active development for production release. Not recommended for deployment at this stage. 

This is a preliminary public release of code for the viral enrichment branch of the pipeline used in:

> Michael A Moso*, George Taiaroa*, Eike Steinig*, Madiyar Zhanduisenov, Grace Butel-Simoes, Ivana Savic, Mona L Taouk, Socheata Chea, Jean Moselen, Jacinta O’Keefe, Jacqueline Prestedge, Georgina L Pollock, Mohammad Khan, Katherine Soloczynskyj, Janath Fernando, Genevieve E Martin, Leon Caly, Ian G Barr, Thomas Tran, Julian Druce, Chuan K Lim, Deborah A Williamson - **Non-SARS-CoV-2 respiratory viral detection and whole genome sequencing from COVID-19 rapid antigen test devices: a laboratory evaluation study** - Lancet Microbe (2024) -[10.1016/S2666-5247(23)00375-0](https://doi.org/10.1016/S2666-5247(23)00375-0)

## Pipeline Testing

```
# Check for errors during development - this will print the startup and completion
# messages to the console and exit the pipeline execution gracefully if not errors
# were found:

nextflow run cerebro/ -profile test_dev

# Check for input checking with minimal database configurations for quality control
# with the human reference database index and 

nextflow run cerebro/ -profile db,db_ont,test_io
nextflow run cerebro/ -profile db,db_sr,test_io

```