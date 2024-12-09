# Cerebro production 

Note that production operations currently work only in local configurations and not web-configurations due to required access to `CerebroFS` (currently not distributable, see also security considerations).

Cerebro operates a (fully distributed) system of watchers and towers, which are registered with the stack:

* Watchers watch a directory path for read files (FASTQ) and upload detected files to `CerebroFS`
* Towers perpetually request read files from a staging area and pull the reads from `CerebroFS` into the configured sub-pipelines of `Cerebro`


## Example production environment

We will now setup an example of the production data ingestion and pipeline configuration to run the `panviral-enrichment` pipeline of `Cerebro` in production.

### Watchers

Open a new terminal and login (assuming localhost default configuration):

```bash
export CEREBRO_API_TOKEN=$(cerebro-client login -e your-email -p your-password)
```

Start watching a new directory:

```bash
cerebro-watcher --team CNS watch --path test_watcher --name TestWatcher --location Home --format fastq-pe
```

Note the `--team` designation before the `watch` command-line task - the team must exist and the user must be a member of the team to suscessfully create a new watcher:

> Here, we create a new directory `test_watcher` and register the watcher with its name (`TestWatcher`), location (`Home`) and the format of the read files to gather and upload (`fastq-pe` for paired-end FASTQ), which defaults to the file glob: `*_{R1,R2}.fastq.gz`. Default globs for `--format` can be overwritten with `--glob`.

### Towers

