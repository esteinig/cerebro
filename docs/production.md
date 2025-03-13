# Cerebro production 

Cerebro operates a (fully distributed) system of watchers and towers, which are registered with the stack:

* Watchers watch a directory path for read files (FASTQ) and upload detected files to `CerebroFS`
* Towers perpetually request read files from a staging area and pull the reads from `CerebroFS` into the configured sub-pipelines of `Cerebro`

Note that production operations currently work only in local configurations due to required access to `CerebroFS` (see security considerations).

## Example production environment

We will now setup an example of the production data ingestion and pipeline configuration to run the production version of the `panviral-enrichment` pipeline with `Cerebro`.

### Watchers

Open a new terminal and login - we are assuming the `localhost-insecure` stack deployment configuration which has a default admin user profile:

```bash
export CEREBRO_API_TOKEN=$(cerebro-client login -e admin -p admin)
```

Start watching a new directory `test_watcher` in the current directory:

```bash
cerebro-watcher --team CNS watch --path test_watcher --name TestWatcher --location Home --format fastq-pe
```

Note the `--team` designation before the `watch` command-line task - the team must exist and the user must be a member of the team to suscessfully create a new watcher. Watchers are team specific so that members of the team with the correct permissions have access to view collected files and launch pipelines.

> Here, we create a new directory `test_watcher` and register the watcher with its name (`TestWatcher`), location (`Home`) and the format of the read files to gather and upload (`fastq-pe` for paired-end FASTQ), which defaults to the file glob: `*_{R1,R2}.fastq.gz`. Default glob formats specified with `--format` can be overwritten with custom `--glob` formats.

### Towers


Open a new terminal and login:


```bash
export CEREBRO_API_TOKEN=$(cerebro-client login -e admin -p admin)
```

Deploy a new tower

```bash
cerebro-tower --team CNS watch --path test_watcher --name TestWatcher --location Home --format fastq-pe
```

it is important that this occurs on the machine where your local stack is running because of the `CerebroFS` components that need to be available to the tower.
