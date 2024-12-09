# Cerebero stack

Install the current version:

```bash 
mamba create -n cerebro -c esteinig cerebro
```

Activate environment:

```bash
mamba activate cerebro
```

## Localhost configuration

We will use the `localhost` configuration template and setup a production stack in the local environment using `--interactive` mode,
which prompts for a set of important admin variables like usernames and passwords:

```bash
cerebro stack deploy --name cerebro-local --outdir cerebro-local --config localhost --interactive
```

All prompts can be entered on the command-line:

```bash
cerebro stack deploy --help
```

**Cerebro**

* `--cerebro-admin-email`: Cerebro admin email to login as admin
* `--cerebro-admin-password`: Cerebro admin password to login as admin
* `--cerebro-admin-name`: Cerebro admin full name for user profile

**MongoDB database configuration**

* `--db-root-username`: MongoDB root username
* `--db-root-password`: MongoDB root password
* `--db-admin-username`: MongoDB admin username
* `--db-admin-password`: MongoDB admin password

**CerebroFS file system configuration**

* `--fs-primary`: Path to primary data storage
* `--fs-secondary`: Path to secondary data storage (backup)

This will create the output directory `cerebro-local` with all required files and configurations for launching the stack.

## Stack management

Change directory and launch the stack:

```bash
cd cerebro-local && docker compose up
```

# Cerebro application

You can now login using your `Cerebro` admin email and password on: `http://app.cerebro.localhost`

## Team creation and member assignment



## Database and project collections


