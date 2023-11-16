# Development

!!! note

    This section is for developers - experience with `docker compose`, `rustc` and `git` for conventional commits and semantic versions is required!

We have tried to make development as easy as possible using a local hot-reload development configuration of the stack. Minimal configuration is required to get started and nearly all stack components (see exceptions) link directly into a repository clone with the ability to quickly change between branches, re-deploy changes and minimize container re-builds. Please see the style guide to observe conventional commits, issue tracking and semantic versioning of pull requests. 

## Styleguide

All commits are made with `cocogitto` in the usual pattern `cog commit {task} "{msg}" {scope}` e.g. `cog commit chore "fix spelling mistake" docs`. Note that apull requests implemented commits that should not appear in the CHANGELOG should be tasked as `chore` as these are excluded from CHANGELOG. Ideally, development commits should always implement `chore` and then squashed into a single commit on merge into `dev` with a meaningful message. Breaking changes indicated as such with `cog` will  will trigger changes to major version, conventional commits tagged as `fix` or `feat` to minor version and any other changes to patch version.

1. Development occurs on the `dev` branch - all pull requests are merged into `dev` before release on `main`
2. Open an issue for your development branch using the templates in the repository and create a new branch linked to the issue. Branch syntax follows conventional commit styles with `{task}/{branch-name}`, e.g. `chore/docs-spelling-11`.
3. When completed, open a pull request of your branch into `dev` and request review
4. Commits from your branch should be squashed and merged as a single conventional commit style message e.g. `chore(docs): fixed to documentation`. This ensures that only the pull request commit merged into `dev` appears in the changelog for new releases. In exceptional cases development commits may be included without squashing but development commits must adhere to conventional commit styles

Release commits on `main` by the repository owners pull changes from `dev` via a `release/{version}` branch. Release actions trigger `cog bump --auto` changelog and version bumps, including compilation of Linux binaris of `cerebro` attached to the release. 

## Requirements

* `Docker` and `docker compose`
* `Rust >= v1.72`, `rust-analyzer` with IDE
* `git` and `cocogitto >= v1.2.1` 

## Basic Deployment

Deploy a development version locally with hot-reloading of application with `npm run dev` and `cargo watch` for the server backend (`cerebro stack run-server`). 

!!! note

    Deployment of the stack via `docker compose up` should be monitored in development mode as informative messages from compiler and components will be logged via `docker compose`. 

    Stack deployments should *always* be downed from the deployment repository, this can be done after upping the stack in the foreground (to monitor logs) by running `docker compose down` in the deployment directory after terminating the running stack. 
    Some arcane reverse proxy errors may occurr otherwise, usually noticeable by a `Bad Gateway` message when navigating to the application domain - in these cases, the containers and network used in the deployment should be removed manually.

Once the stack is deployed, hot-reloads of the server backend via `cargo watch` can be triggered manually by any change to a `<deployment>/cerebro/.trigger` file in the deployment repository, which can be made easily through
your developer environment or for example `echo "trigger me build" >> <deployment>/cerebro/.trigger`. This is to avoid somewhat long compilation times (~ 30 seconds) for the full `cerebro` binary in debug mode (and longer in release mode). 
Note that performance during development in compiled binaries for debugging is usually inferior to release binaries.

Other development branches or specific revisions can be checked out manually in the deployment repository e.g. `cd <deployment>/cerebro && git checkout <branch>`.

!!! warning

    Never deploy a development configured stack in production! Development mode is meant for local stack deployment and testing and is not safe in production.

If you want to deploy a "dev" branch on the web for staging changes before release to production, use the default mode during deployment  and checkout the desired 
state with `cerebro stack deploy --branch` or `cerebro stack deploy --revision` - this will configure and deploy the stack for production on the specific branch or revision,
but requires container re-builds for updates.

!!! warning

    Modification to files in `templates/stack` requires re-deployment of the stack.

## Advanced Deployment

Multiple deployments with unique output directory names can be run simulataneously and require only slight adjustments to the deployment domains.



## Example 

In the following example we deploy an existing development branch from an example issue gixing some documentation (`docs/dev-branch-test`). 

```bash
# get latest cerebro binary to initiate the deployment
curl https://github.com/esteinig/cerebro/releases/download/latest/cerebro-latest-Linux_x86_64.tar.xz -o - | tar -xzO > cerebro


# deploy a local dev configuration 
cerebro stack deploy --outdir cerebro_dev --config dev-local --dev --branch docs/dev-branch-test

# or manually: cerebro stack deploy --outdir cerebro_dev --config dev-local --dev && cd cerebro_dev/cerebro && git checkout docs/dev-branch-test

# enter deployment and up the stack for monitoring rustc compilation
# cerebro is re-compiled in debug mode on first deployment
cd dev_local/ && docker compose up

# changes to the app ui are reflected immediately due to node dev server 
# changes to the rust codebase require rebuild of the app e.g. when changing api code

# [TBD] fmt and clippy hooks are executed for rust changes 

# down the containter stack and remove the stopped containers - this is recommended 
# as a local reverse proxy can return docker network errors and db connections can fail
docker compose down && docker container rm $(docker ps --format "{{.Names}}" | grep "$CEREBRO_PREFIX"-)

# build flag triggers rebuild of changed components in the docker compose file
docker compose up --build

# checkout current feat branch and commit changes
git checkout feat/new_feature
cog commit chore "minor changes to auth api" auth
git push origin feat/new_feature

# create new pr with conventional commit style 
# to merge feat/new_feature into dev

# squash working commits on pr into dev to
# the primary conventional commit
#
# docs/readme branch -> dev
# readme(docs): dev section draft
#
# dev pr squash commits are published in main release changelog
```

## Git progression for release

```bash
# checkout latest main
git checkout main 

# pull latest dev into main and run auto changelog
git fetch dev 

# checkout auto changelog for version release
cog changelog --auto

# get update version 
$RELEASE_VERSION=$(cog bump --dry-run --auto) 

# checkout version release branch
git checkout -b release/$RELEASE_VERSION

# cog commit release chore to version
git add . && cog commit chore "$RELEASE_VERSION" release

# push to remote to trigger release and build actions
git push origin release/$RELEASE_VERSION

# create pr of `release/$RELEASE_VERSION` into 'main'
# release with latest changelog and binaries is created on ci/cd

```