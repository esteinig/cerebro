# Development

Development setup: 

* `Docker` and `docker compose`
* `rust >= v1.72` with `cargo-watch`, `cargo fmt`, `cargo clippy` 
* `cocogitto >= v1.2.1` for convenvional commits and semantic version releases via `cog commit {feat} "{msg}" {scope}` command-line drop-in for `git commit`

Commits must use `cog commit` or manual conventional commit style - these can be squashed or rebased for a single merge commit into `dev` from the conventional commit named branch for example `feat/new_feature` or `chore/docs`. Pull requests only into `dev` never `main`. Release commits on `main` pull changes from `dev` via a `release/{version}` branc. Release action triggers the `cog bump --auto` changelog and version bump with cross compilation pipelines that attach Linux binaries on release publication. 

## Setup dev environment

```bash
# get latest cerebro bin for deployment
curl https://github.com/esteinig/cerebro/releases/download/latest/cerebro-latest-Linux_x86_64.tar.xz -o - | tar -xzO > cerebro

# clone repo to track state
git clone https://github.com/esteinig/cerebro.git@latest

# docker prefix for this deployment
CEREBRO_PREFIX="dev_local" 

# link deployment to the git repo and
# deploy on `dev_local` docker prefix
# with a `dev.cerebro.localhost` domain
# via the contained traefik reverse proxy
cerebro deploy \
    --dev ${pwd}/cerebro 
    --config local.toml \
    --docker-prefix $CEREBRO_PREFIX \
    --outdir dev_local/ \
    --http  # traefik rev proxy for http://{api,app}.dev.cerebro.localhost/
    --http-domain "dev.cerebro.localhost"

# enter deployment and up the stack
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