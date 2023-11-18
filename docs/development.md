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

## Basic deployment

Deploy a development version locally with hot-reloading of application with `npm run dev` and `cargo watch` for the server backend (`cerebro stack run-server`). 

!!! note

    Deployment of the stack via `docker compose up` should be monitored in development mode as informative messages from compiler and components will be logged via `docker compose`. 

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

    Changes to `templates/stack` files and to `src/stack` core deployment functions that alter deployment configuration are currently not considered for hot-reload in a development stack. Modifications must be tested manually and the stack must be re-deployed from scratch.

## Advanced deployment

Multiple stack deployments with unique output directory names can be run simultaneously and require only slight adjustments to the deployment subdomains:


## Known issues

Overall it is usually recommended to down the stack so that `docker` network configurations for the `traefik` deployment are properly removed and do not interfere with subsequent stack launches.

1. Stack deployments should *always* be downed from the deployment repository, this can be done - after upping the stack in the foreground to monitor logs - by running `docker compose down` in the deployment directory after terminating the running stack. Some arcane reverse proxy errors may occurr with `traefik` otherwise, usually noticeable by a `Bad Gateway` or `404` message when navigating to the application domain in your browser - in these cases, the stack should be downed **or** containers and network removed manually if all else fails. Container prefix for the stack is the `<deployment>` directory name.
2. Occasionally the server may not be able to connect to the databases in development mode for some reason - you should see a checkmark and a `Connected to {MongoDB, Redis session, Redis one-time} database` message when starting the server in the compose file on container with tag: `<depoyment>_cerebro-api-1`. Down the stack and up it again, usually this mitigates any issues with internal networking ([tracked in this issue]).
3. After setting up a stack from scratch and navigating to the interface in development mode, the `vite` development server will optimize some components **after login** - sometimes this will result in a failed login attempt without error message. You can monitor this in the `<depoyment>_cerebro-app-1` log output and repeat the login to continue as usual.
 

## Deployment example

In the following example we deploy an existing development branch from a fictional example issue that fixes some documentation (with an imaginary issue opened on branch `docs/dev-branch-test`):

```bash
# get latest cerebro binary to initiate the deployment
curl https://github.com/esteinig/cerebro/releases/download/latest/cerebro-latest-Linux_x86_64.tar.xz -o - | tar -xzO > cerebro

# deploy a local dev configuration 
cerebro stack deploy --outdir cerebro_dev --config http-local --branch docs/dev-branch-test --dev

# or manually checkout the desired branch:
# cerebro stack deploy --outdir cerebro_dev --config http-local --dev && \
#   cd cerebro_dev/cerebro && git checkout docs/dev-branch-test

# up the stack for monitoring compiler and stack operations
# cerebro is re-compiled in debug mode on first deployment
cd cerebro_dev && docker compose up

# open `cerebro_dev/cerebro` subdirectory in your ide and
# make changes to the app or rust code ...
# changes to the app are reflected immediately in the web interface
# changes to the server are triggered on changes to `cerebro_dev/cerebro/.trigger`

# `ctrl + c` to exit the forground compose stack
# `docker compose down` to down the compose stack
```

### Commit example

You can commit changes to the branch in another terminal during development:

```bash
# commit some changes to the documentation 
cd cerebro_dev/cerebro
git add docs/development.md
# use conventional commits with cocogitto
cog commit chore "add deployment section" docs
git push origin docs/dev-branch-test

# open pr for this branch and review 
# squash merge changes with conventional commit style title and pr link
# for example:
# docs(stack): add development documentation section (#10)
```