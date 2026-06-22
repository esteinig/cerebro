# Documentation site

The Cerebro documentation is built with [mdBook](https://github.com/rust-lang/mdBook),
the Rust-native documentation tool, and published to GitHub Pages from the
`feat/production` branch.

## Structure

- `book.toml` — book configuration: title, theme (light default + optional dark), and
  the source directory (`src = "docs"`, so the existing `docs/` tree is the book).
- `docs/SUMMARY.md` — the table of contents / navigation. Every page must be listed
  here to appear in the book.
- `docs/index.md` — the Cerebro landing page (prefix chapter).
- `docs/cerebro-fs/` — the cerebro-fs subsystem documentation (overview, architecture,
  operations runbooks, reference, and production-readiness governance).
- `.github/workflows/docs.yml` — builds and deploys the site to GitHub Pages.
- The build output (`book/`) is git-ignored.

New subsystems are added as `docs/<subsystem>/` plus an entry under the **Subsystems**
part in `docs/SUMMARY.md`, mirroring the cerebro-fs layout.

## Install mdBook

```bash
cargo install mdbook          # from crates.io
# or download a release binary from https://github.com/rust-lang/mdBook/releases
```

## Preview and build locally

```bash
mdbook serve --open    # live-reloading preview at http://localhost:3000
mdbook build           # render the static site into ./book
mdbook test            # optional: check Rust code samples in the docs
```

## Themes

mdBook ships a built-in theme picker. The book defaults to the **light** theme with
**navy** as the dark option (`default-theme` / `preferred-dark-theme` in `book.toml`);
readers can switch via the paintbrush icon. Custom styling can be layered in with
`additional-css` if needed.

## Publishing to GitHub Pages

The `docs` workflow runs on every push to `feat/production` that changes `docs/`,
`book.toml`, or the workflow, builds the book with mdBook, and deploys `./book` to
GitHub Pages.

One-time repository setup:

1. **Settings → Pages → Build and deployment → Source:** select **GitHub Actions**.
2. In `book.toml`, set `site-url` to your Pages base path so assets resolve. For a
   project site this is the repository name, e.g. `site-url = "/cerebro/"`; for a
   user/org root site or a custom domain use `/`.
3. Push to `feat/production`; the site publishes at
   `https://<org-or-user>.github.io/<repo>/`.

To deploy by hand instead of via Actions, `mdbook build` and publish the contents of
`./book` to your `gh-pages` branch (for example with `ghp-import -np book`).
