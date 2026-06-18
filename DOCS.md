# Documentation site

The Cerebro documentation is built with [MkDocs](https://www.mkdocs.org/) and the
[Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) theme, and published
to GitHub Pages from the `feat/production` branch.

## Structure

- `mkdocs.yml` — site configuration: theme, clinical light/dark palette, and navigation.
- `docs/index.md` — the Cerebro landing page (purpose and quick start).
- `docs/cerebro-fs/` — the cerebro-fs subsystem documentation (overview, architecture,
  operations runbooks, reference, and production-readiness governance).
- `docs/assets/extra.css` — small, optional styling refinements.
- `.github/workflows/docs.yml` — builds and deploys the site to GitHub Pages.

New subsystems are added as `docs/<subsystem>/` and a corresponding entry under
**Subsystems** in `mkdocs.yml`, mirroring the cerebro-fs layout (overview →
architecture → operations → reference → governance).

## Preview locally

```bash
pip install "mkdocs-material>=9.5"
mkdocs serve            # live preview at http://127.0.0.1:8000
mkdocs build --strict   # optional: fail on broken links / nav before pushing
```

## Publishing

The `docs` workflow runs on every push to `feat/production` that changes `docs/`,
`mkdocs.yml`, or the workflow itself, builds the site, and deploys it to GitHub Pages.

One-time repository setup:

1. **Settings → Pages → Build and deployment → Source:** select **GitHub Actions**.
2. Set `site_url` in `mkdocs.yml` to the published URL
   (`https://<org-or-user>.github.io/<repo>/`) so canonical links and the sitemap are
   correct.
3. Push to `feat/production`; the site publishes at the URL above.

The theme fetches the Inter and JetBrains Mono fonts from Google Fonts. For
fully self-contained / privacy-sensitive deployments, set `font: false` in
`mkdocs.yml` to fall back to system fonts.
