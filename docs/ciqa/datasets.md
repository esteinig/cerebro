# ciqa — datasets (training & QC)

CIQA works with two kinds of dataset that share one storage shape and one evaluator but serve
different purposes: **training sets**, for interactive model and prompt development, and
**quality-control (QC) sets**, for automated, reproducible regression testing against a baseline.
This page explains both and how they relate. For the scoring maths in context, see
[regression-testing](regression-testing.md).

## The shared shape

Both kinds of dataset are built from the same unit: a **record** plus a **`PrefetchData` blob**.
The record holds the metadata and (for QC) the reference truth; the `PrefetchData` — the tiered,
filtered evidence for one sample ([prefetch & filtering](prefetch-filtering.md)) — is stored in
**GridFS** and referenced by id. This mirrors how the training prefetch records already work, so
a QC dataset "functions the same as the training one" — same persistence, same retrieval, same
`PrefetchData` shape — and the two can share tooling.

```
   dataset record (Mongo)            PrefetchData (GridFS)
   ├─ name / version / description    ├─ primary / secondary / target
   ├─ sample / sample_type            ├─ contamination vectors
   ├─ reference_result   (QC)         └─ MetaGpConfig
   ├─ reference_candidates (QC)
   └─ blob id ──────────────────────► (retrieves the PrefetchData)
```

## Training sets

A **training set** is a collection of prefetch items used to develop and tune the generative step
interactively — trying prompts, models, and filter settings and judging the result by hand.

### The workflow

1. **Assemble.** Prefetch items are collected into a named training collection; each item carries
   its `PrefetchData` (in GridFS) and metadata, and may carry a **preselect** hint (a seeded
   candidate to start from).
2. **Session.** A training **session** presents the items for labelling; a reviewer records a
   **decision** per item.
3. **Evaluate.** The session reduces the decisions to **sensitivity / specificity** so a reviewer
   can see, immediately, how a configuration performs on the set.

### Decisions and the evaluator

Each item's decision is one of:

| `Decision` | Meaning |
|---|---|
| `TP` | True positive — predicted the right pathogen on a positive sample |
| `TN` | True negative — correctly negative |
| `FP` | False positive — called positive (or the wrong pathogen) when negative/absent |
| `FN` | False negative — missed a true positive |
| `Excluded` | Excluded from the metric (e.g. limit-of-detection exclusion) |

From the counts, the evaluator computes:

```
sensitivity = TP / (TP + FN)
specificity = TN / (TN + FP)
```

(with PPV and NPV alongside). This is the **same evaluator** that scores QC and regression — one
definition of TP/FP/TN/FN, so a number from a training session and a number from a QC run mean
exactly the same thing. The preselect hint only seeds the reviewer's starting selection; it does
not change the maths.

## QC sets

A **QC set** is the automated counterpart: a versioned collection of prefetch items **with
reference truth**, used to score META-GPT's diagnoses and to anchor a **baseline** for regression.

### What a QC set adds over a training set

- **Reference truth on every item.** Each record carries the sample's `reference_result`
  (`Positive` / `Negative`) and `reference_candidates` (the true species), plus an optional
  limit-of-detection exclusion flag — the truth the evaluator scores against.
- **Versioning.** A QC dataset carries a **version**. Truth changes (added samples, corrected
  labels) are made by **bumping the version**, never by mutating a version in place, because a
  baseline is tied to a specific dataset version. This keeps every historical comparison
  reproducible.
- **A baseline.** A QC dataset version can have a promoted **baseline** — the known-good
  statistics (and the model manifest that produced them, and the pass thresholds) that a future
  run is judged against (see [regression-testing](regression-testing.md)).

### Building a QC set

Reference-truthed prefetch items are registered as a QC dataset and indexed/staged through Cerebro
Fs (the `upload` command stages the underlying read sets with a `ciqa-<run>` prefix and registers
the dataset record). Thereafter the dataset is addressed by name and version through the `ciqa/*`
API.

```
ciqa/dataset      register / list / delete QC dataset records (PrefetchData in GridFS)
ciqa/baseline     create / list immutable, versioned baselines
ciqa/baseline/promote   mark a baseline active for a dataset version  (operator-gated, audited)
ciqa/regression   submit a run manifest → evaluate vs the active baseline → store + return a report
```

Deleting a QC dataset is a **consequential, audited, operator-gated** action — it is reference
truth other results depend on.

## Training vs QC — when to use which

| | Training set | QC set |
|---|---|---|
| Purpose | Interactive development of model/prompt/filters | Automated regression gating |
| Truth | Reviewer's live judgement | Recorded reference truth on every item |
| Lifecycle | Mutable, iterative | Versioned; truth changes bump the version |
| Output | Immediate sens/spec for a reviewer | A baseline + regression reports |
| Storage shape | Record + `PrefetchData` in GridFS | The same, plus truth + a baseline |
| Evaluator | The one shared evaluator | The one shared evaluator |

Use a training set to find a configuration you trust; promote its performance into a QC baseline;
then use the QC set to keep that performance honest as things change.

## Tested here vs your environment

The evaluator (decisions → sensitivity/specificity) is **pure, unit-tested logic**, cross-checked
so the training and QC numbers are identical for the same inputs. The dataset persistence and the
`ciqa/*` endpoints are **environment-validated** against a running stack, including the audit
entries for promotion and deletion; see [`validation.md`](validation.md).

> **Commands & endpoints.** For copy-pasteable CLI and curl for every step, see [Commands & API](commands.md).
