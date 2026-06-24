# CIQA regression-testing runbook

Operator and QA procedures for measuring META-GPT's **sensitivity / specificity** from its
diagnoses and **regression-testing** a model or configuration change against a frozen **baseline**.
This is the change-control instrument for the generative step. For the datasets it runs on, see
[datasets](datasets.md); for producing the diagnoses it scores, see
[META-GPT diagnosis](meta-gpt-diagnosis.md).

> **Audience.** A method developer or QA operator with a registered QC dataset (reference truth),
> a set of META-GPT outputs and their run manifest, and — for promotion — the authority to change
> the baseline a clinical change is judged against.

> **Golden rule.** A regression run is not done when a command exits `0`. It is done when the
> **Verify** step confirms the report's verdict against the baseline. Every scenario ends in one.
> And: **CIQA never promotes a baseline for you.** A better-looking result is a candidate, not a
> decision.

---

## 1. What regression testing decides

For a run of META-GPT over a QC dataset, regression answers one question: **did this model /
configuration regress relative to the baseline?** It decides with two complementary signals:

- **Absolute floor.** Sensitivity and specificity must each be at or above the baseline's
  **thresholds**. A breach is a regression on that axis — a hard gate.
- **Paired significance (McNemar).** The run's per-sample correctness is compared to the
  baseline's with a **McNemar** test, and — when several configurations are compared at once — with
  a multiple-comparison adjustment (Bonferroni / Holm / Benjamini-Hochberg). A significant adverse
  shift is a regression even if still above the floor; a non-significant result is *not* proof of
  no change (see [RT-4](#rt-4--underpowered-or-unbalanced-panel)).

The output is a **regression report** (current vs baseline statistics, the deltas, the McNemar
p-value, the threshold and regression verdicts, and the per-sample decisions) and an **exit-code
gate** for CI/Nextflow. **Nothing about the baseline is mutated.**

## 2. The inputs you need

| Input | From | Notes |
|---|---|---|
| META-GPT outputs (`{sample}.model.json`) | the diagnose step | [`meta-gpt-diagnosis.md`](meta-gpt-diagnosis.md) |
| Run manifest | the diagnose step | carries model id / quantization / params / config hash |
| QC dataset (truth) + version | the `ciqa/*` API | [`datasets.md`](datasets.md) |
| Baseline (or `active`) | the `ciqa/*` API | tied to a dataset version |

The manifest — not a parsed directory name — attributes the run to a model/configuration.

## 3. Triage — find your scenario

| Situation | Go to |
|---|---|
| First time on this dataset version; no baseline yet | [RT-1](#rt-1--establish-a-baseline) |
| Routine check of a model/config change against the baseline | [RT-2](#rt-2--routine-regression-check) |
| The report says `regressed = true` | [RT-3](#rt-3--a-regressed-result) |
| McNemar is non-significant but you are not reassured | [RT-4](#rt-4--underpowered-or-unbalanced-panel) |
| Report errors that the dataset version does not match the baseline | [RT-5](#rt-5--dataset-version-mismatch) |
| A new configuration is better and you want to adopt it | [RT-6](#rt-6--promote-a-new-baseline) |

---

## RT-1 — Establish a baseline

You have a model/configuration you trust and a QC dataset version with truth, but no baseline yet.

1. Produce diagnoses over the dataset (prefetch → diagnose), giving outputs + a run manifest.
2. Score them against the dataset truth:

   ```bash
   cerebro-ciqa review --diagnostics diagnose/ --plate plate.json --output review.json
   ```

3. When the statistics are acceptable, **create and promote** a baseline for this dataset version,
   recording the manifest, the statistics, and the pass thresholds:

   ```bash
   # create the baseline record (via the API/CLI), then promote it for the dataset version
   #   ciqa/baseline           ← create from (manifest, statistics, thresholds)
   #   ciqa/baseline/promote   ← operator-gated, audited
   ```

**Verify.** The baseline appears for the dataset version, is marked promoted, and an audit entry
records who promoted it. A subsequent [RT-2](#rt-2--routine-regression-check) of the same run
against it reports `regressed = false`, `passed_threshold = true`.

---

## RT-2 — Routine regression check

You changed a model, quantization, threshold, or prompt and want to know if it regressed.

1. Produce diagnoses for the change (prefetch → diagnose) → outputs + manifest.
2. Run regression against the active baseline:

   ```bash
   cerebro-ciqa regression \
     --manifest diagnose/run.manifest.json \
     --baseline active \
     --report regression.report.json
   ```

3. Read the report: `current` vs `baseline` statistics, `delta_sensitivity` / `delta_specificity`,
   `mcnemar_p_value`, `passed_threshold`, `regressed`, and the per-sample decisions.

**Verify.** `regressed = false` **and** `passed_threshold = true` → the change is clear to proceed
under the gate. The exit code reflects the verdict for CI/Nextflow. If either fails, go to
[RT-3](#rt-3--a-regressed-result). The baseline is unchanged regardless.

---

## RT-3 — A regressed result

The report says `regressed = true` (floor breached, or a significant adverse shift).

1. **Do not promote anything.** The current baseline stands.
2. Localise it: read the **per-sample decisions** to see which samples flipped (e.g. new
   false-negatives on a sample type). The deltas and the discordant pairs in the McNemar table show
   whether the loss is broad or concentrated.
3. Decide the cause: a genuine model/config regression, a prefetch change (re-check the prefetch is
   built the same way — [`prefetch-filtering.md`](prefetch-filtering.md)), or a truth-set issue
   (re-check the dataset version's truth — [`datasets.md`](datasets.md)).
4. Fix and re-run [RT-2](#rt-2--routine-regression-check), or abandon the change.

**Verify.** Either a corrected run reports `regressed = false`, or the change is recorded as
rejected. The report is retained as the audit of the decision.

---

## RT-4 — Underpowered or unbalanced panel

McNemar is non-significant, or the panel is small, and you are not reassured.

1. Read the **panel size and discordant counts** the report states. A small panel, or few
   discordant pairs, has low statistical power: a non-significant McNemar is *not* evidence of no
   change.
2. Lean on the **absolute floor** as the primary gate — it does not depend on power — and on the
   per-sample decisions.
3. If the decision is high-stakes, **enlarge the QC dataset version** (more reference-truthed
   samples) and re-baseline before gating on significance ([`datasets.md`](datasets.md)).

**Verify.** The report's stated panel size and discordant counts are adequate for the claim you
are making, or the decision is explicitly made on the floor + per-sample evidence rather than on
the p-value.

---

## RT-5 — Dataset version mismatch

Regression errors (or warns) that the run's dataset version does not match the baseline's.

1. A baseline is valid only against the **dataset version it was measured on**. Comparing across
   versions is meaningless and is refused/flagged by design.
2. Either run against the baseline's dataset version, or establish a baseline for the current
   version ([RT-1](#rt-1--establish-a-baseline)).

**Verify.** The run and the baseline reference the same dataset version, and the report no longer
flags a mismatch.

---

## RT-6 — Promote a new baseline

A new configuration is genuinely better and you intend to adopt it as the standard.

1. Confirm it against the current baseline ([RT-2](#rt-2--routine-regression-check)) and inspect
   the report — better on the floor, no adverse significant shift, adequate panel.
2. **Promote** its statistics as the new baseline for the dataset version — an explicit,
   **operator-gated, audited** action; four-eyes is recommended for a clinical baseline.

   ```bash
   #   ciqa/baseline           ← create from the trusted run's (manifest, statistics, thresholds)
   #   ciqa/baseline/promote   ← gated; recorded in the audit chain
   ```

3. Keep the prior baseline — baselines are immutable and versioned, so historical comparisons
   remain reproducible.

**Verify.** The new baseline is marked active for the dataset version, the prior baseline still
exists, and the audit chain records the promotion (who, what, when). Future
[RT-2](#rt-2--routine-regression-check) runs are now judged against the new baseline.

---

## Tested here vs your environment

The evaluator (decisions → sensitivity/specificity), the McNemar test and its adjustment, and the
regression comparison are **pure, unit-tested logic** — the regulatory core, checked without a
stack and cross-checked so training and QC numbers agree. The `ciqa/*` persistence, the audit
entries for promotion, and the end-to-end report retrieval are **environment-validated** against a
running stack; see [`validation.md`](validation.md). Statistical power is a property of your panel,
not of the code — the report gives you the numbers to judge it.

> **Commands & endpoints.** For copy-pasteable CLI and curl for every step, see [Commands & API](commands.md).
