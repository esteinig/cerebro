//! The single, shared diagnostic evaluator.
//!
//! Sensitivity/specificity scoring previously existed twice — `TrainingSessionRecord::evaluate`
//! (training) and CIQA's `compute_statistics` — with the same counting and the same sens/spec
//! maths duplicated. This module is the **one** home for:
//!
//! - the canonical [`Decision`] (TP/TN/FP/FN/Excluded),
//! - the per-sample decision rule for the common case ([`classify`]), and
//! - the counting → sensitivity/specificity/PPV/NPV statistics ([`DiagnosticStatistics`]).
//!
//! Training delegates its evaluator here (its rule *is* the common rule, so its numbers are
//! preserved); CIQA keeps its richer clinical comparison (orthogonal-test matching, genus
//! dereference) but feeds the resulting decisions into the same counting/statistics, so both
//! report identical numbers from one implementation. The new regression path (Stage 3) uses
//! [`classify`]/[`evaluate`] directly.
//!
//! Everything here is pure and unit-tested. Statistics are stored as **fractions in `0.0..=1.0`**
//! (the form CIQA already uses internally); callers that present percentages multiply by 100.

use serde::{Deserialize, Serialize};

use crate::api::cerebro::schema::TestResult;

/// Canonical decision label. (Moved here from `training::model`; training re-exports it.)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "UPPERCASE")]
pub enum Decision {
    TP,
    TN,
    FP,
    FN,
    Excluded,
}

/// Normalize a candidate species string to the canonical `s__Genus species` truth form:
/// split on whitespace, drop everything from the first `_` in each token, rejoin, prepend `s__`.
///
/// This is the single normaliser used across training and regression so predicted candidates and
/// reference truth are compared on the same footing. (Moved here from `training::model`.)
pub fn normalize_candidate(raw: &str) -> String {
    let core = raw
        .split_whitespace()
        .map(|t| t.split_once('_').map(|(head, _)| head).unwrap_or(t))
        .collect::<Vec<_>>()
        .join(" ");
    format!("s__{}", core)
}

/// One sample's inputs to the evaluator. Reference truth + a prediction (from a human label or
/// from META-GPT, via `metagpt::predicted_from_result`).
#[derive(Debug, Clone)]
pub struct EvaluatedSample {
    pub sample_id: String,
    pub reference_result: Option<TestResult>,
    /// Normalised `s__...` reference truth candidates.
    pub reference_candidates: Option<Vec<String>>,
    pub predicted_result: Option<TestResult>,
    /// Raw predicted candidates; normalised via [`normalize_candidate`] at comparison time.
    pub predicted_candidates: Option<Vec<String>>,
    pub exclude_lod: Option<bool>,
}

/// The per-sample outcome the evaluator emits (and the report renders).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DecisionRow {
    pub sample_id: String,
    pub reference_result: Option<TestResult>,
    pub predicted_result: Option<TestResult>,
    pub decision: Decision,
    /// `Some(true/false)` only for a positive/positive comparison (candidate overlap).
    pub matched_candidate: Option<bool>,
    pub exclude_reason: Option<String>,
}

/// Counts + sensitivity/specificity/PPV/NPV (fractions) + the per-sample rows.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiagnosticStatistics {
    pub sensitivity: f64,
    pub specificity: f64,
    pub ppv: f64,
    pub npv: f64,
    pub total: usize,
    pub tp: usize,
    pub tn: usize,
    pub fp: usize,
    pub fn_: usize,
    pub rows: Vec<DecisionRow>,
}

impl DiagnosticStatistics {
    /// Fold decision rows into counts and statistics. `total` excludes `Excluded` rows — the
    /// same convention as both prior evaluators (training's `tp+tn+fp+fn`, CIQA's filtered len).
    pub fn from_rows(rows: Vec<DecisionRow>) -> Self {
        let (mut tp, mut tn, mut fp, mut fn_) = (0usize, 0usize, 0usize, 0usize);
        for row in &rows {
            match row.decision {
                Decision::TP => tp += 1,
                Decision::TN => tn += 1,
                Decision::FP => fp += 1,
                Decision::FN => fn_ += 1,
                Decision::Excluded => {}
            }
        }
        let frac = |num: usize, den: usize| if den > 0 { num as f64 / den as f64 } else { 0.0 };
        Self {
            sensitivity: frac(tp, tp + fn_),
            specificity: frac(tn, tn + fp),
            ppv: frac(tp, tp + fp),
            npv: frac(tn, tn + fn_),
            total: tp + tn + fp + fn_,
            tp,
            tn,
            fp,
            fn_,
            rows,
        }
    }
}

/// Classify one sample. Reproduces the training rule exactly for the cases training can produce,
/// and defines the `predicted_result == None` cases (META-GPT `Tumor`/`Unknown`) as "not an
/// infectious call" — `FN` against a positive reference, `TN` against a negative one — so such
/// samples stay in the denominator rather than being silently dropped.
pub fn classify(sample: &EvaluatedSample) -> DecisionRow {
    let mut exclude_reason: Option<String> = None;
    let mut matched_candidate: Option<bool> = None;

    let mut decision = match sample.reference_result {
        None => {
            exclude_reason = Some("missing reference_result".into());
            Decision::Excluded
        }
        Some(TestResult::Positive) => match sample.predicted_result {
            Some(TestResult::Positive) => {
                match (&sample.predicted_candidates, &sample.reference_candidates) {
                    (Some(cands), Some(ref_cands)) => {
                        let overlap = cands
                            .iter()
                            .any(|c| ref_cands.iter().any(|rc| *rc == normalize_candidate(c)));
                        matched_candidate = Some(overlap);
                        if overlap {
                            Decision::TP
                        } else {
                            Decision::FP
                        }
                    }
                    _ => {
                        exclude_reason = Some("positive/positive with missing candidates".into());
                        Decision::Excluded
                    }
                }
            }
            Some(TestResult::Negative) => Decision::FN,
            None => Decision::FN,
        },
        Some(TestResult::Negative) => match sample.predicted_result {
            Some(TestResult::Positive) => Decision::FP,
            Some(TestResult::Negative) => Decision::TN,
            None => Decision::TN,
        },
    };

    // LOD exclusion overrides any decision (matches training).
    if sample.exclude_lod.is_some_and(|x| x) {
        exclude_reason = Some("limit of detection".into());
        decision = Decision::Excluded;
    }

    DecisionRow {
        sample_id: sample.sample_id.clone(),
        reference_result: sample.reference_result.clone(),
        predicted_result: sample.predicted_result.clone(),
        decision,
        matched_candidate,
        exclude_reason,
    }
}

/// Evaluate a set of samples: classify each, then fold into statistics. The one implementation.
pub fn evaluate(samples: &[EvaluatedSample]) -> DiagnosticStatistics {
    DiagnosticStatistics::from_rows(samples.iter().map(classify).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(
        id: &str,
        reference: Option<TestResult>,
        ref_cands: Option<&[&str]>,
        predicted: Option<TestResult>,
        pred_cands: Option<&[&str]>,
        lod: Option<bool>,
    ) -> EvaluatedSample {
        EvaluatedSample {
            sample_id: id.to_string(),
            reference_result: reference,
            reference_candidates: ref_cands.map(|c| c.iter().map(|s| s.to_string()).collect()),
            predicted_result: predicted,
            predicted_candidates: pred_cands.map(|c| c.iter().map(|s| s.to_string()).collect()),
            exclude_lod: lod,
        }
    }

    #[test]
    fn normalize_candidate_matches_training_scheme() {
        assert_eq!(normalize_candidate("Escherichia coli"), "s__Escherichia coli");
        assert_eq!(normalize_candidate("Escherichia_X coli_Y"), "s__Escherichia coli");
    }

    #[test]
    fn positive_positive_overlap_is_tp_else_fp() {
        let tp = classify(&sample("a", Some(TestResult::Positive), Some(&["s__Escherichia coli"]),
            Some(TestResult::Positive), Some(&["Escherichia coli"]), None));
        assert_eq!(tp.decision, Decision::TP);
        assert_eq!(tp.matched_candidate, Some(true));

        let fp = classify(&sample("b", Some(TestResult::Positive), Some(&["s__Escherichia coli"]),
            Some(TestResult::Positive), Some(&["Staphylococcus aureus"]), None));
        assert_eq!(fp.decision, Decision::FP);
        assert_eq!(fp.matched_candidate, Some(false));
    }

    #[test]
    fn missing_candidates_on_positive_positive_excludes() {
        let r = classify(&sample("a", Some(TestResult::Positive), None,
            Some(TestResult::Positive), None, None));
        assert_eq!(r.decision, Decision::Excluded);
    }

    #[test]
    fn fn_fp_tn_arms() {
        assert_eq!(classify(&sample("a", Some(TestResult::Positive), None, Some(TestResult::Negative), None, None)).decision, Decision::FN);
        assert_eq!(classify(&sample("b", Some(TestResult::Negative), None, Some(TestResult::Positive), None, None)).decision, Decision::FP);
        assert_eq!(classify(&sample("c", Some(TestResult::Negative), None, Some(TestResult::Negative), None, None)).decision, Decision::TN);
    }

    #[test]
    fn none_prediction_is_fn_for_positive_tn_for_negative() {
        assert_eq!(classify(&sample("a", Some(TestResult::Positive), None, None, None, None)).decision, Decision::FN);
        assert_eq!(classify(&sample("b", Some(TestResult::Negative), None, None, None, None)).decision, Decision::TN);
    }

    #[test]
    fn missing_reference_and_lod_exclude() {
        assert_eq!(classify(&sample("a", None, None, Some(TestResult::Positive), None, None)).decision, Decision::Excluded);
        let lod = classify(&sample("b", Some(TestResult::Positive), Some(&["s__Escherichia coli"]),
            Some(TestResult::Positive), Some(&["Escherichia coli"]), Some(true)));
        assert_eq!(lod.decision, Decision::Excluded);
        assert_eq!(lod.exclude_reason.as_deref(), Some("limit of detection"));
    }

    #[test]
    fn statistics_are_fractions_and_exclude_excluded() {
        let samples = vec![
            sample("tp", Some(TestResult::Positive), Some(&["s__E coli"]), Some(TestResult::Positive), Some(&["E coli"]), None),
            sample("fn", Some(TestResult::Positive), None, Some(TestResult::Negative), None, None),
            sample("tn", Some(TestResult::Negative), None, Some(TestResult::Negative), None, None),
            sample("fp", Some(TestResult::Negative), None, Some(TestResult::Positive), None, None),
            sample("ex", None, None, Some(TestResult::Positive), None, None),
        ];
        let stats = evaluate(&samples);
        assert_eq!((stats.tp, stats.fn_, stats.tn, stats.fp), (1, 1, 1, 1));
        assert_eq!(stats.total, 4); // excluded not counted
        assert!((stats.sensitivity - 0.5).abs() < 1e-9);
        assert!((stats.specificity - 0.5).abs() < 1e-9);
        assert!((stats.ppv - 0.5).abs() < 1e-9);
        assert!((stats.npv - 0.5).abs() < 1e-9);
    }
}
