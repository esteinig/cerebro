//! Regression of a run's diagnostic statistics against a versioned baseline.
//!
//! Pure core: pair the current run's per-sample decisions with the baseline's, compute the
//! McNemar discordant counts, apply the absolute sensitivity/specificity floor, and assemble a
//! [`RegressionReport`]. The McNemar **p-value** is supplied by the caller (CIQA's `stats.rs`
//! retains the `statrs`-backed `mcnemar_test`, per the stage plan); this module owns the pairing
//! and the gate so they are unit-testable without a stats dependency.
//!
//! Report-first: this never mutates a baseline. Baseline promotion is a separate, audited action.

use serde::{Deserialize, Serialize};

use crate::api::diagnostics::eval::{Decision, DecisionRow, DiagnosticStatistics};

/// Correctness of a decision for paired (McNemar) comparison: TP/TN correct, FP/FN incorrect,
/// Excluded contributes no pair. Mirrors CIQA's `is_correct`.
pub fn is_correct(decision: Decision) -> Option<bool> {
    match decision {
        Decision::TP | Decision::TN => Some(true),
        Decision::FP | Decision::FN => Some(false),
        Decision::Excluded => None,
    }
}

/// McNemar discordant counts over samples present in **both** sets (joined by `sample_id`):
/// returns `(b, c, total)` where `b` = baseline-correct & current-incorrect, `c` =
/// baseline-incorrect & current-correct, `total` = number of paired (both-classified) samples.
/// Feed `(b, c, total)` to `mcnemar_test` for the p-value. Mirrors CIQA's `compute_bc_counts`.
pub fn discordant_counts(baseline: &[DecisionRow], current: &[DecisionRow]) -> (u32, u32, u32) {
    use std::collections::HashMap;
    let base: HashMap<&str, bool> = baseline
        .iter()
        .filter_map(|r| is_correct(r.decision).map(|ok| (r.sample_id.as_str(), ok)))
        .collect();
    let curr: HashMap<&str, bool> = current
        .iter()
        .filter_map(|r| is_correct(r.decision).map(|ok| (r.sample_id.as_str(), ok)))
        .collect();

    let (mut b, mut c, mut total) = (0u32, 0u32, 0u32);
    for (sample_id, &base_ok) in &base {
        if let Some(&curr_ok) = curr.get(sample_id) {
            total += 1;
            match (base_ok, curr_ok) {
                (true, false) => b += 1,
                (false, true) => c += 1,
                _ => {}
            }
        }
    }
    (b, c, total)
}

/// A regression evaluation comparing a run's META-GPT statistics to a baseline. Report-first:
/// it records the comparison and the gate; it never modifies the baseline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegressionReport {
    pub run_id: String,
    pub baseline_id: String,
    pub dataset: String,
    pub current: DiagnosticStatistics,
    pub baseline: DiagnosticStatistics,
    pub delta_sensitivity: f64,
    pub delta_specificity: f64,
    pub mcnemar_p_value: f64,
    pub mcnemar_b: u32,
    pub mcnemar_c: u32,
    pub mcnemar_total: u32,
    /// Absolute floor met: `current.sensitivity >= min` and `current.specificity >= min`.
    pub passed_threshold: bool,
    /// Floor breached OR a statistically significant adverse paired shift.
    pub regressed: bool,
    #[serde(default)]
    pub per_sample: Vec<DecisionRow>,
    pub created: String,
}

impl RegressionReport {
    /// Assemble the report and decide the gate. Thresholds are **fractions** (`0.0..=1.0`),
    /// matching [`DiagnosticStatistics`]. A run regresses when the absolute floor is breached, or
    /// when the paired McNemar test is significant (`p < alpha`) *and* the shift is adverse (a
    /// drop in sensitivity or specificity). Statistics-only (non-mutating).
    #[allow(clippy::too_many_arguments)]
    pub fn assemble(
        run_id: &str,
        baseline_id: &str,
        dataset: &str,
        current: DiagnosticStatistics,
        baseline: DiagnosticStatistics,
        mcnemar: (u32, u32, u32),
        mcnemar_p_value: f64,
        min_sensitivity: f64,
        min_specificity: f64,
        alpha: f64,
        created: &str,
    ) -> Self {
        let delta_sensitivity = current.sensitivity - baseline.sensitivity;
        let delta_specificity = current.specificity - baseline.specificity;

        let passed_threshold =
            current.sensitivity >= min_sensitivity && current.specificity >= min_specificity;

        let adverse = delta_sensitivity < 0.0 || delta_specificity < 0.0;
        let significant = mcnemar_p_value < alpha;
        let regressed = !passed_threshold || (significant && adverse);

        let per_sample = current.rows.clone();
        let (b, c, total) = mcnemar;

        Self {
            run_id: run_id.to_string(),
            baseline_id: baseline_id.to_string(),
            dataset: dataset.to_string(),
            current,
            baseline,
            delta_sensitivity,
            delta_specificity,
            mcnemar_p_value,
            mcnemar_b: b,
            mcnemar_c: c,
            mcnemar_total: total,
            passed_threshold,
            regressed,
            per_sample,
            created: created.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::api::cerebro::schema::TestResult;

    fn row(id: &str, decision: Decision) -> DecisionRow {
        DecisionRow {
            sample_id: id.to_string(),
            reference_result: Some(TestResult::Positive),
            predicted_result: Some(TestResult::Positive),
            decision,
            matched_candidate: None,
            exclude_reason: None,
        }
    }

    #[test]
    fn is_correct_mapping() {
        assert_eq!(is_correct(Decision::TP), Some(true));
        assert_eq!(is_correct(Decision::TN), Some(true));
        assert_eq!(is_correct(Decision::FP), Some(false));
        assert_eq!(is_correct(Decision::FN), Some(false));
        assert_eq!(is_correct(Decision::Excluded), None);
    }

    #[test]
    fn discordant_counts_pair_by_sample_id() {
        // s1: base correct, curr incorrect -> b; s2: base incorrect, curr correct -> c;
        // s3: both correct -> neither; s4: only in baseline -> unpaired; s5: excluded -> no pair.
        let baseline = vec![row("s1", Decision::TP), row("s2", Decision::FP), row("s3", Decision::TN), row("s4", Decision::TP), row("s5", Decision::Excluded)];
        let current = vec![row("s1", Decision::FN), row("s2", Decision::TN), row("s3", Decision::TN), row("s5", Decision::TP)];
        let (b, c, total) = discordant_counts(&baseline, &current);
        assert_eq!((b, c, total), (1, 1, 3)); // s1,s2,s3 paired; s5 excluded in baseline -> unpaired
    }

    fn stats(sens: f64, spec: f64) -> DiagnosticStatistics {
        // Minimal stats object with the sens/spec under test; counts/rows not used by assemble's gate.
        DiagnosticStatistics { sensitivity: sens, specificity: spec, ppv: 0.0, npv: 0.0, total: 0, tp: 0, tn: 0, fp: 0, fn_: 0, rows: vec![] }
    }

    #[test]
    fn floor_breach_regresses_regardless_of_significance() {
        let r = RegressionReport::assemble("run", "base", "ds", stats(0.70, 0.95), stats(0.90, 0.95), (0, 0, 10), 1.0, 0.80, 0.80, 0.05, "t");
        assert!(!r.passed_threshold); // 0.70 < 0.80
        assert!(r.regressed);
    }

    #[test]
    fn significant_adverse_shift_regresses_even_above_floor() {
        // Above floor on both axes, but a significant adverse drop in sensitivity.
        let r = RegressionReport::assemble("run", "base", "ds", stats(0.85, 0.95), stats(0.95, 0.95), (8, 1, 20), 0.01, 0.80, 0.80, 0.05, "t");
        assert!(r.passed_threshold);
        assert!(r.delta_sensitivity < 0.0);
        assert!(r.regressed); // significant (p<0.05) and adverse
    }

    #[test]
    fn matching_run_does_not_regress() {
        let r = RegressionReport::assemble("run", "base", "ds", stats(0.95, 0.95), stats(0.95, 0.95), (0, 0, 20), 1.0, 0.80, 0.80, 0.05, "t");
        assert!(r.passed_threshold);
        assert!(!r.regressed);
    }
}
