use crate::api::{
    cerebro::schema::{SampleType, TestResult},
    diagnostics::eval::{evaluate as evaluate_samples, EvaluatedSample},
    training::{
        response::TrainingPrefetchData,
        schema::{CreateTrainingPrefetch, CreateTrainingSession, TrainingRecord},
    },
};
use chrono::{SecondsFormat, Utc};
use rand::{rng, seq::SliceRandom};
use serde::{Deserialize, Serialize};
use uuid::Uuid;

// `Decision` and `normalize_candidate` are now defined once in the shared evaluator
// (api::diagnostics::eval). Re-export them so existing `training::model::{Decision,
// normalize_candidate}` references keep working.
pub use crate::api::diagnostics::eval::{normalize_candidate, Decision};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrainingPrefetchRecord {
    // Unique identifier for retrieving prefetch data from GridFs
    pub id: String,
    /// Name of a training collection
    pub collection: String,
    /// Description of the training collection
    pub description: String,
    /// Human readable label
    pub name: String,
    /// Preselect reference organism
    pub preselect: Option<bool>,
}
impl TrainingPrefetchRecord {
    pub fn from_request(req: &CreateTrainingPrefetch) -> Self {
        Self {
            id: req.id.clone(),
            collection: req.collection.clone(),
            description: req.description.clone(),
            name: req.name.clone(),
            preselect: req.preselect,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrainingSessionRecord {
    /// Unique identifier for the training session
    pub id: String,
    /// Collection name for this training session
    pub collection: String,
    /// User display name
    pub user_name: String,
    /// User unique identifier
    pub user_id: String,
    /// Timestamp when training started
    pub started: String,
    /// Timestamp when training completed (if any)
    pub completed: Option<String>,
    /// Training result data
    pub result: Option<TrainingResult>,
    /// Last training record updated identifier
    pub last_updated: Option<String>,
    /// Number of records in the session
    pub records: Vec<TrainingRecord>,
}

impl TrainingSessionRecord {
    pub fn from_request(
        session: &CreateTrainingSession,
        records: Vec<TrainingPrefetchData>,
        user_id: &str,
        user_name: &str,
    ) -> Self {
        let mut recs: Vec<TrainingRecord> = records
            .into_iter()
            .map(|r| {
                // preselect enabled only if explicitly true
                let preselect_enabled = r.preselect == Some(true);

                let first_candidate = r.prefetch.config.candidates.as_ref().and_then(|v| {
                    // First try to find the first "s__" candidate
                    v.iter()
                        .find(|s| {
                            let s = s.trim();
                            !s.is_empty() && s.starts_with("s__")
                        })
                        // Fallback to the first non-empty candidate
                        .or_else(|| v.iter().find(|s| !s.trim().is_empty()))
                        .cloned()
                });

                // seed user-facing selection if preselect enabled and candidate exists
                let (result, candidates) = if preselect_enabled {
                    if let Some(c) = first_candidate.clone() {
                        // candidate must be trimmed from genus/species prefix that comes with the reference species designation
                        let c_trim = c
                            .trim_start_matches("s__")
                            .trim_start_matches("g__")
                            .to_string();
                        (
                            r.prefetch
                                .config
                                .test_result
                                .clone()
                                .unwrap_or(TestResult::Positive),
                            Some(vec![c_trim]),
                        )
                    } else {
                        (
                            r.prefetch
                                .config
                                .test_result
                                .clone()
                                .unwrap_or(TestResult::Negative),
                            None,
                        )
                    }
                } else {
                    (
                        r.prefetch
                            .config
                            .test_result
                            .clone()
                            .unwrap_or(TestResult::Negative),
                        None,
                    )
                };

                TrainingRecord {
                    id: Uuid::new_v4().to_string(),
                    data_id: r.id,
                    result: result,
                    candidates: candidates,
                    sample_name: Some(r.prefetch.config.sample),
                    sample_type: Some(r.prefetch.config.sample_type),
                    reference_result: r.prefetch.config.test_result,
                    reference_candidates: r.prefetch.config.candidates,
                    exclude_lod: r.prefetch.config.exclude_lod,
                }
            })
            .collect();

        if session.shuffle {
            recs.shuffle(&mut rng());
        }

        Self {
            id: Uuid::new_v4().to_string(),
            collection: session.collection.clone(),
            user_name: user_name.to_string(),
            user_id: user_id.to_string(),
            started: Utc::now().to_rfc3339_opts(SecondsFormat::Secs, true),
            completed: None,
            result: None,
            last_updated: None,
            records: recs,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrainingResult {
    pub sensitivity: f64,
    pub specificity: f64,
    pub total: usize,
    pub true_positive: usize,
    pub true_negative: usize,
    pub false_positive: usize,
    pub false_negative: usize,
    pub data: Vec<TrainingResultRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrainingResultRecord {
    pub record_id: String,
    pub data_id: String,
    pub sample_name: Option<String>,
    pub sample_type: Option<SampleType>,
    pub result: TestResult,
    pub reference_result: Option<TestResult>,
    pub candidates: Option<String>,
    pub reference_candidates: Option<String>,
    /// Only meaningful when both ref/result are Positive
    pub matched_any_candidate: Option<bool>,
    pub decision: Decision,
    pub exclude_reason: Option<String>,
}

impl TrainingSessionRecord {
    pub fn evaluate(&self) -> TrainingResult {
        // Delegate the decision rule and the sens/spec maths to the single shared evaluator.
        // Training records always carry a predicted `result` (so every sample goes through the
        // Some(result) arms), and the candidate-overlap rule is the shared `classify` rule, so
        // the decisions, counts, and statistics are identical to the previous in-place version.
        // Only the repackaging into the training-facing `TrainingResult`/`TrainingResultRecord`
        // (percentages, joined candidate strings) lives here.
        let samples: Vec<EvaluatedSample> = self
            .records
            .iter()
            .map(|r| EvaluatedSample {
                sample_id: r.id.clone(),
                reference_result: r.reference_result.clone(),
                reference_candidates: r.reference_candidates.clone(),
                predicted_result: Some(r.result.clone()),
                predicted_candidates: r.candidates.clone(),
                exclude_lod: r.exclude_lod,
            })
            .collect();

        let stats = evaluate_samples(&samples);

        // `stats.rows` is 1:1 with `self.records` in the same order.
        let data: Vec<TrainingResultRecord> = self
            .records
            .iter()
            .zip(stats.rows.iter())
            .map(|(r, row)| TrainingResultRecord {
                record_id: r.id.clone(),
                data_id: r.data_id.clone(),
                sample_name: r.sample_name.clone(),
                sample_type: r.sample_type.clone(),
                result: r.result.clone(),
                reference_result: r.reference_result.clone(),
                candidates: r.candidates.as_ref().map(|c| c.join(";")),
                reference_candidates: r.reference_candidates.as_ref().map(|c| c.join(";")),
                matched_any_candidate: row.matched_candidate,
                decision: row.decision,
                exclude_reason: row.exclude_reason.clone(),
            })
            .collect();

        TrainingResult {
            // Training presents sensitivity/specificity as percentages; the shared evaluator
            // returns fractions, so scale by 100 to preserve the existing UI numbers.
            sensitivity: stats.sensitivity * 100.,
            specificity: stats.specificity * 100.,
            total: stats.total,
            true_positive: stats.tp,
            true_negative: stats.tn,
            false_positive: stats.fp,
            false_negative: stats.fn_,
            data,
        }
    }
}
