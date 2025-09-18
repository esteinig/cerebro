use chrono::{SecondsFormat, Utc};
use serde::{Deserialize, Serialize};
use uuid::Uuid;
use rand::{seq::SliceRandom, rng};
use crate::api::{cerebro::schema::{SampleType, TestResult}, training::{response::TrainingPrefetchData, schema::{CreateTrainingPrefetch, CreateTrainingSession, TrainingRecord}}};

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
}
impl TrainingPrefetchRecord {
    pub fn from_request(req: &CreateTrainingPrefetch) -> Self {
        Self {
            id: req.id.clone(),
            collection: req.collection.clone(),
            description: req.description.clone(),
            name: req.name.clone()
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
    pub fn from_request(session: &CreateTrainingSession, records: Vec<TrainingPrefetchData>, user_id: &str, user_name: &str) -> Self {

        let mut recs: Vec<TrainingRecord> = records
            .into_iter()
            .map(|r| TrainingRecord {
                id: Uuid::new_v4().to_string(),
                data_id: r.id,
                result: TestResult::Negative,
                candidates: None,
                sample_name: Some(r.prefetch.config.sample),
                sample_type: Some(r.prefetch.config.sample_type),
                reference_result: r.prefetch.config.test_result,
                reference_candidates: r.prefetch.config.candidates
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "UPPERCASE")]
pub enum Decision {
    TP,
    TN,
    FP,
    FN,
    Excluded,
}

/// Normalize a candidate:
/// - split on whitespace
/// - for each token, drop everything from the first '_' to the end
/// - rejoin with single spaces
/// - prepend "s__"
pub fn normalize_candidate(raw: &str) -> String {
    let core = raw
        .split_whitespace()
        .map(|t| t.split_once('_').map(|(head, _)| head).unwrap_or(t))
        .collect::<Vec<_>>()
        .join(" ");
    format!("s__{}", core)
}

impl TrainingSessionRecord {
    pub fn evaluate(&self) -> TrainingResult {
        let mut tp = 0usize;
        let mut tn = 0usize;
        let mut fp = 0usize;
        let mut fn_ = 0usize;
        let mut rows: Vec<TrainingResultRecord> = Vec::with_capacity(self.records.len());

        for r in &self.records {
            let mut exclude_reason: Option<String> = None;
            let mut matched_any_candidate: Option<bool> = None;

            let decision: Decision = match r.reference_result {
                None => {
                    exclude_reason = Some("missing reference_result".into());
                    Decision::Excluded
                }
                Some(TestResult::Positive) => match r.result {
                    TestResult::Positive => match (&r.candidates, &r.reference_candidates) {
                        (Some(cands), Some(ref_cands)) => {
                            let overlap = cands.iter().any(|c| ref_cands.iter().any(|rc| *rc == normalize_candidate(c)));
                            matched_any_candidate = Some(overlap);
                            if overlap { Decision::TP } else { Decision::FP }
                        }
                        _ => {
                            exclude_reason = Some(
                                "positive/positive with missing candidates".into(),
                            );
                            Decision::Excluded
                        }
                    },
                    TestResult::Negative => Decision::FN,
                },
                Some(TestResult::Negative) => match r.result {
                    TestResult::Positive => Decision::FP,
                    TestResult::Negative => Decision::TN,
                },
            };

            match decision {
                Decision::TP => tp += 1,
                Decision::TN => tn += 1,
                Decision::FP => fp += 1,
                Decision::FN => fn_ += 1,
                Decision::Excluded => {}
            }

            rows.push(TrainingResultRecord {
                record_id: r.id.clone(),
                data_id: r.data_id.clone(),
                sample_name: r.sample_name.clone(),
                sample_type: r.sample_type.clone(),
                result: r.result.clone(),
                reference_result: r.reference_result.clone(),
                candidates: r.candidates.as_ref().map(|c| c.join(";")),
                reference_candidates: r.reference_candidates.as_ref().map(|c| c.join(";")),
                matched_any_candidate,
                decision,
                exclude_reason,
            });
        }

        let total = tp + tn + fp + fn_;
        let sensitivity = if tp + fn_ > 0 {
            (tp as f64 / (tp + fn_) as f64)*100.
        } else {
            0.0
        };
        let specificity = if tn + fp > 0 {
            (tn as f64 / (tn + fp) as f64)*100.
        } else {
            0.0
        };

        TrainingResult {
            sensitivity,
            specificity,
            total,
            true_positive: tp,
            true_negative: tn,
            false_positive: fp,
            false_negative: fn_,
            data: rows,
        }
    }
}
