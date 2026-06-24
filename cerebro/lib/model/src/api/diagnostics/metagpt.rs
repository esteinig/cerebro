//! Reading META-GPT diagnostic output without depending on the GPU-heavy `meta-gpt` crate.
//!
//! `cerebro-model` sits below `meta-gpt` in the layering, so it cannot (and should not) depend on
//! it. Instead this module mirrors `meta_gpt::gpt::DiagnosticResult` as a serde DTO that reads the
//! exact on-disk `{sample}.model.json` shape, and maps it to the evaluator's prediction. The
//! on-disk format is unchanged — this only reads it.
//!
//! Keep [`DiagnosisDto`]/[`DiagnosticResultDto`] in lockstep with `meta_gpt::gpt` (a round-trip
//! test guards the JSON shape). The mapping reproduces CIQA's `SampleReview::from_gpt`: the
//! `diagnosis` → `TestResult` rule, and the pathogen string as a (raw) predicted candidate that
//! the evaluator normalises with `normalize_candidate` (the single normaliser — see eval.rs).

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::api::cerebro::schema::TestResult;
use crate::api::diagnostics::eval::EvaluatedSample;
use crate::api::cerebro::model::ModelError;

/// Mirror of `meta_gpt::gpt::Diagnosis` (default serde, i.e. verbatim variant names).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum DiagnosisDto {
    Infectious,
    InfectiousReview,
    NonInfectious,
    NonInfectiousReview,
    Tumor,
    Unknown,
}

/// Mirror of `meta_gpt::gpt::DiagnosticResult`. New/absent fields are defaulted so a forward- or
/// backward-rev META-GPT output still reads.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiagnosticResultDto {
    pub diagnosis: DiagnosisDto,
    #[serde(default)]
    pub candidates: Vec<String>,
    #[serde(default)]
    pub pathogen: Option<String>,
}

impl DiagnosticResultDto {
    pub fn from_json<P: AsRef<Path>>(path: P) -> Result<Self, ModelError> {
        let data = std::fs::read_to_string(path)?;
        serde_json::from_str(&data).map_err(ModelError::JsonDeserialization)
    }
}

/// Map a META-GPT result to a `(predicted_result, predicted_candidates)` prediction.
///
/// Mirrors `SampleReview::from_gpt`: `Infectious`/`InfectiousReview` → `Positive`,
/// `NonInfectious`/`NonInfectiousReview` → `Negative`, `Tumor`/`Unknown` → `None`. The pathogen
/// string is returned **raw** (un-normalised) as the single predicted candidate; the evaluator
/// applies `normalize_candidate` so prediction and reference truth are compared on one footing.
pub fn predicted_from_result(
    result: &DiagnosticResultDto,
) -> (Option<TestResult>, Option<Vec<String>>) {
    let predicted_result = match result.diagnosis {
        DiagnosisDto::Infectious | DiagnosisDto::InfectiousReview => Some(TestResult::Positive),
        DiagnosisDto::NonInfectious | DiagnosisDto::NonInfectiousReview => {
            Some(TestResult::Negative)
        }
        DiagnosisDto::Tumor | DiagnosisDto::Unknown => None,
    };

    let predicted_candidates = result
        .pathogen
        .as_ref()
        .map(|p| vec![p.clone()])
        .filter(|v: &Vec<String>| !v.is_empty());

    (predicted_result, predicted_candidates)
}

/// Build an [`EvaluatedSample`] by joining a META-GPT result to a sample's reference truth.
pub fn evaluated_sample(
    sample_id: &str,
    reference_result: Option<TestResult>,
    reference_candidates: Option<Vec<String>>,
    exclude_lod: Option<bool>,
    result: &DiagnosticResultDto,
) -> EvaluatedSample {
    let (predicted_result, predicted_candidates) = predicted_from_result(result);
    EvaluatedSample {
        sample_id: sample_id.to_string(),
        reference_result,
        reference_candidates,
        predicted_result,
        predicted_candidates,
        exclude_lod,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dto_round_trips_metagpt_json_shape() {
        // The exact on-disk shape written by meta-gpt (verbatim enum variant names).
        let json = r#"{"diagnosis":"InfectiousReview","candidates":["Escherichia coli"],"pathogen":"Escherichia coli"}"#;
        let dto: DiagnosticResultDto = serde_json::from_str(json).unwrap();
        assert_eq!(dto.diagnosis, DiagnosisDto::InfectiousReview);
        assert_eq!(dto.pathogen.as_deref(), Some("Escherichia coli"));
        // re-serialise and re-parse: stable
        let again: DiagnosticResultDto =
            serde_json::from_str(&serde_json::to_string(&dto).unwrap()).unwrap();
        assert_eq!(again.diagnosis, DiagnosisDto::InfectiousReview);
    }

    #[test]
    fn diagnosis_maps_to_result() {
        let mk = |d: DiagnosisDto, p: Option<&str>| DiagnosticResultDto {
            diagnosis: d,
            candidates: vec![],
            pathogen: p.map(String::from),
        };
        assert_eq!(predicted_from_result(&mk(DiagnosisDto::Infectious, Some("E coli"))).0, Some(TestResult::Positive));
        assert_eq!(predicted_from_result(&mk(DiagnosisDto::InfectiousReview, None)).0, Some(TestResult::Positive));
        assert_eq!(predicted_from_result(&mk(DiagnosisDto::NonInfectious, None)).0, Some(TestResult::Negative));
        assert_eq!(predicted_from_result(&mk(DiagnosisDto::NonInfectiousReview, None)).0, Some(TestResult::Negative));
        assert_eq!(predicted_from_result(&mk(DiagnosisDto::Tumor, None)).0, None);
        assert_eq!(predicted_from_result(&mk(DiagnosisDto::Unknown, None)).0, None);
        assert_eq!(predicted_from_result(&mk(DiagnosisDto::Infectious, Some("E coli"))).1, Some(vec!["E coli".to_string()]));
        assert_eq!(predicted_from_result(&mk(DiagnosisDto::Infectious, None)).1, None);
    }
}
