use cerebro_model::api::cerebro::schema::{MetaGpConfig, PrefetchData, TestResult};
use cerebro_pipeline::taxa::taxon::Taxon;
use serde::Serialize;
use cerebro_pipeline::taxa::taxon::LineageOperations;

#[derive(Debug, Clone, Serialize)]
pub struct DomainTally {
    pub prokaryotes: usize,
    pub eukaryotes: usize,
    pub viruses: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct CountsByCategory {
    pub primary: DomainTally,
    pub secondary: DomainTally,
    pub target: DomainTally,
}

#[derive(Debug, Clone, Serialize)]
pub struct PerSampleSummary {
    pub sample: String,
    pub counts: CountsByCategory,
    pub contamination_counts: CountsByCategory,
    pub positive: bool,
    pub positive_match: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct OverallSummary {
    pub total_samples: usize,
    pub total_positive: usize,
    pub total_positive_with_candidate_match: usize,
    pub per_sample: Vec<PerSampleSummary>,
}

pub fn tally_domain(v: &[Taxon]) -> DomainTally {
    let mut p = 0usize;
    let mut e = 0usize;
    let mut vrs = 0usize;

    for t in v {
        let dom = t.lineage.get_domain().unwrap_or_default();
        match dom.as_str() {
            "Bacteria" | "Archaea" => p += 1,
            "Eukaryota" => e += 1,
            "Viruses" => vrs += 1,
            _ => {}
        }
    }
    DomainTally { prokaryotes: p, eukaryotes: e, viruses: vrs }
}

pub fn counts_by_category(d: &PrefetchData) -> CountsByCategory {
    CountsByCategory {
        primary:   tally_domain(&d.primary),
        secondary: tally_domain(&d.secondary),
        target:    tally_domain(&d.target),
    }
}

pub fn counts_by_category_contam(d: &PrefetchData) -> CountsByCategory {
    CountsByCategory {
        primary:   tally_domain(&d.primary_contamination),
        secondary: tally_domain(&d.secondary_contamination),
        target:    tally_domain(&d.target_contamination),
    }
}

pub fn positive_candidate_match(d: &PrefetchData) -> (bool, bool) {
    let positive = matches!(d.config.test_result, Some(TestResult::Positive));
    if !positive {
        return (false, false);
    }

    let labels_in_any = |taxa: &[Taxon]| -> std::collections::HashSet<String> {
        taxa.iter()
            .flat_map(|t| t.lineage.get_labels().into_iter().map(|s| s.to_string()))
            .collect()
    };

    let set_primary   = labels_in_any(&d.primary);
    let set_secondary = labels_in_any(&d.secondary);
    let set_target    = labels_in_any(&d.target);

    let mut matched = false;
    match &d.config.candidates {
        None => {
            log::warn!(
                "Positive sample '{}' has no candidates",
                d.config.sample
            );
        }
        Some(cands) if cands.is_empty() => {
            log::warn!(
                "Positive sample '{}' has empty candidate list",
                d.config.sample
            );
        }
        Some(cands) => {
            'outer: for c in cands {
                if set_primary.contains(c) || set_secondary.contains(c) || set_target.contains(c) {
                    matched = true;
                    break 'outer;
                }
            }
        }
    }

    (true, matched)
}


#[derive(Debug, Clone, Serialize)]
pub enum PrefetchStatus {
    NotDetected
}

#[derive(Debug, Clone, Serialize)]
pub struct MissedDetectionRow {
    pub sample: String,
    pub reference: String,                 
    pub status: PrefetchStatus,            
}

pub fn lineage_label_set(taxa: &[Taxon]) -> std::collections::HashSet<String> {
    taxa.iter()
        .flat_map(|t| t.lineage.get_labels().into_iter().map(|s| s.to_string()))
        .collect()
}

/// True if sample is Positive AND none of the candidate labels match any lineage label
/// in primary|secondary|target. Also warns if candidates missing/empty.
pub fn is_missed_detection(d: &PrefetchData) -> bool {
    if !matches!(d.config.test_result, Some(TestResult::Positive)) {
        return false;
    }
    let cands = match &d.config.candidates {
        None => {
            log::warn!(
                "Positive sample '{}' has no candidates provided",
                d.config.sample
            );
            return true;
        }
        Some(v) if v.is_empty() => {
            log::warn!(
                "Positive sample '{}' has empty candidate list",
                d.config.sample
            );
            return true;
        }
        Some(v) => v,
    };

    let s_primary   = lineage_label_set(&d.primary);
    let s_secondary = lineage_label_set(&d.secondary);
    let s_target    = lineage_label_set(&d.target);

    !cands.iter().any(|c|
        s_primary.contains(c) || s_secondary.contains(c) || s_target.contains(c)
    )
}

pub fn reference_names_from_config(cfg: &MetaGpConfig) -> String {
    match &cfg.candidates {
        Some(v) if !v.is_empty() => v.join(";"),
        _ => String::new(),
    }
}
