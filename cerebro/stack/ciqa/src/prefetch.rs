use cerebro_model::api::cerebro::schema::{MetaGpConfig, PrefetchData, TestResult};
use cerebro_pipeline::taxa::taxon::Taxon;
use serde::Serialize;
use std::path::{Path, PathBuf};
use cerebro_pipeline::taxa::taxon::LineageOperations;
use crate::utils::write_tsv;

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

/// One row of the taxon evidence table.
///
/// Taxa are emitted in long format: one row per profile-evidence record, so a
/// taxon detected by several tools/modes yields several rows. Taxa with no
/// profile evidence still produce a single row (profile columns left empty)
/// so alignment/assembly-only detections are not lost. `n_alignment` and
/// `n_assembly` report the number of alignment and assembly evidence records
/// attached to the taxon.
#[derive(Debug, Clone, Serialize)]
pub struct TaxonEvidenceRow {
    /// Source prefetch file stem
    pub file: String,
    /// Tier the taxon was drawn from: primary | secondary | target
    pub tier: String,
    /// Whether the taxon came from the prevalence-contamination set of that tier
    pub contamination: bool,
    pub taxid: String,
    pub rank: String,
    pub name: String,
    /// Superkingdom/domain parsed from the lineage (d__)
    pub domain: String,
    /// Full GTDB-like lineage string
    pub lineage: String,
    /// Profile tool for this evidence record (empty when the taxon has no profile evidence)
    pub tool: String,
    /// Abundance mode for this evidence record
    pub mode: String,
    pub reads: Option<u64>,
    pub rpm: Option<f64>,
    pub contigs: Option<u64>,
    pub bases: Option<u64>,
    pub bpm: Option<f64>,
    pub abundance: Option<f64>,
    pub n_alignment: usize,
    pub n_assembly: usize,
}

fn enum_to_str<T: Serialize>(value: &T) -> Result<String, serde_json::Error> {
    Ok(serde_json::to_string(value)?.trim_matches('"').to_string())
}

/// Append the evidence rows for a single taxon (one per profile record, or a
/// single profile-less row) to `rows`.
fn push_taxon_rows(
    rows: &mut Vec<TaxonEvidenceRow>,
    file: &str,
    tier: &str,
    contamination: bool,
    taxon: &Taxon,
) -> Result<(), serde_json::Error> {

    let rank = enum_to_str(&taxon.rank)?;
    let domain = taxon.lineage.get_domain().unwrap_or_default();
    let n_alignment = taxon.evidence.alignment.len();
    let n_assembly = taxon.evidence.assembly.len();

    if taxon.evidence.profile.is_empty() {
        rows.push(TaxonEvidenceRow {
            file: file.to_string(),
            tier: tier.to_string(),
            contamination,
            taxid: taxon.taxid.clone(),
            rank,
            name: taxon.name.clone(),
            domain,
            lineage: taxon.lineage.clone(),
            tool: String::new(),
            mode: String::new(),
            reads: None,
            rpm: None,
            contigs: None,
            bases: None,
            bpm: None,
            abundance: None,
            n_alignment,
            n_assembly,
        });
    } else {
        for rec in &taxon.evidence.profile {
            rows.push(TaxonEvidenceRow {
                file: file.to_string(),
                tier: tier.to_string(),
                contamination,
                taxid: taxon.taxid.clone(),
                rank: rank.clone(),
                name: taxon.name.clone(),
                domain: domain.clone(),
                lineage: taxon.lineage.clone(),
                tool: enum_to_str(&rec.tool)?,
                mode: enum_to_str(&rec.mode)?,
                reads: Some(rec.reads),
                rpm: Some(rec.rpm),
                contigs: Some(rec.contigs),
                bases: Some(rec.bases),
                bpm: Some(rec.bpm),
                abundance: Some(rec.abundance),
                n_alignment,
                n_assembly,
            });
        }
    }

    Ok(())
}

/// Parse one or more `PrefetchData` JSON files and write a TSV of taxon identity
/// (taxid, rank, name, domain, lineage) and all associated evidence, across all
/// six tiers (primary/secondary/target and their contamination sets).
pub fn write_taxon_evidence_table(files: &[PathBuf], output: &Path) -> anyhow::Result<()> {

    let mut rows: Vec<TaxonEvidenceRow> = Vec::new();

    for file in files {
        let data = PrefetchData::from_json(file)?;
        let stem = file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or_default()
            .to_string();

        let tiers: [(&str, bool, &Vec<Taxon>); 6] = [
            ("primary",   false, &data.primary),
            ("secondary", false, &data.secondary),
            ("target",    false, &data.target),
            ("primary",   true,  &data.primary_contamination),
            ("secondary", true,  &data.secondary_contamination),
            ("target",    true,  &data.target_contamination),
        ];

        for (tier, contamination, taxa) in tiers {
            for taxon in taxa {
                push_taxon_rows(&mut rows, &stem, tier, contamination, taxon)?;
            }
        }
    }

    write_tsv(&rows, output, true)?;

    Ok(())
}
