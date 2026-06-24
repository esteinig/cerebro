//! Prevalence-based contamination identification.
//!
//! Reagent and environmental contaminants show up across many samples. Prevalence filtering
//! flags a taxon as a contaminant when it is "present" (above an RPM floor) in at least a
//! `threshold` fraction of the samples considered.
//!
//! This is the taxon-only core of that computation. It previously lived in the API server
//! (`apply_prevalence_contamination_filter`, which operated over `Cerebro` documents); lifting
//! it here — the base layer, where `Taxon` is defined and below `cerebro-model` — lets both the
//! server *and* the self-contained (local) prefetch builder identify contaminants with the
//! **same** logic. The server keeps a thin adapter that maps `Cerebro -> taxa` and delegates to
//! this function, so the two paths cannot drift.
//!
//! The function is pure (no I/O, no database, no network) and is unit-tested below.

use std::collections::{HashMap, HashSet};

use crate::taxa::taxon::Taxon;

/// Identify contaminant taxids by cross-sample prevalence.
///
/// A taxon counts as **present** in a sample when any of its profile records has
/// `rpm >= min_rpm`. A taxon is a **contaminant** when it is present in at least `threshold`
/// (a fraction in `0.0..=1.0`) of the supplied samples. When `taxid_subset` is non-empty, only
/// those taxids are considered (all others are ignored).
///
/// Returns the contaminant taxids. Order is unspecified (callers collect into a set).
///
/// Behaviour is identical to the previous server-resident implementation; see the module docs.
/// Each element of `per_sample_taxa` is one sample's (already aggregated) taxa — for a stack
/// query these are the collection's `Cerebro` documents' taxa; for a local run they are the
/// run's per-library model taxa, grouped by nucleic-acid tag by the caller.
pub fn prevalence_contaminant_taxids(
    per_sample_taxa: &[Vec<Taxon>],
    min_rpm: f64,
    threshold: f64,
    taxid_subset: &[String],
) -> Vec<String> {
    // No samples -> no contaminants (avoids a divide-by-zero on `total`).
    if per_sample_taxa.is_empty() {
        return Vec::new();
    }

    // Borrowed lookup set; empty means "consider every taxon".
    let subset: HashSet<&str> = taxid_subset.iter().map(String::as_str).collect();

    let total = per_sample_taxa.len() as f64;
    let mut prevalence: HashMap<String, usize> = HashMap::new();

    for taxa in per_sample_taxa {
        // Count each taxon at most once per sample.
        let mut found: HashSet<&str> = HashSet::new();
        for taxon in taxa {
            if !subset.is_empty() && !subset.contains(taxon.taxid.as_str()) {
                continue;
            }
            // "Present" = any read-profile record at or above the RPM floor.
            if taxon.evidence.profile.iter().any(|record| record.rpm >= min_rpm) {
                found.insert(taxon.taxid.as_str());
            }
        }
        for taxid in found {
            *prevalence.entry(taxid.to_string()).or_insert(0) += 1;
        }
    }

    prevalence
        .into_iter()
        .filter_map(|(taxid, count)| {
            if (count as f64) / total >= threshold {
                Some(taxid)
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modules::pathogen::{AbundanceMode, ProfileRecord, ProfileTool};
    use crate::taxa::taxon::{Taxon, TaxonEvidence};
    use taxonomy::TaxRank;

    /// Minimal taxon carrying a single profile record at the given RPM.
    fn taxon(taxid: &str, rpm: f64) -> Taxon {
        Taxon {
            taxid: taxid.to_string(),
            rank: TaxRank::Species,
            name: format!("taxon_{taxid}"),
            lineage: String::new(),
            evidence: TaxonEvidence {
                alignment: Vec::new(),
                assembly: Vec::new(),
                profile: vec![ProfileRecord {
                    id: "lib".to_string(),
                    tool: ProfileTool::Kraken2,
                    mode: AbundanceMode::Profile,
                    reads: 1,
                    rpm,
                    contigs: 0,
                    bases: 0,
                    bpm: 0.0,
                    abundance: 0.0,
                }],
            },
        }
    }

    fn sorted(mut v: Vec<String>) -> Vec<String> {
        v.sort();
        v
    }

    #[test]
    fn empty_input_yields_no_contaminants() {
        assert!(prevalence_contaminant_taxids(&[], 0.0, 0.5, &[]).is_empty());
    }

    #[test]
    fn taxon_in_every_sample_is_contaminant_at_any_reachable_threshold() {
        let samples = vec![vec![taxon("1", 10.0)], vec![taxon("1", 10.0)], vec![taxon("1", 10.0)]];
        assert_eq!(prevalence_contaminant_taxids(&samples, 1.0, 1.0, &[]), vec!["1".to_string()]);
    }

    #[test]
    fn threshold_is_inclusive_and_fractional() {
        // taxid "1" in 1/4 (0.25), taxid "9" in 1/4 (0.25), taxid "2" in 3/4 (0.75)
        let samples = vec![
            vec![taxon("1", 10.0), taxon("2", 10.0)],
            vec![taxon("2", 10.0)],
            vec![taxon("2", 10.0)],
            vec![taxon("9", 10.0)],
        ];
        // 0.6 -> only "2"; 0.25 -> "1", "2", "9" (all >= 0.25, inclusive)
        assert_eq!(prevalence_contaminant_taxids(&samples, 1.0, 0.6, &[]), vec!["2".to_string()]);
        assert_eq!(
            sorted(prevalence_contaminant_taxids(&samples, 1.0, 0.25, &[])),
            vec!["1".to_string(), "2".to_string(), "9".to_string()]
        );
    }

    #[test]
    fn rpm_floor_excludes_low_abundance_presence() {
        // "1" is present in both samples but below the floor in one; with floor 5.0 it is
        // present in only 1 of 2 -> 0.5, excluded at threshold 0.6.
        let samples = vec![vec![taxon("1", 10.0)], vec![taxon("1", 1.0)]];
        assert!(prevalence_contaminant_taxids(&samples, 5.0, 0.6, &[]).is_empty());
        // With no floor it is present in both -> contaminant.
        assert_eq!(prevalence_contaminant_taxids(&samples, 0.0, 0.6, &[]), vec!["1".to_string()]);
    }

    #[test]
    fn subset_restricts_which_taxa_are_considered() {
        let samples = vec![vec![taxon("1", 10.0), taxon("2", 10.0)], vec![taxon("1", 10.0), taxon("2", 10.0)]];
        // Restrict to "1" only: "2" is never considered even though it is also prevalent.
        assert_eq!(prevalence_contaminant_taxids(&samples, 1.0, 1.0, &["1".to_string()]), vec!["1".to_string()]);
    }

    #[test]
    fn taxon_counted_once_per_sample_even_with_multiple_records() {
        // Two profile records in one sample must not double-count prevalence.
        let mut t = taxon("1", 10.0);
        t.evidence.profile.push(ProfileRecord {
            id: "lib2".to_string(),
            tool: ProfileTool::Sylph,
            mode: AbundanceMode::Profile,
            reads: 1,
            rpm: 10.0,
            contigs: 0,
            bases: 0,
            bpm: 0.0,
            abundance: 0.0,
        });
        let samples = vec![vec![t], vec![taxon("9", 10.0)]];
        // "1" present in 1 of 2 -> 0.5, excluded at 0.6 (not 1.0 from double counting).
        assert!(prevalence_contaminant_taxids(&samples, 1.0, 0.6, &[]).is_empty());
    }
}
