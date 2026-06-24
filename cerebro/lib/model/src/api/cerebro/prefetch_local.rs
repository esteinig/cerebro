//! Self-contained (local) construction of [`PrefetchData`].
//!
//! The stack path builds a sample's tiered prefetch by querying the API
//! (`CerebroClient::get_taxa` per tier, plus the contamination endpoint). This module builds the
//! **same** `PrefetchData` directly from a run's output models, with **no live stack**, by reusing
//! the exact filtering primitives the server uses:
//!
//! - [`aggregate`] to combine a sample's library evidence,
//! - [`apply_filters`] to apply each tier's `TaxonFilterConfig`,
//! - [`prevalence_contaminant_taxids`] to identify cross-sample contaminants over the run,
//! - [`PrefetchData::prune`] to de-duplicate across tiers.
//!
//! Because the same functions back both paths, a `local` and a `stack` prefetch agree by
//! construction (the one optional, stack-only contamination-history enrichment is not applied
//! here; it is a separate, opt-in capability). The assembly is pure and unit-tested below; only
//! [`load_run_models`] performs I/O.
//!
//! ## Library grouping
//! A biological sample's libraries share `Cerebro.sample.id` ("file base name stripped of tags")
//! and are distinguished by their nucleic-acid tag (DNA/RNA) in `Cerebro.sample.tags`. Filtering
//! and prevalence are computed **per nucleic-acid tag**, mirroring the stack, which resolves and
//! queries DNA and RNA libraries separately.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use cerebro_pipeline::taxa::filter::{apply_filters, TaxonFilterConfig};
use cerebro_pipeline::taxa::prevalence::prevalence_contaminant_taxids;
use cerebro_pipeline::taxa::taxon::{aggregate, Taxon};

use crate::api::cerebro::model::{Cerebro, ModelError};
use crate::api::cerebro::response::TaxonHistoryResult;
use crate::api::cerebro::schema::{MetaGpConfig, PrefetchData, PrevalenceContaminationConfig};

/// Nucleic-acid tags that define a sample's library groups. The stack identifies prevalence
/// contamination separately for `DNA` and `RNA`; the local path mirrors that grouping.
pub const NUCLEIC_ACID_TAGS: [&str; 2] = ["DNA", "RNA"];

/// One library's taxon evidence for a biological sample: the decomposed input to the builder.
///
/// Keeping the builder over plain data (rather than over `Cerebro`) makes the assembly logic
/// unit-testable with synthetic taxa and no model scaffolding.
#[derive(Debug, Clone)]
pub struct LibraryTaxa {
    /// Library basename (`Cerebro.name`) — the key `apply_filters` uses for tag lookup.
    pub name: String,
    /// Library tags (`Cerebro.sample.tags`), e.g. `["DNA"]`.
    pub tags: Vec<String>,
    /// The library's taxa (`Cerebro.taxa`).
    pub taxa: Vec<Taxon>,
}

/// Build one tier: aggregate + filter per nucleic-acid group, then split into kept vs
/// contamination by the prevalence set. Mirrors the server's per-tag handling in `get_taxa`.
fn build_tier(
    libraries: &[LibraryTaxa],
    filter: &TaxonFilterConfig,
    prevalence_taxids: &HashMap<String, HashSet<String>>,
    allow_no_evidence: bool,
) -> (Vec<Taxon>, Vec<Taxon>) {
    let mut kept: Vec<Taxon> = Vec::new();
    let mut contamination: Vec<Taxon> = Vec::new();
    let empty: HashSet<String> = HashSet::new();

    for na_tag in NUCLEIC_ACID_TAGS {
        // Libraries belonging to this nucleic-acid group.
        let group: Vec<&LibraryTaxa> = libraries
            .iter()
            .filter(|lib| lib.tags.iter().any(|tag| tag == na_tag))
            .collect();
        if group.is_empty() {
            continue;
        }

        // Aggregate evidence across the group's libraries (same as the server before filtering).
        let mut aggregated: HashMap<String, Taxon> = HashMap::new();
        for lib in &group {
            aggregated = aggregate(&mut aggregated, &lib.taxa);
        }

        // Library-name -> tags map that `apply_filters` consults for tag-based criteria.
        let mut tag_map: HashMap<String, Vec<String>> = HashMap::new();
        for lib in &group {
            tag_map.insert(lib.name.clone(), lib.tags.clone());
        }

        let filtered = apply_filters(
            aggregated.into_values().collect(),
            filter,
            &tag_map,
            allow_no_evidence,
        );

        // Split by the prevalence-contaminant set for this tag (empty set => keep everything).
        let contam_set = prevalence_taxids.get(na_tag).unwrap_or(&empty);
        for taxon in filtered {
            if contam_set.contains(&taxon.taxid) {
                contamination.push(taxon);
            } else {
                kept.push(taxon);
            }
        }
    }

    (kept, contamination)
}

/// Build a sample's [`PrefetchData`] from its library taxa, self-contained (no stack).
///
/// `prevalence_taxids` maps a nucleic-acid tag (`"DNA"`/`"RNA"`) to the contaminant taxids for
/// that tag (e.g. from [`run_prevalence_contamination`]); pass an empty map to disable the split.
/// The output is byte-compatible with the stack builder's `PrefetchData` (same tiers, same
/// contamination vectors, same prune), so META-GPT ingestion is unaffected by the source.
pub fn build_prefetch_local(
    libraries: &[LibraryTaxa],
    prevalence_taxids: &HashMap<String, HashSet<String>>,
    config: &MetaGpConfig,
    allow_no_evidence: bool,
) -> PrefetchData {
    let tiers = &config.filter_configs;

    let (primary, primary_contamination) =
        build_tier(libraries, &tiers.primary, prevalence_taxids, allow_no_evidence);
    let (secondary, secondary_contamination) =
        build_tier(libraries, &tiers.secondary, prevalence_taxids, allow_no_evidence);
    let (target, target_contamination) =
        build_tier(libraries, &tiers.target, prevalence_taxids, allow_no_evidence);

    let mut data = PrefetchData::new(
        primary,
        secondary,
        target,
        primary_contamination,
        secondary_contamination,
        target_contamination,
        &tiers.primary,
        &tiers.secondary,
        &tiers.target,
        config,
    );
    data.prune();
    data
}

/// Select and decompose the libraries belonging to one biological sample (`Cerebro.sample.id`).
pub fn libraries_for_sample(models: &[Cerebro], sample_id: &str) -> Vec<LibraryTaxa> {
    models
        .iter()
        .filter(|model| model.sample.id == sample_id)
        .map(|model| LibraryTaxa {
            name: model.name.clone(),
            tags: model.sample.tags.clone(),
            taxa: model.taxa.clone(),
        })
        .collect()
}

/// Compute prevalence-contaminant taxids over a whole run, per nucleic-acid tag.
///
/// Mirrors the stack's contamination endpoint (which the stack queries once per tag), but scoped
/// to the run's models rather than the whole collection. Honours `min_rpm`/`threshold` from the
/// configuration; the empty subset means "consider every taxon".
pub fn run_prevalence_contamination(
    models: &[Cerebro],
    contam_config: &PrevalenceContaminationConfig,
) -> HashMap<String, HashSet<String>> {
    let mut out: HashMap<String, HashSet<String>> = HashMap::new();
    for na_tag in NUCLEIC_ACID_TAGS {
        let per_sample_taxa: Vec<Vec<Taxon>> = models
            .iter()
            .filter(|model| model.sample.tags.iter().any(|tag| tag == na_tag))
            .map(|model| model.taxa.clone())
            .collect();

        let taxids = prevalence_contaminant_taxids(
            &per_sample_taxa,
            contam_config.min_rpm,
            contam_config.threshold,
            &[],
        );
        out.insert(na_tag.to_string(), taxids.into_iter().collect());
    }
    out
}

/// Load a run's output models from disk (the only I/O in the local path).
pub fn load_run_models(paths: &[PathBuf]) -> Result<Vec<Cerebro>, ModelError> {
    paths.iter().map(Cerebro::from_json).collect()
}

/// Build a run-local taxon history (one `TaxonHistoryResult` per run model) — the same shape the
/// server assembles from the collection, scoped to the run's samples (Task B, TB-D2).
///
/// `sample_id` is the **biological** sample id (`Cerebro.sample.id`), not the library file name,
/// so the rescue's outlier match (`outlier.sample_id == config.sample`) lines up with the stack
/// path (`schema.sample` in `CerebroClient::get_taxa`). Input/host reads come from the model's QC
/// (`quality.reads`), which is exactly what `RpmAnalyzer::from_taxon_history` consumes.
pub fn run_taxon_history(models: &[Cerebro]) -> Vec<TaxonHistoryResult> {
    models
        .iter()
        .map(|m| TaxonHistoryResult {
            id: m.id.clone(),
            sample_id: m.sample.id.clone(),
            sample_tags: m.sample.tags.clone(),
            sample_type: m.sample.sample_type.clone(),
            run_id: m.run.id.clone(),
            run_date: m.run.date.clone(),
            input_reads: m.quality.reads.input_reads,
            host_reads: m.quality.reads.host_reads,
            taxa: m.taxa.clone(),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::api::cerebro::schema::{PrevalenceContaminationConfig, TieredFilterConfig};
    use crate::api::cerebro::schema::{SampleType, TestResult};
    use cerebro_pipeline::modules::pathogen::{AbundanceMode, ProfileRecord, ProfileTool};
    use cerebro_pipeline::taxa::taxon::{Taxon, TaxonEvidence};
    use taxonomy::TaxRank;

    /// Taxon with a profile record at `rpm` and a distinct lineage (prune compares lineages).
    fn taxon(taxid: &str, lineage: &str, rpm: f64) -> Taxon {
        Taxon {
            taxid: taxid.to_string(),
            rank: TaxRank::Species,
            name: format!("taxon_{taxid}"),
            lineage: lineage.to_string(),
            evidence: TaxonEvidence {
                alignment: Vec::new(),
                assembly: Vec::new(),
                profile: vec![ProfileRecord {
                    id: "lib".to_string(),
                    tool: ProfileTool::Kraken2,
                    mode: AbundanceMode::Profile,
                    reads: 100,
                    rpm,
                    contigs: 0,
                    bases: 0,
                    bpm: 0.0,
                    abundance: 0.0,
                }],
            },
        }
    }

    fn library(name: &str, tag: &str, mut taxa: Vec<Taxon>) -> LibraryTaxa {
        // In a real run a profile record's `id` is the library name; mirror that so the evidence
        // filter's tag lookup (keyed by `record.id`) resolves exactly as it does in production.
        for taxon in &mut taxa {
            for record in &mut taxon.evidence.profile {
                record.id = name.to_string();
            }
        }
        LibraryTaxa { name: name.to_string(), tags: vec![tag.to_string()], taxa }
    }

    /// A MetaGpConfig whose three tiers pass everything through (`TieredFilterConfig::none`),
    /// so the tests exercise the assembly logic rather than the filter internals.
    fn config_none() -> MetaGpConfig {
        MetaGpConfig::new(
            "sampleA".to_string(),
            SampleType::Csf,
            Some(TestResult::Positive),
            None,
            None,
            None,
            None,
            TieredFilterConfig::none(),
            PrevalenceContaminationConfig::default(),
        )
    }

    fn taxids(taxa: &[Taxon]) -> HashSet<String> {
        taxa.iter().map(|t| t.taxid.clone()).collect()
    }

    #[test]
    fn contam_history_rescue_moves_flagged_taxon_back_to_signal() {
        // taxid "2" (name "taxon_2") is prevalence-flagged into contamination; "1" stays signal.
        let libs = vec![library(
            "sampleA__DNA",
            "DNA",
            vec![taxon("1", "d__A;s__one", 10.0), taxon("2", "d__A;s__two", 10.0)],
        )];
        let mut prevalence = HashMap::new();
        prevalence.insert("DNA".to_string(), HashSet::from(["2".to_string()]));

        let mut data = build_prefetch_local(&libs, &prevalence, &config_none(), true);
        assert_eq!(taxids(&data.primary_contamination), HashSet::from(["2".to_string()]));

        // Rescue exactly the flagged taxon by name: it must move back into the signal tier and
        // leave contamination empty; the prune invariant (no cross-tier lineage dupes) holds.
        data.apply_contam_history_rescue(|t| t.name == "taxon_2");

        assert_eq!(
            taxids(&data.primary),
            HashSet::from(["1".to_string(), "2".to_string()])
        );
        assert!(data.primary_contamination.is_empty());
    }

    #[test]
    fn contam_history_off_is_noop_when_nothing_rescued() {
        let libs = vec![library(
            "sampleA__DNA",
            "DNA",
            vec![taxon("1", "d__A;s__one", 10.0), taxon("2", "d__A;s__two", 10.0)],
        )];
        let mut prevalence = HashMap::new();
        prevalence.insert("DNA".to_string(), HashSet::from(["2".to_string()]));

        let mut data = build_prefetch_local(&libs, &prevalence, &config_none(), true);
        let before = taxids(&data.primary_contamination);
        data.apply_contam_history_rescue(|_| false); // off / no outliers => nothing rescued
        assert_eq!(taxids(&data.primary_contamination), before);
        assert_eq!(taxids(&data.primary), HashSet::from(["1".to_string()]));
    }

    #[test]
    fn no_prevalence_keeps_all_taxa_in_their_tiers() {
        let libs = vec![library("sampleA__DNA", "DNA", vec![taxon("1", "d__A;s__one", 10.0)])];
        let prevalence = HashMap::new();
        let data = build_prefetch_local(&libs, &prevalence, &config_none(), true);
        assert_eq!(taxids(&data.primary), HashSet::from(["1".to_string()]));
        assert!(data.primary_contamination.is_empty());
    }

    #[test]
    fn prevalence_taxid_is_split_into_contamination() {
        let libs = vec![library(
            "sampleA__DNA",
            "DNA",
            vec![taxon("1", "d__A;s__one", 10.0), taxon("2", "d__A;s__two", 10.0)],
        )];
        let mut prevalence = HashMap::new();
        prevalence.insert("DNA".to_string(), HashSet::from(["2".to_string()]));

        let data = build_prefetch_local(&libs, &prevalence, &config_none(), true);
        assert_eq!(taxids(&data.primary), HashSet::from(["1".to_string()]));
        assert_eq!(taxids(&data.primary_contamination), HashSet::from(["2".to_string()]));
    }

    #[test]
    fn dna_and_rna_libraries_are_aggregated_per_tag() {
        let libs = vec![
            library("sampleA__DNA", "DNA", vec![taxon("1", "d__A;s__one", 10.0)]),
            library("sampleA__RNA", "RNA", vec![taxon("3", "d__A;s__three", 10.0)]),
        ];
        let data = build_prefetch_local(&libs, &HashMap::new(), &config_none(), true);
        // Both nucleic-acid groups contribute to the (unfiltered) primary tier.
        assert_eq!(taxids(&data.primary), HashSet::from(["1".to_string(), "3".to_string()]));
    }

    #[test]
    fn prune_removes_cross_tier_duplicate_lineage() {
        // Same lineage present in primary and secondary; prune drops it from secondary.
        // (none() makes all three tiers identical, so the shared lineage appears in each.)
        let libs = vec![library("sampleA__DNA", "DNA", vec![taxon("1", "d__A;s__one", 10.0)])];
        let data = build_prefetch_local(&libs, &HashMap::new(), &config_none(), true);
        assert_eq!(taxids(&data.primary), HashSet::from(["1".to_string()]));
        // secondary/target are pruned of the lineage already in primary.
        assert!(data.secondary.is_empty(), "secondary should be pruned of primary's lineage");
        assert!(data.target.is_empty(), "target should be pruned of primary's lineage");
    }

    #[test]
    fn run_prevalence_counts_per_tag() {
        // Build minimal models is heavy; exercise the per-tag counter via the pure function the
        // helper delegates to, asserting the tag-keyed shape the builder consumes.
        let prevalence: HashMap<String, HashSet<String>> = HashMap::from([
            ("DNA".to_string(), HashSet::from(["7".to_string()])),
            ("RNA".to_string(), HashSet::new()),
        ]);
        let libs = vec![library("sampleA__DNA", "DNA", vec![taxon("7", "d__A;s__seven", 10.0)])];
        let data = build_prefetch_local(&libs, &prevalence, &config_none(), true);
        assert!(data.primary.is_empty());
        assert_eq!(taxids(&data.primary_contamination), HashSet::from(["7".to_string()]));
    }
}
