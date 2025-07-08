use std::collections::HashMap;
use std::collections::HashSet;
use std::path::Path;

use serde::{Deserialize, Serialize};
use taxonomy::TaxRank;

use crate::error::WorkflowError;
use crate::modules::pathogen::PathogenDetectionRank;
use crate::modules::pathogen::ProfileTool;
use crate::modules::pathogen::AbundanceMode;
use crate::utils::read_tsv;
use crate::utils::write_tsv;

use super::taxon::collapse_taxa;
use super::taxon::LineageOperations;
use super::taxon::Taxon;


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TaxonFilterConfig {
    pub rank: Option<PathogenDetectionRank>,        // Filter by specific taxonomic rank
    pub domains: Vec<String>,                       // Filter by domain names
    pub tools: Vec<ProfileTool>,                    // Filter by specific detection tools
    pub modes: Vec<AbundanceMode>,                  // Filter by detection modes (Sequence/Profile)
    pub min_bases: u64,
    pub max_bases: Option<u64>,
    pub min_bpm: f64,                             
    pub min_reads: u64,                             // Minimum read count for inclusion
    pub min_rpm: f64,                               // Minimum RPM (Reads per million) for inclusion
    pub max_rpm: Option<f64>,
    pub min_abundance: f64,                         // Minimum abundance for inclusion
    pub ntc_ratio: Option<f64>,                     // NTC ratio threshold
    pub lineage: Option<Vec<LineageFilterConfig>>,  // Lineage filter configuration if specified
    pub targets: Option<Vec<String>>,               // Subset taxa to these lineage components
    pub collapse_variants: bool,                    // Collapse species variants by name (GTDB, e.g. Haemophilus influenzae_A or Haemophilus influenzae_H) - sums evidence and adjusts taxon name
    pub ignore_taxstr: Option<Vec<String>>,         // Remove any of these matches 
}

impl TaxonFilterConfig {
    pub fn target_set(&self) -> Option<HashSet<&str>> {
        // Build a HashSet of target names for efficient lookup.
        if let Some(targets) = &self.targets { Some(HashSet::from_iter(targets.into_iter().map(|t| t.as_str()))) } else { None }
    }
    /// Checks if the given taxon passes *all* defined lineage filter configurations.
    pub fn passes_filters(&self, taxon: &Taxon) -> bool {
        if let Some(lineage_filters) = &self.lineage {
            for filter in lineage_filters {
                if !filter.passes_alignment_filters(taxon) {
                    return false;
                }
                // (Additional criteria can be added here.)
            }
        }
        true
    }
}

impl Default for TaxonFilterConfig {
    fn default() -> Self {
        Self {
            rank: Some(PathogenDetectionRank::Species),
            domains: Vec::new(),
            tools: Vec::new(),
            modes: Vec::new(),
            min_bases: 0,
            max_bases: None,
            min_bpm: 0.0,
            min_reads: 0,
            min_rpm: 0.0,
            max_rpm: None,
            min_abundance: 0.0,
            ntc_ratio: None,
            lineage: None,
            targets: None,
            collapse_variants: false,
            ignore_taxstr: None
        }
    }
}

impl TaxonFilterConfig {
    pub fn validation() -> Self {
        Self {
            rank: Some(PathogenDetectionRank::Species),
            domains: Vec::new(),
            tools: vec![ProfileTool::Bracken, ProfileTool::Metabuli, ProfileTool::Ganon2, ProfileTool::Blast, ProfileTool::Vircov],
            modes: vec![AbundanceMode::Mixed],
            min_bases: 200,
            max_bases: None,
            min_bpm: 0.0,
            min_reads: 0,
            min_rpm: 5.0,
            max_rpm: None,
            min_abundance: 0.0,
            ntc_ratio: Some(3.0),
            lineage: Some(vec![
                LineageFilterConfig::validation_bacteria(),
                LineageFilterConfig::validation_viruses(),
                LineageFilterConfig::validation_eukaryota()
            ]),
            targets: None,
            collapse_variants: false,
            ignore_taxstr: None
        }
    }
    pub fn gp_above_threshold(taxstr: Option<Vec<String>>) -> Self {
        Self {
            rank: Some(PathogenDetectionRank::Species),
            domains: Vec::new(),
            tools: vec![ProfileTool::Bracken, ProfileTool::Metabuli, ProfileTool::Ganon2, ProfileTool::Blast, ProfileTool::Vircov],
            modes: vec![AbundanceMode::Mixed],
            min_bases: 0,
            max_bases: None,
            min_bpm: 0.0,
            min_reads: 0,
            min_rpm: 0.0,
            max_rpm: None,
            min_abundance: 0.0,
            ntc_ratio: Some(1.0),
            lineage: Some(vec![
                LineageFilterConfig::gp_bacteria_above_threshold(),
                LineageFilterConfig::gp_viruses_above_threshold(),
                LineageFilterConfig::gp_eukaryota_above_threshold()
            ]),
            targets: None,
            collapse_variants: false,
            ignore_taxstr: taxstr
        }
    }
    pub fn gp_below_threshold(taxstr: Option<Vec<String>>) -> Self {
        Self {
            rank: Some(PathogenDetectionRank::Species),
            domains: Vec::new(),
            tools: vec![ProfileTool::Bracken, ProfileTool::Metabuli, ProfileTool::Ganon2, ProfileTool::Blast, ProfileTool::Vircov],
            modes: vec![AbundanceMode::Mixed],
            min_bases: 0,
            max_bases: None,
            min_bpm: 0.0,
            min_reads: 0,
            min_rpm: 0.0,
            max_rpm: None,
            min_abundance: 0.0,
            ntc_ratio: Some(1.0),
            lineage: Some(vec![
                LineageFilterConfig::gp_bacteria_below_threshold(),
                LineageFilterConfig::gp_viruses_below_threshold(),
                LineageFilterConfig::gp_eukaryota_below_threshold(),
            ]),
            targets: None,
            collapse_variants: false,
            ignore_taxstr: taxstr
        }
    }
    pub fn gp_target_threshold(taxstr: Option<Vec<String>>) -> Self {
        Self {
            rank: Some(PathogenDetectionRank::Species),
            domains: vec![],
            tools: vec![ProfileTool::Bracken, ProfileTool::Metabuli, ProfileTool::Ganon2, ProfileTool::Blast, ProfileTool::Vircov],
            modes: vec![AbundanceMode::Mixed],
            min_bases: 0,
            max_bases: None,
            min_bpm: 0.0,
            min_reads: 0,
            min_rpm: 0.0,
            max_rpm: None,
            min_abundance: 0.0,
            ntc_ratio: None,
            lineage: Some(vec![
                LineageFilterConfig::gp_viruses_target_threshold(),
            ]),
            targets: Some(TargetList::gp_vertebrate_viruses().to_vec()),
            collapse_variants: false,
            ignore_taxstr: taxstr
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetListRecord {
    pub target_name: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetList {
    pub targets: Vec<TargetListRecord>,
}

impl TargetList {
    /// Combines multiple target lists into one.
    /// Duplicates are removed based on the target_name field.
    pub fn combine(lists: &[TargetList]) -> Self {
        let mut combined: Vec<TargetListRecord> = Vec::new();
        for list in lists {
            combined.extend(list.targets.clone());
        }
        let mut seen = HashSet::new();
        let deduped: Vec<TargetListRecord> = combined
            .into_iter()
            .filter(|record| seen.insert(record.target_name.clone()))
            .collect();
        Self { targets: deduped }
    }

    /// Loads a target list from a TSV file at the given path.
    pub fn from_tsv(path: &Path) -> Result<Self, WorkflowError> {
        Ok(Self {
            targets: read_tsv(path, false, true)?,
        })
    }

    /// Writes the target list to a TSV file at the given path.
    pub fn to_tsv(&self, path: &Path) -> Result<(), WorkflowError> {
        write_tsv(&self.targets, path, true)
    }

    /// Converts the target list to a vector of target names.
    pub fn to_vec(&self) -> Vec<String> {
        self.targets
            .clone()
            .into_iter()
            .map(|r| r.target_name)
            .collect()
    }

    /// Parses TSV data from a string.
    /// Assumes the first line is a header and each subsequent line is a target name.
    pub fn from_tsv_str(tsv_data: &str) -> Result<Self, WorkflowError> {
        let mut targets = Vec::new();
        // Split lines and skip the header line.
        for (i, line) in tsv_data.lines().enumerate() {
            if i == 0 {
                continue; // Skip header
            }
            let trimmed = line.trim();
            if !trimmed.is_empty() {
                targets.push(TargetListRecord {
                    target_name: trimmed.to_string(),
                });
            }
        } 
        // Deduplicate the target records based on target_name.
        let mut seen = HashSet::new();
        targets.retain(|record| seen.insert(record.target_name.clone()));
        Ok(Self { targets })
    }

    /// Loads the target list from an embedded TSV file.
    /// The file "templates/ictv_vertebrate_targets.tsv" is embedded at compile time.
    pub fn gp_vertebrate_viruses() -> Self {
        let tsv_data = include_str!("../../templates/ictv_vertebrate_targets.tsv");
        Self::from_tsv_str(tsv_data).unwrap_or_else(|e| {
            log::error!("Error parsing TSV data: {:?}", e);
            Self { targets: Vec::new() }
        })
    }
    /// Loads the target list from an embedded TSV file.
    /// The file "templates/ictv_prokaryote_targets.tsv" is embedded at compile time.
    pub fn gp_prokaryote_viruses() -> Self {
        let tsv_data = include_str!("../../templates/ictv_prokaryote_targets.tsv");
        Self::from_tsv_str(tsv_data).unwrap_or_else(|e| {
            log::error!("Error parsing TSV data: {:?}", e);
            Self { targets: Vec::new() }
        })
    }
    /// Loads the target list from an embedded TSV file.
    /// The file "templates/cns_prokaryote_targets.tsv" is embedded at compile time.
    pub fn gp_cns_bacteria() -> Self {
        let tsv_data = include_str!("../../templates/cns_prokaryote_targets.tsv");
        Self::from_tsv_str(tsv_data).unwrap_or_else(|e| {
            log::error!("Error parsing TSV data: {:?}", e);
            Self { targets: Vec::new() }
        })
    }
}


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct LineageFilterConfig {
    pub lineages: Vec<String>,
    pub tags: Vec<String>,
    pub min_alignment_tools: Option<usize>,
    pub min_alignment_rpm: Option<f64>,
    pub min_alignment_regions: Option<u64>,
    pub min_alignment_regions_coverage: Option<f64>,
    pub min_kmer_tools: Option<usize>,
    pub min_kmer_rpm: Option<f64>,
    pub min_assembly_tools: Option<usize>,
}
impl LineageFilterConfig {
     /// Checks whether the given taxon meets the alignment evidence criteria.
     pub fn passes_alignment_filters(&self, taxon: &Taxon) -> bool {
        // (1) If a minimum number of alignment tools is specified and is greater than 0,
        // then there must be at least one alignment record.
        if let Some(min_tools) = self.min_alignment_tools {
            if min_tools > 0 && taxon.evidence.alignment.is_empty() {
                return false;
            }
        }

        // (2) If a requirement is defined for the minimum number of alignment regions,
        // check that at least one alignment record in the taxon meets both the regions
        // and (optionally) the coverage threshold.
        if let Some(min_regions) = self.min_alignment_regions {
            // Use the specified coverage, or 0.0 if not provided.
            let min_cov = self.min_alignment_regions_coverage.unwrap_or(0.0);
            let valid_alignment = taxon.evidence.alignment.iter().any(|record| {
                record.scan_regions >= min_regions && record.scan_coverage >= min_cov
            });
            if !valid_alignment {
                return false;
            }
        }
        true
    }
    pub fn validation_viruses() -> Self {
        Self {
            lineages: vec!["d__Viruses".to_string()],
            tags: vec!["DNA".to_string(), "RNA".to_string()],
            min_alignment_tools: Some(1),
            min_alignment_rpm: Some(10.0),
            min_alignment_regions: None,
            min_alignment_regions_coverage: None,
            min_kmer_tools: Some(1),
            min_kmer_rpm: Some(10.0),
            min_assembly_tools: None
        }
    }
    pub fn validation_bacteria() -> Self {
        Self {
            lineages: vec!["d__Bacteria".to_string(), "d__Archaea".to_string()],
            tags: vec!["DNA".to_string()],
            min_alignment_tools: None,
            min_alignment_rpm: None,
            min_alignment_regions: None,
            min_alignment_regions_coverage: None,
            min_kmer_tools: Some(3),
            min_kmer_rpm: Some(10.0) ,
            min_assembly_tools: None
        }
    }
    pub fn validation_eukaryota() -> Self {
        Self {
            lineages: vec!["d__Eukaryota".to_string()],
            tags: vec!["DNA".to_string()],
            min_alignment_tools: None,
            min_alignment_rpm: None,
            min_alignment_regions: None,
            min_alignment_regions_coverage: None,
            min_kmer_tools: Some(3),
            min_kmer_rpm: Some(10.0) ,
            min_assembly_tools: Some(1)
        }
    }

    pub fn gp_viruses_above_threshold() -> Self {
        Self {
            lineages: vec!["d__Viruses".to_string()],
            tags: vec!["DNA".to_string(), "RNA".to_string()],
            min_alignment_tools: Some(1),
            min_alignment_rpm: Some(10.0),
            min_alignment_regions: None,
            min_alignment_regions_coverage: None,
            min_kmer_tools: Some(2),
            min_kmer_rpm: Some(10.0),
            min_assembly_tools: None
        }
    }
    pub fn gp_bacteria_above_threshold() -> Self {
        Self {
            lineages: vec!["d__Bacteria".to_string(), "d__Archaea".to_string()],
            tags: vec!["DNA".to_string()],
            min_alignment_tools: None,
            min_alignment_rpm: None,
            min_alignment_regions: None,
            min_alignment_regions_coverage: None,
            min_kmer_tools: Some(3),
            min_kmer_rpm: Some(10.0) ,
            min_assembly_tools: None
        }
    }
    pub fn gp_eukaryota_above_threshold() -> Self {
        Self {
            lineages: vec!["d__Eukaryota".to_string()],
            tags: vec!["DNA".to_string()],
            min_alignment_tools: None,
            min_alignment_rpm: None,
            min_alignment_regions: None,
            min_alignment_regions_coverage: None,
            min_kmer_tools: Some(3),
            min_kmer_rpm: Some(10.0) ,
            min_assembly_tools: Some(1)
        }
    }

    pub fn gp_viruses_below_threshold() -> Self {
        Self {
            lineages: vec!["d__Viruses".to_string()],
            tags: vec!["DNA".to_string(), "RNA".to_string()],
            min_alignment_tools: None,
            min_alignment_rpm: None,
            min_alignment_regions: None,
            min_alignment_regions_coverage: None,
            min_kmer_tools: Some(2),
            min_kmer_rpm: Some(3.0),
            min_assembly_tools: None
        }
    }
    pub fn gp_bacteria_below_threshold() -> Self {
        Self {
            lineages: vec!["d__Bacteria".to_string(), "d__Archaea".to_string()],
            tags: vec!["DNA".to_string()],
            min_alignment_tools: None,
            min_alignment_rpm: None,
            min_alignment_regions: None,
            min_alignment_regions_coverage: None,
            min_kmer_tools: Some(3),
            min_kmer_rpm: Some(3.0),
            min_assembly_tools: None
        }
    }
    pub fn gp_eukaryota_below_threshold() -> Self {
        Self {
            lineages: vec!["d__Eukaryota".to_string()],
            tags: vec!["DNA".to_string()],
            min_alignment_tools: None,
            min_alignment_rpm: None,
            min_alignment_regions: None,
            min_alignment_regions_coverage: None,
            min_kmer_tools: Some(3),
            min_kmer_rpm: Some(5.0) ,
            min_assembly_tools: Some(1)
        }
    }

    pub fn gp_viruses_target_threshold() -> Self {
        Self {
            lineages: vec!["d__Viruses".to_string()],
            tags: vec!["DNA".to_string(), "RNA".to_string()],
            min_alignment_tools: Some(1),
            min_alignment_rpm: Some(0.0),
            min_alignment_regions: Some(2),
            min_alignment_regions_coverage: Some(20.0),
            min_kmer_tools: Some(1),
            min_kmer_rpm: Some(0.0),
            min_assembly_tools: None
        }
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrevalenceContaminationConfig {
    pub threshold: f64,
    pub min_rpm: f64,
    pub sample_type: Option<String>,
    pub collapse_variants: bool
}
impl PrevalenceContaminationConfig {
    pub fn validation() -> Self {
        Self {
            threshold: 0.60,
            min_rpm: 0.0,
            sample_type: None,
            collapse_variants: false
        }
    }
    pub fn gp_default() -> Self {
        Self {
            threshold: 0.60,
            min_rpm: 0.0,
            sample_type: None,
            collapse_variants: false
        }
    }
}

type SampleId = String;
type Tag = String;


pub fn apply_filters(mut taxa: Vec<Taxon>, filter_config: &TaxonFilterConfig, sample_tags: &HashMap<SampleId, Vec<Tag>>, allow_no_evidence: bool) -> Vec<Taxon> {

    // Filter by taxonomic rank
    if let Some(rank) = &filter_config.rank {
        let tax_rank: TaxRank = (*rank).clone().into(); // Convert PathogenDetectionRank to TaxRank --> important because the taxon abstraction uses TaxRank from taxonomy crate (might need to change this)
        taxa = taxa
            .into_iter()
            .filter(|taxon| taxon.rank == tax_rank)
            .collect();
    }



    // Filter by domain
    if !filter_config.domains.is_empty() {
        taxa = taxa
            .into_iter()
            .filter(|taxon| {
                match taxon.lineage.get_domain() {
                    Some(taxon_domain) => filter_config.domains.contains(&taxon_domain),
                    None => false
                }
            })
            .collect();
    }

    // Collapse taxon variant names
    if filter_config.collapse_variants {
        taxa = collapse_taxa(taxa).expect("FAILED TO COLLAPSE TAXA"); // TODO
    }

    // Apply target filter if specified - this is really slow at the moment because of the String checks?
    if let Some(target_set) = &filter_config.target_set() {
        taxa = taxa.into_iter()
            .filter(|taxon| target_set.contains(taxon.name.as_str()))
            .collect();
    }

    // Apply taxstr filter to remove matches in name
    if let Some(taxstr) = &filter_config.ignore_taxstr {
        taxa = taxa.into_iter()
            .filter(|taxon| {
                for s in taxstr {
                    if taxon.name.contains(s) {
                        return false
                    }
                }
                true
            })
            .collect();
    }

    // Filter by detection tools, modes, and thresholds
    taxa = apply_evidence_filters(taxa, filter_config, sample_tags, allow_no_evidence);

    // If lineage filters are specified, apply the filters across retained taxa:
    if let Some(lineage_filters) = &filter_config.lineage {
        taxa = apply_lineage_filters(taxa, lineage_filters, sample_tags);
    }


    taxa
}


pub fn apply_evidence_filters(
    taxa: Vec<Taxon>,
    filter_config: &TaxonFilterConfig,
    sample_tags: &HashMap<SampleId, Vec<Tag>>,
    allow_no_evidence: bool
) -> Vec<Taxon> {

    // Are the evidence-aggregated Taxa only from NTC or ENV controls (no sample library present)
    let only_ntc = sample_tags.clone().into_values().all(|tags| tags.contains(&"NTC".to_string()) || tags.contains(&"ENV".to_string()));

    taxa
    .into_iter()
    .map(|mut taxon| {

        // Step 1: Collect NTC RPMs for matching tools, modes, and tags summing NTC RPM if multiple NTC samples are provided
        let mut ntc_rpms: HashMap<(ProfileTool, AbundanceMode, String), f64> = HashMap::new();

        taxon
            .evidence
            .profile
            .iter()
            .filter(|record| {
                if let Some(tags) = sample_tags.get(&record.id) {
                    tags.contains(&"NTC".to_string()) || tags.contains(&"ENV".to_string()) // Find records with "NTC" tag
                } else {
                    false
                }
            })
            .for_each(|record| {
                if let Some(tags) = sample_tags.get(&record.id) {
                    // Extract nucleic acid tag (DNA/RNA) for NTC record if present
                    if let Some(tag) = tags.iter().find(|tag| tag == &&"DNA".to_string() || tag == &&"RNA".to_string()) {
                        let key = (record.tool.clone(), record.mode.clone(), tag.clone());
                        // Sum RPM values for matching keys
                        ntc_rpms.entry(key)
                            .and_modify(|e| *e += record.rpm) // If key exists, add to existing value
                            .or_insert(record.rpm); // Otherwise, insert the current RPM value
                    }
                }
            });
        
        // Step 2: Filter profile records based on ntc_ratio
        taxon.evidence.profile = taxon
            .evidence
            .profile
            .into_iter()
            .filter(|record| {
                if let Some(tags) = sample_tags.get(&record.id) {
                    if tags.contains(&"NTC".to_string()) || tags.contains(&"ENV".to_string()) {

                        if only_ntc {
                            return true // If we have aggregated NTC evidence only, return the evidence
                        } else {
                            return false // If we have aggregated sample evidence, exclude the NTC evidence from returned taxa
                        }
                    }
                    let nucleic_acid_tag = tags
                        .iter()
                        .find(|tag| tag == &&"DNA".to_string() || tag == &&"RNA".to_string());

                    if let (Some(ratio_threshold), Some(tag)) = (filter_config.ntc_ratio, nucleic_acid_tag) {
                        let key = (record.tool.clone(), record.mode.clone(), tag.clone());
                        if let Some(ntc_rpm) = ntc_rpms.get(&key) {
                            let retain = (ntc_rpm / record.rpm) <= ratio_threshold; // Retain evidence if the NTC RPM / Library RPM ratio is larger than the provided threshold
                            
                            if retain {
                                if taxon.name == "Aquabacterium sp001770725" {
                                    log::info!(
                                        "Retaining: taxon.name = {}, record.rpm = {:.2}, ntc_rpm = {:.2}, ratio = {:.2}, threshold = {}, retained = {} tool={:?} mode={:?} id={}",
                                        taxon.name,
                                        record.rpm,
                                        ntc_rpm,
                                        (ntc_rpm / record.rpm),
                                        ratio_threshold,
                                        retain,
                                        record.tool,
                                        record.mode,
                                        record.id
                                    );
                                }
                                
                            }

                            return retain; 
                            
                        }
                        
                    }
                }
                true // Keep if no NTC comparison or ratio is specified
            })
            .collect();

        // Step 3: Apply remaining filters for tools, modes, and thresholds -
        // must happen after NTC comparison filter otherwise we would exclude 
        // evidence records from NTC
        taxon.evidence.profile = taxon
            .evidence
            .profile
            .into_iter()
            .filter(|record| {
                (filter_config.tools.is_empty() || filter_config.tools.contains(&record.tool))
                &&
                (
                    // Either the modes filter is empty...
                    filter_config.modes.is_empty()
                    // ...or the record's mode is directly contained in the modes filter...
                    || filter_config.modes.contains(&record.mode)
                    // ...or the modes filter contains Mixed and the record matches one of the allowed tool/mode pairs.
                    || (
                        filter_config.modes.contains(&AbundanceMode::Mixed)
                        &&
                        (
                            (record.tool == ProfileTool::Bracken   && record.mode == AbundanceMode::Profile)
                            || (record.tool == ProfileTool::Kraken2   && record.mode == AbundanceMode::Sequence)
                            || (record.tool == ProfileTool::Ganon2    && record.mode == AbundanceMode::Sequence)
                            || (record.tool == ProfileTool::Sylph     && record.mode == AbundanceMode::Sequence)
                            || (record.tool == ProfileTool::Kmcp      && record.mode == AbundanceMode::Sequence)
                            || (record.tool == ProfileTool::Metabuli  && record.mode == AbundanceMode::Sequence)
                            || (record.tool == ProfileTool::Blast     && record.mode == AbundanceMode::Bases)
                            || (record.tool == ProfileTool::Vircov    && record.mode == AbundanceMode::Sequence)
                        )
                    )
                )
                &&
                if record.mode != AbundanceMode::Bases { 
                    // For non-Bases records, check min_reads and min_rpm, and optionally enforce max_rpm
                    record.reads >= filter_config.min_reads && 
                    record.rpm >= filter_config.min_rpm &&
                    filter_config.max_rpm.map_or(true, |max| record.rpm <= max)
                } else {
                    // For Bases records, check min_bases and min_bpm, and optionally enforce max_bases
                    record.bases >= filter_config.min_bases && 
                    record.bpm >= filter_config.min_bpm &&
                    filter_config.max_bases.map_or(true, |max| record.bases <= max)
                }
                &&
                record.abundance >= filter_config.min_abundance
            })
            .collect();
            taxon
        })
        .filter(|taxon| if allow_no_evidence { true } else { !taxon.evidence.profile.is_empty() }) // Remove taxa with no evidence left (or retain if allow_no_evidence is true)
        .collect()

}

pub fn apply_lineage_filters(
    taxa: Vec<Taxon>,
    lineage_filters: &Vec<LineageFilterConfig>,
    sample_tags: &HashMap<String, Vec<String>> // keys: record.id, values: tags for the sample
) -> Vec<Taxon> {
    taxa.into_iter()
        .filter(|taxon| {
            let lineage_str = &taxon.lineage;

            // Step 0: If any of the lineage filters specify required tags,
            // filter evidence records by matching sample tags.
            // For each taxon, create a filtered evidence profile that only retains
            // records whose associated sample tags include one or more required tags.
            // We combine all required tags from filters that match this taxon’s lineage.
            let required_tags: Vec<String> = lineage_filters.iter()
                .filter(|lf| {
                    // Only consider filters that list tag requirements and whose lineages match.
                    !lf.tags.is_empty() && lf.lineages.iter().any(|l| lineage_str.contains(l))
                })
                .flat_map(|lf| lf.tags.iter().cloned())
                .collect();

            // If there are required tags, filter the profile records accordingly.
            // If none are required, keep all records.
            let filtered_profile: Vec<_> = if !required_tags.is_empty() {
                taxon.evidence.profile.iter().filter(|record| {
                    if let Some(record_tags) = sample_tags.get(&record.id) {
                        // Retain record if at least one of its tags is in the required set.
                        record_tags.iter().any(|tag| required_tags.contains(tag))
                    } else {
                        // No sample tags available for this record: exclude it.
                        false
                    }
                }).cloned().collect()
            } else {
                taxon.evidence.profile.clone()
            };

            // If the filtering removed all evidence, we fail this taxon.
            if filtered_profile.is_empty() {
                return false;
            }

            // Step 1: Apply remaining thresholds per lineage filter,
            // using the filtered evidence records.
            let mut passes_all = true;
            for lf in lineage_filters {
                if lf.lineages.iter().any(|l| lineage_str.contains(l)) {

                    // Alignment filter for Vircov records
                    if let Some(min_alignment_tools) = lf.min_alignment_tools {
                        let min_alignment_rpm = lf.min_alignment_rpm.unwrap_or(0.0);
                        
                        let valid_alignment_count = taxon.evidence.alignment.iter().filter(|record| {
                            if record.scan_coverage < lf.min_alignment_regions_coverage.unwrap_or(0.0) {
                                record.scan_regions >= lf.min_alignment_regions.unwrap_or(0)
                            } else {
                                true
                            }
                        }).count();

                        let alignment_count = filtered_profile.iter().filter(|record| {
                            record.tool == ProfileTool::Vircov && record.rpm >= min_alignment_rpm
                        }).count();
                        
                        if alignment_count < min_alignment_tools || valid_alignment_count < 1 {
                            passes_all = false;
                            break;
                        }
                    }

                    // K-mer filter for non-Vircov/Blast records
                    if let Some(min_kmer_tools) = lf.min_kmer_tools {
                        let min_kmer_rpm = lf.min_kmer_rpm.unwrap_or(0.0);
                        let kmer_count = filtered_profile.iter().filter(|record| {
                            match record.tool {
                                ProfileTool::Vircov | ProfileTool::Blast => false,
                                _ => record.rpm >= min_kmer_rpm,
                            }
                        }).count();
                        if kmer_count < min_kmer_tools {
                            passes_all = false;
                            break;
                        }
                    }

                    // Assembly filter for Blast records
                    if let Some(min_assembly_tools) = lf.min_assembly_tools {
                        let assembly_count = filtered_profile.iter().filter(|record| {
                            record.tool == ProfileTool::Blast
                        }).count();
                        if assembly_count < min_assembly_tools {
                            passes_all = false;
                            break;
                        }
                    }
                }
            }
            passes_all
        })
        .collect()
}
