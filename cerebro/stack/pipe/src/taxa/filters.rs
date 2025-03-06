use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use taxonomy::TaxRank;

use crate::modules::pathogen::PathogenDetectionRank;
use crate::modules::pathogen::ProfileTool;
use crate::modules::pathogen::AbundanceMode;

use super::taxon::LineageOperations;
use super::taxon::Taxon;


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TaxonFilterConfig {
    pub rank: Option<PathogenDetectionRank>,      // Filter by specific taxonomic rank
    pub domains: Vec<String>,                     // Filter by domain names
    pub tools: Vec<ProfileTool>,                  // Filter by specific detection tools
    pub modes: Vec<AbundanceMode>,                // Filter by detection modes (Sequence/Profile)
    pub min_bases: u64,
    pub min_bpm: f64,                             
    pub min_reads: u64,                           // Minimum read count for inclusion
    pub min_rpm: f64,                             // Minimum RPM (Reads per million) for inclusion
    pub min_abundance: f64,                       // Minimum abundance for inclusion
    pub ntc_ratio: Option<f64>,                   // NTC ratio threshold
}

impl Default for TaxonFilterConfig {
    fn default() -> Self {
        Self {
            rank: Some(PathogenDetectionRank::Species),
            domains: Vec::new(),
            tools: Vec::new(),
            modes: Vec::new(),
            min_bases: 0,
            min_bpm: 0.0,
            min_reads: 0,
            min_rpm: 0.0,
            min_abundance: 0.0,
            ntc_ratio: None
        }
    }
}

type SampleId = String;
type Tag = String;


pub fn apply_filters(mut taxa: Vec<Taxon>, filter_config: &TaxonFilterConfig, sample_tags: &HashMap<SampleId, Vec<Tag>>) -> Vec<Taxon> {

    let only_ntc = sample_tags.clone().into_values().all(|tags| tags.contains(&"NTC".to_string()));

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

    // Filter by detection tools, modes, and thresholds
    taxa = taxa
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
                    tags.contains(&"NTC".to_string()) // Find records with "NTC" tag
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
                    if tags.contains(&"NTC".to_string()) {
                        if only_ntc {
                            return true // If we have aggregated NTC evidence only, return the evidence
                        } else {
                            return false // If we have aggregated sample evidence, exclude the evidence from returned taxa
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
                                log::debug!(
                                    "Retaining: taxon.name = {}, record.rpm = {:.2}, ntc_rpm = {:.2}, ratio = {:.2}, threshold = {}, retained = {}",
                                    taxon.name,
                                    record.rpm,
                                    ntc_rpm,
                                    (ntc_rpm / record.rpm),
                                    ratio_threshold,
                                    retain
                                );
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
                    && (filter_config.modes.is_empty() || filter_config.modes.contains(&record.mode))
                    && if record.mode != AbundanceMode::Bases { record.reads >= filter_config.min_reads } else { record.bases >= filter_config.min_bases }
                    && if record.mode != AbundanceMode::Bases { record.rpm >= filter_config.min_rpm } else { record.bpm >= filter_config.min_bpm }
                    && record.abundance >= filter_config.min_abundance
            })
            .collect();

         taxon
     })
     .filter(|taxon| !taxon.evidence.profile.is_empty()) // Remove taxa with no evidence left
     .collect();

    taxa
}
