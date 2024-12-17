use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use taxonomy::TaxRank;

use crate::modules::pathogen::PathogenDetectionRank;
use crate::modules::pathogen::ProfileTool;
use crate::modules::pathogen::AbundanceMode;

use super::taxon::Taxon;


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TaxonFilterConfig {
    pub rank: Option<PathogenDetectionRank>,      // Filter by specific taxonomic rank
    pub domains: Vec<String>,                     // Filter by domain names
    pub tools: Vec<ProfileTool>,        // Filter by specific detection tools
    pub modes: Vec<AbundanceMode>,        // Filter by detection modes (Sequence/Profile)
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
            min_reads: 0,
            min_rpm: 0.0,
            min_abundance: 0.0,
            ntc_ratio: None
        }
    }
}

type SampleId = String;
type Tag = String;

// IMPORTANT: change from PathogenDetectionRank to TaxRank
pub fn apply_filters(mut taxa: Vec<Taxon>, filter_config: &TaxonFilterConfig, sample_tags: &HashMap<SampleId, Vec<Tag>>) -> Vec<Taxon> {

    // Filter by taxonomic rank
    if let Some(rank) = &filter_config.rank {
        let tax_rank: TaxRank = (*rank).clone().into(); // Convert PathogenDetectionRank to TaxRank
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
                if let Some(domain_name) = &taxon.level.domain {
                    filter_config.domains.contains(domain_name)
                } else {
                    false // Exclude taxa without a domain_name
                }
            })
            .collect();
    }

     // Filter by detection tools, modes, and thresholds
     taxa = taxa
     .into_iter()
     .map(|mut taxon| {

        // Step 1: Collect NTC RPMs for matching tools, modes, and tags
        let ntc_rpms: HashMap<(ProfileTool, AbundanceMode, String), f64> = taxon
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
            .filter_map(|record| {
                // Extract nucleic acid tag (DNA/RNA) if present
                if let Some(tags) = sample_tags.get(&record.id) {
                    let nucleic_acid_tag = tags
                        .iter()
                        .find(|tag| tag == &&"DNA".to_string() || tag == &&"RNA".to_string());
                    if let Some(tag) = nucleic_acid_tag {
                        Some(((record.tool.clone(), record.mode.clone(), tag.clone()), record.rpm))
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect();

        log::info!("{} {:#?}", taxon.name, ntc_rpms);
        
        // Step 2: Filter profile records based on ntc_ratio
        taxon.evidence.profile = taxon
            .evidence
            .profile
            .into_iter()
            .filter(|record| {
                if let Some(tags) = sample_tags.get(&record.id) {
                    let nucleic_acid_tag = tags
                        .iter()
                        .find(|tag| tag == &&"DNA".to_string() || tag == &&"RNA".to_string());
                    if let (Some(ratio), Some(tag)) = (filter_config.ntc_ratio, nucleic_acid_tag) {
                        let key = (record.tool.clone(), record.mode.clone(), tag.clone());
                        if let Some(ntc_rpm) = ntc_rpms.get(&key) {
                            return record.rpm <= ratio * ntc_rpm; // Compare RPMs
                        }
                    }
                }
                true // Keep if no NTC comparison or ratio is specified
            })
            .collect();
    
        // Step 3: Apply remaining filters for tools, modes, and thresholds
        taxon.evidence.profile = taxon
            .evidence
            .profile
            .into_iter()
            .filter(|record| {
                (filter_config.tools.is_empty() || filter_config.tools.contains(&record.tool))
                    && (filter_config.modes.is_empty() || filter_config.modes.contains(&record.mode))
                    && record.reads >= filter_config.min_reads
                    && record.rpm >= filter_config.min_rpm
                    && record.abundance >= filter_config.min_abundance
            })
            .collect();

         taxon
     })
     .filter(|taxon| !taxon.evidence.profile.is_empty()) // Remove taxa with no evidence left
     .collect();

    taxa
}
