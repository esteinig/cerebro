use serde::{Deserialize, Serialize};
use taxonomy::TaxRank;

use crate::modules::pathogen::PathogenDetectionRank;
use crate::modules::pathogen::PathogenDetectionTool;
use crate::modules::pathogen::PathogenDetectionMode;

use super::taxon::Taxon;


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TaxonFilterConfig {
    pub rank: Option<PathogenDetectionRank>,      // Filter by specific taxonomic rank
    pub domains: Vec<String>,                     // Filter by domain names
    pub tools: Vec<PathogenDetectionTool>,        // Filter by specific detection tools
    pub modes: Vec<PathogenDetectionMode>,        // Filter by detection modes (Sequence/Profile)
    pub min_reads: u64,                           // Minimum read count for inclusion
    pub min_rpm: f64,                             // Minimum RPM (Reads per million) for inclusion
    pub min_abundance: f64,                       // Minimum abundance for inclusion
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
        }
    }
}

// IMPORTANT: change from PathogenDetectionRank to TaxRank
pub fn apply_filters(mut taxa: Vec<Taxon>, filter_config: &TaxonFilterConfig) -> Vec<Taxon> {

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
                if let Some(domain_name) = &taxon.level.domain_name {
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
            taxon.evidence.records = taxon
                .evidence
                .records
                .into_iter()
                .filter(|record| {
                    record.results.iter().any(|result| {
                        (filter_config.tools.is_empty() || filter_config.tools.contains(&result.tool))
                            && (filter_config.modes.is_empty() || filter_config.modes.contains(&result.mode))
                            && result.reads >= filter_config.min_reads
                            && result.rpm >= filter_config.min_rpm
                            && result.abundance >= filter_config.min_abundance
                    })
                })
                .collect();
            taxon
        })
        .filter(|taxon| !taxon.evidence.records.is_empty()) // Remove taxa with no evidence left
        .collect();

    taxa
}
