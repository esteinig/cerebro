use std::collections::HashMap;
use serde::{Deserialize, Serialize};

use cerebro_pipe::taxa::taxon::Taxon;

#[derive(Deserialize)]
pub struct TaxaSummaryMongoPipeline {
    pub cerebro_id: String,
    pub workflow_id: String,
    pub workflow_name: String,
    pub run_id: String,
    pub run_date: String,
    pub sample_id: String,
    pub sample_tags: Vec<String>,
    pub sample_type: String,
    pub sample_group: String,
    pub taxa: HashMap<String, Taxon>
}


// // TaxonOverview + TaxaSummary
// #[derive(Deserialize, Serialize, Debug)]
// pub struct TaxonSummaryOverview {
//     pub cerebro_id: String,
//     pub workflow_id: String,
//     pub workflow_name: String,
//     pub run_id: String,
//     pub run_date: String,
//     pub sample_id: String,
//     pub sample_tag: String,
//     pub sample_type: String,
//     pub sample_group: String,
//     pub taxid: String,
//     pub domain: Option<String>,
//     pub genus: Option<String>,
//     pub name: String,                // used to later map back the tags
//     // pub rpm: f64,                 // total rpm summed from k-mer and alignment evidence
//     // pub rpm_kmer: f64,            // total rpm summed from k-mer evidence
//     // pub rpm_alignment: f64,       // total rpm summed from  alignment evidence
//     // pub contigs: u64,             // total assembled and identified contig evidence
//     // pub contigs_bases: u64,
//     pub kmer: bool,
//     pub alignment: bool,
//     pub assembly: bool
// }
// impl TaxonSummaryOverview {
//     pub fn from_taxon_overview(taxa_summary: &TaxaSummaryMongoPipeline, taxon_overview: &TaxonOverview) -> Self {
//         Self {
//             cerebro_id: taxa_summary.cerebro_id.clone(),
//             workflow_id: taxa_summary.workflow_id.clone(),
//             workflow_name: taxa_summary.workflow_name.clone(),
//             run_id: taxa_summary.run_id.clone(),
//             run_date: taxa_summary.run_date.clone(),
//             sample_id: taxa_summary.sample_id.clone(),
//             sample_tag: taxa_summary.sample_tags.join("-"),
//             sample_type: taxa_summary.sample_type.clone(),
//             sample_group: taxa_summary.sample_group.clone(),
//             taxid: taxon_overview.taxid.clone(),
//             domain: taxon_overview.domain.clone(),
//             genus: taxon_overview.genus.clone(),
//             name: taxon_overview.name.clone(),            
//             // rpm: taxon_overview.rpm,                 
//             // rpm_kmer: taxon_overview.rpm_kmer,            
//             // rpm_alignment: taxon_overview.rpm_alignment,  
//             // contigs: taxon_overview.contigs,             
//             // contigs_bases: taxon_overview.contigs_bases,
//             kmer: taxon_overview.kmer,
//             alignment: taxon_overview.alignment,
//             assembly: taxon_overview.assembly
//         }
//     }
// }


// #[derive(Deserialize)]
// pub struct TaxaSummaryDataResponse {
//     pub status: String,
//     pub message: String,
//     pub data: TaxaSummaryData
// }
// #[derive(Deserialize)]
// pub struct TaxaSummaryData {
//     pub taxa_summary: Vec<TaxonSummaryOverview>,
//     pub csv: String,
// }

