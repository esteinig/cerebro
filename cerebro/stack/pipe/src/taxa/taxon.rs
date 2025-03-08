use anyhow::Result;
use fancy_regex::Regex;
use vircov::vircov::VircovRecord;
use serde::{Deserialize, Serialize};
use std::{path::PathBuf, collections::HashMap};
use taxonomy::{Taxonomy, GeneralTaxonomy, TaxRank};

use crate::modules::assembly::ContigRecord;
use crate::modules::pathogen::ProfileRecord;
use crate::error::WorkflowError;


pub trait TaxonExtraction {
    fn get_taxa(&self, taxonomy_directory: &PathBuf, strict: bool) -> Result<Vec<Taxon>, WorkflowError>;
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Taxon {
    pub taxid: String,
    pub rank: TaxRank,
    pub name: String,
    pub lineage: TaxonLineage,
    pub evidence: TaxonEvidence
}
impl Taxon {
    pub fn from_taxid(taxid: String, taxonomy: &GeneralTaxonomy) -> Result<Self, WorkflowError> {

        let taxid_str = taxid.as_str();

        let rank = taxonomy.rank(taxid_str)
            .map_err(|err|WorkflowError::TaxRankNotAvailable(err, taxid.clone()))?;

        let name = taxonomy.name(taxid_str)
            .map_err(|err|WorkflowError::TaxNameNotAvailable(err, taxid.clone()))?.to_string();

        Ok(Self { 
            taxid: taxid.clone(), 
            rank,
            name,
            lineage: TaxonLineage::from_taxid(taxid_str, taxonomy)?,
            evidence: TaxonEvidence::new() // to be filled

        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TaxonEvidence {
    pub alignment: Vec<VircovRecord>,
    pub assembly: Vec<ContigRecord>,
    pub profile: Vec<ProfileRecord>,
}

impl TaxonEvidence {
    pub fn new() -> Self {
        Self { 
            alignment: Vec::new(),
            assembly: Vec::new(),
            profile: Vec::new()
        }
    }
}

// Lineage in GTDB liek format

pub trait LineageOperations: AsRef<str> {
    fn from_taxid(taxid: &str, taxonomy: &GeneralTaxonomy) -> Result<String, WorkflowError> {
        
        match taxonomy.rank(taxid) {
            Ok(_) => {

                let lineage = vec![
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Superkingdom)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Phylum)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Class)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Order)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Family)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Genus)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Species)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                ]
                .into_iter()
                .map(String::from)
                .collect::<Vec<_>>();

                Ok(lineage_to_gtdb_str(lineage)?)
            },
            Err(err) => {
                log::warn!("Error from taxonomy: {}", err.to_string());
                log::warn!("Rank could not be determined for taxid {taxid} - using unclassified lineage (d__;p__;c__;o__;f__;g__;s__)");
                Ok(String::from("d__;p__;c__;o__;f__;g__;s__"))
            }
        }
    }

    /// Helper method: Returns the taxon name for a given prefix (e.g. "d__", "p__", etc.)
    fn get_taxon_by_prefix(&self, prefix: &str) -> Option<String> {
        self.as_ref()
            .split(";")
            .find(|s| s.trim_start().starts_with(prefix))
            .map(|s| s.trim_start_matches(prefix).to_string())
    }

    fn get_domain(&self) -> Option<String> {
        self.get_taxon_by_prefix("d__")
    }

    fn get_phylum(&self) -> Option<String> {
        self.get_taxon_by_prefix("p__")
    }

    fn get_class(&self) -> Option<String> {
        self.get_taxon_by_prefix("c__")
    }

    fn get_order(&self) -> Option<String> {
        self.get_taxon_by_prefix("o__")
    }

    fn get_family(&self) -> Option<String> {
        self.get_taxon_by_prefix("f__")
    }

    fn get_genus(&self) -> Option<String> {
        self.get_taxon_by_prefix("g__")
    }

    fn get_species(&self) -> Option<String> {
        self.get_taxon_by_prefix("s__")
    }

    fn get_labels(&self) -> Vec<&str> {
        self.as_ref().split(";").collect()
    }
}


pub fn lineage_to_gtdb_str(lineage: Vec<String>) -> Result<String, WorkflowError> {
    if lineage.len() >= 7 {
        let base_lineage = format!(
            "d__{};p__{};c__{};o__{};f__{};g__{};s__{}",
            lineage[0], lineage[1], lineage[2], lineage[3], lineage[4], lineage[5], lineage[6]
        );
        Ok(base_lineage)
    } else {
        Err(WorkflowError::LineageStringTooShort(lineage.join(", ")))
    }
}

pub type TaxonLineage = String;

impl LineageOperations for TaxonLineage {}

// Utility function wrapping the search result handling with name and taxid
pub fn get_taxid_name_from_search(search_result: Option<(&str, f32)>, taxonomy: &GeneralTaxonomy) -> Result<(Option<String>, Option<String>), WorkflowError> {
     match search_result {
        Some((taxid, _)) => {
            let name = taxonomy.name(taxid).map_err(
                |err|WorkflowError::TaxNameNotAvailable(err, taxid.to_owned())
            )?;
            Ok((Some(taxid.to_owned()), Some(name.to_owned())))

        },
        None => Ok((None, None))
    }
}

/* 
=======================
MODULE AND TOOL STRUCTS
=======================
*/

// Utility function to aggregate multiple HashMaps with Taxon
pub fn aggregate(parent_taxa: &mut HashMap<String, Taxon>, taxa: &Vec<Taxon>) -> HashMap<String, Taxon> {
    
    for (taxid, taxon) in taxa.iter().map(|t| (t.taxid.as_str(), t)) {
        // ... check if each taxon exists
        match parent_taxa.get_mut(taxid) {
            Some(tax) => {
                // If so, update the evidence of the taxon in the modules...
                for record in &taxon.evidence.profile {
                    tax.evidence.profile.push(record.clone())
                }
                // If so, update the evidence of the taxon in the modules...
                for record in &taxon.evidence.alignment {
                    tax.evidence.alignment.push(record.clone())
                }
                // If so, update the evidence of the taxon in the modules...
                for record in &taxon.evidence.assembly {
                    tax.evidence.assembly.push(record.clone())
                }
            },
            None => {
                // otherwise insert the taxon from result with its own evidence
                parent_taxa.insert(taxid.to_owned(), taxon.to_owned());
            }
        }
    }  
    parent_taxa.to_owned()
}

pub struct TaxonThresholdConfig {
    min_rpm: f64,
    min_rpm_kmer: f64,
    min_rpm_alignment: f64,
    min_rpm_remap: f64,
    min_contigs: u64,
    min_bases: u64
}
impl TaxonThresholdConfig {
    // pub fn from_args(args: &PipelineTaxaArgs) -> Self {
    //     Self {
    //         min_rpm: args.min_rpm,
    //         min_rpm_kmer: args.min_rpm_kmer,
    //         min_rpm_alignment: args.min_rpm_alignment,
    //         min_rpm_remap: args.min_rpm_remap,
    //         min_contigs: args.min_contigs,
    //         min_bases: args.min_bases
    //     }
    // }
}
impl Default for TaxonThresholdConfig {
    fn default() -> Self {
        Self {
            min_rpm: 0.0,
            min_rpm_kmer: 0.0,
            min_rpm_alignment: 0.0,
            min_rpm_remap: 0.0,
            min_contigs: 0,
            min_bases: 0
        }
    }
}

pub fn taxa_summary(samples: Vec<PathBuf>, output: &PathBuf, sep: char, header: bool, extract: bool, filter_config: Option<PathBuf>, threshold_config: &TaxonThresholdConfig) -> Result<(), WorkflowError> {

    // let mut wtr = csv::WriterBuilder::new()
    //     .delimiter(sep as u8)
    //     .has_headers(header)
    //     .from_path(output)
    //     .map_err(|err| WorkflowError::CreateTaxaOutputFile(err))?;

    // let filter_config = match filter_config {
    //     Some(file) => TaxonFilterConfig::from_path(&file)?,
    //     None => TaxonFilterConfig::default()
    // };

    // for file in samples {
    //     let sample = WorkflowSample::read_json(&file).expect(&format!("Failed to parse sample file: {}", file.display()));
    //     let taxa: Vec<_> = apply_filters(sample.taxa.clone().into_values().collect(), &filter_config);
    //     let taxa_overview: Vec<_> = taxa.iter().map(|taxon| { TaxonOverview::from(taxon) }).collect();
    //     let taxa_overview_summary: Vec<_> = taxa_overview.iter().map(|taxon_overview|{
    //         TaxonSampleOverview::from_taxon_overview(&sample, &taxon_overview, extract).expect("Failed to construct taxon sample overview!") // Fix to error bubbles!
    //     }).filter(|taxon_sample_overview| taxon_sample_overview.pass(&threshold_config)).collect();

    //     if taxa_overview_summary.is_empty() {
    //         log::warn!("No taxa for sample {} - sample is excluded from results table!", sample.id);
    //     }

    //     for record in taxa_overview_summary {
    //         wtr.serialize(record).map_err(|err| WorkflowError::SerializeTaxaTableRow(err))?;
    //     }
    // }

    Ok(())
}


// #[derive(Debug, Clone, Serialize, Deserialize)]
// // A struct representing a high level overview
// // of the aggregated taxon evidence
// pub struct TaxonOverview {
//     pub taxid: String,
//     pub name: String,                    // used to later map back the tags
//     pub domain: Option<String>,
//     pub genus: Option<String>,
//     pub evidence: Vec<ProfileRecord>,
//     pub kmer: bool,
//     pub alignment: bool,
//     pub assembly: bool,
//     pub sample_names: Vec<String>        // the evidence record associated sample names as processed in the pipeline (matching `Cerebro.name`)
// }
// impl TaxonOverview {
//     pub fn from(taxon: &Taxon) -> Self {

//         let (kmer, alignment, assembly) = (false, false, false);

//         let mut results = Vec::new();
//         for record in &taxon.evidence.profile {
//             results.push(record.to_owned())
//         }

//         // Unique record identifiers (sample names)
//         let mut names = Vec::new();
        
//         for record in &taxon.evidence.profile {
//             names.push(record.id.to_owned())
//         };

//         let names = names.into_iter().unique().collect();

//         Self {
//             taxid: taxon.taxid.to_owned(),
//             name: taxon.name.to_owned(),
//             genus: taxon.level.genus.to_owned(),
//             domain: taxon.level.domain.to_owned(),
//             evidence: results,
//             kmer,
//             alignment,
//             assembly,
//             sample_names: names
//         }
//     }
// }

// Utility function to extract the biological sample identifier and library tags [strict]
fn get_sample_regex_matches(file_name: &String)-> Result<(String, Vec<String>), WorkflowError> {
    
    let mut sample_id = String::new();
    let sample_id_regex = Regex::new("^([^_]+)(?=__)").map_err(WorkflowError::SampleIdRegex)?;
    for caps in sample_id_regex.captures_iter(&file_name) {
        let sample_name = caps.map_err(WorkflowError::SampleIdRegexCapture)?.get(0);

        sample_id = match sample_name {
            Some(name) => name.as_str().to_string(),
            None => return Err(WorkflowError::SampleIdRegexCaptureMatch)
        };
        break; // always use the first capture
    }

    let mut library_tags: Vec<String> = Vec::new();
    let sample_id_regex = Regex::new("(?<=__)([A-Za-z0-9]+)").map_err(WorkflowError::SampleTagRegex)?;
    for caps in sample_id_regex.captures_iter(&file_name) {
        let tags = caps.map_err(WorkflowError::SampleTagRegexCapture)?.get(0);
        let tag = match tags {
            Some(name) => name.as_str().to_string(),
            None => return Err(WorkflowError::SampleTagRegexCaptureMatch)
        };
        library_tags.push(tag)
    }
    Ok((sample_id, library_tags))
}
