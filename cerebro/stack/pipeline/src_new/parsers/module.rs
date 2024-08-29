use std::fs::File;
use std::path::PathBuf;
use std::io::{BufReader, BufRead};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use taxonomy::{Taxonomy, GeneralTaxonomy};

use crate::error::WorkflowError;
use crate::output::WorkflowOutputs;
use crate::quality::Phage;
use crate::record::VircovRecord;
use crate::quality::{Fastp, Ercc, Scrubby, Nanoq};
use crate::taxon::{Taxon, TaxonomyWarning, aggregate};
use crate::tools::{Kraken2Uniq, Blast, Diamond, Vircov, VircovScanRemap};
use crate::virus::AnnotationOptions;


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityControlModule {
    pub id: String,
    pub nanoq: Option<Nanoq>,
    pub fastp: Option<Fastp>,
    pub nanoq_scan: Option<Nanoq>,
    pub fastp_scan: Option<Fastp>,
    pub ercc: Option<Ercc>,
    pub phage: Option<Vec<Phage>>,
    pub host_background: Option<Scrubby>,
    pub virus_background: Option<Scrubby>
}
impl QualityControlModule {
    pub fn from(
        id: String,
        files: &WorkflowOutputs
    ) -> Result<Self, WorkflowError>{
        
        let nanoq_scan_results = match &files.nanoq_scan {
            Some(file) => {
                log::info!("Results detected for Nanoq read scanning");
                let input = File::open(file)?;
                let buffered = BufReader::new(input);
                let results: Nanoq = serde_json::from_reader(buffered).map_err(WorkflowError::ParseNanoq)?;
                Some(results)
            },
            None => {
                log::info!("No results detected for Nanoq read scanning");
                None
            }
        };

        let nanoq_results = match &files.nanoq {
            Some(file) => {
                log::info!("Results detected for Nanoq read quality control");
                let input = File::open(file)?;
                let buffered = BufReader::new(input);
                let results: Nanoq = serde_json::from_reader(buffered).map_err(WorkflowError::ParseNanoq)?;
                Some(results)
            },
            None => {
                log::info!("No results detected for Nanoq read quality control");
                None
            }
        };

        let fastp_scan_results = match &files.fastp_scan {
            Some(file) => {
                log::info!("Results detected for Fastp read scanning");
                let input = File::open(file)?;
                let buffered = BufReader::new(input);
                let results: Fastp = serde_json::from_reader(buffered).map_err(WorkflowError::ParseFastp)?;
                Some(results)
            },
            None => {
                log::info!("No results detected for Fastp read scanning");
                None
            }
        };

        let fastp_results = match &files.fastp {
            Some(file) => {
                log::info!("Results detected for Fastp read quality control");
                let input = File::open(file)?;
                let buffered = BufReader::new(input);
                let results: Fastp = serde_json::from_reader(buffered).map_err(WorkflowError::ParseFastp)?;
                Some(results)
            },
            None => {
                log::info!("No results detected for Fastp read quality control");
                None
            }
        };

        let ercc_results = match (&files.ercc_vircov, &files.ercc_scrubby) {
            (Some(vircov_file), Some(scrubby_file)) => {
                log::info!("Results detected for ERCC Control");

                let ercc_vircov_reader = BufReader::new(File::open(vircov_file)?);
                let mut vircov_records: Vec<VircovRecord> = vec![];
                for line in ercc_vircov_reader.lines() {
                    vircov_records.push(VircovRecord::from_str(line?, id.to_string(), "ercc".to_string(), &None,&None)?)
                }
                let buffered = BufReader::new(File::open(scrubby_file)?);
                let scrubby_results: Scrubby = serde_json::from_reader(buffered).map_err(WorkflowError::ParseScrubby)?;

                Some(Ercc::from(&vircov_records, &scrubby_results))
            },
            _ => {
                log::info!("No results detected for ERCC Control");
                None
            }
        };

        let phage_results = match (&files.phage_vircov, &files.phage_scrubby) {
            (Some(vircov_file), Some(scrubby_file)) => {
                log::info!("Results detected for phage control");

                let vircov_reader = BufReader::new(File::open(vircov_file)?);
                let mut vircov_records: Vec<VircovRecord> = vec![];
                for line in vircov_reader.lines() {
                    vircov_records.push(VircovRecord::from_str(line?, id.to_string(), "phage".to_string(), &None,&None)?)
                }

                let buffered = BufReader::new( File::open(scrubby_file)?);
                let _: Scrubby = serde_json::from_reader(buffered).map_err(WorkflowError::ParseScrubby)?;

                Some(vircov_records.iter().map(|record| Phage::from(record)).collect::<Vec<Phage>>())
            },
            _ => {
                log::info!("No results detected for ERCC Control");
                None
            }
        };

        let host_results = match &files.host_scrubby {
            Some(file) => {
                log::info!("Results detected for host depletion");
                let buffered = BufReader::new(File::open(file)?);
                let results: Scrubby = serde_json::from_reader(buffered).map_err(WorkflowError::ParseScrubby)?;
                Some(results)
            },
            None => {
                log::info!("No results detected for host background depletion");
                None
            }
        };

        let virus_results = match &files.virus_scrubby {
            Some(file) => {
                log::info!("Results detected for virus background depletion");
                let input = File::open(file)?;
                let buffered = BufReader::new(input);
                let results: Scrubby = serde_json::from_reader(buffered).map_err(WorkflowError::ParseScrubby)?;
                Some(results)
            },
            None => {
                log::info!("No results detected for virus background depletion");
                None
            }
        };


        Ok(Self { 
            id: id.to_owned(), 
            nanoq: nanoq_results, 
            nanoq_scan: nanoq_scan_results,
            fastp: fastp_results, 
            fastp_scan: fastp_scan_results,
            ercc: ercc_results, 
            phage: phage_results, 
            host_background: host_results, 
            virus_background: virus_results
        })

    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerModule {
    pub taxa: HashMap<String, Taxon>,
    pub kraken2uniq: Option<Vec<Kraken2Uniq>>,
    pub taxonomy_warnings: Vec<TaxonomyWarning>
}
impl KmerModule {
    pub fn new() -> Self {
        Self {
            taxa: HashMap::new(),
            kraken2uniq: None,
            taxonomy_warnings: Vec::new()
        }
    }
    pub fn from(id: String, kraken_reports: Option<Vec<PathBuf>>, taxonomy: &GeneralTaxonomy) -> Result<Self, WorkflowError> {
        // Multiple Kraken2Uniq results due to using multiple databases
        let kraken2uniq = match kraken_reports {
            Some(_kraken_reports) => {
                let mut reports = Vec::new();
                for report in _kraken_reports {
                    let mut data = Kraken2Uniq::from(id.to_owned(), &report)?;
                    log::info!("Results detected for Kraken2Uniq with database `{}`", data.db);
                    // Create the taxon instances from the records
                    data.get_taxa(taxonomy)?;
                    reports.push(data)
                }
                Some(reports)
            },
            None => {
                log::info!("No results detected for Kraken2Uniq");
                None
            }
        };
        KmerModule::aggregate(&mut Self { taxa: HashMap::new(), kraken2uniq, taxonomy_warnings: Vec::new() })
    }
    /// Aggregates taxa from all tools in this module
    fn aggregate(&mut self) -> Result<Self, WorkflowError> {
        // Raise error if this module is already aggregated and function is called another time,
        // otherwise evidence would get updated again. This is because we use vectors for keeping
        // evidence instead of HashSets, as some float values may occurr in data records which
        // do not implement Hash/Eq/PartialEq
        if !self.taxa.is_empty() {
            return Err(WorkflowError::TaxAggregate("KmerModule".to_string()))
        }
    
        if let Some(reports) = &self.kraken2uniq {
            for result in reports {
                if result.taxa.is_empty() { log::warn!("Attempting to aggregate taxa from Kraken2Uniq in KmerModule, but no taxa found (DB: {})", result.db) }
                // For each Kraken2Uniq result (e.g. from different databases) ...
                self.taxa = aggregate(&mut self.taxa, &result.taxa);
                self.taxonomy_warnings.append(&mut result.taxonomy_warnings.clone())
            }
        }
        Ok(self.clone())
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentModule {
    pub taxa: HashMap<String, Taxon>,
    pub alignment: Option<Vec<VircovScanRemap>>,
    pub taxonomy_warnings: Vec<TaxonomyWarning>  // bubbled up taxonomy warnings
}
impl AlignmentModule {
    pub fn new() -> Self {
        Self {
            taxa: HashMap::new(),
            alignment: None,
            taxonomy_warnings: Vec::new()
        }
    }
    pub fn from(id: String, vircov_reports: Option<Vec<PathBuf>>, vircov_scanremap_reports: Option<Vec<PathBuf>>, taxonomy: &GeneralTaxonomy, annotation_options: AnnotationOptions) -> Result<Self, WorkflowError> {

        let vircov_standard_data = AlignmentModule::from_vircov(id.clone(), vircov_reports, taxonomy, &annotation_options)?;
        let vircov_scanremap_data = AlignmentModule::from_vircov_remap(id.clone(), vircov_scanremap_reports, taxonomy, &annotation_options)?;

        let alignment = match (vircov_standard_data, vircov_scanremap_data) {
            (Some(standard_alignments), Some(scanremap_alignments)) => Some([standard_alignments, scanremap_alignments].concat()),
            (Some(standard_alignments), None) => Some(standard_alignments),
            (None, Some(scanremap_alignments)) => Some(scanremap_alignments),
            (None, None) => None
        };

        AlignmentModule::aggregate(&mut Self { taxa: HashMap::new(), alignment, taxonomy_warnings: Vec::new() })
    }

    pub fn from_vircov(id: String, vircov_reports: Option<Vec<PathBuf>>, taxonomy: &GeneralTaxonomy, annotation_options: &AnnotationOptions) -> Result<Option<Vec<VircovScanRemap>>, WorkflowError> {
        // Multiple Vircov results due to using multiple databases
        match vircov_reports {
            Some(_reports) => {
                let mut reports = Vec::new();
                for report in _reports {
                    let mut data = VircovScanRemap::from_vircov(
                        id.to_owned(), Vircov::from(id.to_owned(), 
                        &report, &annotation_options, true
                        )?, &annotation_options
                    )?;
                    log::info!("Results detected for alignment with database `{}`  ", data.db);
                    // Create the taxon instanced from records
                    data.get_taxa(taxonomy)?;
                    reports.push(data)
                }
                Ok(Some(reports))
            },
            None => {
                log::info!("No results detected for alignment");
                Ok(None)
            }
        }
    }

    pub fn from_vircov_remap(id: String, vircov_scanremap_reports: Option<Vec<PathBuf>>, taxonomy: &GeneralTaxonomy, annotation_options: &AnnotationOptions) -> Result<Option<Vec<VircovScanRemap>>, WorkflowError> {
        // Multiple Vircov results due to using multiple databases
        match vircov_scanremap_reports {
            Some(_reports) => {
                let mut reports = Vec::new();
                for report in _reports {
                    let mut data = VircovScanRemap::from(id.to_owned(), &report, &annotation_options, true)?;
                    log::info!("Results detected for alignment with database `{}`  ", data.db);
                    // Create the taxon instanced from records
                    data.get_taxa(taxonomy)?;
                    reports.push(data)
                }
                Ok(Some(reports))
            },
            None => {
                log::info!("No results detected for alignment");
                Ok(None)
            }
        }
    }
    /// Aggregates taxa from all tools in this module
    fn aggregate(&mut self) -> Result<Self, WorkflowError> {
        // Raise error if this module is already aggregated and function is called another time,
        // otherwise evidence would get updated again. This is because we use vectors for keeping
        // evidence instead of HashSets, as some float values may occurr in data records which
        // do not implement Hash/Eq/PartialEq
        if !self.taxa.is_empty() {
            return Err(WorkflowError::TaxAggregate("AlignmentModule".to_string()))
        }
        if let Some(reports) = &self.alignment {
            for result in reports {
                if result.taxa.is_empty() { log::warn!("Attempting to aggregate taxa from Vircov in AlignmentModule, but no taxa found (DB: {})", result.db) }
                // For each Vircov result (e.g. from different databases) ...
                self.taxa = aggregate(&mut self.taxa, &result.taxa);     
                self.taxonomy_warnings.append(&mut result.taxonomy_warnings.clone())
            }
        }
        Ok(self.clone())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AssemblyTools {
    Blast(Blast),
    Diamond(Diamond)
}

// Struct to hold data and produce summary data
// for the quality control module
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyModule {
    pub taxa: HashMap<String, Taxon>,
    pub assembly: Option<Vec<AssemblyTools>>,
    pub taxonomy_warnings: Vec<TaxonomyWarning>  // bubbled up taxonomy warnings
}
impl AssemblyModule {

    pub fn new() -> Self {
        Self {
            taxa: HashMap::new(),
            assembly: None,
            taxonomy_warnings: Vec::new()
        }
    }
    pub fn from(id: String, blastn_reports: Option<Vec<PathBuf>>, diamond_reports: Option<Vec<PathBuf>>, taxonomy: &GeneralTaxonomy) -> Result<Self, WorkflowError> {
        // Multiple Blastn results if using multiple databases
        let mut assembly: Vec<AssemblyTools> = Vec::new();
        match blastn_reports {
            Some(_reports) => {
                for report in _reports {
                    let mut data = Blast::from(id.to_owned(), &report, taxonomy)?;
                    log::info!("Results detected for BLASTn with database `{}`", data.db);
                    // Create the taxon instanced from records
                    data.get_taxa(taxonomy)?;
                    assembly.push(AssemblyTools::Blast(data))
                }
            },
            None => log::info!("No results detected for BLASTn")
        };
        match diamond_reports {
            Some(_reports) => {
                for report in _reports {
                    let mut data = Diamond::from(id.to_owned(), &report, taxonomy)?;
                    log::info!("Results detected for Diamond with database `{}`", data.db);
                    // Create the taxon instanced from records
                    data.get_taxa(taxonomy)?;
                    assembly.push(AssemblyTools::Diamond(data))
                }
            },
            None => log::info!("No results detected for Diamond")
        };

        AssemblyModule::aggregate(&mut Self { taxa: HashMap::new(), assembly: match assembly.is_empty() { true => None, false => Some(assembly) }, taxonomy_warnings: Vec::new() })
    }
    /// Aggregates taxa from all tools in this module
    fn aggregate(&mut self) -> Result<Self, WorkflowError> {

        // Raise error if this module is already aggregated and function is called another time,
        // otherwise evidence would get updated again. This is because we use vectors for keeping
        // evidence instead of HashSets, as some float values may occurr in data records which
        // do not implement Hash/Eq/PartialEq

        if !self.taxa.is_empty() {
            return Err(WorkflowError::TaxAggregate("AlignmentModule".to_string()))
        }
        if let Some(reports) = &self.assembly {
            for result in reports {
                match result {
                    AssemblyTools::Blast(result_data) => {
                        if result_data.taxa.is_empty() { log::warn!("Attempting to aggregate taxa from Blast in AssemblyModule, but no taxa found (DB: {})", result_data.db) }
                        // For each Vircov result (e.g. from different databases) ...
                        self.taxa = aggregate(&mut self.taxa, &result_data.taxa); 
                        self.taxonomy_warnings.append(&mut result_data.taxonomy_warnings.clone())    
                    },
                    AssemblyTools::Diamond(result_data) => {
                        if result_data.taxa.is_empty() { log::warn!("Attempting to aggregate taxa from Diamond in AssemblyModule, but no taxa found (DB: {})", result_data.db) }
                        // For each Vircov result (e.g. from different databases) ...
                        self.taxa = aggregate(&mut self.taxa, &result_data.taxa); 
                        self.taxonomy_warnings.append(&mut result_data.taxonomy_warnings.clone())    
                    },
                }
                
            }
        }
        Ok(self.clone())
    }
}



// Compute the least common ancestor from a vector of unique taxids
pub fn get_lca_taxid(taxids_unique: Vec<String>, taxonomy: &GeneralTaxonomy) -> Result<(Vec<String>, String), WorkflowError> {

    let mut failed_taxids = Vec::new();

    match taxids_unique.len() {
        0 => return Err(WorkflowError::BlastLcaRecordExtraction),
        1 => {
            // If there is just one taxid, simply return it
            let taxid_str = taxids_unique[0].as_str();
            // Check that the taxid exists in the taxonomy, this is a safety check for later taxonomy operations 
            match taxonomy.rank(taxid_str) {
                Err(_) => {
                    failed_taxids.push(taxid_str.to_owned());
                    return Ok((failed_taxids, taxid_str.to_owned()))
                },
                Ok(_) => return Ok((failed_taxids, taxid_str.to_owned()))
            }
        },
        _ => {


            // LCA by getting lineages in ordered form
            let mut lineages: Vec<Vec<&str>> = vec![];
            for taxid in taxids_unique.iter() {
                let taxid_str = taxid.as_str();
                match taxonomy.lineage(taxid_str) {
                    // If we can't get the lineage for this taxid, 
                    // add it to the list of failed taxids and continue
                    Err(_) => failed_taxids.push(taxid_str.to_owned()),
                    Ok(lineage) => lineages.push(lineage)
                }
            };
            // Transform lineages into sets
            let lineage_sets: Vec<HashSet<&str>> = lineages.clone().into_iter().map(|x| HashSet::from_iter(x)).collect();
            // Take the set intersection of lineages
            let common_lineages: HashSet<&str> = if let Some((first, rest)) = lineage_sets.split_first() {
                rest.iter().fold(first.clone(), |acc, i| {
                    acc.intersection(i).copied().collect()
                })
            } else {
                HashSet::new()
            };
            // For each taxid in order from the first lineage  check if its contained in the set intersection 
            // of lineages and return it, otherwise assign root (1)
            let mut lca_taxid = "1";
            if let Some(first) = lineages.first() {
                for taxid in first {
                    if common_lineages.contains(taxid) {
                        lca_taxid = taxid;
                        break;
                    }
                }
            }
            Ok((failed_taxids, lca_taxid.to_string()))
        }
    }
}
