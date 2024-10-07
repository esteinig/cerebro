
use fancy_regex::Regex;
use tabled::Tabled;
use anyhow::Result;
use std::hash::Hash;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use taxonomy::{Taxonomy, GeneralTaxonomy, TaxRank};
use std::{path::PathBuf, fs::File, io::BufReader, collections::HashMap};
use crate::taxon::filters::{apply_filters, TaxonFilterConfig};
use crate::terminal::PipelineTaxaArgs;
use crate::{
    error::WorkflowError, record::{BlastLcaRecord, Kraken2UniqRecord, VircovScanRemapRecord}, sample::WorkflowSample, utils::get_colored_string
};


fn _get_annotations(taxon: &Taxon, annotations: &Vec<TaxonAnnotations>) -> Option<String> {

    let annotation_strings: Vec<String> = annotations.iter().map(|data| {
        
        let mut target_detected = false;
        for annotation in data.annotations.clone().into_iter() {
            if let Some(value) = annotation.taxonomy.species_taxid { 
                if taxon.level.species_taxid == Some(value) {
                    target_detected = true
                }
            }
            if let Some(value) = annotation.taxonomy.species_name { 
                if taxon.level.species_name == Some(value) {
                    target_detected = true
                }
            };
            if let Some(value) = annotation.taxonomy.genus_taxid { 
                if taxon.level.genus_taxid == Some(value) {
                    target_detected = true
                }
            };
            if let Some(value) = annotation.taxonomy.genus_name { 
                if taxon.level.genus_name == Some(value) {
                    target_detected = true
                }
            };
        };
        if target_detected {
            Some(get_colored_string(&data.tag, &data.color))
        } else {
            None
        }
    }).flatten().collect();

    if annotation_strings.is_empty() {
        None
    } else {
        Some(annotation_strings.join("|"))
    }
    
}


// Struct for highlighting taxa in the annotation column
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonAnnotations {
    pub tag: String,
    pub version: String,
    pub taxonomy: String,
    pub date: String,
    pub color: String,
    pub description: String,
    pub annotations: Vec<TaxonAnnotation>
}
impl TaxonAnnotations {
    pub fn from(json: &PathBuf) -> Result<Self, WorkflowError> {
        serde_json::from_reader(BufReader::new(File::open(json)?)).map_err(WorkflowError::JsonSerialization)
    }
}

// Struct for highlighting taxa in the annotation column
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonAnnotation {
    pub description: String,
    pub taxonomy: TaxonMatch
}

// Struct for highlighting taxa in the annotation column
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonMatch {
    pub genus_taxid: Option<String>,
    pub genus_name: Option<String>,
    pub species_taxid: Option<String>,
    pub species_name: Option<String>
}



// Struct to hold data and produce summary data
// for the quality control module
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Taxon {
    pub taxid: String,
    pub rank: TaxRank,
    pub name: String,
    pub lineage: Vec<String>,
    pub level: TaxonLevel,
    pub evidence: TaxonEvidence
}
impl Taxon {
    /// Creates a taxon from a taxid
    pub fn from_taxid(taxid: String, taxonomy: &GeneralTaxonomy, ncbi_domain: bool) -> Result<Self, WorkflowError> {

        let taxid_str = taxid.as_str();

        let rank = taxonomy.rank(taxid_str)
            .map_err(|err|WorkflowError::TaxRankNotAvailable(err, taxid.clone()))?;
        let name = taxonomy.name(taxid_str)
            .map_err(|err|WorkflowError::TaxNameNotAvailable(err, taxid.clone()))?.to_string();
        let lineage = taxonomy.lineage(taxid_str)
            .map_err(|err| WorkflowError::TaxLineageNotAvailable(err, taxid.clone()))?
            .iter().map(|x| x.to_string()).collect();
                                        

        Ok(Self { 
            taxid: taxid.clone(), 
            rank,
            name,
            lineage,
            level: TaxonLevel::new(taxid_str, taxonomy, ncbi_domain)?,
            evidence: TaxonEvidence::new() // to be filled

        })
    }
    /// Get evidence filtered for a specific sample id (from aggregated evidence)
    pub fn get_evidence(&self, id: &str) -> TaxonEvidence {
        TaxonEvidence {
            kmer: self.evidence.kmer.clone().into_iter().filter(|record| record.id == id).collect(),
            alignment: self.evidence.alignment.clone().into_iter().filter(|record| record.id == id).collect(),
            assembly: self.evidence.assembly.clone().into_iter().filter(|record| record.id == id).collect(),
        }
    }

    /// Deduplicate evidence:
    pub fn deduplicate(&mut self) -> Self {

        self.evidence.kmer = self.evidence.kmer.clone().into_iter().dedup().collect();
        self.evidence.assembly = self.evidence.assembly.clone().into_iter().dedup().collect();
        self.evidence.alignment = self.evidence.alignment.clone().into_iter().dedup().collect();

        self.clone()

    }
    /// Update taxon evidence with new evidence object
    pub fn update_evidence(&mut self, evidence: &mut TaxonEvidence) -> Self {
        for ev in &evidence.kmer {
            if !(self.evidence.kmer.contains(&ev)) {
                self.evidence.kmer.push(ev.clone());
            }
        }
        for ev in &evidence.alignment {
            if !(self.evidence.alignment.contains(&ev)) {
                self.evidence.alignment.push(ev.clone());
            }
        }
        for ev in &evidence.assembly {
            if !(self.evidence.assembly.contains(&ev)) {
                self.evidence.assembly.push(ev.clone());
            }
        }
        self.clone()
    }
    /// Update taxon evidence with new evidence object
    pub fn update_evidence_sample_id(&self, sample_id: &str) -> Self {

        let mut cc = self.clone();

        cc.evidence.kmer = self.evidence.kmer.iter().map(|ev| {
            let mut evidence_update = ev.clone();
            evidence_update.id = sample_id.to_string();
            evidence_update
        }).collect::<Vec<Kraken2UniqRecord>>();

        cc.evidence.alignment = self.evidence.alignment.iter().map(|ev| {
            let mut evidence_update = ev.clone();
            evidence_update.id = sample_id.to_string();
            evidence_update
        }).collect::<Vec<VircovScanRemapRecord>>();

        cc.evidence.assembly = self.evidence.assembly.iter().map(|ev| {
            let mut evidence_update = ev.clone();
            evidence_update.id = sample_id.to_string();
            evidence_update
        }).collect::<Vec<BlastLcaRecord>>();

        cc
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TaxonEvidence {
    pub kmer: Vec<Kraken2UniqRecord>,
    pub alignment: Vec<VircovScanRemapRecord>,
    pub assembly: Vec<BlastLcaRecord> // standardised from BLASTN and DIAMOND
}

impl TaxonEvidence {
    /// Creates a new `TaxonEvidence` with empty evidence fields.
    pub fn new() -> Self {
        Self { 
            kmer: Vec::new(), 
            alignment: Vec::new(), 
            assembly: Vec::new() 
        }
    }
    /// Creates a single `TaxonEvidence` instance from a collection of `TaxonEvidence`
    /// instances by flattening their contents into the respecive evidence fields.
    pub fn from(evidence: &Vec<Self>) -> Self {
        Self {
            kmer: evidence.iter().flat_map(|x| x.kmer.clone()).collect(),
            alignment: evidence.iter().flat_map(|x| x.alignment.clone()).collect(),
            assembly: evidence.iter().flat_map(|x| x.assembly.clone()).collect(),
        }
    }
    pub fn test_case() -> Self {
        Self {
            kmer: Vec::new(),
            alignment: Vec::new(),
            assembly: Vec::new()

        }
    }
    // Get sum of reads per million and sum of assembled bases for all evidence 
    pub fn get_evidence_sums(&self) -> (f64, f64, f64, f64) {
        let kmer_sum = self.kmer.iter().map(|x| x.rpm).sum();
        let alignment_sum = self.alignment.iter().map(|x| x.scan_rpm).sum();

        (
            kmer_sum + alignment_sum,
            self.assembly.iter().map(|x| x.bpm).sum(),
            kmer_sum,
            alignment_sum
        )
    }
    // Get sum of reads per million and sum of assembled bases for all evidence 
    pub fn get_evidence_medians(&self) -> (f64, f64, f64, f64) {
        let mut kmer_rpm: Vec<f64> = self.kmer.iter().map(|x| x.rpm).collect();
        let mut alignment_rpm: Vec<f64> = self.alignment.iter().map(|x| x.scan_rpm).collect();
        let assembly_bpm: Vec<f64> = self.assembly.iter().map(|x| x.bpm).collect();
        
        let kmer_median = match kmer_rpm.len() > 0 {
            true => statistical::median(&kmer_rpm),
            false => 0.
        };
        let alignment_median = match alignment_rpm.len() > 0 {
            true => statistical::median(&alignment_rpm),
            false => 0.
        };
        let assembly_median = match assembly_bpm.len() > 0 {
            true => statistical::median(&assembly_bpm),
            false => 0.
        };

        let mut total_rpm = Vec::new();
        total_rpm.append(&mut kmer_rpm);
        total_rpm.append(&mut alignment_rpm);

        let total_rpm_median = match total_rpm.len() > 0 {
            true => statistical::median(&total_rpm),
            false => 0.
        };
        

        (   
            total_rpm_median,
            assembly_median,
            kmer_median,
            alignment_median
        )
    }

    // Get sum of reads per million and sum of assembled bases for all evidence 
    pub fn get_evidence_means(&self) -> (f64, f64, f64, f64) {
        let mut kmer_rpm: Vec<f64> = self.kmer.iter().map(|x| x.rpm).collect();
        let mut alignment_rpm: Vec<f64> = self.alignment.iter().map(|x| x.scan_rpm).collect();
        let assembly_bpm: Vec<f64> = self.assembly.iter().map(|x| x.bpm).collect();
        
        let kmer_mean = match kmer_rpm.len() > 0 {
            true => statistical::mean(&kmer_rpm),
            false => 0.
        };
        let alignment_mean = match alignment_rpm.len() > 0 {
            true => statistical::mean(&alignment_rpm),
            false => 0.
        };
        let assembly_mean = match assembly_bpm.len() > 0 {
            true => statistical::mean(&assembly_bpm),
            false => 0.
        };

        let mut total_rpm = Vec::new();
        total_rpm.append(&mut kmer_rpm);
        total_rpm.append(&mut alignment_rpm);

        let total_rpm_mean = match total_rpm.len() > 0 {
            true => statistical::mean(&total_rpm),
            false => 0.
        };
        

        (   
            total_rpm_mean,
            assembly_mean,
            kmer_mean,
            alignment_mean
        )
    }

    // Get all sample identifiers from the evidence records
    pub fn get_sample_identifiers(&self, unique: bool) -> Vec<String> {
        let mut all_identifiers = Vec::new();

        all_identifiers.append(&mut self.kmer.iter().map(|x| x.id.clone()).collect());
        all_identifiers.append(&mut self.alignment.iter().map(|x| x.id.clone()).collect());
        all_identifiers.append(&mut self.assembly.iter().map(|x| x.id.clone()).collect());

        if unique {
            all_identifiers.into_iter().unique().collect()
        } else {
            all_identifiers
        }

    }
    // Evidence across modules
    pub fn kmer_module(&self) -> bool {
        !(self.kmer.is_empty())  // extend here when new tools in kmer module become available
    }
    pub fn alignment_module(&self) -> bool {
        !(self.alignment.is_empty())    // extend here when new tools in alignment module become available
    }
    pub fn assembly_module(&self) -> bool {
        !(self.assembly.is_empty())     // extend here when new tools in assembly module become available
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TaxonLevel {
    pub domain_taxid: Option<String>,
    pub domain_name: Option<String>,
    pub genus_taxid: Option<String>,
    pub genus_name: Option<String>,
    pub species_taxid: Option<String>,
    pub species_name: Option<String>
}
impl TaxonLevel {
    /// Creates a new instance of `TaxonLevel` given a taxonomic identifier with its corresponding
    /// taxonomy and a flag, whether the domain rank conforms to NCBI definition of domain ("Superkingdom")
    /// or whether the rank conforms to the traditional definition ("Domain"). A taxon level struct is a
    /// convenience abstraction of the taxon lineeage at domain, genus and species levels, so that later
    /// a taxonomy does not have to be consulted to aggegate or filter taxa at these levels.
    /// 
    /// 
    /// # Example
    /// 
    /// ```
    /// let taxonomy = taxonomy::ncbi::load(taxonomy)?;
    /// let taxon_level = TaxonLevel::new("9606", &taxonomy, true)
    /// ```
    pub fn new(taxid: &str, taxonomy: &GeneralTaxonomy, ncbi_domain: bool) -> Result<Self, WorkflowError> {

        // Species abstraction
        let species_search = taxonomy.parent_at_rank(
            taxid, TaxRank::Species
        ).map_err(|err|WorkflowError::TaxRankNotFound(err, TaxRank::Species.to_string(), taxid.to_owned()))?;

        let (species_taxid, species_name) = get_taxid_name_from_search(species_search, taxonomy)?;
        
        // Genus abstraction
        let genus_search = taxonomy.parent_at_rank(
            taxid, TaxRank::Genus
        ).map_err(|err|WorkflowError::TaxRankNotFound(err, TaxRank::Genus.to_string(), taxid.to_owned()))?;

        let (genus_taxid, genus_name) = get_taxid_name_from_search(genus_search, taxonomy)?;
        
        // Domain abstraction
        let domain_rank = match ncbi_domain {
            true => { TaxRank::Superkingdom },
            false => { TaxRank::Domain }
        };

        let domain_taxid_search: Option<(&str, f32)> = taxonomy.parent_at_rank(
            taxid, domain_rank
        ).map_err(|err|WorkflowError::TaxRankNotFound(err, domain_rank.to_string(), taxid.to_owned()))?;

        let (domain_taxid, domain_name) = match domain_taxid_search {
            Some((domain_taxid, _)) => {
                let domain_name = taxonomy.name(domain_taxid).map_err(
                    |err|WorkflowError::TaxNameNotAvailable(err, taxid.to_owned())
                )?;
                (Some(domain_taxid.to_owned()), Some(domain_name.to_owned()))

            },
            None => (None, None)
        };


        Ok(Self { 
            domain_taxid,
            domain_name,
            genus_taxid,
            genus_name,
            species_taxid,
            species_name
        })

    }
    pub fn test_case() -> Self {
        Self {  
            domain_taxid: Some("1".into()),
            domain_name: Some("Eukaryota".into()),
            genus_taxid: Some("10".into()),
            genus_name: Some("Homo".into()),
            species_taxid: Some("9606".into()),
            species_name: Some("Home sapiens".into())

        }
    }
}

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

// Struct for a warning if a specific taxon or lineage
// for a taxon could not be found in the taxonomy;
// bubbles up from taxon assignments and LCA operations
// at tool hierarchical levels
#[derive(Tabled, Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct TaxonomyWarning {
    #[tabled(rename = "ID")]
    pub id: String,
    #[tabled(rename = "Database")]
    pub db: String,
    #[tabled(rename = "Tool")]
    pub tool: String,
    #[tabled(rename = "Function")]
    pub function: String,
    #[tabled(rename = "Details")]
    pub details: String
}
impl TaxonomyWarning {
    pub fn new(id: &str, tool: &str, db: &str, function: &str, details: String) -> Self {
        Self { id: id.to_string(), tool: tool.to_string(), db: db.to_string(), function: function.to_string(), details }
    }
}

/* 
=======================
MODULE AND TOOL STRUCTS
=======================
*/

// Utility function to aggregate multiple HashMaps with Taxon
pub fn aggregate(parent_taxa: &mut HashMap<String, Taxon>, taxa: &HashMap<String, Taxon>) -> HashMap<String, Taxon> {
    for (taxid, taxon) in taxa {
        // ... check if each taxon exists
        match parent_taxa.get_mut(taxid) {
            Some(tax) => {
                // If so, update the evidence of the taxon in the modules...
                for record in &taxon.evidence.kmer {
                    tax.evidence.kmer.push(record.clone())
                }
                for record in &taxon.evidence.alignment {
                    tax.evidence.alignment.push(record.clone())
                }
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
    pub fn from_args(args: &PipelineTaxaArgs) -> Self {
        Self {
            min_rpm: args.min_rpm,
            min_rpm_kmer: args.min_rpm_kmer,
            min_rpm_alignment: args.min_rpm_alignment,
            min_rpm_remap: args.min_rpm_remap,
            min_contigs: args.min_contigs,
            min_bases: args.min_bases
        }
    }
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

    let mut wtr = csv::WriterBuilder::new()
        .delimiter(sep as u8)
        .has_headers(header)
        .from_path(output)
        .map_err(|err| WorkflowError::CreateTaxaOutputFile(err))?;

    let filter_config = match filter_config {
        Some(file) => TaxonFilterConfig::from_path(&file)?,
        None => TaxonFilterConfig::default()
    };

    for file in samples {
        let sample = WorkflowSample::read_json(&file).expect(&format!("Failed to parse sample file: {}", file.display()));
        let taxa: Vec<_> = apply_filters(sample.taxa.clone().into_values().collect(), &filter_config);
        let taxa_overview: Vec<_> = taxa.iter().map(|taxon| { TaxonOverview::from(taxon) }).collect();
        let taxa_overview_summary: Vec<_> = taxa_overview.iter().map(|taxon_overview|{
            TaxonSampleOverview::from_taxon_overview(&sample, &taxon_overview, extract).expect("Failed to construct taxon sample overview!") // Fix to error bubbles!
        }).filter(|taxon_sample_overview| taxon_sample_overview.pass(&threshold_config)).collect();

        if taxa_overview_summary.is_empty() {
            log::warn!("No taxa for sample {} - sample is excluded from results table!", sample.id);
        }

        for record in taxa_overview_summary {
            wtr.serialize(record).map_err(|err| WorkflowError::SerializeTaxaTableRow(err))?;
        }
    }

    Ok(())
}


#[derive(Debug, Clone, Serialize, Deserialize)]
// A struct representing a high level overview
// of the aggregated taxon evidence
pub struct TaxonOverview {
    pub taxid: String,
    pub domain: Option<String>,
    pub genus: Option<String>,
    pub name: String,             // used to later map back the tags
    pub rpm: f64,                 // total rpm summed from k-mer and alignment evidence
    pub rpm_kmer: f64,            // total rpm summed from k-mer evidence
    pub rpm_alignment: f64,       // total rpm summed from  alignment evidence
    pub rpm_remap: f64,           // total rpm summed from  alignment evidence
    pub contigs: u64,             // total assembled and identified contig evidence
    pub contigs_bases: u64,
    pub kmer: bool,
    pub alignment: bool,
    pub assembly: bool,
    pub names: Vec<String>        // the evidence record associated sample names as processed in the pipeline (matching `Cerebro.name`)
}
impl TaxonOverview {
    pub fn from(taxon: &Taxon) -> Self {

        let (kmer, alignment, assembly) = (
            !taxon.evidence.kmer.is_empty(),
            !taxon.evidence.alignment.is_empty(),
            !taxon.evidence.assembly.is_empty()
        );

        // Unique record identifiers (sample names)
        let mut names = Vec::new();
        for record in taxon.evidence.kmer.iter() {
            names.push(record.id.to_owned())
        };
        for record in taxon.evidence.alignment.iter() {
            names.push(record.id.to_owned())
        };
        for record in taxon.evidence.assembly.iter() {
            names.push(record.id.to_owned())
        }
        let names = names.into_iter().unique().collect();

        // Summary values
        let rpm = match kmer || alignment {
            true => taxon.evidence.alignment.iter().map(|e| e.scan_rpm).sum::<f64>() +
                taxon.evidence.kmer.iter().map(|e| e.rpm).sum::<f64>(),
            false => 0.0
        };

        let rpm_kmer = match kmer {
            true => taxon.evidence.kmer.iter().map(|e| e.rpm).sum::<f64>(),
            false => 0.0
        };

        let rpm_alignment = match alignment {
            true => taxon.evidence.alignment.iter().map(|e| e.scan_rpm).sum::<f64>(),
            false => 0.0
        };

        let rpm_remap = match alignment {
            true => taxon.evidence.alignment.iter().map(|e| e.remap_rpm).sum::<f64>(),
            false => 0.0
        };

        let (contigs, contigs_bases) = match assembly {
            true => (taxon.evidence.assembly.len() as u64,taxon.evidence.assembly.iter().map(|e| e.length).sum()),
            false => (0, 0)
        };

        Self {
            taxid: taxon.taxid.to_owned(),
            name: taxon.name.to_owned(),
            genus: taxon.level.genus_name.to_owned(),
            domain: taxon.level.domain_name.to_owned(),
            rpm,
            rpm_kmer,
            rpm_alignment,
            rpm_remap,
            contigs,
            contigs_bases,
            kmer,
            alignment,
            assembly,
            names
        }
    }
}

// DUPLICATED

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


// TaxonOverview + TaxaSummary
#[derive(Deserialize, Serialize, Debug)]
pub struct TaxonSampleOverview {
    pub id: String,
    pub sample_id: Option<String>,
    pub sample_tag: Option<String>,
    pub taxid: String,
    pub domain: Option<String>,
    pub genus: Option<String>,
    pub name: String,             // used to later map back the tags
    pub rpm: f64,                 // total rpm summed from k-mer and alignment evidence
    pub rpm_kmer: f64,            // total rpm summed from k-mer evidence
    pub rpm_alignment: f64,       // total rpm summed from  alignment evidence
    pub rpm_remap: f64,           // total rpm summed from  alignment evidence
    pub contigs: u64,             // total assembled and identified contig evidence
    pub contigs_bases: u64,
    pub kmer: bool,
    pub alignment: bool,
    pub assembly: bool
}
impl TaxonSampleOverview {
    pub fn from_taxon_overview(sample: &WorkflowSample, taxon_overview: &TaxonOverview, tags: bool) -> Result<Self, WorkflowError> {
        
        let (sample_id, sample_tag) = match tags {
            true => {
                let (sample_id, sample_tags) = get_sample_regex_matches(&sample.id)?;
                (Some(sample_id), Some(sample_tags.join("-")))
            }
            false => (None, None)
        };

        Ok(Self {
            id: sample.id.clone(),
            sample_id,
            sample_tag,
            taxid: taxon_overview.taxid.clone(),
            domain: taxon_overview.domain.clone(),
            genus: taxon_overview.genus.clone(),
            name: taxon_overview.name.clone(),            
            rpm: taxon_overview.rpm,                 
            rpm_kmer: taxon_overview.rpm_kmer,            
            rpm_alignment: taxon_overview.rpm_alignment,  
            rpm_remap: taxon_overview.rpm_remap,  
            contigs: taxon_overview.contigs,             
            contigs_bases: taxon_overview.contigs_bases,
            kmer: taxon_overview.kmer,
            alignment: taxon_overview.alignment,
            assembly: taxon_overview.assembly
        })
    }
    pub fn pass(&self, threshold_config: &TaxonThresholdConfig) -> bool {
        self.rpm >= threshold_config.min_rpm &&
        self.rpm_kmer >= threshold_config.min_rpm_kmer &&
        self.rpm_alignment >= threshold_config.min_rpm_alignment &&
        self.rpm_remap >= threshold_config.min_rpm_remap &&
        self.contigs >= threshold_config.min_contigs &&
        self.contigs_bases >= threshold_config.min_bases
    }
}
