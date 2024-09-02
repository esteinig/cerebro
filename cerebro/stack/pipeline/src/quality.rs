use std::path::PathBuf;
use serde::{Deserialize, Serialize, Serializer};

use crate::{record::VircovRecord, module::QualityControlModule, error::WorkflowError};

// NANOQ 
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Nanoq {
    reads: u64,
    bases: u64,
    n50: u64,
    longest: u64,
    shortest: u64,
    mean_length: u64,
    median_length: u64,
    mean_quality: Option<f64>,
    median_quality: Option<f64>,
    filtered: u64
}

// FASTP 

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastpSummary {
    #[serde(rename = "fastp_version")] 
    pub version: String,
    #[serde(rename = "before_filtering")] 
    pub before: FastpReadSummary,
    #[serde(rename = "after_filtering")] 
    pub after: FastpReadSummary
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastpReadSummary {
    #[serde(rename = "total_reads")] 
    pub reads: u64,
    #[serde(rename = "total_bases")] 
    pub bases: u64,
    #[serde(rename = "q20_rate")] 
    pub q20: f64,
    #[serde(rename = "q30_rate")] 
    pub q30: f64,
    #[serde(rename = "read1_mean_length")] 
    pub mean_length_r1: u64,
    #[serde(rename = "read2_mean_length")] 
    pub mean_length_r2: u64,
    #[serde(rename = "gc_content")] 
    pub gc: f64
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastpFilter {
    #[serde(rename = "passed_filter_reads")] 
    pub pass_filter: u64,
    #[serde(rename = "low_quality_reads")] 
    pub low_quality: u64,
    #[serde(rename = "low_complexity_reads")] 
    pub low_complexity:  Option<u64>,
    #[serde(rename = "too_many_N_reads")] 
    pub min_missing: u64,
    #[serde(rename = "too_short_reads")] 
    pub min_length: u64,
    #[serde(rename = "too_long_reads")] 
    pub max_length: u64
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastpDuplication {
    pub rate: f64
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastpAdapters {
    #[serde(rename = "adapter_trimmed_reads")] 
    pub total_reads: u64,
    #[serde(rename = "adapter_trimmed_bases")] 
    pub total_bases: u64,
    #[serde(rename = "read1_adapter_sequence")] 
    pub sequence_r1: String,
    #[serde(rename = "read2_adapter_sequence")] 
    pub sequence_r2: String
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Fastp {
    pub summary: FastpSummary,
    #[serde(rename = "filtering_result")]
    pub filter: FastpFilter,
    // May not be run with duplication option
    pub duplication: Option<FastpDuplication>,
    #[serde(rename = "adapter_cutting")]
    pub adapter: Option<FastpAdapters>
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ercc {
    pub constructs: usize,
    pub constructs_aligned: usize,
    pub reads: u64,
    pub alignments: u64,
    pub depletion: Scrubby,
    pub records: Vec<VircovRecord>,
} 
impl Ercc {
    pub fn from(records: &Vec<VircovRecord>, depletion: &Scrubby) -> Self {

        let constructs = records.len().clone();  // requires zero-output option in Vircov
        
        let mut constructs_aligned = 0;
        let mut reads = 0;
        let mut alignments = 0;
        for rec in records {
            if rec.reads > 0 {
                constructs_aligned += 1
            }
            reads += rec.reads;
            alignments += rec.alignments;
        }
        let reads = reads*2; // always paired
        Self {
            constructs, constructs_aligned, reads, alignments, depletion: depletion.clone(), records: records.clone()
        }
    }
}

// Total phage for now - split into DNA/RNA phage controls
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Phage {
    pub id: String,
    pub reads: u64,
    pub alignments: u64,
    pub coverage: f64,
    pub record: VircovRecord
} 
impl Phage {
    pub fn from(record: &VircovRecord) -> Self {
        
        Self {
            id: record.id.clone(),
            reads: record.reads*2,
            alignments: record.alignments,
            coverage: record.coverage,
            record: record.clone()
        }
    }
}



/* 
=======
Scrubby
=======
*/

// Struct to hold the read depletion for
// each reference/database
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Scrubby {
    pub version: String,
    pub schema_version: String,
    pub summary: ScrubbySummary,
    pub settings: ScrubbySettings,
    pub pipeline: Vec<ScrubbyReferenceSummary>,
}

// Struct to hold a compiled run summary
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScrubbySummary {
    pub total: u64,
    pub depleted: u64,
    pub extracted: u64,
}

// Struct to hold the provided settings
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScrubbySettings {
    #[serde(alias = "kraken_taxa")]
    pub taxa: Vec<String>,
    #[serde(alias = "kraken_taxa_direct")]
    pub taxa_direct: Vec<String>,
    pub min_len: u64,
    pub min_cov: f64,
    pub min_mapq: u8,
    pub extract: bool,
}

// A summary struct to hold counts
// for the read depletion/extraction
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileSummary {
    pub total: u64,
    pub depleted: u64,
    pub extracted: u64,
    pub input_file: PathBuf,
    pub output_file: PathBuf,
}

// Struct to hold the read depletion for
// each reference/database
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScrubbyReferenceSummary {
    pub index: usize,
    pub tool: Option<String>,
    pub name: String,
    pub path: PathBuf,
    pub total: u64,
    pub depleted: u64,
    pub extracted: u64,
    pub files: Vec<FileSummary>,
}

fn round_two<S>(x: &Option<f64>, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{   
    match x {
        Some(value) =>  s.serialize_str(&format!("{:.2}", value)),
        None => s.serialize_str("")
    }
}

fn round_nine<S>(x: &Option<f64>, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{   
    match x {
        Some(value) =>  s.serialize_str(&format!("{:.9}", value)),
        None => s.serialize_str("")
    }
}

pub struct ModelConfig {
    sample_type: Option<String>,
    sample_group: Option<String>,
    ercc_input_mass: Option<f64>,
    library_tag: Option<String>,
    run_id: Option<String>,
    run_date: Option<String>,
    workflow_name: Option<String>,
    workflow_id: Option<String>,
    workflow_date: Option<String>,
    dna_phage_id: String,
    rna_phage_id: String,
    seq_phage_id: String,
}

impl ModelConfig {
    pub fn new(
        sample_type: Option<String>,
        sample_group: Option<String>,
        ercc_input_mass: Option<f64>,
        library_tag: Option<String>,
        run_id: Option<String>,
        run_date: Option<String>,
        workflow_name: Option<String>,
        workflow_id: Option<String>,
        workflow_date: Option<String>,
        dna_phage_id: Option<String>,
        rna_phage_id: Option<String>,
        seq_phage_id: Option<String>,
    ) -> Self {
        
        Self {
            sample_type,
            sample_group,
            ercc_input_mass,
            library_tag,
            run_id,
            run_date,
            workflow_name,
            workflow_id,
            workflow_date,
            dna_phage_id: match dna_phage_id { Some(id) => id, None => String::from("T4-Monash") },
            rna_phage_id: match rna_phage_id { Some(id) => id, None => String::from("T2-Monash") },
            seq_phage_id: match seq_phage_id { Some(id) => id, None => String::from("PhiX") },
        }
    }
}
impl Default for ModelConfig {
    fn default() -> Self { 
        Self {
            sample_type: None,
            sample_group: None,
            ercc_input_mass: None,
            library_tag: None,
            run_id: None,
            run_date: None,
            workflow_name: None,
            workflow_id: None,
            workflow_date: None,
            dna_phage_id: String::from("T4-Monash"),
            rna_phage_id: String::from("T2-Monash"),
            seq_phage_id: String::from("PhiX"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityControlSummary {
    pub id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub model_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run_date: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub library_tag: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sample_group: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sample_type: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub workflow_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub workflow_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub workflow_date: Option<String>,

    pub total_reads: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_bases: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub total_biomass: Option<f64>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub deduplicated_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub deduplicated_percent: Option<f64>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub ercc_constructs: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ercc_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub ercc_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub ercc_input_mass: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_nine")]
    pub ercc_mass_per_read: Option<f64>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub adapter_trimmed_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub adapter_trimmed_percent: Option<f64>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub low_complexity_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub low_complexity_percent: Option<f64>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub mean_length_r1: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mean_length_r2: Option<u64>,

    pub qc_reads: u64,
    #[serde(serialize_with = "round_two")]
    pub qc_percent: Option<f64>,
    pub qc_bases: Option<u64>,

    // Fastp options
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub qc_bases_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qc_missing_bases_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub qc_missing_bases_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qc_min_length_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub qc_min_length_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qc_low_quality_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub qc_low_quality_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub q20_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub q30_percent: Option<f64>,

    pub host_reads: Option<u64>,
    #[serde(serialize_with = "round_two")]
    pub host_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub host_biomass: Option<f64>,

    pub other_reads: Option<u64>,
    #[serde(serialize_with = "round_two")]
    pub other_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub other_biomass: Option<f64>,

    // ONT options
    #[serde(skip_serializing_if = "Option::is_none")]
    pub longest_read: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub shortest_read: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub median_read_length: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mean_read_length: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_length_n50: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub median_read_q: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mean_read_q: Option<u64>,

    // DNA phage control
    pub dna_phage_id: Option<String>,
    pub dna_phage_reads: Option<u64>,
    #[serde(serialize_with = "round_two")]
    pub dna_phage_percent: Option<f64>,
    #[serde(serialize_with = "round_two")]
    pub dna_phage_coverage_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub dna_phage_biomass: Option<f64>,
    pub rna_phage_id: Option<String>,
    pub rna_phage_reads: Option<u64>,
    #[serde(serialize_with = "round_two")]
    pub rna_phage_percent: Option<f64>,
    #[serde(serialize_with = "round_two")]
    pub rna_phage_coverage_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub rna_phage_biomass: Option<f64>,
    // Sequencing phage control
    pub seq_phage_id: Option<String>,
    pub seq_phage_reads: Option<u64>,
    #[serde(serialize_with = "round_two")]
    pub seq_phage_percent: Option<f64>,
    #[serde(serialize_with = "round_two")]
    pub seq_phage_coverage_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(serialize_with = "round_two")]
    pub seq_phage_biomass: Option<f64>,

    // Output reads after quality control
    pub output_reads: u64,
    #[serde(serialize_with = "round_two")]
    pub output_percent: Option<f64>,

    // Threshold important quality control
    #[serde(skip_serializing_if = "Option::is_none")]
    pub control_status_dna_extraction: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub control_status_rna_extraction: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub control_status_library: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub control_status_sequencing: Option<bool>,
    
}
impl QualityControlSummary {
    pub fn from(
        quality_control_module: &QualityControlModule, 
        model_id: Option<String>,
        ercc_input_mass: Option<f64>,
        model_config: Option<ModelConfig>,
    ) -> Result<Self, WorkflowError> {

        let model_cfg = match model_config {
            Some(model_config) => model_config,
            None => ModelConfig::default()
        };

        let ercc_input_mass = match ercc_input_mass {
            Some(manual_ercc_input_mass) => {
                log::debug!("ERCC input mass specified manually, overwriting any sample sheet specifications");
                Some(manual_ercc_input_mass)
            },
            None => {
                log::debug!("ERCC input mass not specified manually, using sample sheet specifications if available");
                model_cfg.ercc_input_mass
            }
        };

        // Illumina derivation of the quality contrrol module for the database model and output tables
        let summary = match (&quality_control_module.fastp, &quality_control_module.nanoq)  {
            (Some(fastp), None) => {
                Self::from_illumina(fastp, quality_control_module, model_id, ercc_input_mass, model_cfg)?
            },
            (None, Some(nanoq)) => {
                Self::from_ont(nanoq, quality_control_module, model_id, model_cfg)?
            },
            _ => return Err(WorkflowError::QualityControlNotConducted)
        };

       

        
        Ok(summary)
    }
    fn from_ont(
        nanoq: &Nanoq,
        quality_control_module: &QualityControlModule, 
        model_id: Option<String>,
        model_cfg: ModelConfig,
    ) -> Result<Self, WorkflowError> {

        let (total_reads, nanoq_scan) = match &quality_control_module.nanoq_scan {
            Some(nanoq_scan) => (nanoq_scan.reads, nanoq_scan.clone()),
            None => return Err(WorkflowError::NanoqScan)
        };

        let (host_reads, other_reads, host_percent, other_percent) = match &quality_control_module.host_background {
            Some(host_depletion) => {
                
                let host_reads = host_depletion.summary.depleted;
                let other_reads = host_depletion.summary.total - host_depletion.summary.depleted;

                let host_percent = (host_reads as f64 / total_reads as f64)*100.;
                let other_percent = (other_reads as f64 / total_reads as f64)*100.;


                (Some(host_reads), Some(other_reads), Some(host_percent), Some(other_percent))
            },
            None => (None, None, None, None)
        };


        let (dna_phage_reads, dna_phage_percent, dna_phage_biomass, dna_phage_coverage_percent) = get_phage_data(
            quality_control_module, &model_cfg.dna_phage_id, total_reads, &None
        );
        let (rna_phage_reads, rna_phage_percent, rna_phage_biomass, rna_phage_coverage_percent) = get_phage_data(
            quality_control_module, &model_cfg.rna_phage_id, total_reads, &None
        );
        let (seq_phage_reads, seq_phage_percent, seq_phage_biomass, seq_phage_coverage_percent) = get_phage_data(
            quality_control_module, &model_cfg.seq_phage_id, total_reads, &None
        );

        let qc_reads = nanoq.reads;
        let (output_reads, output_percent) = compute_output_reads(
            total_reads, qc_reads, other_reads, dna_phage_reads, rna_phage_reads, seq_phage_reads
        );

        Ok(Self {
            id: quality_control_module.id.clone(),
            model_id,
            library_tag: model_cfg.library_tag,
            sample_group: model_cfg.sample_group,
            sample_type: model_cfg.sample_type,
            run_id: model_cfg.run_id,
            run_date: model_cfg.run_date,
            workflow_name: model_cfg.workflow_name,
            workflow_id: model_cfg.workflow_id,
            workflow_date: model_cfg.workflow_date,
            total_reads,
            total_bases: Some(nanoq_scan.bases),
            total_biomass: None,
            deduplicated_reads: None,
            deduplicated_percent: None,
            ercc_constructs: None,
            ercc_reads: None,
            ercc_percent: None,
            ercc_input_mass: None,
            ercc_mass_per_read: None,
            adapter_trimmed_reads: None,
            adapter_trimmed_percent: None,
            low_complexity_reads: None,
            low_complexity_percent: None,
            mean_length_r1: None,
            mean_length_r2: None,
            qc_reads,
            qc_percent: Some((nanoq.reads as f64 / total_reads as f64)*100.),
            qc_bases: Some(nanoq.bases),
            qc_bases_percent: Some((nanoq.bases as f64 / nanoq_scan.bases as f64)*100.),
            qc_missing_bases_reads: None,
            qc_missing_bases_percent: None,
            qc_min_length_reads: None,
            qc_min_length_percent: None,
            qc_low_quality_reads: None,
            qc_low_quality_percent: None,
            q20_percent: None,
            q30_percent: None,
            host_reads,
            host_percent,
            host_biomass: None,
            other_reads,
            other_percent,
            other_biomass: None,
            longest_read: None,
            shortest_read: None,
            median_read_length: None,
            mean_read_length: None,
            read_length_n50: None,
            median_read_q: None,
            mean_read_q: None,
            dna_phage_id: Some(model_cfg.dna_phage_id.into()),
            dna_phage_reads,
            dna_phage_percent,
            dna_phage_coverage_percent,
            dna_phage_biomass,
            rna_phage_id: Some(model_cfg.rna_phage_id.into()),
            rna_phage_reads,
            rna_phage_percent,
            rna_phage_coverage_percent,
            rna_phage_biomass,
            seq_phage_id: Some(model_cfg.seq_phage_id.into()),
            seq_phage_reads,
            seq_phage_percent,
            seq_phage_coverage_percent,
            seq_phage_biomass,
            output_reads,
            output_percent: Some(output_percent),
            control_status_dna_extraction: None,
            control_status_rna_extraction: None,
            control_status_library: None,
            control_status_sequencing: None,
        })
    }
    fn from_illumina(
        fastp: &Fastp,
        quality_control_module: &QualityControlModule, 
        model_id: Option<String>,
        ercc_input_mass: Option<f64>,
        model_cfg: ModelConfig,
    ) -> Result<Self, WorkflowError> {

        let (mut total_reads, deduplicated_reads) = match &quality_control_module.fastp_scan {
            Some(fastp_scan) => {
                // When deduplication is active, the total reads are scanned
                // before any other processing steps. In this case the actual 
                // read quality control run of Fastp (after optional ERCC) is 
                // then the number of deduplicated reads (+ ERCC reads)
                (fastp_scan.summary.before.reads, Some(fastp.summary.before.reads))
            },
            None => (fastp.summary.before.reads, None)
        };

        let (ercc_constructs, ercc_reads, ercc_percent, ercc_mass_per_read, deduplicated_reads) = match &quality_control_module.ercc {
            Some(ercc_module) => {

                let ercc_constructs = ercc_module.constructs_aligned as u64;
                let ercc_reads = ercc_module.reads;

                let deduplicated_reads = match deduplicated_reads {
                    // If deduplication was conducted, and ERCCs removed, 
                    // read quality control is conducted after the ERCC removal, so
                    // we need to add those reads to obtain the total number of
                    // reads after deduplication
                    Some(dedup_reads) => Some(dedup_reads + ercc_reads),
                    None => {
                        // If no deduplication was conducted, and ERCCs were removed,
                        // read quality control is conducted after the ERCC removal, so
                        // we need to add those reads to obtain the total number of reads
                        total_reads += ercc_reads;
                        None
                    }
                };

                let ercc_percent = (ercc_reads as f64 / total_reads as f64)*100.;

                let ercc_mass_per_read = match ercc_input_mass {
                    Some(input_mass) => Some(input_mass as f64 / ercc_reads as f64),
                    None => None
                };           
                
                (Some(ercc_constructs), Some(ercc_reads), Some(ercc_percent), ercc_mass_per_read, deduplicated_reads)  

            }, 
            None => (None, None, None, None, None)
        };

        let (deduplicated_reads, deduplicated_percent) = match deduplicated_reads {
            Some(dedup_reads) => {
                // We only want the deduplicated read count
                let dedup_reads = total_reads - dedup_reads;

                (Some(dedup_reads), Some((dedup_reads as f64 / total_reads as f64)*100.))
            },
            None => (None, None)
        };

        let (host_reads, other_reads, host_percent, other_percent, host_biomass, other_biomass) = match &quality_control_module.host_background {
            Some(host_depletion) => {
                
                let host_reads = host_depletion.summary.depleted;
                let other_reads = host_depletion.summary.total - host_depletion.summary.depleted;

                let host_percent = (host_reads as f64 / total_reads as f64)*100.;
                let other_percent = (other_reads as f64 / total_reads as f64)*100.;

                let (host_biomass, other_biomass) = match ercc_mass_per_read {
                    Some(ercc_read_mass) => {
                        (Some(ercc_read_mass*host_reads as f64), Some(ercc_read_mass*other_reads as f64))
                    },
                    None => (None, None)
                };

                (Some(host_reads), Some(other_reads), Some(host_percent), Some(other_percent), host_biomass, other_biomass)
            },
            None => (None, None, None, None, None, None)
        };

        let (dna_phage_reads, dna_phage_percent, dna_phage_biomass, dna_phage_coverage_percent) = get_phage_data(
            quality_control_module, &model_cfg.dna_phage_id, total_reads, &ercc_mass_per_read
        );
        let (rna_phage_reads, rna_phage_percent, rna_phage_biomass, rna_phage_coverage_percent) = get_phage_data(
            quality_control_module, &model_cfg.rna_phage_id, total_reads, &ercc_mass_per_read
        );
        let (seq_phage_reads, seq_phage_percent, seq_phage_biomass, seq_phage_coverage_percent) = get_phage_data(
            quality_control_module, &model_cfg.seq_phage_id, total_reads, &ercc_mass_per_read
        );

        let (low_complexity_reads, low_complexity_percent) = match fastp.filter.low_complexity {
            Some(low_complexity_reads) => (Some(low_complexity_reads),  Some((low_complexity_reads as f64 / total_reads as f64)*100.)),
            None => (None, None)
        };
        
        let (adapter_trimmed_reads, adapter_trimmed_percent) = match &fastp.adapter {
            Some(adapter_data) => (Some(adapter_data.total_reads), Some((adapter_data.total_reads as f64 / total_reads as f64)*100.)),
            None => (None, None)
        };

        let qc_reads = fastp.filter.min_missing+fastp.filter.min_length+fastp.filter.low_quality+fastp.filter.max_length; 
        let qc_percent = (qc_reads as f64 / total_reads as f64)*100.;

        let qc_missing_bases_reads = fastp.filter.min_missing;
        let qc_missing_bases_percent = (qc_missing_bases_reads as f64 / total_reads as f64)*100.;

        let qc_min_length_reads = fastp.filter.min_length;
        let qc_min_length_percent = (qc_min_length_reads as f64 / total_reads as f64)*100.;

        let qc_low_quality_reads = fastp.filter.low_quality;
        let qc_low_quality_percent = (qc_low_quality_reads as f64 / total_reads as f64)*100.;

        let q20_percent = fastp.summary.after.q20*100.;
        let q30_percent = fastp.summary.after.q30*100.;
        
        let mean_length_r1 = fastp.summary.after.mean_length_r1;
        let mean_length_r2 = fastp.summary.after.mean_length_r2;

        let total_biomass = match (ercc_mass_per_read, ercc_input_mass) {
            (Some(ercc_read_mass), Some(input_mass)) => {
                match deduplicated_reads {
                    Some(deduplicated_reads) => {
                        // If ERCC was deduplicated, then we need to use the total reads - deduplicated reads
                        // for total biomass calculation:
                        Some((ercc_read_mass*(total_reads-deduplicated_reads) as f64)-input_mass)
                    },
                    None => Some((ercc_read_mass*total_reads as f64)-input_mass)
                }
            },
            _ => None
        };

        let (output_reads, output_percent) = compute_output_reads(
            total_reads, qc_reads, other_reads, dna_phage_reads, rna_phage_reads, seq_phage_reads
        );

        Ok(Self {
            id: quality_control_module.id.clone(),
            model_id,
            library_tag: model_cfg.library_tag,
            sample_group: model_cfg.sample_group,
            sample_type: model_cfg.sample_type,
            run_id: model_cfg.run_id,
            run_date: model_cfg.run_date,
            workflow_name: model_cfg.workflow_name,
            workflow_id: model_cfg.workflow_id,
            workflow_date: model_cfg.workflow_date,
            total_reads,
            total_biomass,
            deduplicated_reads,
            deduplicated_percent,
            ercc_constructs,
            ercc_reads,
            ercc_percent,
            ercc_input_mass,
            ercc_mass_per_read,
            adapter_trimmed_reads,
            adapter_trimmed_percent,
            low_complexity_reads,
            low_complexity_percent,
            mean_length_r1: Some(mean_length_r1),
            mean_length_r2: Some(mean_length_r2),
            qc_reads,
            qc_percent: Some(qc_percent),
            qc_bases: None,
            qc_bases_percent: None,
            qc_missing_bases_reads: Some(qc_missing_bases_reads),
            qc_missing_bases_percent: Some(qc_missing_bases_percent),
            qc_min_length_reads: Some(qc_min_length_reads),
            qc_min_length_percent: Some(qc_min_length_percent),
            qc_low_quality_reads: Some(qc_low_quality_reads),
            qc_low_quality_percent: Some(qc_low_quality_percent),
            q20_percent: Some(q20_percent),
            q30_percent: Some(q30_percent),
            host_reads,
            host_percent,
            host_biomass,
            other_reads,
            other_percent,
            other_biomass,
            total_bases: None,
            longest_read: None,
            shortest_read: None,
            median_read_length: None,
            mean_read_length: None,
            read_length_n50: None,
            median_read_q: None,
            mean_read_q: None,
            dna_phage_id: Some(model_cfg.dna_phage_id.into()),
            dna_phage_reads,
            dna_phage_percent,
            dna_phage_coverage_percent,
            dna_phage_biomass,
            rna_phage_id: Some(model_cfg.rna_phage_id.into()),
            rna_phage_reads,
            rna_phage_percent,
            rna_phage_coverage_percent,
            rna_phage_biomass,
            seq_phage_id: Some(model_cfg.seq_phage_id.into()),
            seq_phage_reads,
            seq_phage_percent,
            seq_phage_coverage_percent,
            seq_phage_biomass,
            output_reads,
            output_percent: Some(output_percent),
            control_status_dna_extraction: None,
            control_status_rna_extraction: None,
            control_status_library: None,
            control_status_sequencing: None,
        })

    }
}


fn get_phage_data(quality_control_module: &QualityControlModule, phage_id: &str, total_reads: u64, ercc_mass_per_read: &Option<f64>) -> (Option<u64>, Option<f64>, Option<f64>, Option<f64>) {

    match &quality_control_module.phage {
        Some(phage_data) => {

            // If phage identifier cannot be found returns none 
            let phage = match find_phage_by_ref(&phage_data, &phage_id) {
                Some(phage) => phage, None => return (None, None, None, None)
            };

            let phage_reads = phage.reads;
            let phage_percent = (phage_reads as f64 / total_reads as f64)*100.;
            let phage_coverage = phage.coverage*100.;

            let phage_biomass = match ercc_mass_per_read {
                Some(ercc_read_mass) => Some(ercc_read_mass*phage_reads as f64),
                None => None,
            };
            (Some(phage_reads), Some(phage_percent), phage_biomass, Some(phage_coverage))
        },
        None => (None, None, None, None)
    }

}

fn find_phage_by_ref(vec: &Vec<Phage>, target_id: &str) -> Option<Phage> {
    for phage in vec {
        if &phage.record.reference == target_id {
            return Some(phage.to_owned());
        }
    }
    None
}

fn compute_output_reads(total_reads: u64, qc_reads: u64, other_reads: Option<u64>, dna_phage_reads: Option<u64>, rna_phage_reads: Option<u64>, seq_phage_reads: Option<u64>) -> (u64, f64) {

    match other_reads {
        Some(other_reads) => {
            // Host depletion was active before phage controls - must remove from 'other_reads'
            let out_reads = match dna_phage_reads { Some(r) => other_reads - r, None => other_reads };
            let out_reads = match rna_phage_reads { Some(r) => out_reads - r, None => out_reads };
            let out_reads = match seq_phage_reads { Some(r) => out_reads - r, None => out_reads };
            (out_reads, (out_reads as f64 / total_reads as f64)*100.)
        },
        None => {
            // Host depletion was not active before phage controls - must remove from 'qc_reads'
            let out_reads = match dna_phage_reads { Some(r) => qc_reads - r, None => qc_reads };
            let out_reads = match rna_phage_reads { Some(r) => out_reads - r, None => out_reads };
            let out_reads = match seq_phage_reads { Some(r) => out_reads - r, None => out_reads };
            (out_reads, (out_reads as f64 / total_reads as f64)*100.)
        }
    }
}