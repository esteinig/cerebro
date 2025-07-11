use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use csv::WriterBuilder;
use scrubby::report::ScrubbyReport;
use serde::{Deserialize, Serialize, Serializer};
use vircov::vircov::{VircovRecord, VircovSummary};


use crate::error::WorkflowError;
use crate::modules::pathogen::PathogenDetection;
use crate::nextflow::pathogen::PathogenDetectionOutput;
use crate::nextflow::panviral::PanviralOutput;
use crate::nextflow::quality::QualityControlOutput;
use crate::parsers::fastp::FastpReport;
use crate::parsers::nanoq::NanoqReport;
use crate::taxa::taxon::Taxon;
use crate::tools::scan::ScanReport;
use crate::tools::umi::DeduplicationReport;

// Main module wrapper functions

pub fn write_quality_table(json: &Vec<PathBuf>, reads: &PathBuf, controls: &PathBuf, background: &PathBuf) -> Result<(), WorkflowError> {
    
    let mut data = Vec::new();

    for path in json {
        data.push(QualityControl::from_json(&path)?)
    }

    let mut reads_writer = WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(File::create(reads)?);

    for qc in &data {
        reads_writer.serialize(&qc.reads)?;
    }
    
    reads_writer.flush()?;
    
    let mut controls_writer = WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(controls)?;

    for qc in &data {
        if let Some(ercc) = &qc.controls.ercc {
            for record in &ercc.records {
                controls_writer.serialize(record)?;
            }
        }
        if let Some(organism) = &qc.controls.organism {
            for record in &organism.records {
                controls_writer.serialize(record)?;
            }
        }
    }

    controls_writer.flush()?;
    
    let mut background_writer = WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(background)?;

    for qc in &data {
        if let Some(other) = &qc.background.other {
            for record in &other.records {
                background_writer.serialize(record)?;
            }
        }
    }
    
    background_writer.flush()?;

    Ok(())
}

pub fn write_quality_json(qc_data: &Vec<QualityControl>, output: &PathBuf) -> Result<(), WorkflowError> {
    let writer = BufWriter::new(File::create(output)?);
    serde_json::to_writer_pretty(writer, &qc_data)?;
    Ok(())
}

// Main module

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityControl {
    pub id: String,
    pub reads: ReadQualityControl,
    pub background: Background,
    pub controls: Controls
}
impl QualityControl {
    pub fn from_panviral(output: &PanviralOutput) -> Self {
        Self::from_quality(&output.qc)
    }
    pub fn from_pathogen(output: &PathogenDetectionOutput) -> Self {
        Self::from_quality(&output.qc)
    }
    pub fn from_quality(output: &QualityControlOutput) -> Self {

        let controls = Controls::from(
            &output.id,
            output.controls.clone(), 
            ControlsConfig::default()
        );        

        let background = Background::from(
            &output.id,
            output.host_depletion.clone(),
            output.background_alignment.clone(),
            BackgroundDepletionConfig::default()
        );
        
        let qc: Option<Box<dyn ReadReport>> = if let Some(fastp) = output.reads_qc_fastp.clone() {
            Some(Box::new(fastp))
        } else if let Some(nanoq) = output.reads_qc_nanoq.clone() {
            Some(Box::new(nanoq))
        } else {
            None
        };
    
        let reads = ReadQualityControl::from(
            &output.id,
            Some(output.input_scan.clone()),
            qc,
            Some(controls.clone()),
            Some(background.clone()),
            output.deduplication.clone(),
            Some(output.output_scan.clone())
        );

        Self {
            id: output.id.clone(),
            reads,
            background,
            controls
        }
    }
    pub fn to_json(&self, path: &PathBuf) -> Result<(), WorkflowError> {
        let writer = BufWriter::new(File::create(path)?);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }

    pub fn from_json(path: &PathBuf) -> Result<Self, WorkflowError> {
        let reader = BufReader::new(File::open(path)?);
        let quality_control = serde_json::from_reader(reader)?;
        Ok(quality_control)
    }

    pub fn summary_illumina_pe(&self) -> QualityControlSummary {

        let fastp_map = HashMap::from([
            ("q20_percent".to_string(), self.reads.q20_percent.unwrap()),
            ("q30_percent".to_string(), self.reads.q30_percent.unwrap()),
        ]);

        let builder = QualityControlSummaryBuilder::new(
            &self.id,
            self.reads.input_reads,
            self.reads.output_reads,
            self.reads.ercc_constructs,
            self.controls.organism.clone().unwrap().records,
            fastp_map,
            QcConfig::illumina_pe_sr(),
        );

        builder.build()
    }
}


// For data quality control table output from API

pub struct ModelConfig {
    sample_type: Option<String>,
    sample_group: Option<String>,
    ercc_input_mass: Option<f64>,
    library_tag: Option<String>,
    run_id: Option<String>,
    run_date: Option<String>,
    workflow_name: Option<String>,
    workflow_id: Option<String>,
    workflow_date: Option<String>
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
        workflow_date: Option<String>
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
            workflow_date
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
        }
    }
}

// Serializers

trait IntoOptionF64 {
    fn as_option_f64(&self) -> Option<f64>;
}

impl IntoOptionF64 for f64 {
    fn as_option_f64(&self) -> Option<f64> {
        Some(*self)
    }
}

impl IntoOptionF64 for Option<f64> {
    fn as_option_f64(&self) -> Option<f64> {
        *self
    }
}

fn _round_two<T, S>(x: &T, s: S) -> Result<S::Ok, S::Error>
where
    T: IntoOptionF64,
    S: Serializer,
{
    match x.as_option_f64() {
        Some(value) => s.serialize_str(&format!("{:.2}", value)),
        None => s.serialize_str(""),
    }
}

// Submodules

pub trait ReadReport {
    fn input_reads(&self) -> u64;
    fn input_bases(&self) -> u64;
    fn qc_reads(&self) -> u64;
    fn qc_reads_percent(&self, total: u64) -> f64;
    fn qc_bases(&self) -> u64;
    fn qc_bases_percent(&self, total: u64) -> f64;
    fn input_biomass(&self, ercc: Option<ErccControl>) -> Option<f64>;
    fn host_biomass(&self, ercc: Option<ErccControl>, host: Option<HostBackground>) -> Option<f64>;
    fn control_biomass(&self, ercc: Option<ErccControl>, controls: Option<OrganismControl>) -> Option<f64>;
    fn background_biomass(&self, ercc: Option<ErccControl>, background: Option<OtherBackground>) -> Option<f64>;
    fn output_biomass(&self, ercc: Option<ErccControl>) -> Option<f64>;
    fn deduplicated_reads(&self, report: Option<DeduplicationReport>) -> Option<u64>;
    fn deduplicated_percent(&self, report: Option<DeduplicationReport>) -> Option<f64>;
    fn adapter_trimmed_reads(&self) -> Option<u64>;
    fn adapter_trimmed_percent(&self) -> Option<f64>;
    fn low_complexity_reads(&self) -> Option<u64>;
    fn low_complexity_percent(&self) -> Option<f64>;
    fn mean_read_length_r1(&self) -> Option<u64>;
    fn mean_read_length_r2(&self) -> Option<u64>;
    fn min_length_reads(&self) -> Option<u64>;
    fn min_length_percent(&self) -> Option<f64>;
    fn low_quality_reads(&self) -> Option<u64>;
    fn low_quality_percent(&self) -> Option<f64>;
    fn q20_percent(&self) -> Option<f64>;
    fn q30_percent(&self) -> Option<f64>;
}

impl ReadReport for FastpReport {
    fn input_reads(&self) -> u64 {
        self.summary.before.reads
    }
    fn input_bases(&self) -> u64 {
        self.summary.before.bases
    }
    fn qc_reads(&self) -> u64 {
        self.summary.before.reads - self.summary.after.reads
    }
    fn qc_reads_percent(&self, total: u64) -> f64 {
        if total > 0 {
            ((self.summary.before.reads - self.summary.after.reads) as f64 / total as f64)*100.0
        } else {
            0.0
        }
    }
    fn qc_bases(&self) -> u64 {
        self.summary.after.bases
    }
    fn qc_bases_percent(&self, total: u64) -> f64 {
        if total > 0 {
            (self.summary.after.bases as f64 / total as f64)*100.0
        } else {
            0.0
        }
    }
    fn input_biomass(&self, ercc: Option<ErccControl>) -> Option<f64> {
        ercc.as_ref().map(|ercc| self.summary.before.reads as f64 * ercc.mass_per_read)
    }
    fn host_biomass(&self, ercc: Option<ErccControl>, host: Option<HostBackground>) -> Option<f64> {
        if let Some(ercc) = ercc {
            host.as_ref().map(|host| host.alignments as f64 * ercc.mass_per_read)
        } else {
            None
        }  
    }
    fn control_biomass(&self, ercc: Option<ErccControl>, controls: Option<OrganismControl>) -> Option<f64> {
        if let Some(ercc) = ercc {
            controls.as_ref().map(|organism| organism.alignments as f64 * ercc.mass_per_read)
        } else {
            None
        }  
    }
    fn background_biomass(&self, ercc: Option<ErccControl>, background: Option<OtherBackground>) -> Option<f64> {
        if let Some(ercc) = ercc {
            background.as_ref().map(|background| background.alignments as f64 * ercc.mass_per_read)
        } else {
            None
        }  
    }
    fn output_biomass(&self, ercc: Option<ErccControl>) -> Option<f64> {
        ercc.as_ref().map(|ercc| self.summary.after.reads as f64 * ercc.mass_per_read)
    }
    fn deduplicated_reads(&self, report: Option<DeduplicationReport>) -> Option<u64> {
        report.as_ref().map(|dedup| dedup.deduplicated as u64)
    }
    fn deduplicated_percent(&self,  report: Option<DeduplicationReport>) -> Option<f64> {
        report.as_ref().map(|dedup| dedup.deduplicated_percent)
    }
    fn adapter_trimmed_reads(&self) -> Option<u64> {
        self.adapter.as_ref().map(|adapter| adapter.total_reads)
    }
    fn adapter_trimmed_percent(&self) -> Option<f64> {
        self.adapter.as_ref().map(|adapter| {
            if self.summary.before.reads == 0 {
                0.0
            } else {
                (adapter.total_reads as f64 / self.summary.before.reads as f64) * 100.0
            }
        })
    }
    fn low_complexity_reads(&self) -> Option<u64> {
        self.filter.low_complexity
    }
    fn low_complexity_percent(&self) -> Option<f64> {
        self.filter.low_complexity.map(|reads| {
            if self.summary.before.reads == 0 {
                0.0
            } else {
                (reads as f64 / self.summary.before.reads as f64) * 100.0
            }
        })
    }
    fn mean_read_length_r1(&self) -> Option<u64> {
        Some(self.summary.after.mean_length_r1)
    }
    fn mean_read_length_r2(&self) -> Option<u64> {
        Some(self.summary.after.mean_length_r2)
    }
    fn min_length_reads(&self) -> Option<u64> {
        Some(self.filter.min_length)
    }
    fn min_length_percent(&self) -> Option<f64> {
        Some((self.filter.min_length as f64 / self.summary.before.reads as f64) * 100.0)
    }
    fn low_quality_reads(&self) -> Option<u64> {
        Some(self.filter.low_quality)
    }
    fn low_quality_percent(&self) -> Option<f64> {
        Some((self.filter.low_quality as f64 / self.summary.before.reads as f64) * 100.0)
    }
    fn q20_percent(&self) -> Option<f64> {
        Some(self.summary.after.q20*100.0)
    }
    fn q30_percent(&self) -> Option<f64> {
        Some(self.summary.after.q30*100.0)
    }
}


impl ReadReport for NanoqReport {
    fn input_reads(&self) -> u64 {
        0
    }
    fn input_bases(&self) -> u64 {
        0
    }
    fn qc_reads(&self) -> u64 {
        self.reads
    }
    fn qc_reads_percent(&self, total: u64) -> f64 {
        if total > 0 {
            (self.reads as f64 / total as f64)*100.0
        } else {
            0.0
        }
    }
    fn qc_bases(&self) -> u64 {
        self.bases
    }
    fn qc_bases_percent(&self, total: u64) -> f64 {
        if total > 0 {
            (self.bases as f64 / total as f64)*100.0
        } else {
            0.0
        }
    }
    fn input_biomass(&self, _: Option<ErccControl>) -> Option<f64> {
        None
    }
    fn host_biomass(&self, ercc: Option<ErccControl>, host: Option<HostBackground>) -> Option<f64> {
        if let Some(ercc) = ercc {
            host.as_ref().map(|host| host.alignments as f64 * ercc.mass_per_read)
        } else {
            None
        }  
    }
    fn control_biomass(&self, ercc: Option<ErccControl>, controls: Option<OrganismControl>) -> Option<f64> {
        if let Some(ercc) = ercc {
            controls.as_ref().map(|organism| organism.alignments as f64 * ercc.mass_per_read)
        } else {
            None
        }  
    }
    fn background_biomass(&self, ercc: Option<ErccControl>, background: Option<OtherBackground>) -> Option<f64> {
        if let Some(ercc) = ercc {
            background.as_ref().map(|background| background.alignments as f64 * ercc.mass_per_read)
        } else {
            None
        }  
    }
    fn output_biomass(&self, ercc: Option<ErccControl>) -> Option<f64> {
        ercc.as_ref().map(|ercc| self.reads as f64 * ercc.mass_per_read)
    }
    fn deduplicated_reads(&self, report: Option<DeduplicationReport>) -> Option<u64> {
        report.as_ref().map(|dedup| dedup.deduplicated as u64)
    }
    fn deduplicated_percent(&self,  report: Option<DeduplicationReport>) -> Option<f64> {
        report.as_ref().map(|dedup| dedup.deduplicated_percent)
    }
    fn adapter_trimmed_reads(&self) -> Option<u64> {
        None
    }
    fn adapter_trimmed_percent(&self) -> Option<f64> {
        None
    }
    fn low_complexity_reads(&self) -> Option<u64> {
        None
    }
    fn low_complexity_percent(&self) -> Option<f64> {
        None
    }
    fn mean_read_length_r1(&self) -> Option<u64> {
        Some(self.mean_length)
    }
    fn mean_read_length_r2(&self) -> Option<u64> {
        None
    }
    fn min_length_reads(&self) -> Option<u64> {
        None
    }
    fn min_length_percent(&self) -> Option<f64> {
        None
    }
    fn low_quality_reads(&self) -> Option<u64> {
        None
    }
    fn low_quality_percent(&self) -> Option<f64> {
        None
    }
    fn q20_percent(&self) -> Option<f64> {
        None
    }
    fn q30_percent(&self) -> Option<f64> {
        None
    }
}



#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadQualityControl {
    pub id: String,

    // API
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

    pub input_reads: u64,
    pub input_bases: u64,

    // ERCC
    pub ercc_constructs: Option<u64>,
    pub ercc_reads: Option<u64>,
    pub ercc_reads_percent: Option<f64>,

    // Deduplication
    pub deduplicated_reads: Option<u64>,
    pub deduplicated_percent: Option<f64>,

    // Stages
    pub qc_reads: Option<u64>,
    pub qc_reads_percent: Option<f64>,
    pub qc_bases: Option<u64>,
    pub qc_bases_percent: Option<f64>,
    pub host_reads: Option<u64>,
    pub host_reads_percent: Option<f64>,
    pub control_reads: Option<u64>,
    pub control_reads_percent: Option<f64>,
    pub background_reads: Option<u64>,
    pub background_reads_percent: Option<f64>,
    pub output_reads: u64,
    pub output_reads_percent: f64,
    pub output_bases: u64,
    pub output_bases_percent: f64,
    
    // Biomass
    pub input_biomass: Option<f64>,
    pub host_biomass: Option<f64>,
    pub control_biomass: Option<f64>,
    pub background_biomass: Option<f64>,
    pub output_biomass: Option<f64>,

    // Fastp fields
    pub adapter_trimmed_reads: Option<u64>,
    pub adapter_trimmed_percent: Option<f64>,
    pub low_complexity_reads: Option<u64>,
    pub low_complexity_percent: Option<f64>,
    pub mean_read_length_r1: Option<u64>,
    pub mean_read_length_r2: Option<u64>,
    pub min_length_reads: Option<u64>,
    pub min_length_percent: Option<f64>,
    pub low_quality_reads: Option<u64>,
    pub low_quality_percent: Option<f64>,
    pub q20_percent: Option<f64>,
    pub q30_percent: Option<f64>,
}
impl ReadQualityControl {

    pub fn with_model(
        &mut self,
        model_id: &str, 
        model_config: ModelConfig
    ) -> Result<Self, WorkflowError> {

        self.model_id = Some(model_id.to_string());
        self.run_id = model_config.run_id;
        self.run_date = model_config.run_date;
        self.library_tag = model_config.library_tag;
        self.sample_group = model_config.sample_group;
        self.sample_type = model_config.sample_type;
        self.workflow_name = model_config.workflow_name;
        self.workflow_id = model_config.workflow_id;
        self.workflow_date = model_config.workflow_date;

        Ok(self.clone())
    }
    pub fn from(
        id: &str, 
        input_scan: Option<ScanReport>,
        qc:  Option<Box<dyn ReadReport>>,  // Dynamic dispatch of reports implementing ReadReport trait
        controls: Option<Controls>,
        background: Option<Background>,
        deduplication: Option<DeduplicationReport>,
        output_scan: Option<ScanReport>
    ) -> Self {

        let (input_reads, input_bases) = match input_scan {
            Some(scan_report) => (
                scan_report.reads as u64, 
                scan_report.bases as u64
            ), 
            None => {
                match qc {
                    Some(ref qc_report) => (
                        qc_report.input_reads(), 
                        qc_report.input_bases()
                    ),
                    None => (0, 0)
                }
            }
        };


        let (output_reads, output_reads_percent, output_bases, output_bases_percent) = match output_scan {
            Some(scan_report) => {
                let output_reads_percent = if input_reads == 0 {
                    0.0 // otherwise would be NaN!
                } else {
                    (scan_report.reads as f64 / input_reads as f64)*100.0
                };

                let output_bases_percent = if input_bases == 0 {
                    0.0  // otherwise would be NaN!
                } else {
                    (scan_report.bases as f64 / input_bases as f64)*100.0
                };

                (scan_report.reads as u64, output_reads_percent, scan_report.bases as u64, output_bases_percent)
            }, 
            None => {
                match qc {
                    Some(ref qc_report) => (
                        qc_report.qc_reads(), 
                        qc_report.qc_reads_percent(input_reads), 
                        qc_report.qc_bases(), 
                        qc_report.qc_bases_percent(input_bases)
                    ),
                    None => (0, 0.0, 0, 0.0)
                }
            }
        };

        let ercc = controls.as_ref().map(|c| {
            c.ercc.clone()
        }).flatten();

        let qc_reads = qc.as_ref().map(|r| {
            r.qc_reads()
        });
        let qc_reads_percent = qc.as_ref().map(|r| {
            r.qc_reads_percent(input_reads)
        });
        let qc_bases = qc.as_ref().map(|r| {
            r.qc_bases()
        });
        let qc_bases_percent = qc.as_ref().map(|r| {
            r.qc_bases_percent(input_bases)
        });

        let input_biomass = ercc.as_ref().map(|e| {
            e.biomass(Some(input_reads))
        }).flatten();
        let host_biomass = ercc.as_ref().map(|e| {
            e.biomass(
                background.as_ref().map(|b| b.host_reads()).flatten()
            )
        }).flatten();
        let control_biomass = ercc.as_ref().map(|e| {
            e.biomass(
                controls.as_ref().map(|c| c.organism_reads()).flatten()
            )
        }).flatten();
        let background_biomass = ercc.as_ref().map(|e| {
            e.biomass(
                background.as_ref().map(|b| b.other_reads()).flatten()
            )
        }).flatten();
        let output_biomass = ercc.as_ref().map(|e| {
            e.biomass(Some(output_reads))
        }).flatten();

        let control_reads = controls.as_ref().map(|c| c.organism_reads()).flatten();
        let control_reads_percent = controls.as_ref().map(|c| c.organism_reads_percent(input_reads)).flatten();
        

        Self {
            id: id.to_string(),

            model_id: None,
            run_id: None,
            run_date: None,
            library_tag: None,
            sample_group: None,
            sample_type: None,
            workflow_name: None,
            workflow_id: None,
            workflow_date: None,
            
            input_reads,
            input_bases,

            ercc_constructs: ercc.as_ref().map(|e| e.detected as u64),
            ercc_reads: ercc.as_ref().map(|e| e.alignments),
            ercc_reads_percent: ercc.as_ref().map(|e| (e.alignments as f64 / input_reads as f64)*100.0),

            qc_reads,
            qc_reads_percent,
            qc_bases,
            qc_bases_percent,

            host_reads: background.as_ref().map(|b| b.host_reads()).flatten(),
            host_reads_percent: background.as_ref().map(|b| b.host_reads_percent(input_reads)).flatten(),
            control_reads,
            control_reads_percent,
            background_reads: background.as_ref().map(|b| b.other_reads()).flatten(),
            background_reads_percent: background.as_ref().map(|b| b.other_reads_percent(input_reads)).flatten(),

            output_reads,
            output_reads_percent,
            output_bases,
            output_bases_percent,

            input_biomass,
            host_biomass,
            control_biomass,
            background_biomass,
            output_biomass,

            deduplicated_reads: deduplication.as_ref().map(|dedup| dedup.deduplicated as u64),
            deduplicated_percent: deduplication.as_ref().map(|dedup| (dedup.deduplicated as f64 / input_reads as f64)*100.0),

            adapter_trimmed_reads: qc.as_ref().map(|r| r.adapter_trimmed_reads()).flatten(),
            adapter_trimmed_percent: qc.as_ref().map(|r| r.adapter_trimmed_percent()).flatten(),
            low_complexity_reads: qc.as_ref().map(|r| r.low_complexity_reads()).flatten(),
            low_complexity_percent: qc.as_ref().map(|r| r.low_complexity_percent()).flatten(),
            mean_read_length_r1: qc.as_ref().map(|r| r.mean_read_length_r1()).flatten(),
            mean_read_length_r2: qc.as_ref().map(|r| r.mean_read_length_r2()).flatten(),
            min_length_reads: qc.as_ref().map(|r| r.min_length_reads()).flatten(),
            min_length_percent: qc.as_ref().map(|r| r.min_length_percent()).flatten(),
            low_quality_reads: qc.as_ref().map(|r| r.low_quality_reads()).flatten(),
            low_quality_percent: qc.as_ref().map(|r| r.low_quality_percent()).flatten(),
            q20_percent: qc.as_ref().map(|r| r.q20_percent()).flatten(),
            q30_percent: qc.as_ref().map(|r| r.q30_percent()).flatten(),        
        }
    }
}

// Background depletion

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Background {
    host: Option<HostBackground>,
    other: Option<OtherBackground>,
    config: BackgroundDepletionConfig
}
impl Background {
    pub fn from(id: &str, host: Option<ScrubbyReport>, other: Option<VircovSummary>, config: BackgroundDepletionConfig) -> Self {
        Self {
            host: HostBackground::from_scrubby(id, host),
            other: OtherBackground::from_vircov(id, other, &config),
            config
        }
    }
    pub fn host_reads(&self) -> Option<u64> {
        self.host.as_ref().map(|b| b.alignments)
    }
    pub fn host_reads_percent(&self, total: u64) -> Option<f64> {
        self.host.as_ref().map(|b| (b.alignments as f64 / total as f64)*100.0)
    }
    pub fn other_reads(&self) -> Option<u64> {
        self.other.as_ref().map(|b| b.alignments)
    }
    pub fn other_reads_percent(&self, total: u64) -> Option<f64> {
        self.other.as_ref().map(|b| (b.alignments as f64 / total as f64)*100.0)
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OtherBackground {
    pub id: String,
    pub alignments: u64,
    pub records: Vec<AlignmentRecord>
}

// Requires zero-output option in Vircov
pub trait BackgroundReport {
    fn from_vircov(id: &str, summary: Option<VircovSummary>, config: &BackgroundDepletionConfig) -> Option<Self> where Self: Sized;
}

impl BackgroundReport for OtherBackground {
    fn from_vircov(id: &str, summary: Option<VircovSummary>, config: &BackgroundDepletionConfig) -> Option<Self> {
        match summary { 
            None => return None, 
            Some(summary) => {

                let mut alignments = 0;

                let mut records = Vec::new();
                for record in &summary.records {
                    records.push(AlignmentRecord::from(
                        id, 
                        &record, 
                        RecordClass::from_background(&record, &config)
                    ));
                    alignments += record.scan_alignments;
                }

                Some(Self {
                    id: id.to_string(),
                    alignments,
                    records,
                })
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackgroundDepletionConfig {
    host_id: String,
    control_id: String,
    plasmid_id: String,
    univec_id: String,
    phage_id: String,
    rrna_id: String,
}
impl Default for BackgroundDepletionConfig {
    fn default() -> Self {
        Self {
            host_id: String::from("host::"),
            control_id: String::from("control::"),
            plasmid_id: String::from("plasmid::"),
            univec_id: String::from("univec::"),
            phage_id: String::from("phage::"),
            rrna_id: String::from("rrna::"),
        }
    }
}


pub trait BackgroundDepletionReport {
    fn from_scrubby(id: &str, report: Option<ScrubbyReport>) -> Option<Self> where Self: Sized;
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HostBackground {
    pub id: String,
    pub alignments: u64 
} 
impl BackgroundDepletionReport for HostBackground {
    fn from_scrubby(id: &str, report: Option<ScrubbyReport>) -> Option<Self> {
        if let Some(report) = report {
            Some(Self {
                id: id.to_string(),
                alignments: report.reads_removed
            })
        } else {
            None
        }  
    }
}

// Internal controls

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Controls {
    pub ercc: Option<ErccControl>,
    pub organism: Option<OrganismControl>,
    pub config: ControlsConfig
}
impl Controls {
    pub fn new(ercc: Option<ErccControl>, organism: Option<OrganismControl>, config: ControlsConfig) -> Self {
        Self { ercc, organism, config}
    }
    pub fn from(id: &str, summary: Option<VircovSummary>, config: ControlsConfig) -> Self {
        Self {
            ercc: ErccControl::from_vircov(id, summary.clone(), &config),
            organism: OrganismControl::from_vircov(id, summary.clone(), &config),
            config
        }
    }
    pub fn organism_reads(&self) -> Option<u64> {
        self.organism.as_ref().map(|o| o.alignments)
    }
    pub fn organism_reads_percent(&self, total: u64) -> Option<f64> {
        self.organism.as_ref().map(|o| (o.alignments as f64 / total as f64)*100.0)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ControlsConfig {
    ercc_id: String,
    ercc_mass: f64
}
impl Default for ControlsConfig {
    fn default() -> Self {
        Self {
            ercc_id: String::from("ERCC-"),
            ercc_mass: 25.0
        }
    }
}

pub trait ControlReport {
    fn from_vircov(id: &str, summary: Option<VircovSummary>, config: &ControlsConfig) -> Option<Self> where Self: Sized;
}

// Requires zero-output option in Vircov
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErccControl {
    pub constructs: usize,
    pub detected: usize,
    pub alignments: u64,
    pub reads: u64,
    pub mass_per_read: f64,
    pub records: Vec<AlignmentRecord>,
} 
impl ErccControl {
    pub fn biomass(&self, reads: Option<u64>) -> Option<f64> {
        reads.map(|reads| {
            self.mass_per_read * reads as f64
        })
    }
}
impl ControlReport for ErccControl {
    fn from_vircov(id: &str, summary: Option<VircovSummary>, config: &ControlsConfig) -> Option<Self> {

        match summary { 
            None => return None, 
            Some(summary) => {

                let mut constructs = 0;
                let mut aligned = 0;
                let mut reads = 0;
                let mut alignments = 0;
        
                let mut records = Vec::new();
                for record in &summary.records {
                    if record.reference.starts_with(&config.ercc_id) {
                        if record.scan_alignments > 0 {
                            aligned += 1
                        }
                        constructs += 1;
        
                        reads += record.scan_reads;
                        alignments += record.scan_alignments;
                        records.push(AlignmentRecord::from(id, &record, RecordClass::ErccControl))
                    }
                }
                
                let mass_per_read = if alignments == 0 {
                    0.0
                } else {
                    config.ercc_mass / alignments as f64 // check if need to use reads?
                };
                 
                
                let ercc = Self {
                    constructs, 
                    detected: aligned, 
                    reads, 
                    alignments,
                    mass_per_read,
                    records,
                };
        
                if constructs == 0 {
                    None
                } else {
                    Some(ercc)
                }
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrganismControl {
    pub id: String,
    pub alignments: u64,
    pub records: Vec<AlignmentRecord>
}
impl ControlReport for OrganismControl {
    fn from_vircov(id: &str, summary: Option<VircovSummary>, config: &ControlsConfig) -> Option<Self> {

        match summary { 
            None => return None, 
            Some(summary) => {
                let mut organisms = 0;
                let mut alignments = 0;

                let mut records = Vec::new();
                for record in &summary.records {
                    if !record.reference.starts_with(&config.ercc_id) {
                        records.push(AlignmentRecord::from(id, &record, RecordClass::InternalControl));
                        organisms += 1;
                        alignments += record.scan_alignments;
                    }
                }
                if organisms == 0 {
                    None
                } else {
                    Some(Self { id: id.to_string(), alignments, records })
                }
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RecordClass {
    #[serde(rename = "ercc")]
    ErccControl,
    #[serde(rename = "control")]
    InternalControl,
    #[serde(rename = "host")]
    HostBackground,
    #[serde(rename = "rrna")]
    RibosomalBackground,
    #[serde(rename = "plasmid")]
    PlasmidBackground,
    #[serde(rename = "univec")]
    UnivecBackground,
    #[serde(rename = "phage")]
    PhageBackground,
    #[serde(rename = "other")]
    OtherBackground
}
impl RecordClass {
    pub fn from_background(record: &VircovRecord, config: &BackgroundDepletionConfig) -> Self {
        if record.reference.starts_with(&config.phage_id) {
            return RecordClass::PhageBackground
        } else if record.reference.starts_with(&config.plasmid_id) {
            return RecordClass::PlasmidBackground
        } else if record.reference.starts_with(&config.host_id) {
            return RecordClass::HostBackground
        } else if record.reference.starts_with(&config.control_id) {
            return RecordClass::InternalControl
        } else if record.reference.starts_with(&config.univec_id) {
            return RecordClass::UnivecBackground
        } else if record.reference.starts_with(&config.rrna_id) {
            return RecordClass::RibosomalBackground
        } else {
            return RecordClass::OtherBackground
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentRecord {
    pub id: String,
    pub reference: String,
    pub alignments: u64,
    pub reads: u64,
    pub coverage: f64,
    pub class: RecordClass
}
impl AlignmentRecord {
    pub fn from(id: &str, record: &VircovRecord, class: RecordClass) -> Self {
        Self {
            id: id.to_string(),
            reference: record.reference.clone(),
            alignments: record.scan_alignments,
            reads: record.scan_reads,
            coverage: record.scan_coverage,
            class
        }
    }
}

use std::collections::HashMap;

/// QC status enum
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum QcStatus {
    Fail,
    Pass,
    Ok,
}

/// Numeric threshold for fail/pass/ok checks
#[derive(Debug, Clone)]
pub struct Threshold<T> {
    pub fail_max: T,
    pub pass_max: T,
}

/// Fastp-specific thresholds
#[derive(Debug, Clone)]
pub struct FastpThresholds {
    pub q20_percent: Threshold<f64>,
    pub q30_percent: Threshold<f64>,
}

/// Preset configuration for QC
#[derive(Debug, Clone)]
pub struct QcConfig {
    pub input_reads: Threshold<u64>,
    pub output_reads: Threshold<u64>,
    pub ercc_constructs: Threshold<u64>,
    pub phage_coverage: Threshold<f64>,
    pub fastp: FastpThresholds,
}

impl QcConfig {
    pub fn illumina_pe_sr() -> Self {
        QcConfig {
            input_reads: Threshold { fail_max: 1_000_000, pass_max: 5_000_000 },
            output_reads: Threshold { fail_max: 10_000, pass_max: 100_000 },
            ercc_constructs: Threshold { fail_max: 40, pass_max: 60 },
            phage_coverage: Threshold { fail_max: 0.4, pass_max: 0.6 },
            fastp: FastpThresholds {
                q20_percent: Threshold { fail_max: 80.0, pass_max: 90.0 },
                q30_percent: Threshold { fail_max: 70.0, pass_max: 85.0 },
            },
        }
    }
}


/// Preset configuration for positive controls
#[derive(Debug, Clone)]
pub struct PositiveControlConfig {
    pub bacteria: Vec<String>,
    pub viruses: Vec<String>,
    pub eukaryota: Vec<String>
}

impl PositiveControlConfig {
    pub fn metagp() -> Self {
        Self {
            bacteria: vec![
                String::from("Truepera radiovictrix"),
                String::from("Imtechella halotolerans"),
                String::from("Allobacillus halotolerans"),
            ],
            viruses: vec![
                String::from("Betacoronavirus muris"),
                String::from("Orthopoxvirus vaccinia"),
            ],
            eukaryota: vec![
                String::from("Saccharomyces cerevisiae"),
                String::from("Pneumocystis jirovecii")
            ]
        }
    }
}

/// Trait to extract “detected” names from any source
pub trait DetectionSource {
    fn detected_names(&self) -> Vec<&str>;
}

/// Implement for your two sources:
impl DetectionSource for PathogenDetection {
    fn detected_names(&self) -> Vec<&str> {
        self.profile.iter().map(|rec| rec.name.as_str()).collect()
    }
}
impl DetectionSource for Vec<Taxon> {
    fn detected_names(&self) -> Vec<&str> {
        self.iter().map(|tax| tax.name.as_str()).collect()
    }
}

/// Our summary struct stays the same:
#[derive(Debug, Clone)]
pub struct PositiveControlSummary {
    pub id:        String,
    pub bacteria:  HashMap<String, bool>,
    pub viruses:   HashMap<String, bool>,
    pub eukaryota: HashMap<String, bool>,
}

/// Generic builder that now carries an optional `id`
pub struct PositiveControlSummaryBuilder<'a, S>
where
    S: DetectionSource,
{
    source: &'a S,
    config: &'a PositiveControlConfig,
    id:     String,
}

impl<'a, S> PositiveControlSummaryBuilder<'a, S>
where
    S: DetectionSource,
{
    /// New builder: no id yet
    pub fn new(source: &'a S, config: &'a PositiveControlConfig, id: impl Into<String>) -> Self {
        Self {
            source,
            config,
            id: id.into(),
        }
    }


    /// Build the actual summary
    pub fn build(self) -> PositiveControlSummary {
        let detected = self.source.detected_names();

        let make_map =
            |controls: &Vec<String>, label: &str| -> HashMap<String, bool> {
                controls
                    .iter()
                    .map(|spec| {
                        let seen = detected.contains(&spec.as_str());
                        log::info!(
                            "[{}] {} control `{}` detected: {}",
                            self.id,
                            label,
                            spec,
                            seen
                        );
                        (spec.clone(), seen)
                    })
                    .collect()
            };

        PositiveControlSummary {
            id:        self.id.to_owned(),
            bacteria:  make_map(&self.config.bacteria,  "bacteria"),
            viruses:   make_map(&self.config.viruses,   "virus"),
            eukaryota: make_map(&self.config.eukaryota, "eukaryote"),
        }
    }
}


pub fn write_positive_control_summaries<P: AsRef<Path>>(
    summaries: &[PositiveControlSummary],
    config: &PositiveControlConfig,
    output_path: P,
) -> std::io::Result<()> {
    // Open the output file
    let file = File::create(output_path)?;
    let mut w = BufWriter::new(file);

    // Prepare sorted column order
    let mut bacteria   = config.bacteria.clone();
    let mut eukaryota  = config.eukaryota.clone();
    let mut viruses     = config.viruses.clone();
    bacteria.sort();
    eukaryota.sort();
    viruses.sort();

    // Write header
    write!(w, "id")?;
    for name in bacteria.iter().chain(eukaryota.iter()).chain(viruses.iter()) {
        write!(w, "\t{}", name)?;
    }
    w.write_all(b"\n")?;

    // Write each row
    for summary in summaries {
        write!(w, "{}", summary.id)?;
        // Bacteria columns
        for name in &bacteria {
            let detected = summary
                .bacteria
                .get(name)
                .copied()
                .unwrap_or(false);
            write!(w, "\t{}", detected)?;
        }
        // Eukaryota columns
        for name in &eukaryota {
            let detected = summary
                .eukaryota
                .get(name)
                .copied()
                .unwrap_or(false);
            write!(w, "\t{}", detected)?;
        }
        // Virus columns
        for name in &viruses {
            let detected = summary
                .viruses
                .get(name)
                .copied()
                .unwrap_or(false);
            write!(w, "\t{}", detected)?;
        }
        w.write_all(b"\n")?;
    }

    w.flush()?;
    Ok(())
}


/// Summary results for each category
#[derive(Debug, Serialize, Deserialize)]
pub struct QualityControlSummary {
    pub id: String,
    pub input_reads: QcStatus,
    pub output_reads: QcStatus,
    pub ercc_constructs: QcStatus,
    pub phage_coverage: Option<QcStatus>,
    pub fastp_status: QcStatus,
    pub fastp_details: HashMap<String, QcStatus>,
}
impl QualityControlSummary {
    pub fn to_json(&self, path: &PathBuf) -> Result<(), WorkflowError> {
        let writer = BufWriter::new(File::create(path)?);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }
    pub fn from_json(path: &PathBuf) -> Result<Self, WorkflowError> {
        let reader = BufReader::new(File::open(path)?);
        let quality_control = serde_json::from_reader(reader)?;
        Ok(quality_control)
    }
}

/// Builder for QualityControlSummary
pub struct QualityControlSummaryBuilder<'a> {
    id: &'a str,
    input_reads: u64,
    output_reads: u64,
    ercc_constructs: Option<u64>,
    phage_records: Vec<AlignmentRecord>,
    fastp_fields: HashMap<String, f64>,
    config: QcConfig,
}

impl<'a> QualityControlSummaryBuilder<'a> {
    pub fn new(
        id: &'a str,
        input_reads: u64,
        output_reads: u64,
        ercc_constructs: Option<u64>,
        phage_records: Vec<AlignmentRecord>,
        fastp_fields: HashMap<String, f64>,
        config: QcConfig,
    ) -> Self {
        Self {
            id,
            input_reads,
            output_reads,
            ercc_constructs,
            phage_records,
            fastp_fields,
            config,
        }
    }

    fn status_from_threshold_u64(&self, value: u64, th: &Threshold<u64>) -> QcStatus {
        if value <= th.fail_max {
            QcStatus::Fail
        } else if value <= th.pass_max {
            QcStatus::Pass
        } else {
            QcStatus::Ok
        }
    }

    fn status_from_threshold_f64(&self, value: f64, th: &Threshold<f64>) -> QcStatus {
        if value < th.fail_max {
            QcStatus::Fail
        } else if value < th.pass_max {
            QcStatus::Pass
        } else {
            QcStatus::Ok
        }
    }

    /// Build the summary
    pub fn build(self) -> QualityControlSummary {

        // input_reads
        let in_status = self.status_from_threshold_u64(
            self.input_reads,
            &self.config.input_reads,
        );

        // output_reads
        let out_status = self.status_from_threshold_u64(
            self.output_reads,
            &self.config.output_reads,
        );

        // ercc_constructs
        let ercc_status = self.ercc_constructs
            .map(|c| self.status_from_threshold_u64(c, &self.config.ercc_constructs))
            .unwrap_or(QcStatus::Fail);

        // phage coverage – pass if either T4-DNA or MS2-RNA passes, only fail if both fail
        let phage_status = {
            // collect statuses for both phages
            let statuses: Vec<QcStatus> = self.phage_records.iter()
                .filter(|r| matches!(r.reference.as_str(), "T4-DNA" | "MS2-RNA"))
                .map(|r| {
                    self.status_from_threshold_f64(
                        r.coverage,
                        &self.config.phage_coverage,
                    )
                })
                .collect();

            match statuses.as_slice() {
                [] => None,                                // no phage records at all
                _ if statuses.iter().any(|s| *s == QcStatus::Ok)   => Some(QcStatus::Ok),
                _ if statuses.iter().any(|s| *s == QcStatus::Pass) => Some(QcStatus::Pass),
                _                                                  => Some(QcStatus::Fail),
            }
        };

        // fastp per-field statuses
        let mut details = HashMap::new();
        for (k, v) in &self.fastp_fields {
            let th = match k.as_str() {
                "q20_percent" => &self.config.fastp.q20_percent,
                "q30_percent" => &self.config.fastp.q30_percent,
                _ => continue,
            };
            details.insert(k.clone(), self.status_from_threshold_f64(*v, th));
        }

        // combined: if any fail => Fail, all ok => Ok, else Pass
        
        let fastp_status = if details.values().any(|s| *s == QcStatus::Fail) {
            QcStatus::Fail
        } else if details.values().all(|s| *s == QcStatus::Ok) {
            QcStatus::Ok
        } else {
            QcStatus::Pass
        };

        QualityControlSummary {
            id: self.id.to_string(),
            input_reads: in_status,
            output_reads: out_status,
            ercc_constructs: ercc_status,
            phage_coverage: phage_status,
            fastp_status,
            fastp_details: details,
        }
    }
}
