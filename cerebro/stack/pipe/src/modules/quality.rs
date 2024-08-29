use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;

use csv::WriterBuilder;
use scrubby::report::ScrubbyReport;
use serde::{Deserialize, Serialize, Serializer};
use vircov::vircov::{VircovRecord, VircovSummary};

use crate::error::WorkflowError;
use crate::parsers::fastp::FastpReport;
use crate::nextflow::panviral::PanviralOutput;

// Main module wrapper functions

pub fn write_quality_tsv(qc_data: &Vec<QualityControl>, reads: &PathBuf, controls: &PathBuf) -> Result<(), WorkflowError> {
    
    let mut reads_writer = WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(File::create(reads)?);

    for qc in qc_data {
        reads_writer.serialize(&qc.reads)?;
    }
    
    reads_writer.flush()?;
    
    let mut alignment_writer = WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(File::create(controls)?);

    for qc in qc_data {
        if let Some(controls) = &qc.controls {
            if let Some(ercc) = &controls.ercc {
                for record in &ercc.records {
                    alignment_writer.serialize(record)?;
                }
            }
            if let Some(organism) = &controls.organism {
                for record in &organism.records {
                    alignment_writer.serialize(record)?;
                }
            }
        }
    }
    
    alignment_writer.flush()?;

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
    pub reads: Reads,
    pub background: Option<Background>,
    pub controls: Option<InternalControls>
}
impl QualityControl {
    pub fn from_panviral(output: &PanviralOutput) -> Self {

        let controls = InternalControls::from(
            &output.id,
            output.controls.clone(), 
            InternalControlConfig::default()
        );        

        let background = Background::from(
            &output.id,
            output.host.clone(),
        );

        let reads = Reads::from(
            &output.id,
            output.reads.clone(),
            controls.as_ref().and_then(
                |c| c.ercc.clone()
            )
        );

        Self {
            id: output.id.clone(),
            reads,
            background,
            controls
        }
    }
    pub fn to_json(&self, output: &PathBuf) -> Result<(), WorkflowError> {
        let writer = BufWriter::new(File::create(output)?);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }

    pub fn from_json(input: &PathBuf) -> Result<Self, WorkflowError> {
        let reader = BufReader::new(File::open(input)?);
        let quality_control = serde_json::from_reader(reader)?;
        Ok(quality_control)
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

fn round_two<T, S>(x: &T, s: S) -> Result<S::Ok, S::Error>
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
    fn output_reads(&self) -> u64;
    fn output_reads_percent(&self) -> f64;
    fn output_bases(&self) -> u64;
    fn output_bases_percent(&self) -> f64;
    fn input_biomass(&self, ercc: Option<Ercc>) -> Option<f64>;
    fn output_biomass(&self, ercc: Option<Ercc>) -> Option<f64>;
    fn deduplicated_reads(&self) -> Option<u64>;
    fn deduplicated_percent(&self) -> Option<f64>;
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
    fn output_reads(&self) -> u64 {
        self.summary.after.reads
    }
    fn output_reads_percent(&self) -> f64 {
        if self.summary.before.reads == 0 {
            0.0
        } else {
            (self.summary.after.reads as f64 / self.summary.before.reads as f64)*100.0
        }
    }
    fn output_bases(&self) -> u64 {
        self.summary.after.bases
    }
    fn output_bases_percent(&self) -> f64 {
        if self.summary.before.bases == 0 {
            0.0
        } else {
            (self.summary.after.bases as f64 / self.summary.before.bases as f64)*100.0
        }
    }
    fn input_biomass(&self, ercc: Option<Ercc>) -> Option<f64> {
        ercc.as_ref().map(|ercc| self.summary.before.reads as f64 * ercc.mass_per_read)
    }
    fn output_biomass(&self, ercc: Option<Ercc>) -> Option<f64> {
        ercc.as_ref().map(|ercc| self.summary.after.reads as f64 * ercc.mass_per_read)
    }
    fn deduplicated_reads(&self) -> Option<u64> {
        self.duplication.as_ref().map(|dup| (self.summary.after.reads as f64 * dup.rate) as u64)
    }
    fn deduplicated_percent(&self) -> Option<f64> {
        self.duplication.as_ref().map(|dup| dup.rate * 100.0)
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
        Some(self.summary.after.q20)
    }
    fn q30_percent(&self) -> Option<f64> {
        Some(self.summary.after.q30)
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Reads {
    pub id: String,

    pub input_reads: u64,
    pub input_bases: u64,

    pub qc_reads: u64,
    #[serde(serialize_with = "round_two")]
    pub qc_reads_percent: f64,
    pub qc_bases: u64,
    #[serde(serialize_with = "round_two")]
    pub qc_bases_percent: f64,
    
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ercc_constructs: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ercc_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub ercc_reads_percent: Option<f64>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub organism_controls_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub organism_controls_reads_percent: Option<f64>,
    
    #[serde(skip_serializing_if = "Option::is_none")]
    pub host_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub host_reads_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub other_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub other_reads_percent: Option<f64>,
    
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub input_biomass: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub qc_biomass: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub host_biomass: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub other_biomass: Option<f64>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub deduplicated_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub deduplicated_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub adapter_trimmed_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub adapter_trimmed_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub low_complexity_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub low_complexity_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mean_read_length_r1: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mean_read_length_r2: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_length_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub min_length_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub low_quality_reads: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub low_quality_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub q20_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none", serialize_with = "round_two")]
    pub q30_percent: Option<f64>,
}
impl Reads {
    pub fn from<R: ReadReport>(id: &str, report: R, ercc: Option<Ercc>) -> Self {

        Self {
            id: id.to_string(),
            input_reads: report.input_reads(),
            input_bases: report.input_reads(),
            qc_reads: report.output_reads(),
            qc_reads_percent: report.output_reads_percent(),
            qc_bases: report.output_bases(),
            qc_bases_percent: report.output_bases_percent(),
            input_biomass: report.input_biomass(ercc.clone()),
            qc_biomass: report.output_biomass(ercc.clone()),
            deduplicated_reads: report.deduplicated_reads(),
            deduplicated_percent: report.deduplicated_percent(),
            adapter_trimmed_reads: report.adapter_trimmed_reads(),
            adapter_trimmed_percent: report.adapter_trimmed_percent(),
            low_complexity_reads: report.low_complexity_reads(),
            low_complexity_percent: report.low_complexity_percent(),
            mean_read_length_r1: report.mean_read_length_r1(),
            mean_read_length_r2: report.mean_read_length_r2(),
            min_length_reads: report.min_length_reads(),
            min_length_percent: report.min_length_percent(),
            low_quality_reads: report.low_quality_reads(),
            low_quality_percent: report.low_quality_percent(),
            q20_percent: report.q20_percent(),
            q30_percent: report.q30_percent(),        
        }
    }
}

// Background depletion

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Background {
    host: Option<Host>
}
impl Background {
    pub fn from(id: &str, host: Option<ScrubbyReport>) -> Option<Self> {
        match host {
            Some(ref report) => Some(
                Self {
                    host: Some(Host::from_scrubby(id, report)),
                }
            ),
            None => None
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentRecord {
    pub id: String,
    pub reference: String,
    pub alignments: u64,
    pub reads: u64,
    pub coverage: f64
}
impl AlignmentRecord {
    pub fn from(id: &str, record: &VircovRecord) -> Self {
        Self {
            id: id.to_string(),
            reference: record.reference.clone(),
            alignments: record.scan_alignments,
            reads: record.scan_reads,
            coverage: record.scan_coverage
        }
    }
}

pub trait HostBackgroundReport {
    fn from_scrubby(id: &str, report: &ScrubbyReport) -> Self;
}

// Requires zero-output option in Vircov
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Host {
    pub id: String,
    pub reads: u64 
} 
impl HostBackgroundReport for Host {
    fn from_scrubby(id: &str, report: &ScrubbyReport) -> Self {
        
         Self {
            id: id.to_string(),
            reads: report.reads_removed
        }
    }
}

// Internal controls

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InternalControls {
    pub ercc: Option<Ercc>,
    pub organism: Option<Organism>,
    pub config: InternalControlConfig
}
impl InternalControls {
    pub fn from(id: &str, summary: Option<VircovSummary>, config: InternalControlConfig) -> Option<Self> {
        match summary {
            Some(ref summary) => Some(
                Self {
                    ercc: Ercc::from_vircov(id, summary, &config),
                    organism: Organism::from_vircov(id, summary, &config),
                    config
                }
            ),
            None => None
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InternalControlConfig {
    ercc_id: String,
    ercc_mass: f64
}
impl Default for InternalControlConfig {
    fn default() -> Self {
        Self {
            ercc_id: String::from("ERCC-"),
            ercc_mass: 25.0
        }
    }
}

pub trait InternalControlReport {
    fn from_vircov(id: &str, summary: &VircovSummary, config: &InternalControlConfig) -> Option<Self> where Self: Sized;
}

// Requires zero-output option in Vircov
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ercc {
    pub constructs: usize,
    pub aligned: usize,
    pub alignments: u64,
    pub reads: u64,
    pub mass_per_read: f64,
    pub records: Vec<AlignmentRecord>,
} 
impl InternalControlReport for Ercc {
    fn from_vircov(id: &str, summary: &VircovSummary, config: &InternalControlConfig) -> Option<Self> {

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
                records.push(AlignmentRecord::from(id, &record))
            }
        }

        let mass_per_read = config.ercc_mass / alignments as f64;  // check if need to use reads?
        
        let ercc = Self {
            constructs, 
            aligned, 
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Organism {
    records: Vec<AlignmentRecord>
}
impl InternalControlReport for Organism {
    fn from_vircov(id: &str, summary: &VircovSummary, config: &InternalControlConfig) -> Option<Self> {

        let mut organisms = 0;

        let mut records = Vec::new();
        for record in &summary.records {
            if !record.reference.starts_with(&config.ercc_id) {
                records.push(AlignmentRecord::from(id, &record));
                organisms += 1;
            }
        }
        if organisms == 0 {
            None
        } else {
            Some(Self { records })
        }
    }
}

