use itertools::Itertools;
use memchr::memmem;
use rayon::prelude::*;
use std::{path::PathBuf, collections::{HashSet, HashMap}, str::from_utf8};
use needletail::{parse_fastx_file, parser::write_fastq};
use crate::utils::{get_compression_writer, get_seq_record_identifier, get_file_stem};
use crate::error::WorkflowError;


/// Implements a processing pipeline with initial sequence read processing and collapsing of 
/// unique molecular identifier (UMI) sequences from the NEB UNIQUE DUAL INDEX protocol, which
/// is implemented in the Cerebro metagenomics Illumina sequencing assay. 
/// 
/// This pipeline assumes that the sample sheet for BCLConvert demultiplexing is configured 
/// to the correct threshold cycles. Using the UMI protocol on NextSeq, in this case the setting
/// is: OverrideCycles,Y151;I8U12;I8;Y151
/// 
/// UMI sequence files (I1, I2) are produced by the demultiplexing, where `I1` contains the UMI sequences 
/// lead by the index sequence (I8U12). `I2` can be discarded. Forward and reverse reads contain the UMI 
/// sequence in the read identifier as the last field split by `:`

pub struct Umi { 
    output_format: Option<niffler::compression::Format>,
    compression_level: Option<niffler::compression::Level>
}
impl Umi {
    pub fn new(output_format: Option<niffler::compression::Format>, compression_level: Option<niffler::compression::Level>) -> Self {
        Self {
            output_format,
            compression_level
        }
    }
    // Trim head bases from UMI sequences in I1 read files (NEB Protocol I8U12)
    pub fn trim_umi_index(&self, umi_fastq: &PathBuf, output: &PathBuf, head_trim: usize) -> Result<(), WorkflowError> {
       
        let mut reader = parse_fastx_file(umi_fastq)?;
        let mut writer = get_compression_writer(output, &self.output_format, &self.compression_level)?;
        
        while let Some(record) = reader.next() {
            let rec = record?;
            let rec_len = rec.num_bases();

            let qual = match rec.qual() { 
                Some(qual) => Some(&qual[head_trim..rec_len]),
                None => return Err(WorkflowError::QualityScoresNotPresent)
            };

            write_fastq(rec.id(), &rec.seq()[head_trim..rec_len], qual, &mut writer, rec.line_ending())?;
        }
        Ok(())
    }
    // Append the UMI sequence from the trimmed UMI fastq to the start of R1 reads (NEB Protocol I8U12) for testing Calib
    pub fn prepare_calib(&self, fastq: &PathBuf, output: &PathBuf) -> Result<(), WorkflowError> {       
        
        let mut reader = parse_fastx_file(fastq)?;
        let mut writer = get_compression_writer(output, &self.output_format, &self.compression_level)?;

        while let Some(record) = reader.next() {
            let rec = record?;
            let rec_id = get_seq_record_identifier(&rec)?;

            let mut umi_seq = match rec_id.split(":").collect::<Vec<&str>>().split_last() {
                Some((last_field, _)) => last_field.to_string(),
                None => return Err(WorkflowError::UmiFieldNotFound(rec_id))
            };
            // Calib does not take read quality into consideration, 
            // we can therefore just add some fake quality data:
            let umi_qual = String::from("C").repeat(umi_seq.len());

            umi_seq.push_str(from_utf8(&rec.seq())?);
                    
            let rec_qual = match rec.qual() { 
                Some(qual) => qual,
                None => return Err(WorkflowError::QualityScoresNotPresent)
            };

            let qual = [umi_qual.as_bytes(), rec_qual].concat();

            write_fastq(rec.id(), &umi_seq.as_bytes(), Some(&qual), &mut writer, rec.line_ending())?;
        }
        Ok(())
    }
    /// Cluster sequences naively by their total forward sequence
    pub fn naive_dedup(&self, input: &Vec<PathBuf>, output: &Vec<PathBuf>, forward: &PathBuf, summary: &Option<PathBuf>, prepend_umi: bool) -> Result<Vec<()>, WorkflowError> {

        if !input.len() == output.len() {
            return Err(WorkflowError::InputOutputFileMismatch)
        }

        let mut reader = parse_fastx_file(forward)?;
        let mut clusters: HashMap<String, Vec<NaiveClusterSequence>> = HashMap::new();

        let mut total_reads: usize = 0;
        while let Some(record) = reader.next() {
            let rec = record?;
            let rec_id = get_seq_record_identifier(&rec)?;

            let seq = match prepend_umi {
                true => {
                    let mut umi_seq = match rec_id.split(":").collect::<Vec<&str>>().split_last() {
                        Some((last_field, _)) => last_field.to_string(),
                        None => return Err(WorkflowError::UmiFieldNotFound(rec_id))
                    };
                    umi_seq.push_str(from_utf8(&rec.seq())?);
                    umi_seq
                }
                false => {
                    from_utf8(&rec.seq())?.to_string()
                }
            };
            

            clusters.entry(seq)
                .and_modify(|cluster| cluster.push(NaiveClusterSequence::new(&rec_id, &rec.qual())))
                .or_insert(Vec::from([NaiveClusterSequence::new(&rec_id, &rec.qual())]));

            total_reads += 1;
        }
        
        let mut total_clusters: usize = 0;
        let mut singleton_clusters: usize = 0;
        let mut cluster_records = Vec::new();
        let mut cluster_deduplicated_identifiers = HashSet::new();
        for reads in clusters.values() {
            let cluster_size = reads.len();

            total_clusters += 1;
            
            if cluster_size == 1 {
                singleton_clusters += 1;
            };

            if let Some(_) = summary {
                cluster_records.push(NaiveClusterRecord {
                    id: total_clusters.to_string(),
                    size: cluster_size,
                    reads: reads.iter().map(|c| c.id.to_owned()).collect::<Vec<String>>().join(";")
                });
            }

            let qualities = reads.iter().sorted_by(|a, b| b.mean_quality.partial_cmp(&a.mean_quality).unwrap() ).collect::<Vec<&NaiveClusterSequence>>();
            
            let best_read_mean_quality_forward = match qualities.first() {
                Some(best_quality_record) => best_quality_record.id.clone(),
                None => return Err(WorkflowError::CalibClusterPickNotAvailable)
            };

            cluster_deduplicated_identifiers.insert(best_read_mean_quality_forward);

        }
        log::info!(
            "Detected {} naive clusters ({} singleton clusters) out of {} reads (retaining: {:.2} %)", 
            &total_clusters, &singleton_clusters, &total_reads, (total_clusters  as f32/total_reads as f32)*100.
        );

        if let Some(file) = summary {
            cluster_records.sort_by(|a, b| b.size.cmp(&a.size));
            let mut writer = csv::WriterBuilder::new().delimiter(b'\t').has_headers(true).from_path(&file).unwrap();
            for cluster in cluster_records {
                writer.serialize(&cluster).unwrap();
            }
        }

        log::info!("Removing duplicated reads from input files");
        let result: Result<Vec<()>, WorkflowError> = input.par_iter().enumerate().map(|(idx, fq)| -> Result<(), WorkflowError> {
            let mut reader = parse_fastx_file(fq)?;
            let mut writer = get_compression_writer(&output[idx], &self.output_format, &self.compression_level)?;
            
            while let Some(record) = reader.next() {
                let rec = record?;
                let rec_id = get_seq_record_identifier(&rec)?;
                
                if cluster_deduplicated_identifiers.contains(&rec_id) {
                    rec.write(&mut writer, None)?;
                }
            }
            Ok(())
        }).collect();

        result
    }
    /// Calib cluster assessment
    pub fn calib_dedup(&self, fastq: &Vec<PathBuf>, output: &Vec<PathBuf>, calib_clusters: &PathBuf, summary: Option<PathBuf>, identifiers: Option<PathBuf>) -> Result<Vec<()>, WorkflowError> {

        if fastq.len() != output.len() {
            return Err(WorkflowError::InputOutputFileMismatch)
        }

        let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&calib_clusters).unwrap();
        let mut clusters: HashMap<u64, Vec<CalibRecord>> = HashMap::new();
        let mut total_reads: usize = 0;

        for record in reader.deserialize() {
            let rec: CalibClusterReportRecord = record.map_err(|err| WorkflowError::CalibClusterRecordNotValid(err))?;

            clusters.entry(rec.cluster_id)
                    .and_modify(|cluster| cluster.push(CalibRecord::from(&rec)))
                    .or_insert(Vec::from([CalibRecord::from(&rec)]));

            total_reads += 1;
        }

        let mut total_clusters: usize = 0;
        let mut singleton_clusters: usize = 0;
        let mut cluster_records = Vec::new();
        let mut cluster_deduplicated_identifiers = HashSet::new();

        for (cluster_id, records) in clusters {
            let cluster_size = records.len();

            total_clusters += 1;
            if cluster_size == 1 {
                singleton_clusters += 1;
            };

            cluster_records.push(CalibClusterRecord {
                id: cluster_id.to_string(),
                size: cluster_size,
                reads: records.iter().map(|rec| rec.read_id.to_string()).collect::<Vec<String>>().join(";"),
                read_quality: records.iter().map(|rec| rec.mean_quality_paired.to_string()).collect::<Vec<String>>().join(";")
            });

            let best_read_mean_quality_paired = match records.iter().sorted_by(|a, b| b.mean_quality_paired.partial_cmp(&a.mean_quality_paired).unwrap() ).collect::<Vec<&CalibRecord>>().first() {
                Some(best_quality_record) => best_quality_record.read_id.clone(),
                None => return Err(WorkflowError::CalibClusterPickNotAvailable)
            };

            cluster_deduplicated_identifiers.insert(best_read_mean_quality_paired);
        }

        log::info!(
            "Calib detected {} clusters ({} singleton) out of {} reads", 
            &total_clusters, &singleton_clusters, &total_reads
        );
        log::info!(
            "Deduplicated clusters retaining {} reads ({:.2} %)", 
            &cluster_deduplicated_identifiers.len(), (cluster_deduplicated_identifiers.len() as f32/total_reads as f32)*100.
        );


        if let Some(summary) = summary {
            log::info!("Writing cluster summary for Calib cluster output");
            cluster_records.sort_by(|a, b| b.size.cmp(&a.size));
            let mut writer = csv::WriterBuilder::new().delimiter(b'\t').has_headers(true).from_path(&summary).unwrap();
            for cluster in cluster_records {
                writer.serialize(&cluster).unwrap();
            }
        }

        if let Some(identifiers) = identifiers {
            log::info!("Writing deduplicated read identifiers to file");
            let mut writer = csv::WriterBuilder::new().delimiter(b'\t').has_headers(false).from_path(&identifiers).unwrap();
            for read_id in &cluster_deduplicated_identifiers {
                writer.serialize(&CalibDeduplicatedId { read_id: read_id.to_string() }).unwrap();
            }
        }

        log::info!("Removing duplicated reads from input files");
        
        let result: Result<Vec<()>, WorkflowError> = fastq.par_iter().enumerate().map(|(idx, fq)| -> Result<(), WorkflowError> {
            let mut reader = parse_fastx_file(fq)?;
            let mut writer = get_compression_writer(&output[idx], &self.output_format, &self.compression_level)?;
            
            while let Some(record) = reader.next() {
                let rec = record?;
                let rec_id = get_seq_record_identifier(&rec)?;
                
                if cluster_deduplicated_identifiers.contains(&rec_id) {
                    rec.write(&mut writer, None)?;
                }
            }
            Ok(())
        }).collect();

        result
        
    }
    // Sanity check - is the UMI sequence still in the read or has it been trimmed off
    pub fn check_umi_in_read(&self, fastq: &Vec<PathBuf>, output: &PathBuf) -> Result<(), WorkflowError> {

        let records: Result<Vec<CheckRecord>, WorkflowError> = fastq.par_iter().map(|fq| -> Result<CheckRecord, WorkflowError> {

            let name = get_file_stem(fq)?;
            
            let mut reader = parse_fastx_file(fq)?;

            let mut total_found: usize = 0;
            let mut total_reads: usize = 0;

            while let Some(record) = reader.next() {
                let rec = record?;
    
                let rec_id = get_seq_record_identifier(&rec)?;
                let umi_seq = match rec_id.split(':').last() {
                    Some(umi_seq) => umi_seq.as_bytes(),
                    None => return Err(WorkflowError::RecordIdentifierNotParsed)
                };
                
                let found = memmem::find(&rec.seq(), umi_seq);
    
                if let Some(_) = found {
                    total_found += 1;
    
                }
                total_reads += 1;
            }
            Ok(CheckRecord {
                file_name: name,
                total_reads,
                umi_in_read: total_found,
                umi_in_read_percent: (total_found as f32 / total_reads as f32)*100.
            })
        }).collect();

        match records {
            Ok(mut records) => {
                records.sort_by(|a, b| a.file_name.cmp(&b.file_name));
                let mut writer = csv::WriterBuilder::new().delimiter(b'\t').has_headers(true).from_path(&output).unwrap();
                for cluster in records {
                    writer.serialize(&cluster).unwrap();
                }
                Ok(())
            },
            Err(err) => Err(err)
        }
    }
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct NaiveClusterRecord {
    id: String,
    size: usize,
    reads: String
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct NaiveClusterSequence {
    id: String,
    mean_quality: f32
}
impl NaiveClusterSequence {
    pub fn new(rec_id: &String, rec_qual: &Option<&[u8]>) -> Self {
        Self {
            id: rec_id.to_owned(),
            mean_quality: mean_qual(rec_qual.unwrap())
        }
    }
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct CalibClusterRecord {
    id: String,
    size: usize,
    reads: String,
    read_quality: String
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct CalibDeduplicatedId {
    read_id: String,
}


#[derive(Debug, Clone, serde::Serialize)]
pub struct CheckRecord {
    file_name: String,
    total_reads: usize,
    umi_in_read: usize,
    umi_in_read_percent: f32
}

#[allow(dead_code)]
#[derive(Debug, Clone, serde::Deserialize)]
pub struct CalibClusterReportRecord {
   cluster_id: u64,
   cluster_node: u64,
   read_idx: u64,
   read_id_fwd: String,
   read_seq_fwd: String,
   read_qual_fwd: String,
   read_id_rev: String,
   read_seq_rev: String,
   read_qual_rev: String
}
#[allow(dead_code)]
// Assumes same identifier of forward and reverse reads
#[derive(Debug, Clone, serde::Deserialize)]
pub struct CalibRecord {
   cluster: u64,
   read_id: String,
   mean_quality_forward: f32,
   mean_quality_reverse: f32,
   mean_quality_paired: f32,
}
impl CalibRecord {
    pub fn from(cluster_record: &CalibClusterReportRecord) -> Self {

        let mean_q_fwd: f32 = mean_qual(cluster_record.read_qual_fwd.as_bytes());
        let mean_q_rev = mean_qual(cluster_record.read_qual_rev.as_bytes());
        let mean_q_all = (mean_q_fwd+mean_q_rev) / 2.;

        Self {
            cluster: cluster_record.cluster_id,
            read_id: cluster_record.read_id_fwd[1..].to_string(),  // uses fwd id, cluster ids retain @ from identifier, needs to be stripped
            mean_quality_forward: mean_q_fwd,
            mean_quality_reverse: mean_q_rev,
            mean_quality_paired: mean_q_all
        }
    }
}

fn mean_qual(quality_bytes: &[u8]) -> f32 {
    let mut sum: f32 = 0.0;
    for q in quality_bytes.iter() {
        sum += 10f32.powf((q - 33u8) as f32 / -10f32)
    }
    let mean_error_probability = sum / quality_bytes.len() as f32;
    -10f32 * mean_error_probability.log(10.0)
}