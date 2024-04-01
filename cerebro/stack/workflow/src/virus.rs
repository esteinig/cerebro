// Viral summary tables for enrichment runs

use std::path::PathBuf;
use serde::Deserialize;
use needletail::parse_fastx_file;
use needletail::Sequence;

use crate::record::VircovRecord;
use crate::record::VircovScanRemapRecord;
use crate::utils::get_seq_record_identifier;
use crate::tools::Vircov;
use crate::error::WorkflowError;

#[derive(Debug, Clone, PartialEq)]
pub struct ConsensusRecord {
    pub id: String,
    pub length: u64,
    pub missing: u64,
    pub completeness: f64
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct CoverageRecord {
    pub reference: String,
    pub startpos : u64,
    pub endpos: u64,
    pub numreads: u64,
    pub covbases: u64,
    pub coverage: f64,
    pub meandepth: f64,
    pub meanbaseq: f64,
    pub meanmapq: f64,
    pub description: String
}

pub struct Annotation {
    pub taxid: Option<String>,
    pub taxon: Option<String>,
    pub segment: Option<String>
}
impl Annotation {
    pub fn from(description: &str, options: &AnnotationOptions) -> Self {

        let sep: String = match &options.field_sep {
            None => return Default::default(),
            Some(sep) => sep.to_string()
        };

        let fields: Vec<&str> = description.split(&sep).map(|field| field.trim()).collect();

        let taxid = match &options.field_taxid {
            Some(field_taxid) => {
                let taxid_fields: Vec<&str> = fields.iter().filter(|field| field.starts_with(field_taxid)).map(|x| *x).collect();
                match taxid_fields.first() {
                    Some(f) => Some(f.trim().trim_start_matches(field_taxid).to_string()),
                    _ => None
                }
            },
            None => None
        };

        let taxon = match &options.field_taxon {
            Some(field_taxon) => {
                let taxon_fields: Vec<&str> = fields.iter().filter(|field| field.starts_with(field_taxon)).map(|x| *x).collect();
                match taxon_fields.first() {
                    Some(f) => Some(f.trim().trim_start_matches(field_taxon).to_string()),
                    _ => None
                }
            },
            None => None
        };

        let segment = match &options.field_segment {
            Some(field_segment) => {
                let segment_fields: Vec<&str> = fields.iter().filter(|field| field.starts_with(field_segment)).map(|x| *x).collect();
                match segment_fields.first() {
                    Some(segment) => Some(segment.trim_start_matches(field_segment).to_string()),
                    None => None
                }
            },
            None => None
        };

        Self {
            taxid, taxon, segment
        }       

    }
}
impl Default for Annotation {
    fn default() -> Self {
        Self {
            taxid: None,
            taxon: None,
            segment: None
        }
    }
}

pub struct AnnotationOptions {
    pub field_sep: Option<String>, 
    pub field_taxid: Option<String>, 
    pub field_taxon: Option<String>, 
    pub field_segment: Option<String>, 
    pub field_segment_na: Option<String>
}
impl AnnotationOptions {
    pub fn virosaurus() -> Self {
        Self {
            field_sep: Some(String::from(";")),
            field_taxid: Some(String::from("taxid=")),
            field_taxon: Some(String::from("usual name=")),
            field_segment: Some(String::from("segment=")),
            field_segment_na: Some(String::from("N/A"))
        }
    }
}
impl Default for AnnotationOptions {
    fn default() -> Self {
        Self {
            field_sep: None,
            field_taxid: None,
            field_taxon: None,
            field_segment: None,
            field_segment_na: None
        }
    }
}

#[derive(Debug, Clone)]
pub struct MatchedRecords {
    scan: VircovRecord,
    remap: Option<VircovRecord>
}

pub struct VirusAlignmentSummary {
    _id: String,
    _db: String,
    records: Vec<VircovScanRemapRecord>
}
impl VirusAlignmentSummary {
    pub fn new(
        id: &str, 
        db: &str,
        scan: &PathBuf, 
        remap: &PathBuf, 
        fasta: &Option<PathBuf>, 
        coverage: &Option<PathBuf>,
        annotation_options: &AnnotationOptions
    )-> Result<Self, WorkflowError> {

        let scan_vircov = Vircov::from(
            id.to_string(), 
            &scan,
            &annotation_options,
            false
        ).map_err(|_| WorkflowError::VircovScanParseFailed)?;

        let remap_vircov = Vircov::from(
            id.to_string(), &remap, &annotation_options, false
        ).map_err(|_| WorkflowError::VircovRemapParseFailed)?;

        // Group the records by reference sequence identifier
        let mut reference_records = Vec::new();
        for scan_record in scan_vircov.records {

            let scan_annotation = Annotation::from(&scan_record.description, &annotation_options);

            // We can have multiple results with the same reference  if the reference was segmented. We therefore extract the 
            // segment description from therecord annotations and append the extracted segment

            let matching_records: Vec<&VircovRecord> = remap_vircov.records.iter().filter(|remap_record| {
                let mut scan_id = scan_record.reference.clone();
                let mut remap_id = remap_record.reference.clone();
                if let Some(segment) =  scan_annotation.segment.clone() {
                    scan_id.push_str(&segment);
                    remap_id.push_str(&segment);
                };
                scan_id == remap_id
            }).collect();

            let matched_records = match matching_records.len() {
                0 => MatchedRecords { scan: scan_record, remap: None },
                1 => MatchedRecords { scan: scan_record, remap: Some(matching_records[0].clone()) },
                _ => {
                    log::error!("Found multiple matching records for remap reference with identifier: {}", &scan_record.reference);
                    return Err(WorkflowError::VircovRemapMatchNotIdentifiable(scan_record.reference.clone()))
                }
            };

            reference_records.push(matched_records);

        }
        
        let consensus_records = match fasta {
            Some(fasta) =>  VirusAlignmentSummary::parse_consensus_sequences(fasta)?,
            None => Vec::new()
        };

        let coverage_records = match coverage {
            Some(table) => VirusAlignmentSummary::parse_coverage_data(table)?,
            None => Vec::new()
        };

        let mut summary_records = Vec::new();
        for matched_records in reference_records {
            
            let scan_annotation = Annotation::from(&matched_records.scan.description, &annotation_options);

            let consensus_record = match &matched_records.remap {
                None => None,
                Some(remap_record) => {
                    
                    let matching_records: Vec<&ConsensusRecord> = consensus_records.iter().filter(|record | {
                        let mut consensus_id = record.id.clone();
                        let mut remap_id = remap_record.reference.clone();
                        if let Some(segment) = scan_annotation.segment.clone() {
                            consensus_id.push_str(&segment);
                            remap_id.push_str(&segment);
                        };
                        consensus_id == remap_id
                    }).collect();
        

                    match matching_records.len() {
                        0 => None,
                        1 => Some(matching_records[0].clone()),
                        _ => {
                            log::error!("Found multiple matching records for consensus reference with identifier: {}", &remap_record.reference);
                            return Err(WorkflowError::VircovConsensusSequenceMatchNotIdentifiable(remap_record.reference.clone()))
                        }
                    }
                }
            };

            let coverage_record = match &matched_records.remap {
                None => None,
                Some(remap_record) => {
                    
                    let matching_records: Vec<&CoverageRecord> = coverage_records.iter().filter(|record | {
                        let mut coverage_id = record.reference.clone();
                        let mut remap_id = remap_record.reference.clone();
                        if let Some(segment) = scan_annotation.segment.clone() {
                            coverage_id.push_str(&segment);
                            remap_id.push_str(&segment);
                        };
                        coverage_id == remap_id
                    }).collect();
        

                    match matching_records.len() {
                        0 => None,
                        1 => Some(matching_records[0].clone()),
                        _ => {
                            log::error!("Found multiple matching records for coverage reference with identifier: {}", &remap_record.reference);
                            return Err(WorkflowError::VircovRemapMatchNotIdentifiable(remap_record.reference.clone()))
                        }
                    }
                }
            };

            summary_records.push(
                VircovScanRemapRecord::from(
                    id.to_string(),
                    db.to_string(), 
                    matched_records.scan, 
                    scan_annotation, 
                    matched_records.remap,
                    consensus_record, 
                    coverage_record
                )
            )
        }

        summary_records.sort_by(|a, b| b.scan_reads.cmp(&a.scan_reads));

        Ok(Self {
            _id: id.to_string(),
            _db: db.to_string(),
            records: summary_records
        })
    }
    pub fn write_summary(&self, output: &PathBuf, header: bool) {

        let mut writer = csv::WriterBuilder::new().delimiter(b'\t').has_headers(header).from_path(&output).unwrap();
        for rec in &self.records {
            writer.serialize(rec).unwrap();
        }

    }
    pub fn parse_coverage_data(data: &PathBuf) -> Result<Vec<CoverageRecord>, WorkflowError> {
        let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(&data).unwrap();
        let mut coverage_records = Vec::new();
        for rec in reader.deserialize() {
            coverage_records.push(rec.map_err(|err| WorkflowError::CsvRecordNotParsed(err))?)
        }
        Ok(coverage_records)
    }
    pub fn parse_consensus_sequences(fasta: &PathBuf) -> Result<Vec<ConsensusRecord>, WorkflowError>{
        let parse_result = parse_fastx_file(fasta).ok();

        // File may be empty, return no records
        let mut reader = match parse_result {
            Some(reader) => reader,
            None => return Ok(Vec::new())
        };

        let mut consensus_records = Vec::new();

        while let Some(record) = reader.next() {
            
            let rec = record?;
            let count = rec.sequence().iter().filter(|x| *x == &b'N').count() as u64;
            let rec_id = get_seq_record_identifier(&rec)?;
            let rec_len = rec.num_bases() as u64;
            let completeness  = 100.0 - (count as f64 / rec_len as f64)*100.0;

            consensus_records.push(ConsensusRecord {
                id: rec_id,
                length: rec_len,
                missing: count,
                completeness
            })
        }
        Ok(consensus_records)
    }
}