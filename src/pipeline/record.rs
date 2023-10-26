

use anyhow::Result;
use itertools::Itertools;
use taxonomy::GeneralTaxonomy;
use serde::{Deserialize, Serialize, Serializer};
use crate::tools::modules::virus::{ConsensusRecord, Annotation, CoverageRecord, AnnotationOptions};
use super::{error::WorkflowError, module::get_lca_taxid};


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct VircovRecord {
    pub id: String,
    pub db: String,
    pub tool: String,
    pub taxid: String,
    pub reference: String,
    pub regions: u64,
    pub reads: u64,
    pub alignments: u64,
    pub bases: u64,
    pub length: u64,
    pub coverage: f64,
    pub description: String,
    pub tags: String,
    pub rpm: f64
}

impl VircovRecord {
    pub fn from_str(line: String, id: String, db: String, tax_sep: &Option<String>, tax_field: &Option<String>) -> Result<Self, WorkflowError> {
        let fields: Vec<&str> = line.split('\t').collect();
        let descr = fields[7].trim().to_string();
        
        let record = Self {
            id,
            db: db.to_string(),
            tool: "vircov".to_string(),
            taxid: match (tax_sep, tax_field) { 
                (Some(sep), Some(field)) => {
                    let fields: Vec<&str> = descr.split(sep).into_iter().filter(|x| x.trim().starts_with(field)).collect();
                    match fields.first() {
                        Some(f) => f.trim().trim_start_matches(field).to_string(),
                        _ => return Err(WorkflowError::VircovTaxidFieldMissing(db))
                    }
                },
                _ => "".to_string() // not necessary e.g. in ercc qc 
            }, // parse from description
            reference: fields[0].to_string().trim().to_string(),
            regions: fields[1]
                .parse::<u64>()
                .map_err(|_| WorkflowError::VircovRegionFieldConversion)?,
            reads: fields[2]
                .parse::<u64>()
                .map_err(|_| WorkflowError::VircovReadFieldConversion)?,
            alignments: fields[3]
                .parse::<u64>()
                .map_err(|_| WorkflowError::VircovAlignmentFieldConversion)?,
            bases: fields[4]
                .parse::<u64>()
                .map_err(|_| WorkflowError::VircovBasepairFieldConversion)?,
            length: fields[5]
                .parse::<u64>()
                .map_err(|_| WorkflowError::VircovLengthFieldConversion)?,
            coverage: fields[6]
                .parse::<f64>()
                .map_err(|_| WorkflowError::VircovCoverageFieldConversion)?,
            description: descr,
            tags: fields[8].trim().to_string(),
            rpm: 0.0
        };
        Ok(record)
    }
}


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct VircovScanRemapRecord {
    pub id: String,
    pub db: String,
    pub tool: String,
    pub reference: String,
    pub reference_length: u64,
    pub scan_regions: u64,
    pub scan_reads: u64,
    pub scan_alignments: u64,
    pub scan_bases_covered: u64,
    pub scan_coverage: f64,
    pub remap_regions: Option<u64>,
    pub remap_reads: Option<u64>,
    pub remap_alignments: Option<u64>,
    pub remap_bases_covered: Option<u64>,
    pub remap_coverage: Option<f64>,
    pub remap_mean_depth: Option<f64>,
    pub consensus_length: Option<u64>,
    pub consensus_missing: Option<u64>,
    pub consensus_completeness: Option<f64>,
    pub taxid: Option<String>,
    pub name: Option<String>,
    pub segment: Option<String>,
    pub reference_description: String,
    pub scan_rpm: f64,
    pub remap_rpm: f64,
}

impl VircovScanRemapRecord {
    pub fn from(id: String, db: String, scan_record: VircovRecord, scan_annotation: Annotation, remap_record: Option<VircovRecord>, consensus_record: Option<ConsensusRecord>, coverage_record: Option<CoverageRecord>) -> Self {

        let (
            remap_regions,
            remap_reads,
            remap_alignments,
            remap_bases_covered,
            remap_coverage
        ) = match remap_record { 
            Some(r) =>(
                Some(r.regions),
                Some(r.reads), 
                Some(r.alignments),
                Some(r.bases),
                Some(r.coverage*100.)
            ),
            None => (None, None, None, None, None)
        };

        let (
            consensus_length, consensus_missing, consensus_completeness
        ) = match consensus_record {
            Some(r) => (Some(r.length), Some(r.missing), Some(r.completeness)),
            None => (None, None, None)
        };

        let remap_mean_depth = match coverage_record {
            Some(r) => Some(r.meandepth),
            None => None
        };

        Self {
            id,
            db,
            tool: "vircov_scan_remap".to_string(),
            reference: scan_record.reference,
            reference_length: scan_record.length,
            scan_regions: scan_record.regions,
            scan_reads: scan_record.reads,
            scan_alignments: scan_record.alignments,
            scan_bases_covered: scan_record.bases,
            scan_coverage: scan_record.coverage*100.,
            remap_regions,
            remap_reads,
            remap_alignments,
            remap_bases_covered,
            remap_coverage,
            remap_mean_depth,
            consensus_length,
            consensus_missing,
            consensus_completeness,
            taxid: scan_annotation.taxid,
            name: scan_annotation.taxon,
            segment: scan_annotation.segment,
            reference_description: scan_record.description,
            scan_rpm: 0.0,
            remap_rpm: 0.0
        }

    }    
    pub fn from_str(line: &str, annotation_options: &AnnotationOptions) -> Result<Self, WorkflowError> {

        let fields: Vec<&str> = line.split('\t').collect(); 
        let reference_description = fields[22].to_string();

        // Harmonized annotation options for aggregation upstream
        let annotation = Annotation::from(&reference_description, &annotation_options);
        
        let record = Self {
            id: fields[0].trim().to_string(),
            db: fields[1].trim().to_string(),
            tool: fields[2].trim().to_string(),
           
            reference: fields[3].trim().to_string(),
            reference_length: fields[4].parse::<u64>().map_err(|_| WorkflowError::VircovRegionFieldConversion)?,
            scan_regions: fields[5].parse::<u64>().map_err(|_| WorkflowError::VircovRegionFieldConversion)?,
            scan_reads: fields[6].parse::<u64>().map_err(|_| WorkflowError::VircovRegionFieldConversion)?,
            scan_alignments: fields[7].parse::<u64>().map_err(|_| WorkflowError::VircovRegionFieldConversion)?,
            scan_bases_covered: fields[8].parse::<u64>().map_err(|_| WorkflowError::VircovRegionFieldConversion)?,
            scan_coverage: fields[9].parse::<f64>().map_err(|_| WorkflowError::VircovCoverageFieldConversion)?,
            remap_regions: fields[10].parse::<u64>().ok(),
            remap_reads: fields[11].parse::<u64>().ok(),
            remap_alignments: fields[12].parse::<u64>().ok(),
            remap_bases_covered: fields[13].parse::<u64>().ok(),
            remap_coverage: fields[14].parse::<f64>().ok(),
            remap_mean_depth: fields[15].parse::<f64>().ok(),
            consensus_length: fields[16].parse::<u64>().ok(),
            consensus_missing: fields[17].parse::<u64>().ok(),
            consensus_completeness: fields[18].parse::<f64>().ok(),
            taxid: annotation.taxid,
            name: annotation.taxon,
            segment: annotation.segment,
            reference_description,
            scan_rpm: 0.0,
            remap_rpm: 0.0
        };
        Ok(record)
    }   
}


// Kraken2Uniq report record
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Kraken2UniqRecord {
    pub id: String,
    pub db: String,
    pub tool: String,
    pub percent: String,
    pub reads: u64,
    pub reads_direct: u64,
    pub kmers: u64,
    pub kmers_unique: u64,
    pub tax_level: String,
    pub taxid: String,
    pub taxname: String,
    pub rpm: f64
}

impl Kraken2UniqRecord {
    // Create a record from a parsed line
    pub fn from_str(line: String, id: String, db: String) -> Result<Self, WorkflowError> {

        let fields: Vec<&str> = line.split('\t').collect();

        let record = Self {
            id,
            db,
            tool: "kraken2uniq".to_string(),
            percent: fields[0].to_string().trim().to_string(),
            reads: fields[1]
                .parse::<u64>()
                .map_err(|_| WorkflowError::KrakenReportReadFieldConversion)?,
            reads_direct: fields[2]
                .parse::<u64>()
                .map_err(|_| WorkflowError::KrakenReportDirectReadFieldConversion)?,
            kmers: fields[3]
                .parse::<u64>()
                .map_err(|_| WorkflowError::KrakenReportKmerFieldConversion)?,
            kmers_unique: fields[4]
                .parse::<u64>()
                .map_err(|_| WorkflowError::KrakenReportKmerUniqueFieldConversion)?,
            tax_level: fields[5].trim().to_string(),
            taxid: fields[6].trim().to_string(),
            taxname: fields[7].trim().to_string(),
            rpm: 0.0
        };

        Ok(record)
    }
    pub fn as_row_field(&self) -> String {
        format!("{}: {}::{}::{}::{:.2}", self.db, self.reads, self.kmers, self.kmers_unique, self.kmers as f64 / self.kmers_unique as f64)        
    }
}


// Highest Bitscore of LCA identified taxa
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BlastLcaRecord {
    pub id: String,
    pub db: String,
    pub tool: String,
    pub length: u64,    
    pub alignment: u64,
    pub coverage: f64,         
    pub identity: f64,
    pub evalue: f64,
    pub bitscore: f64,
    pub taxid: String,
    pub title: String,
    pub reference: String,
    pub reference_length: u64,
    pub read_coverage: f64,           // if present in contig name e.g. from spades output
    pub bpm: f64,
}
impl BlastLcaRecord {
    pub fn from_blast_record(id: &str, db: &str, records: Vec<BlastRecord>, taxonomy: &GeneralTaxonomy) -> Result<Self, WorkflowError> {
        
        let (record_length, record_name) = match records.get(0) {
            Some(record) => (record.qlen, record.qid.clone()),
            None => return Err(WorkflowError::BlastLcaRecordExtraction)
        };

        // Extraction of coverage from MetaSpades assembly output
        let record_name_components: Vec<&str> = record_name.split("_").collect();
        let read_coverage = match record_name_components.contains(&"cov") {
            true => match record_name_components.last(){
                Some(cov) => cov.parse::<f64>().map_err(WorkflowError::BlastLcaFloatFieldConversion)?,
                None => return Err(WorkflowError::BlastLcaCovExtraction)
            },
            false => 0.0
        };
        
        // Find the LCA of all BlastRecords for this contig
        let taxids_unique: Vec<String> = records.iter().map(|x| x.taxid.to_string()).unique().collect();

        let (failed_taxids, lca_taxid) = get_lca_taxid(taxids_unique.clone(), taxonomy)?;

        if !failed_taxids.is_empty() {
            let proportion_failed = failed_taxids.len() as f64 / taxids_unique.len() as f64;
            if proportion_failed >= 0.5 {
                log::warn!("Failed to detect a significant proportion of unique taxids from all hits for contig {}: {:?} ({:.2} %) in the reference taxonomy. LCA classification for this contig is ignored.", &record_name, &failed_taxids, proportion_failed*100.0);
                return Err(WorkflowError::BlastLcaTaxidsNotFound)
            }
        }

        let highest_bitscore = match records.into_iter().max_by(|a, b| a.bitscore.total_cmp(&b.bitscore)){
            Some(record) => record,
            None => return Err(WorkflowError::BlastLcaRecordExtractionEvalue)
        };
        
        Ok(Self { 
            id: id.to_string(),
            db: db.to_string(), 
            tool: "blastn".to_string(),
            length: record_length,
            alignment: highest_bitscore.length,
            coverage: (highest_bitscore.length as f64/highest_bitscore.qlen as f64)*100.0,
            identity: highest_bitscore.pident,
            evalue: highest_bitscore.evalue,
            bitscore: highest_bitscore.bitscore,
            taxid: lca_taxid.to_string(),
            read_coverage,
            reference: highest_bitscore.sid,
            reference_length: highest_bitscore.slen,
            title: highest_bitscore.title,
            bpm: 0.
        } )
    }
    pub fn from_diamond_record(id: &str, db: &str, records: Vec<DiamondRecord>, taxonomy: &GeneralTaxonomy) -> Result<Self, WorkflowError> {
        let (record_length, record_name) = match records.get(0) {
            Some(record) => (record.qlen, record.qid.clone()),
            None => return Err(WorkflowError::BlastLcaRecordExtraction)
        };

        // Extraction of coverage from MetaSpades assembly output
        let record_name_components: Vec<&str> = record_name.split("_").collect();
        let read_coverage = match record_name_components.contains(&"cov") {
            true => match record_name_components.last(){
                Some(cov) => cov.parse::<f64>().map_err(WorkflowError::BlastLcaFloatFieldConversion)?,
                None => return Err(WorkflowError::BlastLcaCovExtraction)
            },
            false => 0.0
        };
        
        let mut taxids: Vec<String> = Vec::new();
        for record in &records {
            for taxid in &record.taxids {
                taxids.push(taxid.to_string())
            }
        };
        let taxids_unique: Vec<String> = taxids.iter().unique().map(|x| x.to_string()).collect();

        // Find the LCA of all hits against this contig
        let (failed_taxids, lca_taxid) = get_lca_taxid(taxids_unique.clone(), taxonomy)?;
       
        if !failed_taxids.is_empty() {
            let proportion_failed = failed_taxids.len() as f64 / taxids_unique.len() as f64;
            if proportion_failed >= 0.5 {
                log::warn!("Failed to detect a significant proportion (>= 50%) of unique taxids from all hits for contig {}: {:?} ({:.2} %) in the reference taxonomy. LCA classification for this contig is ignored.", &record_name, &failed_taxids, proportion_failed*100.0);
                return Err(WorkflowError::BlastLcaTaxidsNotFound)
            }
        }

        let highest_bitscore = match records.into_iter().max_by(|a, b| a.bitscore.total_cmp(&b.bitscore)){
            Some(record) => record,
            None => return Err(WorkflowError::BlastLcaRecordExtractionEvalue)
        };
        
        Ok(Self { 
            id: id.to_string(),
            db: db.to_string(), 
            tool: "diamond".to_string(),
            length: record_length,
            alignment: highest_bitscore.length,
            coverage: (highest_bitscore.length as f64/highest_bitscore.qlen as f64)*100.0,
            identity: highest_bitscore.pident,
            evalue: highest_bitscore.evalue,
            bitscore: highest_bitscore.bitscore,
            taxid: lca_taxid.to_string(),
            read_coverage,
            reference: highest_bitscore.sid,
            reference_length: highest_bitscore.slen,
            title: highest_bitscore.title,
            bpm: 0.
        } )
    }
    pub fn as_row_field(&self) -> String {
        format!("{}: {}::{}::{:.2}::{:.2}::{:.2}", self.db, self.length, self.alignment, self.coverage, self.identity, self.read_coverage)        
    }
    
}

// qseqid qlen qstart qend sseqid slen sstart send length nident pident evalue bitscore staxid ssciname stitle
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlastRecord {
    pub id: String,
    pub db: String,
    pub qid: String,
    pub qlen: u64,
    pub qstart: u64,
    pub qend: u64,
    pub sid: String,
    pub slen: u64,
    pub sstart: u64,
    pub send: u64,
    pub length: u64,
    pub nident: u64,
    pub pident: f64,
    pub evalue: f64,
    pub bitscore: f64,
    pub taxid: String,
    pub taxname: String,
    pub title: String
}

impl BlastRecord {
    pub fn from_str(line: String, id: &str, db: &str) -> Result<Self, WorkflowError> {
        let fields: Vec<&str> = line.split('\t').collect();
        
        let record = Self {
            id: id.to_string(),
            db: db.to_string(),
            qid: fields[0].to_string(),
            qlen: fields[1]
                .parse::<u64>()
                .map_err(WorkflowError::BlastQueryLengthFieldConversion)?,
            qstart: fields[2]
                .parse::<u64>()
                .map_err(WorkflowError::BlastQueryStartFieldConversion)?,
            qend: fields[3]
                .parse::<u64>()
                .map_err(WorkflowError::BlastQueryEndFieldConversion)?,
            sid: fields[4].to_string(),
            slen: fields[5]
                .parse::<u64>()
                .map_err(WorkflowError::BlastSubjectLengthFieldConversion)?,
            sstart: fields[6]
                .parse::<u64>()
                .map_err( WorkflowError::BlastSubjectStartFieldConversion)?,
            send: fields[7]
                .parse::<u64>()
                .map_err(WorkflowError::BlastSubjectEndFieldConversion)?,
            length: fields[8]
                .parse::<u64>()
                .map_err(WorkflowError::BlastLengthFieldConversion)?,
            nident: fields[9]
                .parse::<u64>()
                .map_err(WorkflowError::BlastIdentityNumberFieldConversion)?,
            pident: fields[10]
                .parse::<f64>()
                .map_err(WorkflowError::BlastIdentityPercentFieldConversion)?,
            evalue: fields[11]
                .parse::<f64>()
                .map_err(WorkflowError::BlastEvalueFieldConversion)?,
            bitscore: fields[12]
                .parse::<f64>()
                .map_err(WorkflowError::BlastBitscoreFieldConversion)?,
            taxid: fields[13].to_string(),
            taxname: fields[14].to_string(),
            title: fields[15].to_string()
        };
        Ok(record)
    }
}



// qseqid qlen qstart qend sseqid slen sstart send length nident pident evalue bitscore staxids sscinames stitle 
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct DiamondRecord {
    pub id: String,
    pub db: String,
    pub qid: String,
    pub qlen: u64,
    pub qstart: u64,
    pub qend: u64,
    pub sid: String,
    pub slen: u64,
    pub sstart: u64,
    pub send: u64,
    pub length: u64,
    pub nident: u64,
    pub pident: f64,
    pub evalue: f64,
    pub bitscore: f64,
    pub taxids: Vec<String>,
    pub taxnames: Vec<String>,
    pub title: String
}

impl DiamondRecord {
    pub fn from_str(line: String, id: &str, db: &str) -> Result<Self, WorkflowError> {
        let fields: Vec<&str> = line.split('\t').collect();
        
        let record = Self {
            id: id.to_string(),
            db: db.to_string(),
            qid: fields[0].to_string(),
            qlen: fields[1]
                .parse::<u64>()
                .map_err(WorkflowError::BlastQueryLengthFieldConversion)?,
            qstart: fields[2]
                .parse::<u64>()
                .map_err(WorkflowError::BlastQueryStartFieldConversion)?,
            qend: fields[3]
                .parse::<u64>()
                .map_err(WorkflowError::BlastQueryEndFieldConversion)?,
            sid: fields[4].to_string(),
            slen: fields[5]
                .parse::<u64>()
                .map_err(WorkflowError::BlastSubjectLengthFieldConversion)?,
            sstart: fields[6]
                .parse::<u64>()
                .map_err( WorkflowError::BlastSubjectStartFieldConversion)?,
            send: fields[7]
                .parse::<u64>()
                .map_err(WorkflowError::BlastSubjectEndFieldConversion)?,
            length: fields[8]
                .parse::<u64>()
                .map_err(WorkflowError::BlastLengthFieldConversion)?,
            nident: fields[9]
                .parse::<u64>()
                .map_err(WorkflowError::BlastIdentityNumberFieldConversion)?,
            pident: fields[10]
                .parse::<f64>()
                .map_err(WorkflowError::BlastIdentityPercentFieldConversion)?,
            evalue: fields[11]
                .parse::<f64>()
                .map_err(WorkflowError::BlastEvalueFieldConversion)?,
            bitscore: fields[12]
                .parse::<f64>()
                .map_err(WorkflowError::BlastBitscoreFieldConversion)?,
            taxids: fields[13].to_string().split(";").map(String::from).collect(),
            taxnames: fields[13].to_string().split(";").map(String::from).collect(),
            title: fields[15].to_string()
        };
        Ok(record)
    }
}
 