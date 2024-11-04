use std::{collections::HashMap, fs::File, io::{BufReader, BufWriter}, path::Path};

use serde::{Deserialize, Serialize};

use crate::{error::WorkflowError, nextflow::pathogen::{BrackenReportRecord, GanonReportRecord, KmcpReportRecord, KrakenReportRecord, MetabuliReportRecord, PathogenOutput, SylphReport, SylphReportRecord}, utils::{read_tsv, write_tsv}};

use super::quality::QualityControl;

fn compute_rpm(reads: u64, input_reads: u64) -> Option<f64> {
    if input_reads == 0 {
        None  // Cannot compute RPM when input_reads is 0
    } else {
        Some((reads as f64 / input_reads as f64) * 1_000_000.0)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenDetectionFilter {
    pub taxids: Option<Vec<String>>,
    pub names: Option<Vec<String>>,
    pub ranks: Option<Vec<String>>
}
impl PathogenDetectionFilter {
    pub fn new(
        taxids: Option<Vec<String>>,
        names: Option<Vec<String>>,
        ranks: Option<Vec<String>>
    ) -> Result<Self, WorkflowError> {
        Ok(Self {
            taxids, names, ranks
        })
    }
    pub fn from_json(path: &Path) -> Result<Self, WorkflowError> {
        let reader = BufReader::new(File::open(path)?);
        let pathogen = serde_json::from_reader(reader)?;
        Ok(pathogen)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenDetection {
    pub id: String,
    pub records: Vec<PathogenDetectionRecord>
}
impl PathogenDetection {
    pub fn from_pathogen(output: &PathogenOutput, quality: &QualityControl) -> Result<Self, WorkflowError> {

        let input_reads = quality.reads.input_reads;
        let mut detection_map: HashMap<String, PathogenDetectionRecord> = HashMap::new();

        // Process KrakenReport
        if let Some(kraken_report) = &output.profile.kraken2 {
            for record in &kraken_report.records {

                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                    .or_insert(PathogenDetectionRecord::default());
                entry.set_kraken(taxid, record, input_reads); 
            }
        }

        // Process MetabuliReport
        if let Some(metabuli_report) = &output.profile.metabuli {
            for record in &metabuli_report.records {

                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                    .or_insert(PathogenDetectionRecord::default());
                entry.set_metabuli(taxid, record, input_reads); 
            }
        }

        // Process BrackenReport
        if let Some(bracken_report) = &output.profile.bracken {
            for record in &bracken_report.records {
                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                .or_insert(PathogenDetectionRecord::default());
                entry.set_bracken(taxid, record, input_reads); 
            }
        }


        // Process KmcpReport
        if let Some(kmcp_report) = &output.profile.kmcp {
            for record in &kmcp_report.records {
                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                .or_insert(PathogenDetectionRecord::default());
                entry.set_kmcp(taxid, record, input_reads); 
            }
        }
        
        // Process SylphReport
        if let Some(sylph_report) = &output.profile.sylph {
            for record in &sylph_report.records {
                let (taxid, rank) = Self::get_sylph_taxinfo(record)?;
                let entry = detection_map.entry(taxid.to_string())
                    .or_insert(PathogenDetectionRecord::default());
                entry.set_sylph(taxid, rank, record, input_reads);
            }
        }

        // Process GanonReadsReport
        if let Some(ganon_reads_report) = &output.profile.ganon_reads {
            for record in &ganon_reads_report.records {
                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.to_string())
                    .or_insert(PathogenDetectionRecord::default());
                entry.set_ganon_sequence(taxid, record, input_reads);
            }
        }


        // Process GanonAbundanceReport
        if let Some(ganon_reads_report) = &output.profile.ganon_reads {
            for record in &ganon_reads_report.records {
                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.to_string())
                    .or_insert(PathogenDetectionRecord::default());
                entry.set_ganon_sequence(taxid, record, input_reads);
            }
        }

        let records = detection_map
            .into_iter()
            .map(|(_, record)| record)
            .collect();

        Ok(Self { id: output.id.to_string(), records })
    }
    fn get_sylph_taxinfo(record: &SylphReportRecord) -> Result<(String, PathogenDetectionRank), WorkflowError> {

        let tax_str = record.clade_name.trim().to_string();
        let tax_str = tax_str.split("|").last();

        let (taxid, taxlevel) = match tax_str {
            Some(tax_str) => { 
                let taxid_str_split = tax_str.split("__").collect::<Vec<_>>();
                if taxid_str_split.len() != 2 {
                    log::warn!("Failed to recover taxinfo from Sylph output");
                    return Err(WorkflowError::SylphTaxInfoRecoveryFailure(record.clade_name.clone()))
                } else {
                    (taxid_str_split[0], taxid_str_split[1])
                }
            },
            None => { 
                log::warn!("Failed to recover taxid from Sylph output");
                return Err(WorkflowError::SylphTaxInfoRecoveryFailure(record.clade_name.clone()))
            }
        };
        
        Ok((taxid.to_string(), PathogenDetectionRank::from_str(taxlevel)))

    }
    pub fn filter_by_taxonomy(
        &self,
        taxids: Option<Vec<String>>, 
        names: Option<Vec<String>>, 
        ranks: Option<Vec<String>>
    ) -> Vec<PathogenDetectionRecord> {

        let ranks = match ranks {
            Some(ranks) => {
                let transformed_ranks = ranks.into_iter().map(|r| PathogenDetectionRank::from_str(&r)).collect::<Vec<_>>();
                log::info!("Transformed ranks: {:?}", transformed_ranks);
                Some(transformed_ranks)
            },
            None => None
        };

        // Filter records by taxid, name, or rank
        self.records.iter().filter(|record| {
            // Check if the taxid matches
            let taxid_match = if let Some(taxids) = &taxids {
                taxids.contains(&record.taxid)
            } else {
                true
            };

            // Check if any name matches (from kraken, bracken, or metabuli)
            let name_match = if let Some(names) = &names {
                let kraken_name_match = record.kraken_sequence_name.as_ref().map_or(false, |name| names.contains(name));
                let bracken_name_match = record.bracken_profile_name.as_ref().map_or(false, |name| names.contains(name));
                let metabuli_name_match = record.metabuli_sequence_name.as_ref().map_or(false, |name| names.contains(name));
                kraken_name_match || bracken_name_match || metabuli_name_match
            } else {
                true
            };

            // Check if any rank matches (from kraken, bracken, or metabuli)
            let rank_match = if let Some(ranks) = &ranks {

                let kraken_rank_match = record.kraken_sequence_rank.as_ref().map_or(false, |rank| ranks.contains(rank));
                let bracken_rank_match = record.bracken_profile_rank.as_ref().map_or(false, |rank| ranks.contains(rank));
                let metabuli_rank_match = record.metabuli_sequence_rank.as_ref().map_or(false, |rank| ranks.contains(rank));
                kraken_rank_match || bracken_rank_match || metabuli_rank_match
            } else {
                true
            };

            // Record should be included if it matches taxid, name, or rank criteria
            taxid_match && name_match && rank_match
        })
        .cloned()
        .collect()
    }

    pub fn write_records(records: &Vec<PathogenDetectionRecord>, path: &Path) -> Result<(), WorkflowError> {
        write_tsv(&records, path, true)
    }
    pub fn read_records(path: &Path) -> Result<Vec<PathogenDetectionRecord>, WorkflowError> {
        read_tsv(path, false, true)
    }
    pub fn to_tsv(&self, path: &Path) -> Result<(), WorkflowError> {
        write_tsv(&self.records, path, true)
    }
    pub fn from_tsv(&self, path: &Path, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self {
            id: id.to_string(),
            records: read_tsv(path, false, true)?
        })
    }
    pub fn to_json(&self, path: &Path) -> Result<(), WorkflowError> {
        let writer = BufWriter::new(File::create(path)?);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }
    pub fn from_json(path: &Path) -> Result<Self, WorkflowError> {
        let reader = BufReader::new(File::open(path)?);
        let pathogen = serde_json::from_reader(reader)?;
        Ok(pathogen)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum PathogenDetectionRank {
    Superkingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    NoRank,
    Other(String)
}
impl PathogenDetectionRank {
    pub fn from_str(rank_str: &str) -> Self {
        match rank_str {
            "d" | "D" | "d__" | "superkingdom"  => Self::Superkingdom,
            "p" | "P" | "p__" | "phylum"  => Self::Phylum,
            "c" | "C" | "c__" | "class"  => Self::Class,
            "o" | "O" | "o__" | "order"  => Self::Order,
            "f" | "F" | "f__" | "family"  => Self::Family,
            "g" | "G" | "g__" | "genus"  => Self::Genus,
            "s" | "S" | "s__" | "species" => Self::Species,
            "r" | "R" | "root" | "no rank" |  ""  => Self::NoRank,
            _ => Self::Other(rank_str.to_string())
        }
    }
}



#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenDetectionRecord {
    // pub id: String,
    // pub rank: String,
    pub taxid: String,
    // pub input_reads: u64,
    // pub qc_reads: u64,
    pub kraken_sequence_name: Option<String>,
    pub kraken_sequence_rank: Option<PathogenDetectionRank>,
    pub kraken_sequence_reads: Option<u64>,
    pub kraken_sequence_rpm: Option<f64>,
    pub kraken_sequence_abundance: Option<f64>,
    pub metabuli_sequence_name: Option<String>,
    pub metabuli_sequence_rank: Option<PathogenDetectionRank>,
    pub metabuli_sequence_reads: Option<u64>,
    pub metabuli_sequence_rpm: Option<f64>,
    pub metabuli_sequence_abundance: Option<f64>,
    pub ganon_sequence_name: Option<String>,
    pub ganon_sequence_rank: Option<PathogenDetectionRank>,
    pub ganon_sequence_reads: Option<u64>,
    pub ganon_sequence_rpm: Option<f64>,
    pub ganon_sequence_abundance: Option<f64>,
    pub sylph_sequence_name: Option<String>,
    pub sylph_sequence_rank: Option<PathogenDetectionRank>,
    pub sylph_sequence_reads: Option<u64>,
    pub sylph_sequence_rpm: Option<f64>,
    pub sylph_sequence_abundance: Option<f64>,
    pub bracken_profile_name: Option<String>,
    pub bracken_profile_rank: Option<PathogenDetectionRank>,
    pub bracken_profile_reads: Option<u64>,
    pub bracken_profile_rpm: Option<f64>,
    pub bracken_profile_abundance: Option<f64>,
    pub kmcp_profile_name: Option<String>,
    pub kmcp_profile_rank: Option<PathogenDetectionRank>,
    pub kmcp_profile_reads: Option<u64>,
    pub kmcp_profile_rpm: Option<f64>,
    pub kmcp_profile_abundance: Option<f64>,
    pub ganon_profile_name: Option<String>,
    pub ganon_profile_rank: Option<PathogenDetectionRank>,
    pub ganon_profile_reads: Option<u64>,
    pub ganon_profile_rpm: Option<f64>,
    pub ganon_profile_abundance: Option<f64>,
    pub sylph_profile_name: Option<String>,
    pub sylph_profile_rank: Option<PathogenDetectionRank>,
    pub sylph_profile_reads: Option<u64>,
    pub sylph_profile_rpm: Option<f64>,
    pub sylph_profile_abundance: Option<f64>,
}
impl Default for PathogenDetectionRecord {
    fn default() -> Self {
        Self {
            taxid: String::from(""),
            kraken_sequence_reads: None,
            kraken_sequence_name: None,
            kraken_sequence_rank: None,
            kraken_sequence_rpm: None,
            kraken_sequence_abundance: None,
            metabuli_sequence_reads: None,
            metabuli_sequence_rpm: None,
            metabuli_sequence_name: None,
            metabuli_sequence_rank: None,
            metabuli_sequence_abundance: None,
            ganon_sequence_reads: None,
            ganon_sequence_rpm: None,
            ganon_sequence_name: None,
            ganon_sequence_rank: None,
            ganon_sequence_abundance: None,
            sylph_sequence_reads: None,
            sylph_sequence_rpm: None,
            sylph_sequence_name: None,
            sylph_sequence_rank: None,
            sylph_sequence_abundance: None,
            bracken_profile_reads: None,
            bracken_profile_rpm: None,
            bracken_profile_name: None,
            bracken_profile_rank: None,
            bracken_profile_abundance: None,
            kmcp_profile_reads: None,
            kmcp_profile_rpm: None,
            kmcp_profile_name: None,
            kmcp_profile_rank: None,
            kmcp_profile_abundance: None,
            ganon_profile_reads: None,
            ganon_profile_rpm: None,
            ganon_profile_name: None,
            ganon_profile_rank: None,
            ganon_profile_abundance: None,
            sylph_profile_reads: None,
            sylph_profile_rpm: None,
            sylph_profile_name: None,
            sylph_profile_rank: None,
            sylph_profile_abundance: None,
        }
    }
}
impl PathogenDetectionRecord {
    pub fn with_default(taxid: String) -> Self {
        Self {
            taxid,
            ..Default::default()
        }
    }
    pub fn set_kraken(&mut self, taxid: String, record: &KrakenReportRecord, input_reads: u64) -> () {
        
        let sequence_abundance = (record.reads as f64 / input_reads as f64)*100.0;

        self.taxid = taxid;
        self.kraken_sequence_reads = Some(record.reads);
        self.kraken_sequence_rpm = compute_rpm(record.reads, input_reads);
        self.kraken_sequence_name = Some(record.taxname.trim().to_string());
        self.kraken_sequence_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
        self.kraken_sequence_abundance = Some(sequence_abundance);
    }
    pub fn set_bracken(&mut self, taxid: String, record: &BrackenReportRecord, input_reads: u64) -> () {

        let profile_abundance = (record.reads as f64 / input_reads as f64)*100.0;

        self.taxid = taxid;
        self.bracken_profile_reads = Some(record.reads);
        self.bracken_profile_name = Some(record.taxname.trim().to_string());
        self.bracken_profile_rpm = compute_rpm(record.reads, input_reads);
        self.bracken_profile_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));

        self.bracken_profile_abundance = Some(profile_abundance);
    }
    pub fn set_metabuli(&mut self, taxid: String, record: &MetabuliReportRecord, input_reads: u64) -> () {

        let sequence_abundance = (record.reads as f64 / input_reads as f64)*100.0;

        self.taxid = taxid;
        self.metabuli_sequence_reads = Some(record.reads);
        self.metabuli_sequence_rpm = compute_rpm(record.reads, input_reads);
        self.metabuli_sequence_name = Some(record.taxname.trim().to_string());
        self.metabuli_sequence_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
        self.metabuli_sequence_abundance = Some(sequence_abundance);
    }
    pub fn set_kmcp(&mut self, taxid: String, record: &KmcpReportRecord, input_reads: u64) -> () {
        
        self.taxid = taxid;

        let estimated_reads = ((record.abundance/100.0)*input_reads as f64) as u64;
        let taxname = record.taxname_lineage.split("|").last();

        self.kmcp_profile_reads = Some(estimated_reads);
        self.kmcp_profile_rpm = compute_rpm(estimated_reads, input_reads);
        self.kmcp_profile_name = taxname.map(|s| s.to_string());
        self.kmcp_profile_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
        self.kmcp_profile_abundance = Some(record.abundance)
    }
    pub fn set_ganon_sequence(&mut self, taxid: String, record: &GanonReportRecord, input_reads: u64) -> () {
        
        self.taxid = taxid;

        self.ganon_sequence_reads = Some(record.cumulative);
        self.ganon_sequence_rpm = compute_rpm(record.cumulative, input_reads);
        self.ganon_sequence_name = Some(record.taxname.trim().to_string());
        self.ganon_sequence_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
        self.ganon_sequence_abundance = Some(record.cumulative_percent);

    }
    pub fn set_ganon_profile(&mut self, taxid: String, record: &GanonReportRecord, input_reads: u64) -> () {
        
        self.taxid = taxid;

        let estimated_reads = ((record.cumulative_percent/100.0)*input_reads as f64) as u64;

        self.ganon_profile_reads = Some(estimated_reads);
        self.ganon_profile_rpm = compute_rpm(estimated_reads, input_reads);
        self.ganon_profile_name = Some(record.taxname.trim().to_string());
        self.ganon_profile_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
        self.ganon_profile_abundance = Some(record.cumulative_percent);

    }
    pub fn set_sylph(&mut self, taxid: String, tax_level: PathogenDetectionRank, record: &SylphReportRecord, input_reads: u64) -> () {
        
        self.taxid = taxid.clone();

        let estimated_sequence_reads = ((record.sequence_abundance/100.0)*input_reads as f64) as u64;

        self.sylph_sequence_reads = Some(estimated_sequence_reads);
        self.sylph_sequence_rpm = compute_rpm(estimated_sequence_reads, input_reads);
        self.sylph_sequence_name = Some(taxid.to_string());
        self.sylph_sequence_rank = Some(tax_level.clone());
        self.sylph_sequence_abundance = Some(record.sequence_abundance);

        let estimated_abundance_reads = ((record.relative_abundance/100.0)*input_reads as f64) as u64;

        self.sylph_profile_reads = Some(estimated_abundance_reads);
        self.sylph_profile_rpm = compute_rpm(estimated_abundance_reads, input_reads);
        self.sylph_profile_name = Some(taxid.to_string());
        self.sylph_profile_rank = Some(tax_level.clone());
        self.sylph_profile_abundance = Some(record.relative_abundance);
    }
}

pub struct PathogenDetectionToolRecord {
    name: String,
    rank: String,
    reads: u64,
    rpm: f64,
}
impl PathogenDetectionToolRecord {
    pub fn new(name: String, rank: String, reads: u64, rpm: f64) -> Self {
        Self { name, rank, reads, rpm }
    }
}