use std::{collections::HashMap, fs::File, io::{BufReader, BufWriter}, path::{Path, PathBuf}};

use serde::{Deserialize, Serialize};
use taxonomy::{ncbi, GeneralTaxonomy, TaxRank, Taxonomy};

use crate::{error::WorkflowError, nextflow::pathogen::{BrackenReportRecord, GanonReportRecord, KmcpAbundanceReportRecord, KmcpReadsReportRecord, KrakenReportRecord, MetabuliReportRecord, PathogenOutput, SylphReportRecord}, utils::{read_tsv, write_tsv}};

use super::quality::QualityControl;


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenDetectionTableRecord {
    id: String,
    taxid: String,
    rank: Option<TaxRank>,
    name: Option<String>,
    lineage: Option<String>,
    kraken_reads: Option<u64>,
    kraken_rpm: Option<f64>,
    bracken_reads: Option<u64>,
    bracken_rpm: Option<f64>,
    metabuli_reads: Option<u64>,
    metabuli_rpm: Option<f64>,
    ganon_reads: Option<u64>,
    ganon_rpm: Option<f64>,
    kmcp_reads: Option<u64>,
    kmcp_rpm: Option<f64>,
    sylph_reads: Option<u64>,
    sylph_rpm: Option<f64>,
}
impl PathogenDetectionTableRecord {
    pub fn from_record(id: &str, record: &PathogenDetectionRecord, taxonomy: Option<&GeneralTaxonomy>) -> Result<Self, WorkflowError> {
        let taxid = record.taxid.as_str();

        let (rank, name, lineage) = match taxonomy {
            Some(taxonomy) => {
                let rank = taxonomy.rank(record.taxid.as_str())?;
                let name = taxonomy.name(record.taxid.as_str())?.to_string();
                
                let lineage = vec![
                    taxonomy.name(taxonomy.parent_at_rank(taxid, TaxRank::Superkingdom)?.unwrap_or(("", 0.0)).0).unwrap_or(""),
                    taxonomy.name(taxonomy.parent_at_rank(taxid, TaxRank::Phylum)?.unwrap_or(("", 0.0)).0).unwrap_or(""),
                    taxonomy.name(taxonomy.parent_at_rank(taxid, TaxRank::Class)?.unwrap_or(("", 0.0)).0).unwrap_or(""),
                    taxonomy.name(taxonomy.parent_at_rank(taxid, TaxRank::Order)?.unwrap_or(("", 0.0)).0).unwrap_or(""),
                    taxonomy.name(taxonomy.parent_at_rank(taxid, TaxRank::Family)?.unwrap_or(("", 0.0)).0).unwrap_or(""),
                    taxonomy.name(taxonomy.parent_at_rank(taxid, TaxRank::Genus)?.unwrap_or(("", 0.0)).0).unwrap_or(""),
                    taxonomy.name(taxonomy.parent_at_rank(taxid, TaxRank::Species)?.unwrap_or(("", 0.0)).0).unwrap_or(""),
                ].into_iter().map(String::from).collect::<Vec<_>>();

                (Some(rank), Some(name), Some(Self::lineage_to_str(lineage)?))
            },
            None => (None, None, None)
        };

        Ok(Self {
            id: id.to_string(),
            taxid: record.taxid.clone(),
            rank,
            name,
            lineage,
            kraken_reads: record.kraken_sequence_reads,
            kraken_rpm: record.kraken_sequence_rpm,
            bracken_reads: record.bracken_profile_reads,
            bracken_rpm: record.bracken_profile_rpm,
            metabuli_reads: record.metabuli_sequence_reads,
            metabuli_rpm: record.metabuli_sequence_rpm,
            ganon_reads: record.ganon_sequence_reads,
            ganon_rpm: record.ganon_sequence_rpm,
            kmcp_reads: record.kmcp_sequence_reads,
            kmcp_rpm: record.kmcp_sequence_rpm,
            sylph_reads: record.sylph_sequence_reads,
            sylph_rpm: record.sylph_sequence_rpm
        })
    }

    pub fn lineage_to_str(lineage: Vec<String>) -> Result<String, WorkflowError> {
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
}

// Write the evidence table for multiple samples
pub fn write_pathogen_table(json: &Vec<PathBuf>, path: &Path, taxonomy_directory: Option<PathBuf>, filter_json: Option<PathBuf>) -> Result<(), WorkflowError> {

    let taxonomy = match taxonomy_directory {
        Some(dir) => {
            let tax = ncbi::load(&dir)?;
            Some(tax)
        },
        None => None
    };

    let pd_filter = match filter_json {
        Some(path) => Some(PathogenDetectionFilter::from_json(&path)?),
        None => None
    };

    let mut table_records = Vec::new();
    for file in json {
        let pd = PathogenDetection::from_json(file)?;
        let records = match pd_filter {
            Some(ref f) => pd.filter_by_taxonomy(f.taxids.clone(), f.names.clone(), f.ranks.clone()),
            None => pd.records,
        };
        for record in records {
            table_records.push(
                PathogenDetectionTableRecord::from_record(
                    &pd.id, 
                    &record, 
                    taxonomy.as_ref()
                )
            ?)
        }
    }

    write_tsv(&table_records, path, true)?;

    Ok(())
}

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
    pub paired_end: bool,
    pub records: Vec<PathogenDetectionRecord>
}
impl PathogenDetection {
    pub fn from_pathogen(output: &PathogenOutput, quality: &QualityControl, paired_end: bool) -> Result<Self, WorkflowError> {

        let input_reads = quality.reads.input_reads;
        let classifier_reads = quality.reads.output_reads;

        let mut detection_map: HashMap<String, PathogenDetectionRecord> = HashMap::new();

        // Process KrakenReport
        if let Some(kraken_report) = &output.profile.kraken2 {
            for record in &kraken_report.records {

                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                    .or_insert(PathogenDetectionRecord::default());
                entry.set_kraken(taxid, record, input_reads, classifier_reads, paired_end); 
            }
        }

        // Process MetabuliReport
        if let Some(metabuli_report) = &output.profile.metabuli {
            for record in &metabuli_report.records {

                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                    .or_insert(PathogenDetectionRecord::default());
                entry.set_metabuli(taxid, record, input_reads, classifier_reads, paired_end); 
            }
        }

        // Process BrackenReport
        if let Some(bracken_report) = &output.profile.bracken {
            for record in &bracken_report.records {
                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                .or_insert(PathogenDetectionRecord::default());
                entry.set_bracken(taxid, record, input_reads, classifier_reads, paired_end); 
            }
        }


        // Process KmcpAbundanceReport
        if let Some(kmcp_report) = &output.profile.kmcp_abundance {
            for record in &kmcp_report.records {
                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                .or_insert(PathogenDetectionRecord::default());
                entry.set_kmcp_profile(taxid, record); 
            }
        }


        // Process KmcpReadsReport
        if let Some(kmcp_report) = &output.profile.kmcp_reads {

            // Group records by taxonomic identifier and sum the read counts
            let kmcp_taxid_report = kmcp_report.get_taxid_report()?;

            for record in &kmcp_taxid_report.records {
                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.clone())
                .or_insert(PathogenDetectionRecord::default());
                entry.set_kmcp_sequence(taxid, record, input_reads, classifier_reads, paired_end); 
            }
        }
        
        // Process SylphReport
        if let Some(sylph_report) = &output.profile.sylph {
            for record in &sylph_report.records {
                let (taxid, rank) = Self::get_sylph_taxinfo(record)?;

                if rank == PathogenDetectionRank::Strain {
                    continue;
                }

                let entry = detection_map.entry(taxid.to_string())
                    .or_insert(PathogenDetectionRecord::default());
                entry.set_sylph(taxid, rank, record, input_reads, classifier_reads);
            }
        }

        // Process GanonReadsReport
        if let Some(ganon_reads_report) = &output.profile.ganon_reads {
            for record in &ganon_reads_report.records {
                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.to_string())
                    .or_insert(PathogenDetectionRecord::default());
                entry.set_ganon_sequence(taxid, record, input_reads, classifier_reads, paired_end);
            }
        }

        // Process GanonAbundanceReport
        if let Some(ganon_profile_report) = &output.profile.ganon_abundance {
            for record in &ganon_profile_report.records {
                let taxid = record.taxid.trim().to_string();
                let entry = detection_map.entry(taxid.to_string())
                    .or_insert(PathogenDetectionRecord::default());
                entry.set_ganon_profile(taxid, record);
            }
        }

        let records = detection_map
            .into_iter()
            .map(|(_, record)| record)
            .collect();

        Ok(Self { id: output.id.to_string(), paired_end, records })
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
                    (taxid_str_split[1], taxid_str_split[0])
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
    pub fn from_tsv(&self, path: &Path, id: &str, paired_end: bool) -> Result<Self, WorkflowError> {
        Ok(Self {
            id: id.to_string(),
            paired_end,
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
    Strain,
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
            "t" | "T" | "t__" | "strain" => Self::Strain,
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
    pub kmcp_sequence_name: Option<String>,
    pub kmcp_sequence_rank: Option<PathogenDetectionRank>,
    pub kmcp_sequence_reads: Option<u64>,
    pub kmcp_sequence_rpm: Option<f64>,
    pub kmcp_sequence_abundance: Option<f64>,
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
            kmcp_sequence_reads: None,
            kmcp_sequence_rpm: None,
            kmcp_sequence_name: None,
            kmcp_sequence_rank: None,
            kmcp_sequence_abundance: None,
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
    pub fn set_kraken(&mut self, taxid: String, record: &KrakenReportRecord, input_reads: u64, classifier_reads: u64, paired_end: bool) -> () {
        
        let reads = if paired_end { record.reads*2 } else { record.reads };

        let sequence_abundance = (reads as f64 / classifier_reads as f64)*100.0;

        self.taxid = taxid;
        self.kraken_sequence_reads = Some(reads);
        self.kraken_sequence_rpm = compute_rpm(reads, input_reads);
        self.kraken_sequence_name = Some(record.taxname.trim().to_string());
        self.kraken_sequence_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
        self.kraken_sequence_abundance = Some(sequence_abundance);
    }
    pub fn set_bracken(&mut self, taxid: String, record: &BrackenReportRecord, input_reads: u64, classifier_reads: u64, paired_end: bool) -> () {

        let reads = if paired_end { record.reads*2 } else { record.reads };

        let profile_abundance = (reads as f64 / classifier_reads as f64)*100.0;

        self.taxid = taxid;
        self.bracken_profile_reads = Some(reads);
        self.bracken_profile_name = Some(record.taxname.trim().to_string());
        self.bracken_profile_rpm = compute_rpm(reads, input_reads);
        self.bracken_profile_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));

        self.bracken_profile_abundance = Some(profile_abundance);
    }
    pub fn set_metabuli(&mut self, taxid: String, record: &MetabuliReportRecord, input_reads: u64, classifier_reads: u64, paired_end: bool) -> () {

        let reads = if paired_end { record.reads*2 } else { record.reads };

        let sequence_abundance = (reads as f64 / classifier_reads as f64)*100.0;

        self.taxid = taxid;
        self.metabuli_sequence_reads = Some(reads);
        self.metabuli_sequence_rpm = compute_rpm(reads, input_reads);
        self.metabuli_sequence_name = Some(record.taxname.trim().to_string());
        self.metabuli_sequence_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
        self.metabuli_sequence_abundance = Some(sequence_abundance);
    }
    pub fn set_kmcp_sequence(&mut self, taxid: String, record: &KmcpReadsReportRecord, input_reads: u64, classifier_reads: u64, paired_end: bool) -> () {
        
        let reads = if paired_end { record.reads*2 } else { record.reads };
        let sequence_abundance = (reads as f64 / classifier_reads as f64)*100.0;

        self.taxid = taxid;

        self.kmcp_sequence_reads = Some(reads);
        self.kmcp_sequence_rpm = compute_rpm(reads, input_reads);
        self.kmcp_sequence_name = Some(record.taxname.clone());
        self.kmcp_sequence_rank = Some(PathogenDetectionRank::from_str(&record.rank));
        self.kmcp_sequence_abundance = Some(sequence_abundance)
    }
    pub fn set_kmcp_profile(&mut self, taxid: String, record: &KmcpAbundanceReportRecord) -> () {
        
        self.taxid = taxid;

        let taxname = record.taxname_lineage.split("|").last();

        self.kmcp_profile_reads = None;
        self.kmcp_profile_rpm = None;
        self.kmcp_profile_name = taxname.map(|s| s.to_string());
        self.kmcp_profile_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
        self.kmcp_profile_abundance = Some(record.abundance)
    }
    pub fn set_ganon_sequence(&mut self, taxid: String, record: &GanonReportRecord, input_reads: u64, classifier_reads: u64, paired_end: bool) -> () {
        
        let reads = if paired_end { record.cumulative*2 } else { record.cumulative };
        let sequence_abundance = (reads as f64 / classifier_reads as f64)*100.0;

        self.taxid = taxid;

        self.ganon_sequence_reads = Some(reads);
        self.ganon_sequence_rpm = compute_rpm(reads, input_reads);
        self.ganon_sequence_name = Some(record.taxname.trim().to_string());
        self.ganon_sequence_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
        self.ganon_sequence_abundance = Some(sequence_abundance);

    }
    pub fn set_ganon_profile(&mut self, taxid: String, record: &GanonReportRecord) -> () {
        
        self.taxid = taxid;

        // let estimated_reads = ((record.cumulative_percent/100.0)*classifier_reads as f64) as u64;

        self.ganon_profile_reads = None;
        self.ganon_profile_rpm = None;
        self.ganon_profile_name = Some(record.taxname.trim().to_string());
        self.ganon_profile_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
        self.ganon_profile_abundance = Some(record.cumulative_percent);

    }
    pub fn set_sylph(&mut self, taxid: String, tax_level: PathogenDetectionRank, record: &SylphReportRecord, input_reads: u64, classifier_reads: u64) -> () {
        
        self.taxid = taxid.clone();

        let estimated_sequence_reads = ((record.sequence_abundance/100.0)*classifier_reads as f64) as u64;

        self.sylph_sequence_reads = Some(estimated_sequence_reads);
        self.sylph_sequence_rpm = compute_rpm(estimated_sequence_reads, input_reads);
        self.sylph_sequence_name = Some(taxid.to_string());
        self.sylph_sequence_rank = Some(tax_level.clone());
        self.sylph_sequence_abundance = Some(record.sequence_abundance);

        // let estimated_abundance_reads = ((record.relative_abundance/100.0)*classifier_reads as f64) as u64;

        self.sylph_profile_reads = None;
        self.sylph_profile_rpm = None;
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