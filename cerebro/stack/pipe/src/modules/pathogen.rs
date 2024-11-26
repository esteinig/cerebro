use std::{collections::HashMap, fs::File, io::{BufReader, BufWriter}, path::{Path, PathBuf}};

use serde::{Deserialize, Serialize};
use taxonomy::{ncbi, GeneralTaxonomy, TaxRank, Taxonomy};

use crate::{error::WorkflowError, nextflow::pathogen::{PathogenOutput, SylphReportRecord}, taxa::taxon::Taxon, utils::{read_tsv, write_tsv}};

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
    pub fn from_record(
        id: &str,
        record: &PathogenDetectionRecord,
        taxonomy: Option<&GeneralTaxonomy>,
    ) -> Result<Self, WorkflowError> {
        let taxid = record.taxid.as_str();

        // Map lineage and taxonomy information
        let (rank, name, lineage) = match taxonomy {
            Some(taxonomy) => {
                let rank = taxonomy.rank(taxid)?;
                let name = taxonomy.name(taxid)?.to_string();

                let lineage = vec![
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Superkingdom)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Phylum)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Class)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Order)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Family)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Genus)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                    taxonomy
                        .name(taxonomy.parent_at_rank(taxid, TaxRank::Species)?.unwrap_or(("", 0.0)).0)
                        .unwrap_or(""),
                ]
                .into_iter()
                .map(String::from)
                .collect::<Vec<_>>();

                (Some(rank), Some(name), Some(Self::lineage_to_str(lineage)?))
            }
            None => (None, None, None),
        };

        // Extract reads and RPM values from the results
        let mut kraken_reads = None;
        let mut kraken_rpm = None;
        let mut bracken_reads = None;
        let mut bracken_rpm = None;
        let mut metabuli_reads = None;
        let mut metabuli_rpm = None;
        let mut ganon_reads = None;
        let mut ganon_rpm = None;
        let mut kmcp_reads = None;
        let mut kmcp_rpm = None;
        let mut sylph_reads = None;
        let mut sylph_rpm = None;

        for result in &record.results {
            match (result.tool.clone(), result.mode.clone()) {
                (PathogenDetectionTool::Kraken2, PathogenDetectionMode::Sequence) => {
                    kraken_reads = Some(result.reads);
                    kraken_rpm = Some(result.rpm);
                }
                (PathogenDetectionTool::Bracken, PathogenDetectionMode::Profile) => {
                    bracken_reads = Some(result.reads);
                    bracken_rpm = Some(result.rpm);
                }
                (PathogenDetectionTool::Metabuli, PathogenDetectionMode::Sequence) => {
                    metabuli_reads = Some(result.reads);
                    metabuli_rpm = Some(result.rpm);
                }
                (PathogenDetectionTool::Ganon, PathogenDetectionMode::Sequence) => {
                    ganon_reads = Some(result.reads);
                    ganon_rpm = Some(result.rpm);
                }
                (PathogenDetectionTool::Kmcp, PathogenDetectionMode::Sequence) => {
                    kmcp_reads = Some(result.reads);
                    kmcp_rpm = Some(result.rpm);
                }
                (PathogenDetectionTool::Sylph, PathogenDetectionMode::Sequence) => {
                    sylph_reads = Some(result.reads);
                    sylph_rpm = Some(result.rpm);
                }
                _ => {}
            }
        }

        Ok(Self {
            id: id.to_string(),
            taxid: record.taxid.clone(),
            rank,
            name,
            lineage,
            kraken_reads,
            kraken_rpm,
            bracken_reads,
            bracken_rpm,
            metabuli_reads,
            metabuli_rpm,
            ganon_reads,
            ganon_rpm,
            kmcp_reads,
            kmcp_rpm,
            sylph_reads,
            sylph_rpm,
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
    pub fn from_pathogen(
        output: &PathogenOutput,
        quality: &QualityControl,
        paired_end: bool,
    ) -> Result<Self, WorkflowError> {
        let input_reads = quality.reads.input_reads;
        let classifier_reads = quality.reads.output_reads;

        let mut detection_map: HashMap<String, PathogenDetectionRecord> = HashMap::new();

        // Process Kraken2 reports
        if let Some(kraken_report) = &output.profile.kraken2 {
            for record in &kraken_report.records {
                let taxid = record.taxid.trim().to_string();
                let reads = if paired_end { record.reads * 2 } else { record.reads };
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::with_default(&output.id, &taxid));
                entry.add_result(
                    PathogenDetectionTool::Kraken2,
                    PathogenDetectionMode::Sequence,
                    record.taxname.trim().to_string(),
                    record.tax_level.clone(),
                    reads,
                    rpm,
                    abundance,
                );
            }
        }

        // Process Metabuli reports
        if let Some(metabuli_report) = &output.profile.metabuli {
            for record in &metabuli_report.records {
                let taxid = record.taxid.trim().to_string();
                let reads = if paired_end { record.reads * 2 } else { record.reads };
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::with_default(&output.id, &taxid));
                entry.add_result(
                    PathogenDetectionTool::Metabuli,
                    PathogenDetectionMode::Sequence,
                    record.taxname.trim().to_string(),
                    record.tax_level.clone(),
                    reads,
                    rpm,
                    abundance,
                );
            }
        }

        // Process Bracken reports
        if let Some(bracken_report) = &output.profile.bracken {
            for record in &bracken_report.records {
                let taxid = record.taxid.trim().to_string();
                let reads = if paired_end { record.reads * 2 } else { record.reads };
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::with_default(&output.id, &taxid));
                entry.add_result(
                    PathogenDetectionTool::Bracken,
                    PathogenDetectionMode::Profile,
                    record.taxname.trim().to_string(),
                    record.tax_level.clone(),
                    reads,
                    rpm,
                    abundance,
                );
            }
        }

        // Process KMCP reads reports
        if let Some(kmcp_reads) = &output.profile.kmcp_reads {
            let kmcp_taxid_report = kmcp_reads.get_taxid_report()?;
            for record in &kmcp_taxid_report.records {
                let taxid = record.taxid.trim().to_string();
                let reads = if paired_end { record.reads * 2 } else { record.reads };
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::with_default(&output.id, &taxid));
                entry.add_result(
                    PathogenDetectionTool::Kmcp,
                    PathogenDetectionMode::Sequence,
                    record.taxname.clone(),
                    record.rank.clone(),
                    reads,
                    rpm,
                    abundance,
                );
            }
        }

        // Process Sylph reports
        if let Some(sylph_report) = &output.profile.sylph {
            for record in &sylph_report.records {
                let (taxid, rank) = Self::get_sylph_taxinfo(record)?;
                if rank == PathogenDetectionRank::Strain {
                    continue;
                }

                let estimated_reads =
                    ((record.sequence_abundance / 100.0) * classifier_reads as f64) as u64;
                let rpm = compute_rpm(estimated_reads, input_reads).unwrap_or(0.0);

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::with_default(&output.id, &taxid));
                entry.add_result(
                    PathogenDetectionTool::Sylph,
                    PathogenDetectionMode::Sequence,
                    taxid.clone(),
                    format!("{:?}", rank),
                    estimated_reads,
                    rpm,
                    record.sequence_abundance,
                );
            }
        }

        // Process Ganon reads reports
        if let Some(ganon_reads) = &output.profile.ganon_reads {
            for record in &ganon_reads.records {
                let taxid = record.taxid.trim().to_string();
                let reads = if paired_end { record.cumulative * 2 } else { record.cumulative };
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::with_default(&output.id, &taxid));
                entry.add_result(
                    PathogenDetectionTool::Ganon,
                    PathogenDetectionMode::Sequence,
                    record.taxname.trim().to_string(),
                    record.tax_level.clone(),
                    reads,
                    rpm,
                    abundance,
                );
            }
        }

        // Process Ganon abundance reports
        if let Some(ganon_abundance) = &output.profile.ganon_abundance {
            for record in &ganon_abundance.records {
                let taxid = record.taxid.trim().to_string();
                let abundance = record.cumulative_percent;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::with_default(&output.id, &taxid));
                entry.add_result(
                    PathogenDetectionTool::Ganon,
                    PathogenDetectionMode::Profile,
                    record.taxname.trim().to_string(),
                    record.tax_level.clone(),
                    0, // Ganon abundance doesn't report read counts
                    0.0, // No RPM for abundance mode
                    abundance,
                );
            }
        }

        let records = detection_map.into_iter().map(|(_, record)| record).collect();
        Ok(Self {
            id: output.id.clone(),
            paired_end,
            records,
        })
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
        ranks: Option<Vec<String>>,
    ) -> Vec<PathogenDetectionRecord> {
        self.records
            .iter()
            .filter(|record| {
                let taxid_match = taxids.as_ref().map_or(true, |ids| ids.contains(&record.taxid));
                let name_match = names.as_ref().map_or(true, |names| {
                    record.results.iter().any(|res| names.contains(&res.name))
                });
                let rank_match = ranks.as_ref().map_or(true, |ranks| {
                    record
                        .results
                        .iter()
                        .any(|res| ranks.contains(&res.rank))
                });

                taxid_match && name_match && rank_match
            })
            .cloned()
            .collect()
    }
    // Uses the provided taxonomy to create the taxon structs for the database model
    pub fn get_taxa(&self, taxonomy_directory: &PathBuf, strict: bool) -> Result<HashMap<String, Taxon>, WorkflowError> {

        let taxonomy = ncbi::load(taxonomy_directory)?;

        let mut taxa = HashMap::new();
        for record in &self.records {
            let mut taxon = match Taxon::from_taxid(record.taxid.clone(), &taxonomy, true) {
                Ok(taxon) => taxon,
                Err(err) => {
                    log::error!("Failed to find taxid '{}' in provided taxonomy", record.taxid);
                    if strict {
                        return Err(err)
                    } else {
                        log::warn!("Strict mode is not enabled - detection of taxid '{}' will be skipped", record.taxid);
                        continue;
                    }
                }
            };
            taxon.evidence.records.push(record.clone());
            taxa.insert(record.taxid.clone(), taxon);
        }
        Ok(taxa)
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum PathogenDetectionTool {
    Kraken2,
    Metabuli,
    Ganon,
    Kmcp,
    Bracken,
    Sylph,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum PathogenDetectionMode {
    Sequence,
    Profile,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PathogenDetectionResult {
    pub tool: PathogenDetectionTool,
    pub mode: PathogenDetectionMode,
    pub name: String,
    pub rank: String,
    pub reads: u64,
    pub rpm: f64,
    pub abundance: f64,
}

impl PathogenDetectionResult {
    pub fn new(
        tool: PathogenDetectionTool,
        mode: PathogenDetectionMode,
        name: String,
        rank: String,
        reads: u64,
        rpm: f64,
        abundance: f64,
    ) -> Self {
        Self {
            tool,
            mode,
            name,
            rank,
            reads,
            rpm,
            abundance
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PathogenDetectionRecord {
    pub id: String,
    pub taxid: String,
    pub results: Vec<PathogenDetectionResult>, // Results from various tools and modes
}

impl PathogenDetectionRecord {
    pub fn with_default(id: &str, taxid: &str) -> Self {
        Self {
            id: id.to_string(),
            taxid: taxid.to_string(),
            results: Vec::new(),
        }
    }

    pub fn add_result(
        &mut self,
        tool: PathogenDetectionTool,
        mode: PathogenDetectionMode,
        name: String,
        rank: String,
        reads: u64,
        rpm: f64,
        abundance: f64
    ) {
        let result = PathogenDetectionResult::new(tool, mode, name, rank, reads, rpm, abundance);
        self.results.push(result);
    }
}

// #[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
// pub struct PathogenDetectionRecord {
//     pub id: String,
//     pub taxid: String,
//     pub results: Vec<PathogenDetectionResult>,
//     // pub kraken_sequence_name: Option<String>,
//     // pub kraken_sequence_rank: Option<PathogenDetectionRank>,
//     // pub kraken_sequence_reads: Option<u64>,
//     // pub kraken_sequence_rpm: Option<f64>,
//     // pub kraken_sequence_abundance: Option<f64>,
//     // pub metabuli_sequence_name: Option<String>,
//     // pub metabuli_sequence_rank: Option<PathogenDetectionRank>,
//     // pub metabuli_sequence_reads: Option<u64>,
//     // pub metabuli_sequence_rpm: Option<f64>,
//     // pub metabuli_sequence_abundance: Option<f64>,
//     // pub ganon_sequence_name: Option<String>,
//     // pub ganon_sequence_rank: Option<PathogenDetectionRank>,
//     // pub ganon_sequence_reads: Option<u64>,
//     // pub ganon_sequence_rpm: Option<f64>,
//     // pub ganon_sequence_abundance: Option<f64>,
//     // pub kmcp_sequence_name: Option<String>,
//     // pub kmcp_sequence_rank: Option<PathogenDetectionRank>,
//     // pub kmcp_sequence_reads: Option<u64>,
//     // pub kmcp_sequence_rpm: Option<f64>,
//     // pub kmcp_sequence_abundance: Option<f64>,
//     // pub sylph_sequence_name: Option<String>,
//     // pub sylph_sequence_rank: Option<PathogenDetectionRank>,
//     // pub sylph_sequence_reads: Option<u64>,
//     // pub sylph_sequence_rpm: Option<f64>,
//     // pub sylph_sequence_abundance: Option<f64>,
//     // pub bracken_profile_name: Option<String>,
//     // pub bracken_profile_rank: Option<PathogenDetectionRank>,
//     // pub bracken_profile_reads: Option<u64>,
//     // pub bracken_profile_rpm: Option<f64>,
//     // pub bracken_profile_abundance: Option<f64>,
//     // pub kmcp_profile_name: Option<String>,
//     // pub kmcp_profile_rank: Option<PathogenDetectionRank>,
//     // pub kmcp_profile_reads: Option<u64>,
//     // pub kmcp_profile_rpm: Option<f64>,
//     // pub kmcp_profile_abundance: Option<f64>,
//     // pub ganon_profile_name: Option<String>,
//     // pub ganon_profile_rank: Option<PathogenDetectionRank>,
//     // pub ganon_profile_reads: Option<u64>,
//     // pub ganon_profile_rpm: Option<f64>,
//     // pub ganon_profile_abundance: Option<f64>,
//     // pub sylph_profile_name: Option<String>,
//     // pub sylph_profile_rank: Option<PathogenDetectionRank>,
//     // pub sylph_profile_reads: Option<u64>,
//     // pub sylph_profile_rpm: Option<f64>,
//     // pub sylph_profile_abundance: Option<f64>,
// }
// impl Default for PathogenDetectionRecord {
//     fn default() -> Self {
//         Self {
//             id: String::from(""),
//             taxid: String::from(""),
//             results: Vec::new()
//             // kraken_sequence_reads: None,
//             // kraken_sequence_name: None,
//             // kraken_sequence_rank: None,
//             // kraken_sequence_rpm: None,
//             // kraken_sequence_abundance: None,
//             // metabuli_sequence_reads: None,
//             // metabuli_sequence_rpm: None,
//             // metabuli_sequence_name: None,
//             // metabuli_sequence_rank: None,
//             // metabuli_sequence_abundance: None,
//             // ganon_sequence_reads: None,
//             // ganon_sequence_rpm: None,
//             // ganon_sequence_name: None,
//             // ganon_sequence_rank: None,
//             // ganon_sequence_abundance: None,
//             // kmcp_sequence_reads: None,
//             // kmcp_sequence_rpm: None,
//             // kmcp_sequence_name: None,
//             // kmcp_sequence_rank: None,
//             // kmcp_sequence_abundance: None,
//             // sylph_sequence_reads: None,
//             // sylph_sequence_rpm: None,
//             // sylph_sequence_name: None,
//             // sylph_sequence_rank: None,
//             // sylph_sequence_abundance: None,
//             // bracken_profile_reads: None,
//             // bracken_profile_rpm: None,
//             // bracken_profile_name: None,
//             // bracken_profile_rank: None,
//             // bracken_profile_abundance: None,
//             // kmcp_profile_reads: None,
//             // kmcp_profile_rpm: None,
//             // kmcp_profile_name: None,
//             // kmcp_profile_rank: None,
//             // kmcp_profile_abundance: None,
//             // ganon_profile_reads: None,
//             // ganon_profile_rpm: None,
//             // ganon_profile_name: None,
//             // ganon_profile_rank: None,
//             // ganon_profile_abundance: None,
//             // sylph_profile_reads: None,
//             // sylph_profile_rpm: None,
//             // sylph_profile_name: None,
//             // sylph_profile_rank: None,
//             // sylph_profile_abundance: None,
//         }
//     }
// }
// impl PathogenDetectionRecord {
//     pub fn with_default(id: &str, taxid: &str) -> Self {
//         Self {
//             id: id.to_string(),
//             taxid: taxid.to_string(),
//             ..Default::default()
//         }
//     }
//     pub fn set_kraken(&mut self, taxid: String, record: &KrakenReportRecord, input_reads: u64, classifier_reads: u64, paired_end: bool) -> () {
        
//         let reads = if paired_end { record.reads*2 } else { record.reads };

//         let sequence_abundance = (reads as f64 / classifier_reads as f64)*100.0;

//         self.taxid = taxid;
//         self.kraken_sequence_reads = Some(reads);
//         self.kraken_sequence_rpm = compute_rpm(reads, input_reads);
//         self.kraken_sequence_name = Some(record.taxname.trim().to_string());
//         self.kraken_sequence_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
//         self.kraken_sequence_abundance = Some(sequence_abundance);
//     }
//     pub fn set_bracken(&mut self, taxid: String, record: &BrackenReportRecord, input_reads: u64, classifier_reads: u64, paired_end: bool) -> () {

//         let reads = if paired_end { record.reads*2 } else { record.reads };

//         let profile_abundance = (reads as f64 / classifier_reads as f64)*100.0;

//         self.taxid = taxid;
//         self.bracken_profile_reads = Some(reads);
//         self.bracken_profile_name = Some(record.taxname.trim().to_string());
//         self.bracken_profile_rpm = compute_rpm(reads, input_reads);
//         self.bracken_profile_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));

//         self.bracken_profile_abundance = Some(profile_abundance);
//     }
//     pub fn set_metabuli(&mut self, taxid: String, record: &MetabuliReportRecord, input_reads: u64, classifier_reads: u64, paired_end: bool) -> () {

//         let reads = if paired_end { record.reads*2 } else { record.reads };

//         let sequence_abundance = (reads as f64 / classifier_reads as f64)*100.0;

//         self.taxid = taxid;
//         self.metabuli_sequence_reads = Some(reads);
//         self.metabuli_sequence_rpm = compute_rpm(reads, input_reads);
//         self.metabuli_sequence_name = Some(record.taxname.trim().to_string());
//         self.metabuli_sequence_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
//         self.metabuli_sequence_abundance = Some(sequence_abundance);
//     }
//     pub fn set_kmcp_sequence(&mut self, taxid: String, record: &KmcpReadsReportRecord, input_reads: u64, classifier_reads: u64, paired_end: bool) -> () {
        
//         let reads = if paired_end { record.reads*2 } else { record.reads };
//         let sequence_abundance = (reads as f64 / classifier_reads as f64)*100.0;

//         self.taxid = taxid;

//         self.kmcp_sequence_reads = Some(reads);
//         self.kmcp_sequence_rpm = compute_rpm(reads, input_reads);
//         self.kmcp_sequence_name = Some(record.taxname.clone());
//         self.kmcp_sequence_rank = Some(PathogenDetectionRank::from_str(&record.rank));
//         self.kmcp_sequence_abundance = Some(sequence_abundance)
//     }
//     pub fn set_kmcp_profile(&mut self, taxid: String, record: &KmcpAbundanceReportRecord) -> () {
        
//         self.taxid = taxid;

//         let taxname = record.taxname_lineage.split("|").last();

//         self.kmcp_profile_reads = None;
//         self.kmcp_profile_rpm = None;
//         self.kmcp_profile_name = taxname.map(|s| s.to_string());
//         self.kmcp_profile_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
//         self.kmcp_profile_abundance = Some(record.abundance)
//     }
//     pub fn set_ganon_sequence(&mut self, taxid: String, record: &GanonReportRecord, input_reads: u64, classifier_reads: u64, paired_end: bool) -> () {
        
//         let reads = if paired_end { record.cumulative*2 } else { record.cumulative };
//         let sequence_abundance = (reads as f64 / classifier_reads as f64)*100.0;

//         self.taxid = taxid;

//         self.ganon_sequence_reads = Some(reads);
//         self.ganon_sequence_rpm = compute_rpm(reads, input_reads);
//         self.ganon_sequence_name = Some(record.taxname.trim().to_string());
//         self.ganon_sequence_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
//         self.ganon_sequence_abundance = Some(sequence_abundance);

//     }
//     pub fn set_ganon_profile(&mut self, taxid: String, record: &GanonReportRecord) -> () {
        
//         self.taxid = taxid;

//         // let estimated_reads = ((record.cumulative_percent/100.0)*classifier_reads as f64) as u64;

//         self.ganon_profile_reads = None;
//         self.ganon_profile_rpm = None;
//         self.ganon_profile_name = Some(record.taxname.trim().to_string());
//         self.ganon_profile_rank = Some(PathogenDetectionRank::from_str(&record.tax_level));
//         self.ganon_profile_abundance = Some(record.cumulative_percent);

//     }
//     pub fn set_sylph(&mut self, taxid: String, tax_level: PathogenDetectionRank, record: &SylphReportRecord, input_reads: u64, classifier_reads: u64) -> () {
        
//         self.taxid = taxid.clone();

//         let estimated_sequence_reads = ((record.sequence_abundance/100.0)*classifier_reads as f64) as u64;

//         self.sylph_sequence_reads = Some(estimated_sequence_reads);
//         self.sylph_sequence_rpm = compute_rpm(estimated_sequence_reads, input_reads);
//         self.sylph_sequence_name = Some(taxid.to_string());
//         self.sylph_sequence_rank = Some(tax_level.clone());
//         self.sylph_sequence_abundance = Some(record.sequence_abundance);

//         // let estimated_abundance_reads = ((record.relative_abundance/100.0)*classifier_reads as f64) as u64;

//         self.sylph_profile_reads = None;
//         self.sylph_profile_rpm = None;
//         self.sylph_profile_name = Some(taxid.to_string());
//         self.sylph_profile_rank = Some(tax_level.clone());
//         self.sylph_profile_abundance = Some(record.relative_abundance);
//     }
// }
