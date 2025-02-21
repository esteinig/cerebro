use std::{collections::HashMap, fs::File, io::{BufReader, BufWriter}, path::{Path, PathBuf}};

use serde::{Deserialize, Serialize};
use taxonomy::{ncbi, GeneralTaxonomy, TaxRank, Taxonomy};
use vircov::vircov::VircovRecord;

use crate::{error::WorkflowError, nextflow::pathogen::{PathogenDetectionOutput, SylphReportRecord}, taxa::taxon::{Taxon, TaxonExtraction}, utils::{read_tsv, write_tsv}};

use super::{alignment::Alignment, assembly::{MetagenomeAssembly, ContigRecord}, quality::QualityControl};


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenDetectionTableRecord {
    id: String,
    taxid: String,
    rank: Option<TaxRank>,
    name: Option<String>,
    lineage: Option<String>,
    vircov_reads: Option<u64>,
    vircov_rpm: Option<f64>,
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
    blast_contigs: Option<u64>,
    blast_bases: Option<u64>,
    blast_bpm: Option<f64>,
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
                
                match taxonomy.rank(taxid) {
                    Ok(rank) => {
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
                    },
                    Err(_) => (Some(TaxRank::Unspecified), Some(record.name.clone()), None)
                }
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
        let mut vircov_reads = None;
        let mut vircov_rpm = None;
        let mut blast_contigs = None;
        let mut blast_bases = None;
        let mut blast_bpm = None;

        for result in &record.results {
            match (result.tool.clone(), result.mode.clone()) {
                (ProfileTool::Kraken2, AbundanceMode::Sequence) => {
                    kraken_reads = Some(result.reads);
                    kraken_rpm = Some(result.rpm);
                }
                (ProfileTool::Bracken, AbundanceMode::Profile) => {
                    bracken_reads = Some(result.reads);
                    bracken_rpm = Some(result.rpm);
                }
                (ProfileTool::Metabuli, AbundanceMode::Sequence) => {
                    metabuli_reads = Some(result.reads);
                    metabuli_rpm = Some(result.rpm);
                }
                (ProfileTool::Ganon2, AbundanceMode::Sequence) => {
                    ganon_reads = Some(result.reads);
                    ganon_rpm = Some(result.rpm);
                }
                (ProfileTool::Kmcp, AbundanceMode::Sequence) => {
                    kmcp_reads = Some(result.reads);
                    kmcp_rpm = Some(result.rpm);
                }
                (ProfileTool::Sylph, AbundanceMode::Sequence) => {
                    sylph_reads = Some(result.reads);
                    sylph_rpm = Some(result.rpm);
                }
                 // Alignment can be multiple per taxonomic identifier (segments, contigs)
                (ProfileTool::Vircov, AbundanceMode::Sequence) => {
                    vircov_reads = Some(match vircov_reads {
                        Some(existing) => existing + result.reads,
                        None => result.reads,
                    });
                    vircov_rpm = Some(match vircov_rpm {
                        Some(existing) => existing + result.rpm,
                        None => result.rpm,
                    });
                }
                // Assembly can be multiple per taxonomic identifier (contigs)
                (ProfileTool::Blast, AbundanceMode::Bases) => {
                    blast_contigs = Some(match blast_contigs {
                        Some(existing) => existing + 1,
                        None => 1,
                    });
                    blast_bases = Some(match blast_bases {
                        Some(existing) => existing + result.bases,
                        None => result.bases,
                    });
                    blast_bpm = Some(match blast_bpm {
                        Some(existing) => existing + result.bpm,
                        None => result.bpm,
                    });
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
            vircov_reads,
            vircov_rpm,
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
            blast_contigs,
            blast_bases,
            blast_bpm
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
            None => pd.profile,
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

fn compute_bpm(bases: u64, input_bases: u64) -> Option<f64> {
    if input_bases == 0 {
        None  // Cannot compute BPM when input_bases is 0
    } else {
        Some((bases as f64 / input_bases as f64) * 1_000_000.0)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenDetectionFilter {
    pub taxids: Option<Vec<String>>,
    pub names: Option<Vec<String>>,
    pub ranks: Option<Vec<PathogenDetectionRank>>
}
impl PathogenDetectionFilter {
    pub fn new(
        taxids: Option<Vec<String>>,
        names: Option<Vec<String>>,
        ranks: Option<Vec<PathogenDetectionRank>>
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
    pub profile: Vec<PathogenDetectionRecord>,
    pub alignment: Vec<VircovRecord>,
    pub assembly: Vec<ContigRecord>
}
impl PathogenDetection {
    pub fn from_pathogen(
        output: &PathogenDetectionOutput,
        quality: &QualityControl,
        paired_end: bool,
    ) -> Result<Self, WorkflowError> {

        // Submodules for alignment and metagenome assembly run in pathogen detection

        let alignment = Alignment::from_pathogen(&output)?;
        let assembly = MetagenomeAssembly::from_pathogen(&output)?;

        // Relevant quality control data
        let input_reads = quality.reads.input_reads;
        let input_bases = quality.reads.input_bases;
        let classifier_reads = quality.reads.output_reads;
        let classifier_bases = quality.reads.output_bases;

        let mut detection_map: HashMap<String, PathogenDetectionRecord> = HashMap::new();

        for record in &alignment.records {

            let taxid = match &record.taxid {
                Some(taxid) => taxid,
                None => return Err(WorkflowError::PathogenTaxidAnnotationMissing)
            };

            let taxname = match &record.name {
                Some(taxname) => taxname,
                None => match &record.bin {
                    Some(bin) => bin,  // fallback to bin name
                    None => return Err(WorkflowError::PathogenTaxnameAnnotationMissing)
                }
            };

            let remap_reads = match record.remap_alignments {
                Some(reads) => reads,
                None => continue // remapping stage executed but unsuccessful - can occur when very few reads detected in scanning stage
            };

            let taxid = taxid.trim().to_string();
            let name = taxname.trim().to_string();
            let rank = PathogenDetectionRank::Species; // alignment always species-level for now
            let reads = remap_reads;
            let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
            let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

            let entry = detection_map
                .entry(taxid.clone())
                .or_insert_with(|| PathogenDetectionRecord::new(&output.id, &taxid, &name, rank));
            entry.add_result(
                ProfileTool::Vircov,
                AbundanceMode::Sequence,
                reads,
                rpm,
                0,
                0.0,
                abundance,
            );
        }

        if let Some(records) = &output.assembly.blast {

            for record in records {
                let mut taxid = record.taxid.trim().to_string();
                let mut name = record.taxname.trim().to_string();
                
                // Custom BLAST databases don't usually have the 'ssciname' associated 
                // we could use a taxonomy to add this later, but here we use the
                // 'stitle' as taxname as configured in 'Cipher' with sequence description
                // in the format: {taxid}:::{name} 

                if name == "N/A".to_string() || taxid == "N/A".to_string() {
                    let description = record.title.trim().split(":::").collect::<Vec<_>>();
                    taxid = description.first().ok_or(WorkflowError::PathogenTaxidAnnotationMissing)?.to_string();
                    name = description.last().ok_or(WorkflowError::PathogenTaxnameAnnotationMissing)?.to_string();
                }

                let rank = PathogenDetectionRank::from_str(&record.taxrank);
                let bases = record.length;
                let bpm = compute_bpm(record.length, input_bases).unwrap_or(0.0);
                let abundance = (bases as f64 / classifier_bases as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::new(&output.id, &taxid, &name, rank));
                entry.add_result(
                    ProfileTool::Blast,
                    AbundanceMode::Bases,
                    0,
                    0.0,
                    bases,
                    bpm,
                    abundance,
                );
            }
        }

        // Process Kraken2 reports
        if let Some(kraken_report) = &output.profile.kraken2 {
            for record in &kraken_report.records {
                let taxid = record.taxid.trim().to_string();
                let name = record.taxname.trim().to_string();
                let rank = PathogenDetectionRank::from_str(&record.tax_level);
                let reads = if paired_end { record.reads * 2 } else { record.reads };
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::new(&output.id, &taxid, &name, rank));
                entry.add_result(
                    ProfileTool::Kraken2,
                    AbundanceMode::Sequence,
                    reads,
                    rpm,
                    0,
                    0.0,
                    abundance,
                );
            }
        }

        // Process Metabuli reports
        if let Some(metabuli_report) = &output.profile.metabuli {
            for record in &metabuli_report.records {
                let taxid = record.taxid.trim().to_string();
                let name = record.taxname.trim().to_string();
                let rank = PathogenDetectionRank::from_str(&record.tax_level);
                let reads = if paired_end { record.reads * 2 } else { record.reads };
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::new(&output.id, &taxid, &name, rank));
                entry.add_result(
                    ProfileTool::Metabuli,
                    AbundanceMode::Sequence,
                    reads,
                    rpm,
                    0,
                    0.0,
                    abundance,
                );
            }
        }

        // Process Bracken reports
        if let Some(bracken_report) = &output.profile.bracken {
            for record in &bracken_report.records {
                let taxid = record.taxid.trim().to_string();
                let name = record.taxname.trim().to_string();
                let rank = PathogenDetectionRank::from_str(&record.tax_level);
                let reads = if paired_end { record.reads * 2 } else { record.reads };
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::new(&output.id, &taxid, &name, rank));
                entry.add_result(
                    ProfileTool::Bracken,
                    AbundanceMode::Profile,
                    reads,
                    rpm,
                    0,
                    0.0,
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

                let name = taxid.clone();
                let estimated_reads = ((record.sequence_abundance / 100.0) * classifier_reads as f64) as u64;
                let rpm = compute_rpm(estimated_reads, input_reads).unwrap_or(0.0);

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::new(&output.id, &taxid, &name, rank));
                entry.add_result(
                    ProfileTool::Sylph,
                    AbundanceMode::Sequence,
                    estimated_reads,
                    rpm,
                    0,
                    0.0,
                    record.sequence_abundance,
                );
            }
        }

        // Process KMCP reads reports
        if let Some(kmcp_reads) = &output.profile.kmcp_reads {
            let kmcp_taxid_report = kmcp_reads.get_taxid_report()?;
            for record in &kmcp_taxid_report.records {
                let taxid = record.taxid.trim().to_string();
                let name = record.taxname.clone();
                let rank = PathogenDetectionRank::from_str(&record.rank);
                let reads = if paired_end { record.reads * 2 } else { record.reads };
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::new(&output.id, &taxid, &name, rank));
                entry.add_result(
                    ProfileTool::Kmcp,
                    AbundanceMode::Sequence,
                    reads,
                    rpm,
                    0,
                    0.0,
                    abundance,
                );
            }
        }

        // Process KMCP abundance reports
        if let Some(kmcp_abundance) = &output.profile.kmcp_abundance {
            for record in &kmcp_abundance.records {
                let taxid = record.taxid.trim().to_string();
                let name =  record.taxname_lineage.split("|").last().unwrap_or("ERROR");  // TODO: check what we can do to put a placeholder
                let rank = PathogenDetectionRank::from_str(&record.tax_level);
                let reads = 0;
                let rpm = 0.0;
                let abundance = record.abundance;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::new(&output.id, &taxid, &name, rank));
                entry.add_result(
                    ProfileTool::Kmcp,
                    AbundanceMode::Profile,
                    reads,
                    rpm,
                    0,
                    0.0,
                    abundance,
                );
            }
        }

        // Process Ganon reads reports
        if let Some(ganon_reads) = &output.profile.ganon_reads {
            for record in &ganon_reads.records {
                let taxid = record.taxid.trim().to_string();
                let name = record.taxname.trim().to_string();
                let rank = PathogenDetectionRank::from_str(&record.tax_level);
                let reads = record.cumulative; // Ganon respects modern Illumina PE read identifier format whereas Kraken2/Bracken and Metabuli do not!
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = (reads as f64 / classifier_reads as f64) * 100.0;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::new(&output.id, &taxid, &name, rank));
                entry.add_result(
                    ProfileTool::Ganon2,
                    AbundanceMode::Sequence,
                    reads,
                    rpm,
                    0,
                    0.0,
                    abundance,
                );
            }
        }


        // Process Ganon abundance reports
        if let Some(ganon_abundance) = &output.profile.ganon_abundance {
            for record in &ganon_abundance.records {
                let taxid = record.taxid.trim().to_string();
                let name = record.taxname.trim().to_string();
                let rank = PathogenDetectionRank::from_str(&record.tax_level);
                let reads = if paired_end { record.cumulative * 2 } else { record.cumulative };
                let rpm = compute_rpm(reads, input_reads).unwrap_or(0.0);
                let abundance = record.cumulative_percent;

                let entry = detection_map
                    .entry(taxid.clone())
                    .or_insert_with(|| PathogenDetectionRecord::new(&output.id, &taxid, &name, rank));
                entry.add_result(
                    ProfileTool::Ganon2,
                    AbundanceMode::Profile,
                    reads,
                    rpm,
                    0,
                    0.0,
                    abundance,
                );
            }
        }

        let profile = detection_map.into_iter().map(|(_, record)| record).collect();
        
        Ok(Self {
            id: output.id.clone(),
            paired_end,
            profile,
            alignment: alignment.records,
            assembly: assembly.records,
        })
    }

    fn get_sylph_taxinfo(
        record: &SylphReportRecord,
    ) -> Result<(String, PathogenDetectionRank), WorkflowError> {
        let tax_str = record.clade_name.trim().to_string();
        let tax_str = tax_str.split('|').last();

        let (taxid, taxlevel) = match tax_str {
            Some(tax_str) => {
                let taxid_str_split = tax_str.split("__").collect::<Vec<_>>();
                if taxid_str_split.len() != 2 {
                    log::warn!("Failed to recover taxinfo from Sylph output");
                    return Err(WorkflowError::SylphTaxInfoRecoveryFailure(
                        record.clade_name.clone(),
                    ));
                } else {
                    (taxid_str_split[1], taxid_str_split[0])
                }
            }
            None => {
                log::warn!("Failed to recover taxid from Sylph output");
                return Err(WorkflowError::SylphTaxInfoRecoveryFailure(
                    record.clade_name.clone(),
                ));
            }
        };

        Ok((taxid.to_string(), PathogenDetectionRank::from_str(taxlevel)))
    }
    pub fn filter_by_taxonomy(
        &self,
        taxids: Option<Vec<String>>,
        names: Option<Vec<String>>,
        ranks: Option<Vec<PathogenDetectionRank>>,
    ) -> Vec<PathogenDetectionRecord> {
        self.profile
            .iter()
            .filter(|record| {
                let taxid_match = taxids.as_ref().map_or(true, |ids| ids.contains(&record.taxid));
                let name_match = names.as_ref().map_or(true, |names| names.contains(&record.name));
                let rank_match = ranks.as_ref().map_or(true, |ranks| {
                    ranks.contains(&record.rank)
                });

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

impl TaxonExtraction for PathogenDetection {
    fn get_taxa(&self, taxonomy_directory: &PathBuf, strict: bool) -> Result<HashMap<String, Taxon>, WorkflowError> {
        let taxonomy = ncbi::load(taxonomy_directory)?;

        let mut taxa: HashMap<String, Taxon> = HashMap::new();
        for record in &self.profile {

            if let Some(existing_taxon) = taxa.get_mut(&record.taxid) {
                // If the taxid is already present, update its evidence
                for profile_record in &record.results {
                    existing_taxon.evidence.profile.push(profile_record.clone());
                }
            } else {
                // Otherwise, create a new Taxon and insert it into the HashMap
                let mut taxon = match Taxon::from_taxid(record.taxid.clone(), &taxonomy, true) {
                    Err(err) => {
                        log::error!("Failed to find taxid '{}' in provided taxonomy", &record.taxid);
                        if strict {
                            return Err(err);
                        } else {
                            log::warn!(
                                "Strict mode is not enabled - detection of taxid '{}' will be skipped",
                                record.taxid
                            );
                            continue;
                        }
                    },
                    Ok(taxon) => taxon,
                };
                
                for profile_record in &record.results {
                    taxon.evidence.profile.push(profile_record.clone());
                }
                taxa.insert(record.taxid.clone(), taxon);
            }
        }

        for record in &self.alignment {

            let taxid = match &record.taxid {
                Some(taxid) => taxid,
                None => return Err(WorkflowError::PathogenTaxidAnnotationMissing)
            };

            if let Some(existing_taxon) = taxa.get_mut(taxid) {
                // If the taxid is already present, update its evidence
                existing_taxon.evidence.alignment.push(record.clone());
            } else {
                // Otherwise, create a new Taxon and insert it into the HashMap
                let mut taxon = match Taxon::from_taxid(taxid.clone(), &taxonomy, true) {
                    Err(err) => {
                        log::error!("Failed to find taxid '{}' in provided taxonomy", &taxid);
                        if strict {
                            return Err(err);
                        } else {
                            log::warn!(
                                "Strict mode is not enabled - detection of taxid '{}' will be skipped",
                                taxid
                            );
                            continue;
                        }
                    },
                    Ok(taxon) => taxon,
                };
                
                taxon.evidence.alignment.push(record.clone());
                taxa.insert(taxid.clone(), taxon);
            }
        }

        Ok(taxa)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, clap::ValueEnum)]
pub enum PathogenDetectionRank {
    Superkingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Strain,
    Root,
    Unclassified,
    NoRank,
    Other
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
            "r" | "R" | "root" => Self::Root, 
            "unspecified" | "no rank" | "" => Self::NoRank,
            "unclassified" | "-" => Self::Unclassified,
            _ => Self::Other
        }
    }
}
impl PathogenDetectionRank {
    pub fn to_string(&self) -> String {
        match self {
            Self::Superkingdom => "Superkingdom".to_string(),
            Self::Phylum => "Phylum".to_string(),
            Self::Class => "Class".to_string(),
            Self::Order => "Order".to_string(),
            Self::Family => "Family".to_string(),
            Self::Genus => "Genus".to_string(),
            Self::Species => "Species".to_string(),
            Self::Strain => "Strain".to_string(),
            Self::NoRank => "NoRank".to_string(),
            Self::Root => "Root".to_string(),
            Self::Unclassified => "Unclassified".to_string(),
            Self::Other => "Other".to_string(),
        }
    }
}

impl From<TaxRank> for PathogenDetectionRank {
    fn from(tax_rank: TaxRank) -> Self {
        match tax_rank {
            TaxRank::Superkingdom => PathogenDetectionRank::Superkingdom,
            TaxRank::Phylum => PathogenDetectionRank::Phylum,
            TaxRank::Class => PathogenDetectionRank::Class,
            TaxRank::Order => PathogenDetectionRank::Order,
            TaxRank::Family => PathogenDetectionRank::Family,
            TaxRank::Genus => PathogenDetectionRank::Genus,
            TaxRank::Species => PathogenDetectionRank::Species,
            TaxRank::Strain => PathogenDetectionRank::Strain,
            TaxRank::Unspecified => PathogenDetectionRank::NoRank,
            _ => PathogenDetectionRank::Other,
        }
    }
}

impl From<PathogenDetectionRank> for TaxRank {
    fn from(rank: PathogenDetectionRank) -> Self {
        match rank {
            PathogenDetectionRank::Superkingdom => TaxRank::Superkingdom,
            PathogenDetectionRank::Phylum => TaxRank::Phylum,
            PathogenDetectionRank::Class => TaxRank::Class,
            PathogenDetectionRank::Order => TaxRank::Order,
            PathogenDetectionRank::Family => TaxRank::Family,
            PathogenDetectionRank::Genus => TaxRank::Genus,
            PathogenDetectionRank::Species => TaxRank::Species,
            PathogenDetectionRank::Strain => TaxRank::Strain,
            _ => TaxRank::Unspecified,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum ProfileTool {
    Blast,
    Kraken2,
    Metabuli,
    Ganon2,
    Kmcp,
    Bracken,
    Sylph,
    Vircov
}
impl ProfileTool {
    pub fn to_string(&self) -> String {
        match self {
            Self::Blast => "Blast".to_string(),
            Self::Kraken2 => "Kraken2".to_string(),
            Self::Metabuli => "Metabuli".to_string(),
            Self::Ganon2 => "Ganon2".to_string(),
            Self::Kmcp => "Kmcp".to_string(),
            Self::Bracken => "Bracken".to_string(),
            Self::Sylph => "Sylph".to_string(),
            Self::Vircov => "Vircov".to_string(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum AbundanceMode {
    Bases,
    Profile,
    Sequence,
}
impl AbundanceMode {
    pub fn to_string(&self) -> String {
        match self {
            Self::Bases => "Bases".to_string(),
            Self::Profile => "Profile".to_string(),
            Self::Sequence => "Sequence".to_string(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ProfileRecord {
    pub id: String,
    pub tool: ProfileTool,
    pub mode: AbundanceMode,
    pub reads: u64,
    pub rpm: f64,
    pub bases: u64,
    pub bpm: f64,
    pub abundance: f64,
}

impl ProfileRecord {
    pub fn new(
        id: String,
        tool: ProfileTool,
        mode: AbundanceMode,
        reads: u64,
        rpm: f64,
        bases: u64, 
        bpm: f64,
        abundance: f64,
    ) -> Self {
        Self {
            id,
            tool,
            mode,
            reads,
            rpm,
            bases,
            bpm,
            abundance
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct PathogenDetectionRecord {
    pub id: String,
    pub taxid: String,
    pub name: String,
    pub rank: PathogenDetectionRank,
    pub results: Vec<ProfileRecord>,
}

impl PathogenDetectionRecord {
    pub fn new(id: &str, taxid: &str, name: &str, rank: PathogenDetectionRank) -> Self {
        Self {
            id: id.to_string(),
            taxid: taxid.to_string(),
            name: name.to_string(),
            rank,
            results: Vec::new(),
        }
    }

    pub fn add_result(
        &mut self,
        tool: ProfileTool,
        mode: AbundanceMode,
        reads: u64,
        rpm: f64,
        bases: u64,
        bpm: f64,
        abundance: f64
    ) {
        let result = ProfileRecord::new(self.id.clone(), tool, mode, reads, rpm, bases, bpm, abundance);
        self.results.push(result);
    }
}