use anyhow::Result;
use vircov::vircov::VircovRecord;
use serde::{Deserialize, Serialize};
use std::hash::{Hash, Hasher};
use std::{path::PathBuf, collections::HashMap};
use taxonomy::{Taxonomy, GeneralTaxonomy, TaxRank};
use std::collections::hash_map::DefaultHasher;

use crate::modules::assembly::ContigRecord;
use crate::modules::pathogen::{AbundanceMode, ProfileRecord, ProfileTool};
use crate::error::WorkflowError;


pub trait TaxonExtraction {
    fn get_taxa(&self, taxonomy_directory: &PathBuf, strict: bool) -> Result<Vec<Taxon>, WorkflowError>;
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Taxon {
    pub taxid: String,
    pub rank: TaxRank,
    pub name: String,
    pub lineage: TaxonLineage,
    pub evidence: TaxonEvidence
}
impl Taxon {
    pub fn from_taxid(taxid: String, taxonomy: &GeneralTaxonomy) -> Result<Self, WorkflowError> {

        let taxid_str = taxid.as_str();

        let rank = taxonomy.rank(taxid_str)
            .map_err(|err|WorkflowError::TaxRankNotAvailable(err, taxid.clone()))?;

        let name = taxonomy.name(taxid_str)
            .map_err(|err|WorkflowError::TaxNameNotAvailable(err, taxid.clone()))?.to_string();

        Ok(Self { 
            taxid: taxid.clone(), 
            rank,
            name,
            lineage: TaxonLineage::from_taxid(taxid_str, taxonomy)?,
            evidence: TaxonEvidence::new() // to be filled

        })
    }
    pub fn species_data(&self, evidence: bool) -> String {

        let mut str = String::new();

        if let Some(species) = self.lineage.get_species() {
            str.push_str(
                &format!("Species: {species}\n")
            );
        }

        if evidence {
            if self.evidence.profile.is_empty() {
                str.push_str("No evidence available.");
            } else {
                // For each profile record, display key metrics.
                for record in &self.evidence.profile {
                    if record.contigs > 0 {
                        str.push_str(
                            &format!(
                                "{} contigs {} bases | ",
                                record.contigs,
                                record.bases
                            )
                        );
                    } else {
                        str.push_str(
                            &format!(
                                "{:.2} RPM | ",
                                record.rpm
                            )
                        )
                    }   
                }
            }
        }
        str
    }
}

impl std::fmt::Display for Taxon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        let elaborate = false;
        
        if elaborate {
            // Header with basic taxon information in Markdown format.
            writeln!(f, "## Taxon: {} (TaxID: {})", self.name, self.taxid)?;
            writeln!(f, "**Rank:** {}", self.rank)?;
            writeln!(f, "**Lineage:** {}", self.lineage)?;
            writeln!(f)?;

            // --- Profile Evidence Summary ---
            writeln!(f, "### Profile Evidence")?;
            if self.evidence.profile.is_empty() {
                writeln!(f, "_No profile evidence available._")?;
            } else {
                // Markdown table header for profile evidence.
                writeln!(f, "| Tool    | Mode    | Reads | RPM   | Abundance |")?;
                writeln!(f, "|---------|---------|-------|-------|-----------|")?;
                // For each profile record, display key metrics.
                for record in &self.evidence.profile {
                    writeln!(
                        f,
                        "| {} | {} | {} | {:.2} | {:.2} |",
                        record.tool.to_string(),
                        record.mode.to_string(),
                        record.reads,
                        record.rpm,
                        record.abundance
                    )?;
                }
            }
            writeln!(f)?;

            // --- Alignment Evidence Summary ---
            writeln!(f, "### Alignment Evidence")?;
            if self.evidence.alignment.is_empty() {
                writeln!(f, "_No alignment evidence available._")?;
            } else {
                // For alignment evidence, we focus on counts and (if available) aggregated RPM.
                // (Replace the placeholder aggregation with real field extraction as needed.)
                let alignment_count = self.evidence.alignment.len();
                let summed_alignment_rpm: f64 = self.evidence.alignment.iter().map(|_rec| {
                    // Assume each VircovRecord has an 'rpm' field.
                    // For demonstration, we use a placeholder value of 0.0.
                    0.0
                }).sum();
                writeln!(f, "- **Total Alignment Records:** {}", alignment_count)?;
                writeln!(f, "- **Summed RPM:** {:.2}", summed_alignment_rpm)?;
            }
            writeln!(f)?;

            // --- Assembly Evidence Summary ---
            writeln!(f, "### Assembly Evidence")?;
            if self.evidence.assembly.is_empty() {
                writeln!(f, "_No assembly evidence available._")?;
            } else {
                // For assembly evidence, we assume you wish to report total contigs, bases, and an RPM metric.
                // Replace the following placeholder aggregations with actual data extraction.
                let total_contigs: u64 = self.evidence.assembly.iter().map(|_rec| {
                    // Assume each ContigRecord has a field 'contigs'
                    0
                }).sum();
                let total_bases: u64 = self.evidence.assembly.iter().map(|_rec| {
                    // Assume each ContigRecord has a field 'bases'
                    0
                }).sum();
                let summed_assembly_rpm: f64 = self.evidence.assembly.iter().map(|_rec| {
                    // Assume each ContigRecord has an RPM metric.
                    0.0
                }).sum();

                // Markdown table for assembly evidence.
                writeln!(f, "| Metric         | Value   |")?;
                writeln!(f, "|----------------|---------|")?;
                writeln!(f, "| Total Contigs  | {} |", total_contigs)?;
                writeln!(f, "| Total Bases    | {} |", total_bases)?;
                writeln!(f, "| Summed RPM     | {:.2} |", summed_assembly_rpm)?;
            }
            writeln!(f)?;

        } else {

            if let Some(species) = self.lineage.get_species() {
                writeln!(f, "Species: {species}")?;
            }


            // --- Profile Evidence Summary ---
            if self.evidence.profile.is_empty() {
                writeln!(f, "No evidence available.")?;
            } else {
                // For each profile record, display key metrics.
                let mut results = String::new();
                for record in &self.evidence.profile {
                    if record.contigs > 0 {
                        results.push_str(
                            &format!(
                                "{} contigs {} bases |",
                                record.contigs,
                                record.bases
                            )
                        );
                    } else {
                        results.push_str(
                            &format!(
                                "{:.2} RPM |",
                                record.rpm
                            )
                        )
                    }
                    
                }
                writeln!(f, "{}", results)?;
            }
            writeln!(f)?;

        }
        
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct TaxonEvidence {
    pub alignment: Vec<VircovRecord>,
    pub assembly: Vec<ContigRecord>,
    pub profile: Vec<ProfileRecord>,
}

impl TaxonEvidence {
    pub fn new() -> Self {
        Self { 
            alignment: Vec::new(),
            assembly: Vec::new(),
            profile: Vec::new()
        }
    }
}

// Lineage in GTDB liek format

pub trait LineageOperations: AsRef<str> {
    fn from_taxid(taxid: &str, taxonomy: &GeneralTaxonomy) -> Result<String, WorkflowError> {
        
        match taxonomy.rank(taxid) {
            Ok(_) => {

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

                Ok(lineage_to_gtdb_str(lineage)?)
            },
            Err(err) => {
                log::warn!("Error from taxonomy: {}", err.to_string());
                log::warn!("Rank could not be determined for taxid {taxid} - using unclassified lineage (d__;p__;c__;o__;f__;g__;s__)");
                Ok(String::from("d__;p__;c__;o__;f__;g__;s__"))
            }
        }
    }

    /// Helper method: Returns the taxon name for a given prefix (e.g. "d__", "p__", etc.)
    fn get_taxon_by_prefix(&self, prefix: &str) -> Option<String> {
        self.as_ref()
            .split(";")
            .find(|s| s.trim_start().starts_with(prefix))
            .map(|s| s.trim_start_matches(prefix).to_string())
    }

    fn get_domain(&self) -> Option<String> {
        self.get_taxon_by_prefix("d__")
    }

    fn get_phylum(&self) -> Option<String> {
        self.get_taxon_by_prefix("p__")
    }

    fn get_class(&self) -> Option<String> {
        self.get_taxon_by_prefix("c__")
    }

    fn get_order(&self) -> Option<String> {
        self.get_taxon_by_prefix("o__")
    }

    fn get_family(&self) -> Option<String> {
        self.get_taxon_by_prefix("f__")
    }

    fn get_genus(&self) -> Option<String> {
        self.get_taxon_by_prefix("g__")
    }

    fn get_species(&self) -> Option<String> {
        self.get_taxon_by_prefix("s__")
    }

    fn get_labels(&self) -> Vec<&str> {
        self.as_ref().split(";").collect()
    }
}


/// For a given taxon rank, return the corresponding index in the GTDB-style lineage string
/// and the prefix (e.g., "d__" for Superkingdom, "s__" for Species, etc.).
pub fn rank_index_and_prefix(rank: &TaxRank) -> Option<(usize, &'static str)> {
    match rank {
        TaxRank::Superkingdom => Some((0, "d__")),
        TaxRank::Phylum       => Some((1, "p__")),
        TaxRank::Class        => Some((2, "c__")),
        TaxRank::Order        => Some((3, "o__")),
        TaxRank::Family       => Some((4, "f__")),
        TaxRank::Genus        => Some((5, "g__")),
        TaxRank::Species      => Some((6, "s__")),
        _ => None,
    }
}

pub fn lineage_to_gtdb_str(lineage: Vec<String>) -> Result<String, WorkflowError> {
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

pub type TaxonLineage = String;

impl LineageOperations for TaxonLineage {}

// Utility function wrapping the search result handling with name and taxid
pub fn get_taxid_name_from_search(search_result: Option<(&str, f32)>, taxonomy: &GeneralTaxonomy) -> Result<(Option<String>, Option<String>), WorkflowError> {
     match search_result {
        Some((taxid, _)) => {
            let name = taxonomy.name(taxid).map_err(
                |err|WorkflowError::TaxNameNotAvailable(err, taxid.to_owned())
            )?;
            Ok((Some(taxid.to_owned()), Some(name.to_owned())))

        },
        None => Ok((None, None))
    }
}

/* 
=======================
MODULE AND TOOL STRUCTS
=======================
*/

// Utility function to aggregate multiple HashMaps with Taxon
pub fn aggregate(parent_taxa: &mut HashMap<String, Taxon>, taxa: &Vec<Taxon>) -> HashMap<String, Taxon> {
    
    for (taxid, taxon) in taxa.iter().map(|t| (t.taxid.as_str(), t)) {
        // ... check if each taxon exists
        match parent_taxa.get_mut(taxid) {
            Some(tax) => {
                // If so, update the evidence of the taxon in the modules...
                for record in &taxon.evidence.profile {
                    tax.evidence.profile.push(record.clone())
                }
                // If so, update the evidence of the taxon in the modules...
                for record in &taxon.evidence.alignment {
                    tax.evidence.alignment.push(record.clone())
                }
                // If so, update the evidence of the taxon in the modules...
                for record in &taxon.evidence.assembly {
                    tax.evidence.assembly.push(record.clone())
                }
            },
            None => {
                // otherwise insert the taxon from result with its own evidence
                parent_taxa.insert(taxid.to_owned(), taxon.to_owned());
            }
        }
    }  
    parent_taxa.to_owned()
}

/// Helper: Remove trailing underscore and block of uppercase letters if present.
fn get_base_name(name: &str) -> &str {
    if let Some(pos) = name.rfind('_') {
        let suffix = &name[pos + 1..];
        if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_uppercase()) {
            return &name[..pos];
        }
    }
    name
}


/// Collapses variants of taxa with the same base name (e.g. "Haemophilus influenzae_AC" and
/// "Haemophilus influenzae_B" become "Haemophilus influenzae"). In the collapsed Taxon,
/// alignment and assembly records are concatenated and the lineage is updated so that the
/// component for its rank is replaced with the base name. For profile evidence, for each unique
/// (id, tool, mode) combination, the numeric fields are summed.
pub fn collapse_taxa(taxa: Vec<Taxon>) -> Result<Vec<Taxon>, WorkflowError> {

    // Group taxa by their base name.
    let mut groups: HashMap<String, Vec<Taxon>> = HashMap::new();
    for taxon in taxa {
        let base = get_base_name(&taxon.name).to_string();
        groups.entry(base).or_default().push(taxon);
    }

    let mut collapsed = Vec::new();

    // Process each group.
    for (base_name, group) in groups {
        
        // If only a single taxon in group, leave it untouched and continue with groups
        if group.len() == 1 {
            collapsed.push(group[0].clone());
            continue;
        }

        // Choose representative: search the group for any Taxon whose name matches the base name.
        // If not found, use the first taxon in the group.
        let representative = group
            .iter()
            .find(|taxon| taxon.name == base_name)
            .unwrap_or(&group[0]);

        let rank = &representative.rank;

        // Update the lineage: replace the part corresponding to the taxon's rank with the base name.
        let mut new_lineage = representative.lineage.clone();
        if let Some((idx, prefix)) = rank_index_and_prefix(rank) {
            let mut parts: Vec<String> = new_lineage.split(';').map(|s| s.to_string()).collect();
            if parts.len() > idx {
                parts[idx] = format!("{}{}", prefix, base_name);
                new_lineage = parts.join(";");
            } else {
                return Err(WorkflowError::LineageStringTooShort(new_lineage));
            }
        }

        // Initialize accumulators for evidence.
        let mut alignment = Vec::new();
        let mut assembly = Vec::new();
        let mut profile_map: HashMap<(String, ProfileTool, AbundanceMode), ProfileRecord> = HashMap::new();

        // Accumulate evidence from all taxa in the group.
        for taxon in &group {
            // Concatenate alignment and assembly records.
            alignment.extend(taxon.evidence.alignment.iter().cloned());
            assembly.extend(taxon.evidence.assembly.iter().cloned());

            // For profile evidence, sum the numeric fields for each (id, tool, mode) key.
            for record in &taxon.evidence.profile {
                let key = (record.id.clone(), record.tool.clone(), record.mode.clone());
                let entry = profile_map.entry(key).or_insert(ProfileRecord {
                    id: record.id.clone(),
                    tool: record.tool.clone(),
                    mode: record.mode.clone(),
                    reads: 0,
                    rpm: 0.0,
                    contigs: 0,
                    bases: 0,
                    bpm: 0.0,
                    abundance: 0.0,
                });
                entry.reads += record.reads;
                entry.rpm += record.rpm;
                entry.contigs += record.contigs;
                entry.bases += record.bases;
                entry.bpm += record.bpm;
                entry.abundance += record.abundance;
            }
        }

        let profile = profile_map.into_values().collect();

        // Create a new taxid for this collapsed taxon in the format: collapsed-{hash of rank and taxon basename} 
        // this is so the taxids can be found in subsequent requests that use taxids e.g. for prevalence  
        // contamination or taxon history - which also need to call this function on the taxa!

        let mut hasher = DefaultHasher::new();
        format!("collapsed-{}-{}", representative.rank, base_name).hash(&mut hasher);
        
        let new_taxid = hasher.finish();

        // Build the new, collapsed Taxon.
        let new_taxon = Taxon {
            taxid: new_taxid.to_string(), 
            rank: representative.rank.clone(),
            name: base_name,
            lineage: new_lineage,
            evidence: TaxonEvidence {
                alignment,
                assembly,
                profile,
            },
        };

        collapsed.push(new_taxon);
    }

    Ok(collapsed)
}