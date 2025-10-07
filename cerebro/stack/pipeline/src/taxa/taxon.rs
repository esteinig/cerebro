use anyhow::Result;
use itertools::Itertools;
use vircov::vircov::VircovRecord;
use serde::{Deserialize, Serialize};
use std::{path::PathBuf, collections::HashMap};
use taxonomy::{Taxonomy, GeneralTaxonomy, TaxRank};

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use crate::modules::assembly::ContigRecord;
use crate::modules::pathogen::{AbundanceMode, ProfileRecord, ProfileTool};
use crate::error::WorkflowError;


pub trait TaxonExtraction {
    fn get_taxa(&self, taxonomy_directory: &PathBuf, strict: bool, gtdb_break_monophyly: bool) -> Result<Vec<Taxon>, WorkflowError>;
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

        let mut virus_data = Vec::new();
        if let Some(domain) = self.lineage.get_domain() {
            if domain.as_str() == "Viruses" {

                // If there are (viral) alignment records, also list the reference description 
                // for name and traditional abbreviation for clarity on species name
                for record in &self.evidence.alignment {
                    virus_data.push(format!("{}\n", &record.reference_description.trim_start_matches("bin=")));
                }
            }
        }

        let virus_dedup: Vec<_> = virus_data.iter().dedup().collect();

        for virus_descr in virus_dedup {
            str.push_str(virus_descr)
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
    /// Compute a simple profile-based score.
    ///
    /// `base_weight` scales the log10(total_bases) term to the same order as RPM
    pub fn profile_score(&self, base_weight: f64) -> f64 {
        let total_rpm: f64 = self.profile.iter().map(|r| r.rpm).sum();
        let total_bases: f64 = self.profile.iter().map(|r| r.bases as f64).sum();
        let base_score = total_bases.log10().max(0.0) * base_weight;
        total_rpm + base_score
    }
}

// Lineage in GTDB like format
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


/// Collapse variant suffixes to a single "{Genus species}" name.
/// If `gtdb == true`, only collapse taxa in domains Bacteria or Archaea.
/// If `taxonomy` is provided and a taxon with the collapsed name exists,
/// use its taxid instead of a hashed "collapsed-..." id.
pub fn collapse_taxa(
    taxa: Vec<Taxon>,
    gtdb: bool,
    taxonomy: Option<&GeneralTaxonomy>,
) -> Result<Vec<Taxon>, WorkflowError> {
    // choose grouping key per taxon
    let mut groups: HashMap<String, Vec<Taxon>> = HashMap::new();
    for taxon in taxa {
        let applicable = !gtdb || matches!(taxon.lineage.get_domain().as_deref(), Some("Bacteria" | "Archaea"));
        let key = if applicable { collapse_gs(&taxon.name) } else { taxon.name.clone() };
        groups.entry(key).or_default().push(taxon);
    }

    let mut out = Vec::with_capacity(groups.len());
    for (base_name, group) in groups {
        // representative
        let representative = group.iter().find(|t| t.name == base_name).unwrap_or(&group[0]);
        let rank = &representative.rank;

        // lineage update (replace the component at the rank with the collapsed base name)
        let mut new_lineage = representative.lineage.clone();
        if let Some((idx, prefix)) = rank_index_and_prefix(rank) {
            let mut parts: Vec<String> = new_lineage.split(';').map(|s| s.to_string()).collect();
            if parts.len() <= idx {
                return Err(WorkflowError::LineageStringTooShort(new_lineage));
            }
            parts[idx] = format!("{}{}", prefix, base_name);
            new_lineage = parts.join(";");
        }

        // accumulate evidence (concat align/assembly, sum profile by (id,tool,mode))
        let mut alignment = Vec::new();
        let mut assembly = Vec::new();
        let mut profile_map: HashMap<(String, ProfileTool, AbundanceMode), ProfileRecord> = HashMap::new();

        for taxon in &group {
            alignment.extend(taxon.evidence.alignment.iter().cloned());
            assembly.extend(taxon.evidence.assembly.iter().cloned());
            for r in &taxon.evidence.profile {
                let key = (r.id.clone(), r.tool.clone(), r.mode.clone());
                let e = profile_map.entry(key).or_insert(ProfileRecord {
                    id: r.id.clone(),
                    tool: r.tool.clone(),
                    mode: r.mode.clone(),
                    reads: 0,
                    rpm: 0.0,
                    contigs: 0,
                    bases: 0,
                    bpm: 0.0,
                    abundance: 0.0,
                });
                e.reads += r.reads;
                e.rpm += r.rpm;
                e.contigs += r.contigs;
                e.bases += r.bases;
                e.bpm += r.bpm;
                e.abundance += r.abundance;
            }
        }
        let profile = profile_map.into_values().collect();

        // decide taxid
        let changed_name = base_name != representative.name || group.len() > 1;
        let new_taxid = if changed_name {
            if let Some(tax) = taxonomy {
                taxid_for_name_with_rank(tax, &base_name, rank)
                    .unwrap_or_else(|| {
                        if group.len() == 1 { 
                            representative.taxid.clone() 
                        } else { 
                            log::warn!("Failed to find base name in taxonomy for taxon variant collapse - falling back to hashed identifier for: '{base_name}' 
                            (group representative = {}, n = {})", representative.name, group.len());
                            hashed_collapsed_id(rank, &base_name)
                        }
                    })
            } else {
                hashed_collapsed_id(rank, &base_name)
            }
        } else {
            representative.taxid.clone()
        };

        out.push(Taxon {
            taxid: new_taxid,
            rank: representative.rank.clone(),
            name: base_name,
            lineage: new_lineage,
            evidence: TaxonEvidence { alignment, assembly, profile },
        });
    }

    Ok(out)
}

/// Collapse only the species token.
/// "Haemophilus_A influenzae_CD" -> "Haemophilus_A influenzae"
/// "Haemophilus influenzae_CD"   -> "Haemophilus influenzae"
/// "Haemophilus_A"               -> "Haemophilus_A"  // no species present
fn collapse_gs(name: &str) -> String {
    let mut it = name.split_whitespace();
    let g = it.next().unwrap_or("");
    let s = it.next();

    match s {
        Some(species) if !g.is_empty() => format!("{} {}", g, strip_variant_token(species)),
        _ => name.to_string(),
    }
}

/// Remove a trailing "_[A-Z]+" suffix from a token if present.
fn strip_variant_token(tok: &str) -> &str {
    if let Some(pos) = tok.rfind('_') {
        let suf = &tok[pos + 1..];
        if !suf.is_empty() && suf.chars().all(|c| c.is_ascii_uppercase()) {
            return &tok[..pos];
        }
    }
    tok
}

/// Stable hashed id for collapsed taxa when no taxonomy match is found.
fn hashed_collapsed_id(rank: &TaxRank, base_name: &str) -> String {
    let mut h = DefaultHasher::new();
    format!("collapsed-{rank}-{base_name}").hash(&mut h);
    h.finish().to_string()
}

fn taxid_for_name_with_rank(
    tax: &GeneralTaxonomy,
    name: &str,
    want: &TaxRank,
) -> Option<String> {
    let cands = tax.find_all_by_name(name);
    if cands.is_empty() {
        return None;
    }
    // prefer exact-rank match, else first
    if let Some(tid) = cands.iter().find(|&&tid| tax.rank(tid).ok().as_ref() == Some(want)) {
        return Some((*tid).to_string());
    }
    Some(cands[0].to_string())
}