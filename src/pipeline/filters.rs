use taxonomy::TaxRank;
use crate::pipeline::taxon::Taxon;
use serde::{Serialize, Deserialize};


// Taxon rank filters (direct rank)
pub fn filter_by_rank(taxa: Vec<Taxon>, rank: TaxRank) -> Vec<Taxon> {
    taxa.into_iter().filter(|taxon| taxon.rank == rank).collect()
}
pub fn filter_by_rank_name(taxa: Vec<Taxon>, name: String) -> Vec<Taxon> {
    taxa.into_iter().filter(|taxon| taxon.name == name).collect()
}
pub fn filter_by_rank_taxid(taxa: Vec<Taxon>, taxid: String) -> Vec<Taxon> {
    taxa.into_iter().filter(|taxon| taxon.taxid == taxid).collect()
}

// Taxon level filters (parent rank - species, genus or domain) 
//
// This implementation is necessary as the parent levels are always linked (if a name is present, a taxid is present)
pub fn filter_by_parent(taxa: Vec<Taxon>, level: &str, name: Option<String>, taxid: Option<String>) -> Vec<Taxon> {
    match (name.clone(), taxid.clone()) {
        (Some(_), Some(_)) => {
            taxa.into_iter().filter(|taxon| {
                match level.to_lowercase().as_str() {
                    "domain" => taxon.level.domain_name == name && taxon.level.domain_taxid == taxid,
                    "genus" => taxon.level.genus_name == name && taxon.level.genus_taxid == taxid,
                    "species" => taxon.level.species_name == name && taxon.level.species_taxid == taxid,
                    _ => false // TODO: Error implementation
            }}).collect()
        },
        (Some(_), None) => {
            taxa.into_iter().filter(|taxon| {
                match level.to_lowercase().as_str() {
                    "domain" => taxon.level.domain_name == name,
                    "genus" => taxon.level.genus_name == name,
                    "species" => taxon.level.species_name == name,
                    _ => false // TODO: Error implementation
                
            }}).collect()
        },
        (None, Some(_)) => {
            taxa.into_iter().filter(|taxon| {
                match level.to_lowercase().as_str() {
                    "domain" => taxon.level.domain_taxid == taxid,
                    "genus" => taxon.level.genus_taxid == taxid,
                    "species" => taxon.level.species_taxid == taxid,
                    _ => false // TODO: Error implementation
                
            }}).collect()
        },
        (None, None) => {
            taxa.into_iter().filter(|taxon| {
                match level.to_lowercase().as_str() {
                    "domain" => taxon.level.domain_name == None && taxon.level.domain_taxid == None,
                    "genus" => taxon.level.genus_name == None && taxon.level.genus_taxid == None,
                    "species" => taxon.level.species_name == None && taxon.level.species_taxid == None,
                    _ => false // TODO: Error implementation
                
            }}).collect()
        }
    }
}

// Evidence filters - add here when integrating new tool
pub fn filter_by_evidence(taxa: Vec<Taxon>) -> Vec<Taxon> {
    taxa.into_iter().filter(|taxon| !(taxon.evidence.kmer.is_empty()) || !(taxon.evidence.alignment.is_empty()) || !(taxon.evidence.assembly.is_empty())).collect()
}

pub fn filter_by_tags(taxa: Vec<Taxon>, tags: &Vec<String>) -> Vec<Taxon> {
    taxa.into_iter().map(|mut taxon| {

        // NEEDS TO BE CHECKED TO FILTER BY EVIDENCE TAGS - FOR NOW FILTER BY ANY OCCURRENCE OF THE TAGS IN TAXON EVIDENCE

        let mut filter_all_evidence = false;

        for record in &taxon.evidence.kmer {
            if tags.contains(&record.id){
                filter_all_evidence = true;
            }
        }
        for record in &taxon.evidence.alignment {
            if tags.contains(&record.id){
                filter_all_evidence = true;
            }
        }
        for record in &taxon.evidence.assembly {
            if tags.contains(&record.id){
                filter_all_evidence = true;
            }
        }

        if filter_all_evidence {
            taxon.evidence.kmer = Vec::new();
            taxon.evidence.alignment = Vec::new();
            taxon.evidence.assembly = Vec::new();
        }

        // if !taxon.evidence.kmer.clone().into_iter().filter(|record| tags.contains(&record.id)).collect::<Vec<Kraken2UniqRecord>>().is_empty() {
        //     taxon.evidence.kmer = Vec::new();
        // }
        // if !taxon.evidence.alignment.clone().into_iter().filter(|record| tags.contains(&record.id)).collect::<Vec<VircovRecord>>().is_empty() {
        //     taxon.evidence.alignment = Vec::new();
        // }
        // if !taxon.evidence.assembly.clone().into_iter().filter(|record| tags.contains(&record.id)).collect::<Vec<BlastLcaRecord>>().is_empty() {
        //     taxon.evidence.assembly = Vec::new();
        // }
        
        taxon
    }).collect()
}

pub fn filter_by_kraken2uniq(taxa: Vec<Taxon>, min_reads: u64, databases: Vec<String>) -> Vec<Taxon> {

    taxa.into_iter().map(|mut taxon| {
        taxon.evidence.kmer = taxon.evidence.kmer.into_iter().filter(|record| record.reads >= min_reads).collect();
        
        if !databases.is_empty() {
            taxon.evidence.kmer = taxon.evidence.kmer.into_iter().filter(|record| databases.contains(&record.db)).collect();
        }
        
        taxon
    }).collect()
}


pub fn filter_by_vircov_scan(taxa: Vec<Taxon>, min_reads: u64, min_bases: u64, min_regions: u64, min_coverage: f64, ref_length: u64) -> Vec<Taxon> {
    taxa.clone().into_iter().map(|mut taxon| {
        taxon.evidence.alignment = taxon.evidence.alignment.into_iter().filter(|record| {
            record.scan_reads >= min_reads && 
            record.scan_bases_covered >= min_bases &&
            record.scan_regions >= min_regions && 
            record.scan_coverage >= min_coverage && 
            record.reference_length >= ref_length
    }).collect();
        taxon
    }).collect()
}

pub fn filter_by_blast(taxa: Vec<Taxon>, min_contig_length: u64, min_coverage: f64, min_identity: f64) -> Vec<Taxon> {
    taxa.clone().into_iter().map(|mut taxon| {
        taxon.evidence.assembly = taxon.evidence.assembly.into_iter().filter(|record| {
            record.length >= min_contig_length && 
            record.coverage >= min_coverage &&
            record.identity >= min_identity
    }).collect();
        taxon
    }).collect()
}



#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonFilterSettings {
    pub parent_rank: TaxonRankFilterSettings,
    pub parent_rank_exclude: TaxonRankExcludeFilterSettings,
    pub kraken2uniq: Kraken2UniqFilterSettings,
    pub alignment: AlignmentFilterSettings

}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonRankFilterSettings {
    pub domain: Option<String>,
    pub genus: Option<String>,
    pub species: Option<String>
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonRankExcludeFilterSettings {
    pub domain: Option<Vec<String>>,
    pub genus: Option<Vec<String>>,
    pub species: Option<Vec<String>>
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Kraken2UniqFilterSettings {
    pub min_reads: u64 
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentFilterSettings {
    pub min_reads: u64,
    pub min_alignments: u64,
    pub min_bases: u64,
    pub min_regions: u64,
    pub min_coverage: f64
}
impl TaxonFilterSettings {
    pub fn new() -> Self {
        Self {
            parent_rank: TaxonRankFilterSettings { domain: None, genus: None, species: None },
            parent_rank_exclude: TaxonRankExcludeFilterSettings { domain: None, genus: None, species: None },
            kraken2uniq: Kraken2UniqFilterSettings { min_reads: 0 },
            alignment: AlignmentFilterSettings { min_reads: 0, min_alignments: 0, min_bases: 0, min_regions: 0, min_coverage: 0.0 }
        }
    }
}
