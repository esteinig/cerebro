use anyhow::Result;
use taxonomy::{GeneralTaxonomy, Taxonomy, TaxRank};
use std::{path::{Path, PathBuf}, collections::HashMap, io::BufReader, io::{BufRead, BufWriter}, fs::File};
use needletail::{parse_fastx_file, Sequence};
use std::str::from_utf8;
use crate::tools::error::ToolError;


pub struct Annotate {
    acc2taxid: HashMap<String, String>
}

impl Annotate {
    pub fn new(paths: &[PathBuf]) -> Result<Self, ToolError> {
        // Read the accession to taxonomic identifier file - can be very large
        let mut acc2taxid = HashMap::new();
        for path in paths {
            log::info!("Loading acc2taxid file: {:?}", path);
            let map_file = BufReader::new(File::open(path)?);
            for line in map_file.lines(){
                let data = line?;
                let content: Vec<&str> = data.split('\t').collect();
                acc2taxid.insert(content[0].to_owned(), content[2].to_owned());
            }
        }
        Ok(Self { acc2taxid })
    }

    pub fn annotate(
        &self,
        input: &Path,
        output: &Path,
        exclude: &Vec<String>,
        taxonomy: &Option<GeneralTaxonomy>
    ) -> Result<(), ToolError>{
        // Input output of read files includes compression detection
        let mut reader = parse_fastx_file(input)?;
        let file = File::create(&output)?;
        let mut writer = BufWriter::new(file);

        let mut total = 0;
        let mut excluded = 0;
        'seqs: while let Some(record) = reader.next() {
            let rec = record?;
            let mut header = from_utf8(rec.id())?.to_string();
            let rec_id = header.split(' ').next().unwrap_or("").to_string(); // needletail parses the entire header as identifier (including description)
            let rec_id_clean = rec_id.split('.').next().unwrap_or("").to_string();
            
            total +=1;

            for term in exclude {
                if header.contains(term){
                    log::debug!("Excluding record {} due to presence of term: {}", rec_id_clean, term);
                    continue 'seqs;
                }
            }
            match self.acc2taxid.get(&rec_id_clean) {
                Some(taxid) => {

                    let mut found_taxid = taxid.to_owned();

                    if let Some(_taxonomy) = taxonomy {
                        
                        let species_search = _taxonomy.parent_at_rank(
                            found_taxid.as_str(), TaxRank::Species
                        ).unwrap_or(None);  // TODO: update bacterial database current with tax-dump, otherwise this will revert to old identifier

                        let (species_taxid, _) = get_taxid_name_from_search(species_search, _taxonomy)?;
                
                        if let Some(stax) = species_taxid {
                            found_taxid = stax.clone();
                            log::debug!("Found species taxid for accession: {} ({})", rec_id_clean, found_taxid);
                        } else {
                            log::warn!("Could not find species taxid for accession: {} ({})", rec_id_clean, found_taxid);
                        }
                                
                        
                    }

                    header += format!("; taxid={}", found_taxid).as_str();
                    needletail::parser::write_fasta(header.as_bytes(), rec.sequence(), &mut writer, needletail::parser::LineEnding::Unix)?;
                },
                None => {
                    log::warn!("Could not find taxid for accession: {}", rec_id_clean);
                    excluded += 1
                }
            };
        }
        log::info!("Excluded {}/{} sequences without traceable taxonomic identifier", excluded, total);
        Ok(())
    }

}

// Utility function wrapping the search result handling with name and taxid
pub fn get_taxid_name_from_search(search_result: Option<(&str, f32)>, taxonomy: &GeneralTaxonomy) -> Result<(Option<String>, Option<String>), ToolError> {
    match search_result {
       Some((taxid, _)) => {
           let name = taxonomy.name(taxid).map_err(
               |err|ToolError::TaxNameNotAvailable(err, taxid.to_owned())
           )?;
           Ok((Some(taxid.to_owned()), Some(name.to_owned())))

       },
       None => Ok((None, None))
   }
}