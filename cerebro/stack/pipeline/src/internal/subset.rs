use anyhow::Result;
use crate::error::WorkflowError;
use crate::utils::get_compression_writer;
use needletail::parse_fastx_file;
use serde::Serialize;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::{PathBuf, Path};
use std::str::from_utf8;
use std::collections::{HashSet, HashMap};
use std::collections::hash_map::Entry;

/*
=================
Fasta subset data
=================
*/

#[derive(Serialize)]
pub struct SubsetData {
    pub total: u64,
    pub depleted: u64,
    pub retained: u64,
    pub input_file: PathBuf,
    pub output_file: PathBuf,
}

/*
=================
Fasta subsetter
=================
*/


pub struct FastaSubset {
    pub mash_records: Option<Vec<MashScreenRecord>>
}

impl FastaSubset {
    pub fn from_mash(path: &Path, min_identity: &f64, min_shared_hashes: &u32) -> Result<Self, WorkflowError> {

        let reader = BufReader::new(File::open(path)?);
        
        let mut mash_records: Vec<MashScreenRecord> = Vec::new();

        for result in reader.lines() {
            let record: MashScreenRecord = MashScreenRecord::from_str(result?)?;
            if record.identity >= *min_identity && record.shared_hashes >= *min_shared_hashes 
            {   
                mash_records.push(record);
            }
        }

        Ok(Self { mash_records: Some(mash_records)  } )
    }
    pub fn select_grouped(&self, group_index: &Option<usize>, group_by: &Option<String>, group_sep: &String) -> Result<HashSet<String>, WorkflowError> {

        match &self.mash_records {
            Some(mash_data) => {
            
                let grouped = match group_index {
                    Some(group_idx) => {

                        let mut grouped: HashMap<String, Vec<MashScreenRecord>> = HashMap::new();
                        for record in mash_data {
                            // Get the group identifier from the record query identifier
                            let group_id = record.query_id.split(group_sep).nth(*group_idx);
                            // Group the records by unique group identifiers
                            match group_id {
                                Some(value) => {
        
                                    match grouped.entry(value.to_string()) {
                                        Entry::Occupied(mut entry) => {
                                            entry.get_mut().push(record.clone());
                                        }
                                        Entry::Vacant(entry) => {
                                            entry.insert(vec![record.clone()]);
                                        }
                                    }
        
                                },
                                None => return Err(WorkflowError::GroupIndexError(group_idx.to_string(), record.query_id.clone()))
                            }
                        }
                        Some(grouped)
                    },
                    None => {
                        match group_by {
                            Some(group_field) => {

                                let mut grouped: HashMap<String, Vec<MashScreenRecord>> = HashMap::new();
                                for record in mash_data {
                                    // Get the group identifier from the record query identifier
                                    let group_id: Vec<Option<&str>> = record.query_comment.split(group_sep).map(
                                        |x| match x.starts_with(group_field) {
                                            true => x.strip_prefix(group_field),
                                            false => None
                                        }).filter(|x| x.is_none()).collect();
                                    
                                        if group_id.is_empty() {
                                            return Err(WorkflowError::GroupIndexError("no id field".to_string(), record.query_id.clone()))
                                        }
                                        
                                    // Group the records by unique group 
                                    match group_id[0] {
                                        Some(value) => {
                
                                            match grouped.entry(value.to_string()) {
                                                Entry::Occupied(mut entry) => {
                                                    entry.get_mut().push(record.clone());
                                                }
                                                Entry::Vacant(entry) => {
                                                    entry.insert(vec![record.clone()]);
                                                }
                                            }
                
                                        },
                                        None => return Err(WorkflowError::GroupIndexError("no id field value".to_string(), record.query_id.clone()))
                                    }
                                }
                                Some(grouped)
                            },
                            None => None
                        }
                        
                    }
                };

            // Select the best record by a selection criterion and
            // add to a HashSet of query identifiers
            match grouped {
                Some(grouped_data) => {
                    let mut seq_ids: HashSet<String> = HashSet::new();
                    for (group, records) in grouped_data {
                        let record = records.into_iter().max_by_key(|r| r.shared_hashes);
                        match record {
                            Some(selected) => {
                                seq_ids.insert(selected.query_id);
                            },
                            None => return Err(WorkflowError::GroupSelectionError(group))
                        }
                    }
                    Ok(seq_ids)
                },
                None =>  {
                    let mut seq_ids: HashSet<String> = HashSet::new();
                    for record in mash_data {
                        seq_ids.insert(record.query_id.clone());  // all seq ids with shared hashes from mash if not grouped
                    }
                    Ok(seq_ids)
                }
            }

            },
            None => Ok(HashSet::new())
        }

    }
    pub fn get_record_ids(&self, group_index: &Option<usize>, group_by: &Option<String>, group_sep: &String) -> Result<HashSet<String>, WorkflowError> {
        Ok(self.select_grouped(group_index, group_by, group_sep)?)
    }
    pub fn subset(&self, seq_ids: HashSet<String>, fasta: &PathBuf, output: &PathBuf, output_format: &Option<niffler::compression::Format>, compression_level: &niffler::compression::Level) -> Result<(), WorkflowError> {
        // Input output of read files includes compression detection
        let mut reader = parse_fastx_file(fasta)?;
        let mut writer = get_compression_writer(output, output_format, &Some(*compression_level))?;
        
        while let Some(record) = reader.next() {
            let rec = record?;
            let rec_id = from_utf8(rec.id())?.split(' ').next().unwrap_or("").to_string(); // needletail parses the entire header as identifier (including description)
            
            if seq_ids.contains(&rec_id) {
                rec.write(&mut writer, None)?;
            }
        }
        Ok(())
    }

}

/*
=================
Mash screen data
=================
*/
#[derive(Debug, Clone)]
pub struct MashScreenRecord {
    pub identity: f64,
    pub shared_hashes: u32,
    pub total_hashes: u32,
    pub median_multiplicity: u32,
    pub p_value: f64,
    pub query_id: String,
    pub query_comment: String
}

impl MashScreenRecord {
    pub fn from_str(line: String) -> Result<Self, WorkflowError>{
        let fields: Vec<&str> = line.split('\t').collect();
        let hash_data: Vec<&str> = fields[1].split('/').collect();

        let record = Self {
            identity: fields[0].parse::<f64>()?,
            shared_hashes: hash_data[0].parse::<u32>()?,
            total_hashes: hash_data[1].parse::<u32>()?,
            median_multiplicity: fields[2].parse::<u32>()?,
            p_value: fields[3].parse::<f64>()?,
            query_id: fields[4].to_string(),
            query_comment: fields[5].to_string()
        };

        Ok(record)
    }
}
