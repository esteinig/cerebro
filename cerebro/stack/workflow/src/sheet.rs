//! Sample sheet for pipeline


use wax::Glob;
use csv::Writer;
use std::path::Component;
use serde::{Serialize, Deserialize};
use std::{path::{Path, PathBuf}, collections::{HashMap, HashSet}};

use crate::error::WorkflowUtilityError;



#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, PartialOrd)]
pub struct SampleSheetEntry{
    pub sample_id: String,
    pub run_id: String,
    pub run_date: String,
    pub aneuploidy: bool,
    pub comment: Option<String>,
    pub sample_group: Option<String>,
    pub sample_type: Option<String>,
    pub ercc_input: Option<f64>,
    pub forward_path: Option<PathBuf>,  // option because of production sample sheet ionput
    pub reverse_path: Option<PathBuf>
}
impl SampleSheetEntry {
    pub fn new(run_date: &str, run_id: &str, sample_id: &str, sample_group: Option<&str>, sample_type: Option<&str>, ercc_input: Option<f64>, aneuploidy: bool, files: &Vec<PathBuf>) -> Result<Self, WorkflowUtilityError> {
        let (forward_path, reverse_path) = match files.len() {
            1 => (Some(files[0].clone()), None),
            2 => (Some(files[0].clone()), Some(files[1].clone())),
            _ => return Err(WorkflowUtilityError::SampleSheetEntryFiles(format!("{:?}", files)))
        };
        Ok(SampleSheetEntry { 
            sample_id: sample_id.into(), 
            run_date: run_date.into(), 
            run_id: run_id.into(), 
            aneuploidy,
            comment: None,
            sample_group: sample_group.map(str::to_string), 
            sample_type: sample_type.map(str::to_string), 
            ercc_input,
            forward_path, 
            reverse_path 
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, PartialOrd)]
pub struct SampleSheet{
    pub entries: Vec<SampleSheetEntry>
}
impl SampleSheet {
    pub fn new(
        directories: &Vec<PathBuf>, 
        paired_glob: &str, 
        allow_single: bool, 
        run_id: Option<String>, 
        run_date: Option<String>, 
        sample_group: Option<String>, 
        sample_type: Option<String>, 
        ercc_input: Option<f64>, 
        symlinks: bool
    ) -> Result<SampleSheet, WorkflowUtilityError> {
        
        let mut entries = Vec::new();
        for directory in directories {
            let sample_files = get_paired_files(directory, paired_glob, allow_single, symlinks)?;
        
            let entries_directory = sample_files.iter().map(|(sample_id, file_paths)| -> Result<SampleSheetEntry, WorkflowUtilityError> {
                
                let run = match &run_id {
                    Some(id) => id.to_string(),
                    None => {
                        crate::utils::get_file_stem(directory).map_err(
                            |_| WorkflowUtilityError::SampleSheetEntryRunDirectoryPath(format!("Failed to extract path as string: {:?}", &directory))
                        )?
                    }
                };

                let date = match &run_date {
                    Some(date) => date.to_string(),
                    None => chrono::Utc::now().format("%Y-%m-%d").to_string()
                };
    
                SampleSheetEntry::new(
                    &date, 
                    &run, 
                    sample_id, 
                    sample_group.as_deref(), 
                    sample_type.as_deref(),
                    ercc_input,
                    false,
                    file_paths
                )
            }).collect::<Result<Vec<SampleSheetEntry>, WorkflowUtilityError>>()?;

            entries.push(entries_directory);
        } 
        
        let mut entries: Vec<SampleSheetEntry> = entries.into_iter().flatten().collect();

        // Check if all sample identifiers are unique - otherwise spooky things may happen!
        check_unique_sample_id(&entries)?;

        // Sort by sample identifier
        entries.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));

        Ok(Self { entries })
    }
    pub fn from(input: &PathBuf) -> Result<Self, WorkflowUtilityError> {
        let mut reader = csv::ReaderBuilder::new().from_path(input).map_err(|_| WorkflowUtilityError::SampleSheetCsvReader(format!("{:?}", input)))?;
        let mut entries = Vec::new();
        for record in reader.deserialize() {
            entries.push(record?)
        }
        check_unique_sample_id(&entries)?;
        Ok(Self { entries })
    }
    pub fn write(&self, output: &PathBuf) -> Result<(), WorkflowUtilityError> {
        check_unique_sample_id(&self.entries)?; // in case generated without constructor
        let mut writer = Writer::from_path(output).map_err(|_| WorkflowUtilityError::SampleSheetCsvWriter(format!("{:?}", output)))?;
        for entry in &self.entries {
            writer.serialize(entry)?;
        }
        writer.flush()?;
        Ok(())
    }
    pub fn get_run_id(&self, sample_id: &str) -> Option<String> {
        let matches: Vec<&SampleSheetEntry> = self.entries.iter().filter(|record| record.sample_id == sample_id).collect();
        match matches.len() {
            0 => None,
            1 => Some(matches[0].run_id.clone()),
            _ => None // sample_id should be unique
        }
    }     
    pub fn get_run_date(&self, sample_id: &str) -> Option<String> {
        let matches: Vec<&SampleSheetEntry> = self.entries.iter().filter(|record| record.sample_id == sample_id).collect();
        match matches.len() {
            0 => None,
            1 => Some(matches[0].run_date.clone()),
            _ => None // sample_id should be unique
        }
    }  
    pub fn get_sample_group(&self, sample_id: &str) -> Option<String> {
        let matches: Vec<&SampleSheetEntry> = self.entries.iter().filter(|record| record.sample_id == sample_id).collect();
        match matches.len() {
            0 => None,
            1 => match &matches[0].sample_group {
                Some(v) => Some(v.to_owned()),
                None => Some(String::from(""))
            },
            _ => None // sample_id should be unique
        }
    }  
    pub fn get_sample_type(&self, sample_id: &str) -> Option<String> {
        let matches: Vec<&SampleSheetEntry> = self.entries.iter().filter(|record| record.sample_id == sample_id).collect();
        match matches.len() {
            0 => None,
            1 => match &matches[0].sample_group {
                Some(v) => Some(v.to_owned()),
                None => Some(String::from(""))
            },
            _ => None // sample_id should be unique
        }
    }  
    pub fn get_ercc_input(&self, sample_id: &str) -> Option<f64> {
        let matches: Vec<&SampleSheetEntry> = self.entries.iter().filter(|record| record.sample_id == sample_id).collect();
        match matches.len() {
            0 => None,
            1 => matches[0].ercc_input.clone(),
            _ => None // sample_id should be unique
        }
    }  
}

// A helper function to check if sample identifiers are unique
fn check_unique_sample_id(entries: &Vec<SampleSheetEntry>) -> Result<(), WorkflowUtilityError> {
    let n_unique = entries.iter().map(|x| x.sample_id.as_str()).collect::<HashSet<&str>>().len();
    if n_unique != entries.len() {
        return Err(WorkflowUtilityError::SampleSheetIdsNotUnique)
    } 
    Ok(())
}

// A helper function to get paired files from a suitable glob match
fn get_paired_files(directory: &Path, paired_glob: &str, single: bool, symlinks: bool) -> Result<HashMap<String, Vec<PathBuf>>, WorkflowUtilityError> {
   
    let glob = Glob::new(paired_glob).map_err(|_| WorkflowUtilityError::GlobCreate(paired_glob.to_string()))?;

    // Get potentially paired file paths into a HashMap
    let mut paired_files = HashMap::new();
    for entry in glob.walk_with_behavior(directory, match symlinks { true => wax::LinkBehavior::ReadTarget, false => wax::LinkBehavior::ReadFile }) {
        let entry = entry.map_err(|_| WorkflowUtilityError::GlobWalk(format!("{:?}", &directory)))?;
        let file_path = match symlinks {
            true => entry.path().canonicalize()?,
            false => to_lexical_absolute(&entry.path().to_path_buf())?
        };
        let sample_id = entry.matched().get(1).ok_or_else(||WorkflowUtilityError::GlobMatchSampleIdentifier(format!("{:?}", file_path)))?;
        log::debug!("Sample sheet utility - [{:?}] - detected paired-end file: {:?}", sample_id, file_path);
        paired_files.entry(sample_id.to_owned()).or_insert_with(Vec::new).push(file_path.to_path_buf());
    } 

    // Check the entries of the HashMap if single files are not allowed:
    let paired_files = match single {
        false => {
            let mut sorted = HashMap::new();
            for (sample_id, mut file_paths) in paired_files.clone().into_iter() {
                if file_paths.len() != 2 {
                    return Err(WorkflowUtilityError::GlobPairedFiles(sample_id.to_string()))
                }
                // For paired files the entries are unsorted! Sort them here by their names 
                // this will only work for traditional reverse/forward names like R1 and R2
                file_paths.sort();
                sorted.insert(sample_id, file_paths);
            }
            sorted
        },
        true => paired_files
    };
    log::info!("{:?}", paired_files);
    Ok(paired_files)
}


fn to_lexical_absolute(path: &PathBuf) -> std::io::Result<PathBuf> {
    let mut absolute = if path.is_absolute() {
        PathBuf::new()
    } else {
        std::env::current_dir()?
    };
    for component in path.components() {
        match component {
            Component::CurDir => {},
            Component::ParentDir => { absolute.pop(); },
            component @ _ => absolute.push(component.as_os_str()),
        }
    }
    Ok(absolute)
}