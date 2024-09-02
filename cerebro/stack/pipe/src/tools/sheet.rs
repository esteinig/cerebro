//! Sample sheet for pipeline


use wax::Glob;
use csv::{ReaderBuilder, Writer};
use std::path::Component;
use serde::{Serialize, Deserialize};
use std::{path::{Path, PathBuf}, collections::{HashMap, HashSet}};
use crate::{error::WorkflowError, utils::FileComponent};



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
    #[serde(skip_serializing_if = "Option::is_none")]
    pub forward_path: Option<PathBuf>,  // option because of production sample sheet input
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reverse_path: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fastq: Option<PathBuf>
}
impl SampleSheetEntry {
    pub fn new(
        run_date: &str, 
        run_id: &str, 
        sample_id: &str, 
        sample_group: Option<&str>, 
        sample_type: Option<&str>, 
        ercc_input: Option<f64>, 
        aneuploidy: bool, 
        files: &Vec<PathBuf>
    ) -> Result<Self, WorkflowError> {
        let (forward_path, reverse_path, fastq) = match files.len() {
            1 => (None, None, Some(files[0].clone())),
            2 => (Some(files[0].clone()), Some(files[1].clone()), None),
            _ => return Err(WorkflowError::SampleSheetEntryFiles(format!("{:?}", files)))
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
            reverse_path,
            fastq
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
        prefixes: Option<PathBuf>, 
        paired_glob: &str, 
        allow_single: bool, 
        run_id: Option<String>, 
        run_date: Option<String>, 
        sample_group: Option<String>, 
        sample_type: Option<String>, 
        ercc_input: Option<f64>, 
        symlinks: bool,
        recursive: bool,
        not_unique: bool
    ) -> Result<SampleSheet, WorkflowError> {
        

        let mut entries = Vec::new();
        for directory in directories {
            let sample_files = get_files(
                directory, 
                paired_glob, 
                allow_single, 
                symlinks, 
                recursive, 
                prefixes.clone()
            )?;
        
            let entries_directory = sample_files.iter().map(|(sample_id, file_paths)| -> Result<SampleSheetEntry, WorkflowError> {
                
                let run = match &run_id {
                    Some(id) => id.to_string(),
                    None => {
                        crate::utils::get_file_component(directory, FileComponent::FileStem).map_err(
                            |_| WorkflowError::SampleSheetEntryRunDirectoryPath(format!("Failed to extract path as string: {:?}", &directory))
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
            }).collect::<Result<Vec<SampleSheetEntry>, WorkflowError>>()?;

            entries.push(entries_directory);
        } 
        
        let mut entries: Vec<SampleSheetEntry> = entries.into_iter().flatten().collect();

        // Check if all sample identifiers are unique - otherwise spooky things may happen!
        if !not_unique {
            check_unique_sample_id(&entries)?;
        }

        // Sort by sample identifier
        entries.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));

        Ok(Self { entries })
    }
    pub fn from(input: &PathBuf) -> Result<Self, WorkflowError> {
        let mut reader = csv::ReaderBuilder::new().from_path(input).map_err(|_| WorkflowError::SampleSheetCsvReader(format!("{:?}", input)))?;
        let mut entries = Vec::new();
        for record in reader.deserialize() {
            entries.push(record?)
        }
        
        check_unique_sample_id(&entries)?;

        Ok(Self { entries })
    }
    pub fn write(&self, output: &PathBuf) -> Result<(), WorkflowError> {
        let mut writer = Writer::from_path(output).map_err(|_| WorkflowError::SampleSheetCsvWriter(format!("{:?}", output)))?;
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

// A helper function to read the prefix file
fn get_prefixes(path: &PathBuf) -> Result<Vec<String>, WorkflowError> {

    let mut reader = ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(path)?;

    let mut prefixes = Vec::new();
    for record in reader.records() {
        let rec = record?;
        prefixes.push(rec.as_slice().to_string())
    }

    Ok(prefixes)
}

fn check_sample_prefix(prefixes: &Vec<String>, sample_id: &str) -> (bool, String) {
    for prefix in prefixes {
        if sample_id.starts_with(prefix) {
            return (true, prefix.to_owned())
        }
    }
    return (false, String::new())
}

fn sample_prefix_summary(prefixes: &Vec<String>, detected: &Vec<String>) -> Vec<(String, usize)> {
    
    let mut counts = Vec::new();
    for prefix in prefixes {
        let found = detected.iter().filter(|p| p == &prefix).collect::<Vec<_>>();
        let num_found = found.len();
        counts.push((prefix.to_owned(), num_found));
    }
    return counts
}


// A helper function to check if sample identifiers are unique
fn check_unique_sample_id(entries: &Vec<SampleSheetEntry>) -> Result<(), WorkflowError> {
    let n_unique = entries.iter().map(|x| x.sample_id.as_str()).collect::<HashSet<&str>>().len();
    if n_unique != entries.len() {
        return Err(WorkflowError::SampleSheetIdsNotUnique)
    } 
    Ok(())
}

// A helper function to get paired files from a suitable glob match
pub fn get_files(directory: &Path, paired_glob: &str, single: bool, symlinks: bool, recursive: bool, prefixes: Option<PathBuf>) -> Result<HashMap<String, Vec<PathBuf>>, WorkflowError> {
   
    let prefix_subset = match prefixes {
        Some(ref path) => get_prefixes(path)?,
        None => Vec::new()
    };

    let glob = Glob::new(paired_glob).map_err(|_| WorkflowError::GlobCreate(paired_glob.to_string()))?;

    let mut walk_behaviour = wax::WalkBehavior::default();

    walk_behaviour.link = if symlinks { wax::LinkBehavior::ReadTarget } else { wax::LinkBehavior::ReadFile };
    walk_behaviour.depth = if recursive { usize::MAX } else { 1 };

    // Get potentially paired file paths into a HashMap
    let mut paired_files = HashMap::new();
    let mut prefix_detected = Vec::new();
    for entry in glob.walk_with_behavior(directory, walk_behaviour) {
        let entry = entry.map_err(|_| WorkflowError::GlobWalk(format!("{:?}", &directory)))?;
        
        let file_path = match symlinks {
            true => entry.path().canonicalize()?,
            false => to_lexical_absolute(&entry.path().to_path_buf())?
        };
        let sample_id = entry.matched().get(1).ok_or_else(||WorkflowError::GlobMatchSampleIdentifier(format!("{:?}", file_path)))?;
        
        log::debug!("Sample sheet utility - [{:?}] - detected paired-end file: {:?}", sample_id, file_path);

        if prefix_subset.is_empty() {

            paired_files.entry(sample_id.to_owned()).or_insert_with(Vec::new).push(file_path.to_path_buf());
        } else {
            let (include, prefix) = check_sample_prefix(&prefix_subset, sample_id); 
            if include {
                paired_files.entry(sample_id.to_owned())
                    .or_insert_with(Vec::new)
                    .push(file_path.to_path_buf());

                prefix_detected.push(prefix)
            }
        }
    } 

    // Check the entries of the HashMap if single files are not allowed:
    let paired_files = match single {
        false => {
            let mut sorted = HashMap::new();
            for (sample_id, mut file_paths) in paired_files.clone().into_iter() {
                if file_paths.len() != 2 {
                    return Err(WorkflowError::GlobPairedFiles(sample_id.to_string()))
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

    log::info!("{:#?}", paired_files);

    if let Some(_) = prefixes {
        let counts = sample_prefix_summary(&prefix_subset, &prefix_detected);
        for (prefix, count) in counts {
            log::info!("Prefix '{prefix}' => {count}")
        }
    }
    
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