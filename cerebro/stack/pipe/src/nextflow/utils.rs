use std::path::PathBuf;

use crate::error::WorkflowError;

pub fn get_file_by_name(path: &PathBuf, id: &str, extension: &str) -> Result<Option<PathBuf>, WorkflowError> {
    
    let file_path = path.join(format!("{id}{extension}"));
    
    match file_path.exists() && file_path.is_file() { 
        true => Ok(Some(file_path)), 
        false => Ok(None) 
    }
}
pub fn get_files_from_patterns(path: &PathBuf, patterns: &[&str]) -> Result<Option<Vec<PathBuf>>, WorkflowError> {
    
    // Need to have canonical path for GlobWalk
    let full_path = std::fs::canonicalize(path).map_err(|_| WorkflowError::InvalidReferencePath)?;

    let walker = globwalk::GlobWalkerBuilder::from_patterns(
        full_path, patterns
    )
        .max_depth(1)
        .follow_links(false)
        .build().map_err(|_| WorkflowError::GlobWalkBuilder)?
        .into_iter()
        .filter_map(Result::ok);
    
    let mut alignments = Vec::new();
    for file in walker {
        alignments.push(file.path().to_path_buf())
    };
    match alignments.is_empty() {
        true => Ok(None),
        false => Ok(Some(alignments))
    }
}