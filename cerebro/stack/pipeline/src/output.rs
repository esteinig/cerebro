/* 
===============
Workflow Sample
===============
*/

// v0.3.0

use std::path::PathBuf;
use crate::error::WorkflowError;

// Workflow sample aggregate which represents the output of the workflow module 
// for a single processed sample (distinct from biological sample, which is a set of WorkflowSamples)
pub struct WorkflowOutputs {
    pub nanoq: Option<PathBuf>,
    pub nanoq_scan: Option<PathBuf>,
    pub fastp: Option<PathBuf>,
    pub fastp_scan: Option<PathBuf>,
    pub ercc_vircov: Option<PathBuf>,
    pub ercc_scrubby: Option<PathBuf>,
    pub phage_vircov: Option<PathBuf>,
    pub phage_scrubby: Option<PathBuf>,
    pub host_scrubby: Option<PathBuf>,
    pub virus_scrubby: Option<PathBuf>,
    pub alignment: Option<Vec<PathBuf>>,
    pub alignment_remap: Option<Vec<PathBuf>>,
    pub kraken2uniq: Option<Vec<PathBuf>>,
    pub blastn: Option<Vec<PathBuf>>,
    pub diamond: Option<Vec<PathBuf>>,
}
impl WorkflowOutputs {
    pub fn from_results(path: &PathBuf) -> Result<Self, WorkflowError> {
        Ok(Self{
            nanoq: WorkflowOutputs::get_file_by_name(&path, "qc__nanoq__reads")?,
            nanoq_scan: WorkflowOutputs::get_file_by_name(&path, "qc__nanoq__scan")?,
            fastp: WorkflowOutputs::get_file_by_name(&path, "qc__fastp__reads")?,
            fastp_scan: WorkflowOutputs::get_file_by_name(&path, "qc__fastp__scan")?,
            ercc_vircov: WorkflowOutputs::get_file_by_name(&path, "qc__vircov__ercc")?,
            ercc_scrubby: WorkflowOutputs::get_file_by_name(&path, "qc__scrubby__ercc")?,
            phage_vircov: WorkflowOutputs::get_file_by_name(&path, "qc__vircov__phage")?,
            phage_scrubby: WorkflowOutputs::get_file_by_name(&path, "qc__scrubby__phage")?,
            host_scrubby: WorkflowOutputs::get_file_by_name(&path, "qc__scrubby__host")?,
            virus_scrubby: WorkflowOutputs::get_file_by_name(&path, "qc__scrubby__virus_background")?,
            alignment: WorkflowOutputs::get_files_from_patterns(&path, &["align__vircov__*"])?,
            alignment_remap: WorkflowOutputs::get_files_from_patterns(&path, &["align__vircov_remap__*"])?,
            kraken2uniq: WorkflowOutputs::get_files_from_patterns(&path, &["kmer__kraken2uniq__*"])?,
            blastn: WorkflowOutputs::get_files_from_patterns(&path, &["assembly__blastn__*"])?,
            diamond: WorkflowOutputs::get_files_from_patterns(&path, &["assembly__diamond__*"])?,
        })
    }
    pub fn get_file_by_name(path: &PathBuf, filename: &str) -> Result<Option<PathBuf>, WorkflowError> {
        let file_path = path.join(filename);
        match file_path.exists() { true => Ok(Some(file_path)), false => Ok(None) }
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
}