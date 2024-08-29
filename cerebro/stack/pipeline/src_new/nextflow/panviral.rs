use super::utils::get_file_by_name;

pub struct PanviralOutput {
    pub qc: Option<PathBuf>,
    pub host: Option<PathBuf>,
    pub viruses: Option<PathBuf>,
    pub controls: Option<PathBuf>,
    
}
impl PanviralOutput {
    pub fn from_results(path: &PathBuf, id: &str) -> Result<Self, WorkflowError> {
        Ok(Self{
            qc: get_file_by_name(&path, &id, ".qc.json")?,
            host: get_file_by_name(&path, &id, ".host.json")?,
            viruses: get_file_by_name(&path, &id, ".viruses.tsv")?,
            controls: get_file_by_name(&path, &id, ".controls.tsv")?,
        })
    }

}