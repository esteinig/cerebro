use crate::api::files::model::SeaweedFile;

pub struct CiqaDataset {
    pub files: Vec<SeaweedFile>,
    pub description: String,
    pub version: String,
    pub qa_thresholds: CiqaThresholds,
    pub dataset_prefix: String
}

pub struct CiqaThresholds {
    pub sensitivity: f64,
    pub specificty: f64
}
