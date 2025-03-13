use serde::{Serialize, Deserialize};

use crate::api::files::model::FileType;

use super::schema::RegisterWatcherSchema;

#[derive(Debug, Clone, Serialize, Deserialize, clap::ValueEnum)]
pub enum WatcherFormat {
    Fastq,
    FastqPe,
    Iseq,
    Nextseq
}
impl WatcherFormat {
    pub fn default_glob(&self) -> String {
        match self {
            Self::Fastq => String::from("*.fastq.gz"),
            Self::FastqPe => String::from("*_{R1,R2}.fastq.gz"),
            Self::Iseq => String::from("*_{L001_R1_001,L001_R2_001}.fastq.gz"),
            Self::Nextseq => String::from("*_{R1_001,R2_001}.fastq.gz")
        }
    }
    pub fn file_type(&self) -> FileType {
        match self {
            Self::Fastq => FileType::ReadSingle,
            Self::FastqPe => FileType::ReadPaired,
            Self::Iseq =>  FileType::ReadPaired,
            Self::Nextseq => FileType::ReadPaired
        }
    }
}

impl std::fmt::Display for WatcherFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fastq => write!(f, "fastq"),
            Self::FastqPe => write!(f, "fastq-pe"),
            Self::Iseq => write!(f, "iseq"),
            Self::Nextseq => write!(f, "nextseq")
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProductionWatcher {
    pub id: String,
    pub date: String,
    pub name: String,
    pub location: String,
    pub format: WatcherFormat,
    pub glob: String,
    pub last_ping: String,
}
impl ProductionWatcher {
    pub fn from_schema(schema: &RegisterWatcherSchema) -> Self {
        Self {
            id: schema.id.clone(),
            date: schema.date.clone(),
            name: schema.name.clone(),
            location: schema.location.clone(),
            format: schema.format.clone(),
            glob: schema.glob.clone(),
            last_ping: schema.last_ping.clone()
        }
    }
}

