use core::fmt;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::WorkflowError;


/// Enum representing the available classifiers.
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum Classifier {
    #[serde(rename="kraken2")]
    Kraken2,
    #[serde(rename="metabuli")]
    Metabuli,
}


/// Enum representing the available profilers.
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum Profiler {
    #[serde(rename="sylph")]
    Sylph,
    #[serde(rename="kmcp")]
    Kmcp,
    #[serde(rename="bracken")]
    Bracken,

}

/// Enum representing the available aligners.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum Aligner {
    #[serde(rename="bowtie2")]
    Bowtie2,
    #[serde(rename="minimap2")]
    Minimap2,
    #[serde(rename="strobealign")]
    Strobealign,
}


impl Aligner {
    // Used for identification of pre-built-indices
    pub fn short_name(&self) -> &str {
        match self {
            Aligner::Bowtie2 => "bt2",
            Aligner::Minimap2 => "mm2",
            Aligner::Strobealign => "sti",
        }
    }
}
impl fmt::Display for Aligner {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Aligner::Bowtie2 => write!(f, "bowtie2"),
            Aligner::Minimap2 => write!(f, "minimap2"),
            Aligner::Strobealign => write!(f, "strobealign"),
        }
    }
}


impl Classifier {
    // Used for identification of pre-built-indices
    pub fn short_name(&self) -> &str {
        match self {
            Classifier::Kraken2 => "kraken2",
            Classifier::Metabuli => "metabuli",
        }
    }
}
impl fmt::Display for Classifier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Classifier::Kraken2 => write!(f, "kraken2"),
            Classifier::Metabuli => write!(f, "metabuli"),
        }
    }
}


impl Profiler {
    // Used for identification of pre-built-indices
    pub fn short_name(&self) -> &str {
        match self {
            Profiler::Sylph => "sylph",
            Profiler::Kmcp => "kmcp",
            Profiler::Bracken => "bracken",
        }
    }
}
impl fmt::Display for Profiler {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Profiler::Sylph => write!(f, "sylph"),
            Profiler::Bracken => write!(f, "bracken"),
            Profiler::Kmcp  =>  write!(f, "kmcp")
        }
    }
}


struct KrakenReportRecord {}
struct MetabuliReportRecord {}
struct KmcpReportRecord {}
struct BrackenReportRecord {}


pub struct KrakenReport {
    records: Vec<KrakenReportRecord>
}
impl KrakenReport {
    pub fn from_report(path: &Path) -> Result<Self, WorkflowError> {
        Ok(Self { records: Vec:: new() })
    }
}

pub struct MetabuliReport {
    records: Vec<MetabuliReportRecord>
}
impl MetabuliReport {
    pub fn from_report(path: &Path) -> Result<Self, WorkflowError> {
        Ok(Self { records: Vec:: new() })
    }
}

pub struct KmcpReport {
    records: Vec<KmcpReportRecord>
}
impl KmcpReport {
    pub fn from_report(path: &Path) -> Result<Self, WorkflowError> {
        Ok(Self { records: Vec:: new() })
    }
}

pub struct BrackenReport {
    records: Vec<BrackenReportRecord>
}
impl BrackenReport {
    pub fn from_report(path: &Path) -> Result<Self, WorkflowError> {
        Ok(Self { records: Vec:: new() })
    }
}

