use core::fmt;

use serde::{Deserialize, Serialize};


/// Enum representing the available classifiers.
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum Classifier {
    #[serde(rename="kraken2")]
    Kraken2,
    #[serde(rename="metabuli")]
    Metabuli,
    #[serde(rename="sylph")]
    Sylph,
    #[serde(rename="kmcp")]
    Kmcp,
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
            Classifier::Sylph => "sylph",
            Classifier::Kmcp => "kmcp",
        }
    }
}
impl fmt::Display for Classifier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Classifier::Kraken2 => write!(f, "k2"),
            Classifier::Metabuli => write!(f, "mt"),
            Classifier::Sylph => write!(f, "sy"),
            Classifier::Kmcp => write!(f, "kp"),
        }
    }
}