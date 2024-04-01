use thiserror::Error;

#[derive(Error, Debug)]
pub enum ReportError {
    #[error("Failed to open configuration file")]
    TomlConfigFile(#[source] std::io::Error),
    #[error("Failed to parse configuration")]
    TomlConfigParse(#[from] toml::de::Error),
}
