
use cerebro_pipe::error::WorkflowError;
use reqwest::StatusCode;
use thiserror::Error;
use cerebro_model::api::cerebro::model::ModelError;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum HttpClientError {
    #[error(transparent)]
    IoError(#[from] std::io::Error),
    /// Represents failure to get password input from prompt
    #[error("failed to obtain password from prompt")]
    PasswordInput,
    /// Represents failure to get team databases
    #[error("team does not have any databases configured, please contact sysadmin")]
    TeamDatabasesNotFound,
    /// Represents failure to provide a parameter for the insert model function
    #[error("one of `project_id` or `project_name` is required and must be a valid project in the team database ({0})")]
    InsertModelProjectParameter(String),
    /// Represents failure to provide a parameter for the insert model function
    #[error("one of `database_id` or `database_name` is required and must be a valid database for the team ({0})")]
    InsertModelDatabaseParameter(String),
    /// Represents failure to ping server status route
    #[error("failed input/output")]
    IOFailure(#[source] std::io::Error),
    /// Represents failure to ping server status route
    #[error("{0} - {1}")]
    PingServer(StatusCode, String),
    /// Represents failure to process the response of a request made with the client
    #[error("{0}")]
    ResponseFailure(StatusCode),
    /// Represents failure to process the response data of a request made with the client
    #[error("{0} ({1})")]
    DataResponseFailure(StatusCode, String),
    /// Represents all other cases of `reqwest::Error`.
    #[error("failed to make request")]
    ReqwestFailure(#[from] reqwest::Error),
    /// Represents all other cases of `serde_json::Error`.
    #[error("failed to serialize/deserialize data")]
    SerdeFailure(#[from] serde_json::Error),
    /// Represents failure to deserialize a model from file
    #[error("failed to read data model from file")]
    ModelError(#[from] ModelError),
    /// Represents failure to deserialize the taxon filter from file
    #[error("failed to read filter configuration")]
    DeserializeFilter(#[from] WorkflowError),
    /// Represents failure to use a valid sample identifier
    #[error("sample identifier is an empty string")]
    ModelSampleIdentifierEmpty,
    #[error("tower identifier could not be found (--id | --json)")]
    TowerIdentifierNotFound,
    #[error("watcher identifier could not be found (--id | --json)")]
    WatcherIdentifierArgNotFound,
    #[error("failed to send request with 'team' access parameter - did you provide a team name or identifier?")]
    RequireTeamNotConfigured,
    #[error("failed to send request with 'db' access parameter  - did you provide a database name or identifier?")]
    RequireDbNotConfigured,
    #[error("failed to send request with 'project' access parameter  - did you provide a project name or identifier?")]
    RequireProjectNotConfigured,

    #[error("failed to match quality control data identifier ({0}) with pathogen detection data identifier ({0})")]
    IdentifiersNotMatched(String, String)
}