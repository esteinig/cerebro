
use cerebro_workflow::error::WorkflowError;
use thiserror::Error;
use cerebro_model::api::cerebro::model::ModelError;

/*
========================
Custom error definitions
========================
*/

#[derive(Error, Debug)]
pub enum HttpClientError {
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
    PingServer(String, String),
    /// Represents failure to upload model to database
    #[error("{0} - {1}")]
    InsertModelResponseFailure(String, String),
    /// Represents failure to retrieve taxa summary data
    #[error("{0} - {1}")]
    TaxaSummaryDataResponseFailure(String, String),
    /// Represents failure to retrieve file registration response
    #[error("{0} - {1}")]
    RegisterFileResponseFailure(String, String),
    /// Represents failure to retrieve file listing response
    #[error("{0} - {1}")]
    ListFilesResponseFailure(String, String),
    /// Represents failure to obtain user teams
    #[error("{0} - {1}")]
    UserTeamsResponseFailure(String, String),
    /// Represents failure to upload model to database
    #[error("{0} - {1}")]
    LoginUserResponseFailure(String, String),
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
    ModelSampleIdentifierEmpty
}