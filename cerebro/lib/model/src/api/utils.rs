use serde::Deserialize;

// Generic error response
#[derive(Debug, Deserialize)]
pub struct ErrorResponse {
    pub status: String,
    pub message: String,
}
