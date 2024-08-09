use serde::Deserialize;

// Generic error response
#[derive(Debug, Deserialize)]
pub struct ErrorResponse {
    pub status: String,
    pub message: String,
}

pub enum HttpMethod {
    Get,
    Post,
    Delete,
    Put,
}
impl std::fmt::Display for HttpMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Get => write!(f, "GET"),
            Self::Post => write!(f, "POST"),
            Self::Delete => write!(f, "DELETE"),
            Self::Put => write!(f, "PUT"),

        }
    }
}