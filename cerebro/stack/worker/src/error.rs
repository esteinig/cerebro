//! Worker error taxonomy (S3-1).
//!
//! [`WorkerError`] is the internal error type for job runners and context setup.
//! Faktory's `JobRunner` requires `Error: Into<Box<dyn Error>>`/`std::io::Error`,
//! so [`WorkerError::into_io`] gives a clean conversion at the runner boundary.

use std::io;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum WorkerError {
    /// A required client was not configured (missing API/FS settings). The job is
    /// not runnable until the worker is configured; surfaced rather than silently
    /// skipped.
    #[error("{0} client is not configured (check the CEREBRO_* worker env)")]
    ClientNotConfigured(&'static str),

    /// The worker reached a job kind whose behaviour is delivered in a later
    /// package (S3-2 / S3-3). Carried so stub runners are explicit, not silent.
    #[error("job kind '{0}' is not yet implemented (delivered in a later Stage 3 package)")]
    NotImplemented(&'static str),

    /// Malformed job arguments.
    #[error("invalid job args: {0}")]
    InvalidArgs(String),

    /// An error from the Cerebro API client (lifecycle endpoints).
    #[error("api client error: {0}")]
    Api(String),

    /// An error from the Cerebro FS client (capture/verify/restore/physical ops).
    #[error("fs client error: {0}")]
    Fs(String),

    /// Any other failure.
    #[error("{0}")]
    Other(String),
}

impl WorkerError {
    /// Convert to the `std::io::Error` Faktory's `JobRunner` expects.
    pub fn into_io(self) -> io::Error {
        io::Error::new(io::ErrorKind::Other, self.to_string())
    }
}

impl From<WorkerError> for io::Error {
    fn from(err: WorkerError) -> Self {
        err.into_io()
    }
}

/// Convenience for mapping any `Display` error into an `io::Error` at a boundary.
pub fn io_err<E: std::fmt::Display>(e: E) -> io::Error {
    io::Error::new(io::ErrorKind::Other, e.to_string())
}
