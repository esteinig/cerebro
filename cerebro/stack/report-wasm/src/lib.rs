#[cfg(feature = "lib")]
pub mod world;
#[cfg(feature = "lib")]
pub mod compiler;
#[cfg(feature = "lib")]
pub mod report;

#[cfg(feature = "lib")]
pub use compiler::LibraryReportCompiler;
#[cfg(feature = "lib")]
pub use report::{ReportFormat, ReportType, ReportConfig, PathogenDetectionReport, TrainingCompletionReport};