use serde::{Serialize, Deserialize};
use serde_json::Value as JsonValue;
use chrono::{DateTime, Utc};
use mongodb::bson::serde_helpers::chrono_datetime_as_bson_datetime;
use mongodb::bson::serde_helpers::chrono_datetime_as_bson_datetime_optional;
use uuid::Uuid;

/// Primary key type (string UUID works well with existing patterns).
pub type ScheduleId = String;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScheduleJob {
    /// App-level identifier (UUID string).
    pub id: ScheduleId,

    /// Worker type (Faktory "kind"), e.g. "gridfs_process".
    pub kind: String,

    /// Single JSON object with job arguments (will be wrapped as [ { ... } ] for Faktory).
    pub args: JsonValue,

    /// Target Faktory queue.
    pub queue: String,

    /// Next time the job should run.
    #[serde(with = "chrono_datetime_as_bson_datetime")]
    pub run_at: DateTime<Utc>,

    /// If Some(N), the job recurs every N seconds; if None, one-shot.
    pub interval_seconds: Option<i64>,

    /// Default Faktory retry count.
    pub retry: i32,

    /// Faktory reservation TTL (seconds).
    pub reserve_for_seconds: u64,

    /// Whether the scheduler should consider this job.
    pub enabled: bool,

    /// Last time this job was enqueued.
    #[serde(with = "chrono_datetime_as_bson_datetime_optional")]
    pub last_run_at: Option<DateTime<Utc>>,

    /// Last scheduling error (if any).
    pub last_error: Option<String>,

    /// Audit fields.
    #[serde(with = "chrono_datetime_as_bson_datetime")]
    pub created_at: DateTime<Utc>,
    #[serde(with = "chrono_datetime_as_bson_datetime")]
    pub updated_at: DateTime<Utc>,
}

impl ScheduleJob {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        kind: String,
        args: JsonValue,
        queue: Option<String>,
        run_at: DateTime<Utc>,
        interval_seconds: Option<i64>,
        retry: Option<i32>,
        reserve_for_seconds: Option<u64>,
        enabled: Option<bool>,
    ) -> Self {
        let now = Utc::now();
        Self {
            id: uuid::Uuid::new_v4().to_string(),
            kind,
            args,
            queue: queue.unwrap_or_else(|| "default".into()),
            run_at,
            interval_seconds,
            retry: retry.unwrap_or(10),
            reserve_for_seconds: reserve_for_seconds.unwrap_or(3600),
            enabled: enabled.unwrap_or(true),
            last_run_at: None,
            last_error: None,
            created_at: now,
            updated_at: now,
        }
    }
}

#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct JobRunDoc {
    #[serde(default)]
    pub id: Option<uuid::Uuid>,            // optional if you also track a UUID
    pub jid: String,                       // Faktory job id (string)
    pub kind: String,
    pub queue: String,
    pub status: String,                    // "queued" | "running" | "succeeded" | "failed"
    #[serde(default)]
    pub result: Option<serde_json::Value>,
    #[serde(default)]
    pub error: Option<String>,
    pub created_at: chrono::DateTime<chrono::Utc>,
    pub updated_at: chrono::DateTime<chrono::Utc>,
}