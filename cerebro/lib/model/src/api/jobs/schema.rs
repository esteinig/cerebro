use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnqueueJobRequest {
    pub kind: String,
    pub args: serde_json::Value,   // single object; server wraps as args[0]
    pub queue: Option<String>,
    pub retry: Option<i32>,
    pub reserve_for_seconds: Option<u64>,
}

// Matches your response.rs EnqueueJobResponse
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnqueueJobResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScheduleJobRequest {
    pub kind: String,
    pub args: serde_json::Value,
    pub queue: Option<String>,
    pub run_at: chrono::DateTime<chrono::Utc>, // first execution time (UTC)
    pub interval_seconds: Option<i64>,         // None => one-shot
    pub retry: Option<i32>,
    pub reserve_for_seconds: Option<u64>,
    pub enabled: Option<bool>,
}

// Matches your response.rs CreateScheduleResponse
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CreateScheduleResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>, // schedule_id
}