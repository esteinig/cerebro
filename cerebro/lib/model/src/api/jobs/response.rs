use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use bson::serde_helpers::chrono_datetime_as_bson_datetime_optional;
use bson::serde_helpers::chrono_datetime_as_bson_datetime;
use super::model::ScheduleJob;

#[derive(Serialize, Deserialize)]
pub struct EnqueueJobResponse {
    pub status: String,
    pub message: String,
    /// Optionally return some identifier or echo the kind/queue; here we return the job id.
    pub data: Option<String>,
}
impl EnqueueJobResponse {
    pub fn success(id: &str) -> Self {
        Self {
            status: "success".into(),
            message: "Job enqueued".into(),
            data: Some(id.to_string()),
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: "error".into(),
            message: format!("Error enqueueing job: {error_message}"),
            data: None,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct CreateScheduleResponse {
    pub status: String,
    pub message: String,
    /// Return the created schedule id
    pub data: Option<String>,
}
impl CreateScheduleResponse {
    pub fn success(id: &str) -> Self {
        Self {
            status: "success".into(),
            message: "Schedule created".into(),
            data: Some(id.to_string()),
        }
    }
    pub fn conflict(msg: &str) -> Self {
        Self {
            status: "fail".into(),
            message: msg.to_string(),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: "error".into(),
            message: format!("Error creating schedule: {error_message}"),
            data: None,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct JobCompletionResponse {
    pub completed: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub result: Option<serde_json::Value>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    // optional: echo current status for UIs
    #[serde(skip_serializing_if = "Option::is_none")]
    pub status: Option<String>,
}

#[derive(Serialize, Deserialize)]
pub struct ListSchedulesResponse {
    pub status: String,
    pub message: String,
    pub data: Option<Vec<ScheduleJob>>,
}
impl ListSchedulesResponse {
    pub fn success(schedules: Vec<ScheduleJob>) -> Self {
        Self {
            status: "success".into(),
            message: "Schedules found".into(),
            data: Some(schedules),
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: "fail".into(),
            message: "No schedules found".into(),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: "error".into(),
            message: format!("Error in database query: {error_message}"),
            data: None,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct DeleteScheduleResponse {
    pub status: String,
    pub message: String,
    pub data: Option<ScheduleJob>,
}
impl DeleteScheduleResponse {
    pub fn success(removed: ScheduleJob) -> Self {
        Self {
            status: "success".into(),
            message: "Schedule deleted".into(),
            data: Some(removed),
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: "fail".into(),
            message: "Schedule not found".into(),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: "error".into(),
            message: format!("Error deleting schedule: {error_message}"),
            data: None,
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct SimpleScheduleActionResponse {
    pub status: String,
    pub message: String,
    pub data: Option<String>,
}
impl SimpleScheduleActionResponse {
    pub fn success(msg: &str, id: &str) -> Self {
        Self {
            status: "success".into(),
            message: msg.to_string(),
            data: Some(id.to_string()),
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: "fail".into(),
            message: "Schedule not found".into(),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: "error".into(),
            message: format!("Error updating schedule: {error_message}"),
            data: None,
        }
    }
}


// Lightweight public summary
#[derive(Serialize, Deserialize, Clone)]
pub struct ScheduleJobSummary {
    pub id: String,
    pub kind: String,
    pub queue: String,
    #[serde(with = "chrono_datetime_as_bson_datetime")]
    pub run_at: DateTime<Utc>,
    pub interval_seconds: Option<i64>,
    pub enabled: bool,
    #[serde(with = "chrono_datetime_as_bson_datetime_optional")]
    pub last_run_at: Option<DateTime<Utc>>,
    pub last_error: Option<String>,
}
impl From<&ScheduleJob> for ScheduleJobSummary {
    fn from(j: &ScheduleJob) -> Self {
        Self {
            id: j.id.clone(),
            kind: j.kind.clone(),
            queue: j.queue.clone(),
            run_at: j.run_at,
            interval_seconds: j.interval_seconds,
            enabled: j.enabled,
            last_run_at: j.last_run_at,
            last_error: j.last_error.clone(),
        }
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct JobsStatusData {
    pub next: Vec<ScheduleJobSummary>,
    pub previous: Vec<ScheduleJobSummary>,
}

#[derive(Serialize, Deserialize)]
pub struct JobsStatusResponse {
    pub status: String,
    pub message: String,
    pub data: Option<JobsStatusData>,
}
impl JobsStatusResponse {
    pub fn success(data: JobsStatusData) -> Self {
        Self {
            status: "success".into(),
            message: "Job status summary".into(),
            data: Some(data),
        }
    }
    pub fn not_found() -> Self {
        Self {
            status: "fail".into(),
            message: "No job status available".into(),
            data: None,
        }
    }
    pub fn server_error(error_message: String) -> Self {
        Self {
            status: "error".into(),
            message: format!("Error retrieving job status: {error_message}"),
            data: None,
        }
    }
}