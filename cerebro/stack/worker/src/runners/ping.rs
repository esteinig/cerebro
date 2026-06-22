//! `ping` — a real liveness job.
//!
//! Cheap, dependency-free job that proves the producer → Faktory → consumer path
//! end-to-end. Enqueue it (e.g. via the server's `POST /jobs/enqueue` with
//! `{"n": 3}`) to smoke-test a deployment.

use std::io;

use async_trait::async_trait;
use faktory::{Job, JobRunner};
use serde::Deserialize;

use crate::runners::parse_args;
use crate::telemetry::{JobOutcome, Metrics};

#[derive(Debug, Deserialize)]
struct PingArgs {
    #[serde(default = "default_n")]
    n: u64,
}
fn default_n() -> u64 {
    1
}

pub struct Ping {
    metrics: Metrics,
}

impl Ping {
    pub fn new(metrics: Metrics) -> Self {
        Self { metrics }
    }
}

#[async_trait]
impl JobRunner for Ping {
    type Error = io::Error;

    async fn run(&self, job: Job) -> Result<(), Self::Error> {
        let kind = job.kind().to_string();
        self.metrics.record_job(&kind, JobOutcome::Started);

        let args: PingArgs = match parse_args::<PingArgs>(&job) {
            Ok(a) => a,
            Err(e) => {
                self.metrics.record_job(&kind, JobOutcome::Failed);
                return Err(e.into());
            }
        };

        for i in 1..=args.n {
            tracing::info!(i, "pong");
            tokio::time::sleep(std::time::Duration::from_millis(200)).await;
        }

        self.metrics.record_job(&kind, JobOutcome::Succeeded);
        Ok(())
    }
}
