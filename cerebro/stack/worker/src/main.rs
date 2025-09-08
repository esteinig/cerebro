use std::io;
use async_trait::async_trait;
use faktory::{Job, JobRunner, Worker};
use serde::Deserialize;
use tracing::Level;

// =======================
// helpers
// =======================
fn io_err<E: std::fmt::Display>(e: E) -> io::Error {
    io::Error::new(io::ErrorKind::Other, e.to_string())
}

/// Pull the first JSON object from job.args() and deserialize into T.
fn parse_args<T: for<'de> Deserialize<'de>>(job: &Job) -> Result<T, io::Error> {
    let first = job
        .args()
        .get(0)
        .ok_or_else(|| io_err("missing job args (expected single JSON object)"))?;
    serde_json::from_value::<T>(first.clone()).map_err(io_err)
}

// =======================
// 1) generate_report
// =======================
#[derive(Debug, Deserialize)]
struct ReportArgs {
    range: String,                 // e.g. "yesterday"
    #[serde(default)]
    format: Option<String>,        // e.g. "pdf"
}

struct GenerateReport;

#[async_trait]
impl JobRunner for GenerateReport {
    type Error = io::Error;

    async fn run(&self, job: Job) -> Result<(), Self::Error> {
        let args: ReportArgs = parse_args(&job)?;
        let format = args.format.as_deref().unwrap_or("pdf");

        tracing::info!(kind=%job.kind(), range=%args.range, %format, "generate_report started");
        // Simulate heavy work
        tokio::time::sleep(std::time::Duration::from_secs(3)).await;
        tracing::info!("generate_report done");
        Ok(())
    }
}

// =======================
// 2) aggregate_metrics
// =======================
#[derive(Debug, Deserialize)]
struct AggregateArgs {
    region: String,                // e.g. "eu-west"
    #[serde(default = "default_hours")]
    hours: u32,                    // default 24
}
fn default_hours() -> u32 { 24 }

struct AggregateMetrics;

#[async_trait]
impl JobRunner for AggregateMetrics {
    type Error = io::Error;

    async fn run(&self, job: Job) -> Result<(), Self::Error> {
        let args: AggregateArgs = parse_args(&job)?;
        tracing::info!(kind=%job.kind(), region=%args.region, hours=%args.hours, "aggregate_metrics started");

        for step in 1..=3 {
            tokio::time::sleep(std::time::Duration::from_secs(1)).await;
            tracing::info!(%step, "aggregate_metrics progress");
        }

        tracing::info!("aggregate_metrics done");
        Ok(())
    }
}

// =======================
// 3) ping
// =======================
#[derive(Debug, Deserialize)]
struct PingArgs {
    #[serde(default = "default_n")]
    n: u64,
}
fn default_n() -> u64 { 1 }

struct Ping;

#[async_trait]
impl JobRunner for Ping {
    type Error = io::Error;

    async fn run(&self, job: Job) -> Result<(), Self::Error> {
        let args: PingArgs = parse_args(&job)?;
        for i in 1..=args.n {
            tracing::info!(i, "pong");
            tokio::time::sleep(std::time::Duration::from_millis(400)).await;
        }
        Ok(())
    }
}

// =======================
// 4) gridfs_process
//     downloads from Mongo GridFS and processes
// =======================
#[derive(Debug, Deserialize)]
struct GridFsArgs {
    file_id: String,               // hex ObjectId
    #[serde(default = "default_bucket")]
    bucket: String,                // default "fs"
}
fn default_bucket() -> String { "fs".into() }

struct GridFsProcess;

#[async_trait]
impl JobRunner for GridFsProcess {
    type Error = io::Error;

    async fn run(&self, job: Job) -> Result<(), Self::Error> {

        let args: GridFsArgs = parse_args(&job)?;

        // Bucket access

        tracing::info!(bytes=0, file_id=%args.file_id, bucket=%args.bucket, "gridfs_process done");
        Ok(())
    }
}


// =======================
// main: bootstrap worker
// =======================
#[tokio::main]
async fn main() -> io::Result<()> {

    tracing_subscriber::fmt()
        .with_env_filter(
            std::env::var("RUST_LOG").unwrap_or_else(|_| "info,faktory_worker=info".into()),
        )
        .with_level(true)
        .with_target(false)
        .with_max_level(Level::INFO)
        .init();

    // Uses FAKTORY_URL if present (e.g. tcp://:devpass@faktory:7419)
    match std::env::var("FAKTORY_URL") {
        Ok(url) => println!("FAKTORY_URL={}", url),
        Err(_)  => println!("FAKTORY_URL not set; using faktory::Client defaults"),
    }

    let mut w = Worker::builder()
        .register("generate_report", GenerateReport)
        .register("aggregate_metrics", AggregateMetrics)
        .register("ping", Ping)
        .register("gridfs_process", GridFsProcess)
        .connect()
        .await
        .expect("connect to faktory");

    // Match queues the app uses for this worker
    if let Err(e) = w.run(&["default", "backend", "maintenance"]).await {
        eprintln!("worker failed: {e}");
    }

    Ok(())
}