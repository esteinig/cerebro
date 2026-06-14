//! Worker health & metrics HTTP server (S3-1).
//!
//! A tiny actix-web server, spawned alongside the Faktory consumer, exposing:
//! * `GET /health`  — liveness/readiness probe (always `200 ok` once running),
//! * `GET /metrics` — Prometheus exposition of the worker metrics.
//!
//! It runs on the same tokio runtime as the worker. Binding happens synchronously
//! so a bad address surfaces immediately; the served future is then spawned.

use actix_web::{get, web, App, HttpResponse, HttpServer};

use crate::telemetry::Metrics;

#[get("/health")]
async fn health() -> HttpResponse {
    HttpResponse::Ok().content_type("text/plain").body("ok")
}

#[get("/metrics")]
async fn metrics_endpoint(metrics: web::Data<Metrics>) -> HttpResponse {
    HttpResponse::Ok()
        .content_type("text/plain; version=0.0.4")
        .body(metrics.encode())
}

/// Bind and spawn the health/metrics server on `addr` (`host:port`).
pub fn spawn(metrics: Metrics, addr: String) -> std::io::Result<()> {
    let data = web::Data::new(metrics);
    let server = HttpServer::new(move || {
        App::new()
            .app_data(data.clone())
            .service(health)
            .service(metrics_endpoint)
    })
    .bind(&addr)?
    .run();

    tracing::info!(%addr, "worker health/metrics server listening");
    tokio::spawn(server);
    Ok(())
}
