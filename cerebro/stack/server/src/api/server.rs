

use clap::Parser;
use actix_cors::Cors;
use redis::Client as RedisClient;
use actix_web::middleware::Logger;
use actix_web::{web, App, HttpServer};
use mongodb::{bson::doc, Client as MongoClient};

use crate::api::logs::handler::logs_config;
use crate::api::auth::handler::auth_config;
use crate::api::users::handler::user_config;
use crate::api::teams::handler::team_config;
use crate::api::files::handler::files_config;
use crate::api::cerebro::handler::cerebro_config;
use crate::api::watchers::handler::watchers_config;
use crate::api::towers::handler::towers_config;


use crate::terminal::{App as Cli, Commands};

use super::stage::handler::stage_config;



/*
=============================
MAIN ASYNC LAUNCH THROUCH CLI
=============================
*/


/// Main application configuration
fn app_config(cfg: &mut web::ServiceConfig) {
    let json_cfg = web::JsonConfig::default().limit(16 * 1024 * 1024); // 16 MB transfer limit for MongoDB model size limits
    cfg.app_data(json_cfg);
        
}

pub struct AppState {
    pub db: MongoClient,
    pub auth_session: RedisClient,
    pub auth_onetime: RedisClient,
    pub env: cerebro_model::api::config::Config,
}

/// Main controller function for the async server
/// 
/// Launched through the command-line interface client,
/// which is why we are re-parse the arguments - allows
/// for the main asynchronous routine to be run through 
/// the iniutial synchronous main routine.
#[actix_web::main]
pub async fn main() -> std::io::Result<()> {

    let cli = Cli::parse();

    match cli.command {
        Commands::Run ( args ) => {

            // Setup logging for development - 
            // integrate with global verbosity setting
            if std::env::var_os("RUST_LOG").is_none() {
                std::env::set_var("RUST_LOG", "actix_web=info");
            }

            // Database and authentication configuration
            let config = match cerebro_model::api::config::Config::from_toml(&args.config) {
                Ok(config) => {
                    log::info!("Received valid server configuration");
                    config
                },
                Err(err) => {
                    log::error!("Failed to read configuration file: {}", err.to_string());
                    std::process::exit(1);
                }
            };


            // URI validation and client creation
            let mongo_client = match MongoClient::with_uri_str(&config.database.connections.mongodb_uri).await  {
                Ok(client) => {
                    log::info!("Created MongoDB client for main database");
                    client
                }
                Err(e) => {
                    log::error!("Invalid MongoDB URI for main database: {}", e);
                    std::process::exit(1);
                }
            };

            let redis_client_auth_session = match RedisClient::open(config.database.connections.redis_auth_session_uri.to_owned()) {
                Ok(client) => {
                    log::info!("Created Redis client for session database");
                    client
                }
                Err(e) => {
                    log::error!("Invalid Redis URI for session database: {}", e);
                    std::process::exit(1);
                }
            };

            let redis_client_auth_onetime = match RedisClient::open(config.database.connections.redis_auth_onetime_uri.to_owned()) {
                Ok(client) => {
                    log::info!("Created Redis client for one-time database");
                    client
                }
                Err(e) => {
                    log::error!("InvalidRedis URI for one-time database: {}", e);
                    std::process::exit(1);
                }
            };

            // Database connection checks
            match mongo_client.list_database_names(None, None).await {
                Ok(_) => log::info!("Connected to MongoDB main database"),
                Err(err) => {
                    log::error!("Could not connect to MongoDB main database: {}", err.to_string());
                    std::process::exit(1);
                }
            };

            match redis_client_auth_session.get_connection() {
                Ok(_) => log::info!("Connected to Redis session database"),
                Err(err) => {
                    log::error!("Could not connect to Redis session database: {}", err.to_string());
                    std::process::exit(1);
                }
            };

            match redis_client_auth_onetime.get_connection() {
                Ok(_) => log::info!("Connected to Redis one-time database"),
                Err(err) => {
                    log::error!("Could not connect to Redis one-time database: {}", err.to_string());
                    std::process::exit(1);
                }
            };

            HttpServer::new(move || {
                App::new()
                    .wrap(
                        Cors::default()
                        .allowed_origin(&config.security.cors.app_origin_public_url)
                        .allowed_origin(&config.security.cors.app_origin_docker_url)
                        .allowed_methods(vec![
                            "GET", "POST", "PATCH", "DELETE"
                        ])
                        .allowed_headers(vec![
                            actix_web::http::header::AUTHORIZATION,
                            actix_web::http::header::CONTENT_TYPE,
                            actix_web::http::header::ACCEPT,
                        ])
                        .expose_headers(vec![
                            actix_web::http::header::CONTENT_DISPOSITION
                        ])
                        .supports_credentials()
                        .max_age(3600)
                    )
                    .app_data(web::Data::new(AppState {
                        db: mongo_client.clone(),
                        auth_session: redis_client_auth_session.clone(),
                        auth_onetime: redis_client_auth_onetime.clone(),
                        env: config.clone(),
                    }))
                    .configure(app_config)
                    // Email authentication configuration for global activation
                    .configure(|cfg| auth_config(cfg, &config))
                    // Endpoint configurations
                    .configure(user_config)
                    .configure(team_config)
                    .configure(logs_config)
                    .configure(files_config)
                    .configure(towers_config)
                    .configure(watchers_config)
                    .configure(stage_config)
                    // Application functionality configuration for global security
                    .configure(|cfg| cerebro_config(cfg, &config))
                    .wrap(Logger::default())
            })
            .workers(args.threads)
            .bind((args.host, args.port))?
            .run()
            .await
        }
    }
}
