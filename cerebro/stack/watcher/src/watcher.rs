use cerebro_model::api::files::model::WatcherConfig;
use notify::EventKind;
use notify::event::CreateKind;
use notify::{Config, PollWatcher, RecursiveMode, Watcher};
use std::path::{Path, PathBuf};
use std::time::Duration;
use std::thread;
use serde::{Serialize, Deserialize};

use crate::filer::{CerebroClientConfig, WatchFiler};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlackConfig {
    pub channel: String,
    pub token: String
}

pub fn watch_production<P: AsRef<Path>>(watch_path: P, interval: Duration, timeout: Duration, timeout_interval: Duration, slack_config: Option<SlackConfig>) -> anyhow::Result<()> {
    
    let watch_config = WatcherConfig::default();
    let cerebro_client_config = CerebroClientConfig::default();

    let (tx, rx) = std::sync::mpsc::channel();
    let tx_c = tx.clone();

    let mut watcher = PollWatcher::new(
        move |watch_event| {
            tx_c.send(watch_event).unwrap();
        },
        Config::default().with_poll_interval(interval)
    )?;

    // Instantiate watch filer to health check API 
    // and FS prior to starting watcher
    WatchFiler::new(
        PathBuf::new(), PathBuf::from("PLACEHOLDER"), 
        watch_config.clone() , cerebro_client_config.clone(), slack_config.clone()
    )?;

    log::info!("Starting watcher {} @ {} ...", watch_config.name, watch_config.location);

    watcher.watch(watch_path.as_ref(), RecursiveMode::Recursive)?;

    'event: for e in rx {

        match e {
            Ok(event) => {
                match event.kind {
                    // We are using a poller, no check for folders or file types
                    EventKind::Create(CreateKind::Any) => {
                        
                        let (path, is_dir) = match event.paths.first() {
                            Some(path) => {                             
                                (path.to_owned(), path.is_dir())
                            },
                            None => {
                                log::error!("[{}] Could not extract input path from event", watch_config.name);
                                continue 'event;
                            }
                        };

                        if is_dir {

                            // Check that the created directory is one level down from the watch path
                            // otherwise any other subdirectories  will be captured
                            match path.parent() {
                                Some(parent_dir) => {
                                    if parent_dir != watch_path.as_ref() {
                                        continue 'event;
                                    }
                                },
                                None => {
                                    log::error!("[{}@{}] Could not extract parent directory from event", watch_config.name, watch_config.location);
                                    continue 'event;
                                }
                            }

                            // Spawn a thread that handles input directory completion 
                            // and validation to launch the workflow run 
                            let slack_cfg = slack_config.clone();
                            let watch_cfg = watch_config.clone();
                            let cerebro_client_cfg = cerebro_client_config.clone();

                            thread::spawn(move || {
                                log::info!("[{}@{}] Input directory detected: {}", watch_cfg.name, watch_cfg.location, path.display());

                                match watch_event_timeout(path.clone(), timeout_interval, timeout, &watch_cfg.name, &watch_cfg.location) {
                                    Ok(_) => {
                                        
                                        match WatchFiler::new(
                                            PathBuf::new(), path, watch_cfg.clone() , cerebro_client_cfg, slack_cfg,   // TODO BASE PATH
                                        ) {
                                            Ok(mut filer) => {
                                                
                                                if let Err(err) = filer.validate_inputs() {
                                                    log::warn!("[{}@{}] Failed to validate inputs (error: {})", watch_cfg.name, watch_cfg.location, err.to_string())
                                                } else {
                                                    if let Err(err) = filer.upload_and_register_illumina_pe(
                                                        &watch_cfg.cerebro_team_name, 
                                                        &watch_cfg.cerebro_db_name,
                                                        "*_{R1,R2}.fastq.gz",
                                                        false,
                                                        true
                                                    ) {
                                                        log::warn!("[{}@{}] Failed to upload files to Cerebro FS API (error: {})", watch_cfg.name, watch_cfg.location, err.to_string())
                                                    }
                                                }
                                            }, 
                                            Err(err) => {
                                                log::warn!("[{}@{}] Could not initiate watch filer (error: {})", watch_cfg.name, watch_cfg.location, err.to_string())
                                            }
                                        };
                                    }, 
                                    Err(err) => {
                                        log::warn!("[{}@{}] Failed to watch input directory for completion (error: {})", watch_cfg.name, watch_cfg.location, err.to_string());
                                    }
                                };
                            });
                        }
                        // Not joining thread - if the main watch process fails, threads will exit!
                    },
                    _ => {}
                }
            },
            Err(_) => {}
        }
    }

    Ok(())
}


pub fn watch_event_timeout<P: AsRef<Path>>(path: P, interval: Duration, timeout: Duration, watcher_name: &str, watcher_loc: &str) -> notify::Result<()> {

    let path_name = match path.as_ref().file_name() {
        Some(name) => name.to_str().unwrap_or("unknown"), None => "unknown"
    };
    
    log::info!("[{watcher_name}@{watcher_loc}::{path_name}] Watching input directory for changes...");

    let (tx, rx) = std::sync::mpsc::channel();

    let tx_c = tx.clone();
    
    let mut watcher = PollWatcher::new(
        move |watch_event| {
            tx_c.send(watch_event).unwrap();
        },
        Config::default().with_poll_interval(interval)
    )?;

    watcher.watch(path.as_ref(), RecursiveMode::Recursive)?;

    loop {
        match rx.recv_timeout(timeout) {
            Ok(event) => {
                log::info!("[{watcher_name}@{watcher_loc}::{path_name}] Event received, continue polling...");

                // If the event is an error (e.g. timeout directory deleted)
                // this handle will terminate the timeout watcher 
                if let Err(err) = event {
                    log::warn!("[{watcher_name}@{watcher_loc}::{path_name}] Error in poll watcher, terminating polling");
                    return Err(err)
                }

            },
            Err(_) => {
                log::info!("[{watcher_name}@{watcher_loc}::{path_name}] No event received before timeout");
                break
            }
        }
        thread::sleep(Duration::from_secs(1));
    }

    log::info!("[{watcher_name}@{watcher_loc}::{path_name}] Timeout watcher thread completed");
    log::info!("[{watcher_name}@{watcher_loc}::{path_name}] Continue with input checks and notifications");

    Ok(())
}