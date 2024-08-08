use cerebro_client::client::CerebroClient;
use cerebro_fs::client::{FileSystemClient, UploadConfig};
use cerebro_model::api::files::model::{WatcherFormat, WatcherConfig};
use notify::EventKind;
use notify::event::CreateKind;
use notify::{Config, PollWatcher, RecursiveMode, Watcher};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::Duration;
use std::thread;
use serde::{Serialize, Deserialize};

use crate::error::WatcherError;
use crate::slack::{SlackClient, SlackConfig, SlackTools};
use crate::terminal::{App, WatchArgs};
use crate::utils::FileGetter;


#[derive(Clone, Debug)]
pub struct CerebroClientConfig {
    pub api_url: String,
    pub api_token: Option<String>,
    pub api_token_file: Option<PathBuf>,
    pub _danger_invalid_certificate: bool,
    pub fs_url: String,
    pub fs_port: String,
}
impl Default for CerebroClientConfig {
    fn default() -> Self {
        Self {
            api_url: String::from("http://api.cerebro.localhost"),
            api_token: std::env::var("CEREBRO_API_TOKEN").ok(),
            api_token_file: None,
            _danger_invalid_certificate: false,
            fs_url: String::from("http://fs.cerebro.localhost"),
            fs_port: String::from("9333"),
        }
    }
}

#[derive(Clone, Debug)]
pub struct CerebroWatcher {
    pub watcher_config: WatcherConfig,
    pub upload_config: UploadConfig,
    pub api_client: CerebroClient,
    pub fs_client: FileSystemClient,
    pub slack_tools: Option<SlackTools>
}
impl CerebroWatcher {
    pub fn new(
        watcher_config: WatcherConfig, 
        client_config: CerebroClientConfig, 
        upload_config: UploadConfig,
        slack_config: Option<SlackConfig>
    ) -> Result<Self, WatcherError> {
        

        // Setup the Cerebro API and FS clients and ping status

        let api_client = CerebroClient::new(
            &client_config.api_url,
            &client_config.api_token,
            false,
            false,
            &client_config.api_token_file
        )?;

        let fs_client = FileSystemClient::new(
            &api_client, 
            &client_config.fs_url, 
            &client_config.fs_port
        );

        log::info!("Checking status of Cerebro API at {}",  &api_client.url);
        api_client.ping_servers()?;

        log::info!("Checking status of SeaweedFS master at {}",  &fs_client.fs_url);
        fs_client.ping_status()?;

        let slack_tools = match slack_config {
            Some(slack_config) => {
                let slack_tools = SlackTools::from_config(&slack_config);
                log::info!("Sending watcher initialisation message to Slack channel: {}", slack_config.channel);
                slack_tools.client.send(
                    &slack_tools.message.watcher_setup(&watcher_config.name, &watcher_config.location)
                )?;
                Some(slack_tools)
            },
            _ => None
        };

        Ok(Self {
            watcher_config,
            upload_config,
            api_client,
            fs_client,
            slack_tools,
        })
    }
    pub fn watch<P: AsRef<Path>>(&self, path: P, fastq_glob: Option<String>) -> Result<(), WatcherError> {

        let (tx, rx) = std::sync::mpsc::channel();
        let tx_c = tx.clone();

        // Setup a poll watcher so that paths on networked drives can be watched 
        // that are not supported by internal notification systems
        let mut watcher = PollWatcher::new(
            move |watch_event| {
                tx_c.send(watch_event).unwrap();
            },
            Config::default().with_poll_interval(self.watcher_config.interval)
        )?;

        log::info!("Starting watcher {} @ {} ...", self.watcher_config.name, self.watcher_config.location);

        watcher.watch(path.as_ref(), RecursiveMode::Recursive)?;


        // Main event loop for watcher
        'event: for e in rx {

            match e {
                Ok(event) => {
                    match event.kind {
                        // We are using a poller, no check for folders or file types so
                        // we do a manual check on the first element in the returned 
                        // event paths (input path)
                        EventKind::Create(CreateKind::Any) => {
                            
                            let (new_dir, is_dir) = match event.paths.first() {
                                Some(path) => {                             
                                    (path.to_owned(), path.is_dir())
                                },
                                None => {
                                    log::error!("[{}] Could not extract input path from event", self.watcher_config.name);
                                    continue 'event;
                                }
                            };

                            if is_dir {
                                // Check that the created directory is one level down from the watch path
                                // otherwise any other subdirectories will also be captured
                                match new_dir.parent() {
                                    Some(parent_dir) => {
                                        if parent_dir != path.as_ref() {
                                            continue 'event;
                                        }
                                    },
                                    None => {
                                        log::error!("[{}@{}] Could not extract parent directory from event", self.watcher_config.name, self.watcher_config.location);
                                        continue 'event;
                                    }
                                }

                                let watcher = self.clone();
                                let glob = fastq_glob.clone();

                                thread::spawn(move || {

                                    let run_id = match new_dir.file_name() {
                                        Some(name) => name.to_str().unwrap_or("unknown"), 
                                        None => "unknown"
                                    };    

                                    log::info!("[{}@{}::{}] Run directory detected: {}", watcher.watcher_config.name, watcher.watcher_config.location, run_id, new_dir.display());

                                    // We watch the input directory for any changes in the given timeout period with the given timeout interval...
                                    match watch_event_timeout(
                                        new_dir.clone(), 
                                        watcher.watcher_config.timeout_interval, 
                                        watcher.watcher_config.timeout, 
                                        &watcher.watcher_config.name, 
                                        &watcher.watcher_config.location
                                    ) {
                                        Ok(_) => {
                                            // ... if there are no changes to the input directory, we start the file registration and upload 

                                            match watcher.watcher_config.format.get_fastq_files(&new_dir, glob) {
                                                Err(err) => log::error!("Error getting read files: {}", err.to_string()),
                                                Ok(fastq_files) => {

                                                    match watcher.fs_client.upload_samples(&fastq_files, run_id.to_string(), watcher.upload_config, watcher.watcher_config) {
                                                        Err(err) => log::error!("Error uploading read files to Cerebro FS: {}", err.to_string()),
                                                        Ok(()) => log::info!("Uploaded files for run: {run_id}")
                                                    }

                                                }
                                            }

                                        }, 
                                        Err(err) => {
                                            log::warn!("[{}@{}] Failed to watch input directory for completion (error: {})", watcher.watcher_config.name, watcher.watcher_config.location, err.to_string());
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

}

pub fn watch_event_timeout<P: AsRef<Path>>(path: P, interval: Duration, timeout: Duration, watcher_name: &str, watcher_loc: &str) -> notify::Result<()> {

    let run_id = match path.as_ref().file_name() {
        Some(name) => name.to_str().unwrap_or("unknown"), 
        None => "unknown"
    };  
    
    log::info!("[{watcher_name}@{watcher_loc}::{run_id}] Watching input directory for changes...");

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
                log::info!("[{watcher_name}@{watcher_loc}::{run_id}] Event received, continue polling...");

                // If the event is an error (e.g. timeout directory deleted) this handle will terminate the timeout watcher 
                if let Err(err) = event {
                    log::warn!("[{watcher_name}@{watcher_loc}::{run_id}] Error in poll watcher, terminating polling");
                    return Err(err)
                }

            },
            Err(_) => {
                log::info!("[{watcher_name}@{watcher_loc}::{run_id}] No event received before timeout");
                break
            }
        }
        thread::sleep(Duration::from_secs(1));
    }

    log::info!("[{watcher_name}@{watcher_loc}::{run_id}] Timeout watcher thread completed");
    log::info!("[{watcher_name}@{watcher_loc}::{run_id}] Continue with input checks and notifications");

    Ok(())
}