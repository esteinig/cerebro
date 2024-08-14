


use notify::{Config, PollWatcher, RecursiveMode, Watcher};
use cerebro_fs::client::{FileSystemClient, UploadConfig};
use cerebro_model::api::watchers::model::ProductionWatcher;
use cerebro_client::client::CerebroClient;
use notify::event::CreateKind;
use std::path::Path;
use std::time::Duration;
use notify::EventKind;
use std::thread;

use crate::utils::FileGetter;
use crate::error::WatcherError;
use crate::slack::{SlackConfig, SlackTools};

#[derive(Clone, Debug)]
pub struct CerebroWatcher {
    pub config: ProductionWatcher,
    pub upload_config: UploadConfig,
    pub api_client: CerebroClient,
    pub fs_client: FileSystemClient,
    pub slack_tools: Option<SlackTools>
}
impl CerebroWatcher {
    pub fn new(
        config: ProductionWatcher, 
        api_client: CerebroClient, 
        fs_client: FileSystemClient,
        upload_config: UploadConfig,
        slack_config: Option<SlackConfig>
    ) -> Result<Self, WatcherError> {
        
        let slack_tools = match slack_config {
            Some(slack_config) => {
                let slack_tools = SlackTools::from_config(&slack_config);
                log::info!("Sending watcher initialisation message to Slack channel: {}", slack_config.channel);
                slack_tools.client.send(
                    &slack_tools.message.watcher_setup(&config.name, &config.location)
                )?;
                Some(slack_tools)
            },
            _ => None
        };

        Ok(Self {
            config,
            upload_config,
            api_client,
            fs_client,
            slack_tools,
        })
    }
    pub fn watch<P: AsRef<Path>>(&self, path: P, interval: Duration, timeout: Duration, timeout_interval: Duration) -> Result<(), WatcherError> {

        let (tx, rx) = std::sync::mpsc::channel();
        let tx_c = tx.clone();

        // Setup a poll watcher so that paths on networked drives can be watched 
        // that are not supported by internal notification systems
        let mut watcher = PollWatcher::new(
            move |watch_event| {
                tx_c.send(watch_event).unwrap();
            },
            Config::default().with_poll_interval(interval)
        )?;

        let ping_clone = self.clone();
        thread::spawn(move || {
            loop {
                if let Err(err) = ping_clone.api_client.ping_watcher(
                    &ping_clone.config.id, 
                    false
                ) {
                    log::error!("Error in updating watcher activity: {}", err.to_string())
                };
                thread::sleep(std::time::Duration::from_secs(60));
            }
        });

        log::info!("Starting watcher {} @ {} ...", self.config.name, self.config.location);

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
                                    log::error!("[{}] Could not extract input path from event", self.config.name);
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
                                        log::error!("[{}@{}] Could not extract parent directory from event", self.config.name, self.config.location);
                                        continue 'event;
                                    }
                                }

                                let watcher = self.clone();
                                let watcher_config = watcher.config.clone();

                                thread::spawn(move || {

                                    let run_id = match new_dir.file_name() {
                                        Some(name) => name.to_str().unwrap_or("unknown"), 
                                        None => "unknown"
                                    };    

                                    log::info!("[{}@{}::{}] Run directory detected: {}", watcher.config.name, watcher.config.location, run_id, new_dir.display());

                                    // We watch the input directory for any changes in the given timeout period with the given timeout interval...
                                    match watch_event_timeout(
                                        new_dir.clone(), 
                                        timeout_interval, 
                                        timeout, 
                                        &watcher.config.name, 
                                        &watcher.config.location
                                    ) {
                                        Ok(_) => {
                                            // ... if there are no changes to the input directory, we start the file registration and upload 

                                            match watcher.config.format.get_fastq_files(&new_dir, Some(watcher.config.glob)) {
                                                Err(err) => log::error!("Error getting read files: {}", err.to_string()),
                                                Ok(fastq_files) => {

                                                    match watcher.fs_client.upload_files_from_watcher(
                                                        &fastq_files, 
                                                        run_id.to_string(),
                                                        Some(watcher.config.format.file_type()), 
                                                        watcher.upload_config,
                                                        watcher_config
                                                    ) {
                                                        Err(err) => log::error!("Error uploading read files to Cerebro FS: {}", err.to_string()),
                                                        Ok(()) => log::info!("Uploaded files for run: {run_id}")
                                                    }
                                                }
                                            }
                                        }, 
                                        Err(err) => {
                                            log::warn!("[{}@{}] Failed to watch input directory for completion (error: {})", watcher.config.name, watcher.config.location, err.to_string());
                                        }
                                    };
                                });
                            }
                            // Not joining timeout threads - if the main watch process fails, threads will exit
                        },
                        _ => {} // Nothing happens if other events are registered
                    }
                },
                Err(_) => {} // Nothing happens if we get an error in the event channel
            }
        }
        Ok(())
    }

}

pub fn watch_event_timeout<P: AsRef<Path>>(path: P, timeout_interval: Duration, timeout: Duration, watcher_name: &str, watcher_loc: &str) -> notify::Result<()> {

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
        Config::default().with_poll_interval(timeout_interval)
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