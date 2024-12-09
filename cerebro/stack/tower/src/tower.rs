use std::path::PathBuf;
use std::sync::Arc;
use std::time::Duration;
use cerebro_model::api::stage::model::StagedSample;
use tokio::fs::{create_dir_all, remove_dir_all, File};
use tokio::process::{Child, Command};
use tokio::sync::{mpsc, Semaphore};
use tokio::task;
use tokio::io::AsyncWriteExt;
use tokio::time::sleep;

use crate::client::TowerClient;
use crate::error::TowerError;

#[derive(Clone, Debug)]
pub struct NextflowConfig {
    main: PathBuf,
    config: PathBuf,
    workdir: PathBuf,
    databases: PathBuf,
    profile: Vec<String>,

    cleanup: bool,
    execution_directory: Option<PathBuf>
}
impl NextflowConfig {
    pub async fn new(
        main: &PathBuf,
        config: &PathBuf,
        workdir: &PathBuf, 
        databases: &PathBuf,
        profile: &Vec<String>,
        cleanup: bool
    ) -> Result<Self, TowerError> {

        if !databases.exists() {
            log::error!("Failed to detect database directory: {}", databases.display())
        }
        if !main.exists() {
            log::error!("Failed to detect main script for Nextflow: {}", databases.display())
        }
        if !config.exists() {
            log::error!("Failed to detect main configuration for Nextflow: {}", databases.display())
        }
        if !workdir.exists() {
            create_dir_all(&workdir).await?;
        }
        
        Ok(Self { 
            main: main.to_owned(),
            config: config.to_owned(),
            workdir: workdir.to_owned(),
            databases: databases.to_owned(),
            profile: profile.to_owned(),
            execution_directory: None,
            cleanup
        })
    }
    pub async fn create(&mut self, samples: &Vec<StagedSample>) -> Result<PathBuf, TowerError> {

        let uuid = uuid::Uuid::new_v4().to_string();

        let path = self.workdir.join(uuid);
        let stage_path = path.join("stage");

        log::info!("Creating execution directory for Nextflow at: {}", path.display());

        create_dir_all(&path).await?;
        create_dir_all(&stage_path).await?;

        log::info!("Writing staged samples to execution directory (n = {})", samples.len());

        self.to_json(&samples, &stage_path).await?;

        self.execution_directory = Some(path.clone());

        Ok(path)

    }
    pub async fn cleanup(&self) -> Result<(), TowerError> {

        if let Some(ref dir) = self.execution_directory {
            if self.cleanup {
                log::info!("Cleaning up execution directory at: {}", dir.display());
                remove_dir_all(&dir).await?;
            }
        }

        Ok(())
    }
    pub async fn to_json(&self, staged_samples: &Vec<StagedSample>, path: &PathBuf) -> Result<(), TowerError> {

        for sample in staged_samples {
            let file_name = format!("{}.json", sample.id);
            let file_path = path.join(file_name);

            let json_data = serde_json::to_string_pretty(&sample)?;
            let mut file = File::create(file_path).await?;
            file.write_all(json_data.as_bytes()).await?;

            log::info!("Wrote staged sample: {}", sample.id);
        }
        
        Ok(())
    }
    async fn launch(&mut self, samples: &Vec<StagedSample>) -> Result<Child, TowerError> {

        let dir = self.create(samples).await?;
        log::info!("Executing pipeline at: {}", dir.display());

        let command = format!(
            "nextflow run {} -config {} {} -entry production --databaseDirectory {} --executionDirectory {}", 
            self.main.canonicalize()?.display(), 
            self.config.canonicalize()?.display(),
            if self.profile.is_empty() { String::new() } else { format!("-profile {}", self.profile.join(",")) },
            self.databases.canonicalize()?.display(),
            dir.canonicalize()?.display()
        );

        log::info!("Running execution command: '{}'", command);

        let process = Command::new("sh")
            .current_dir(
                dir.canonicalize()?
            )
            .arg("-c")
            .arg(&command)
            .spawn()?;

        Ok(process)
    }

}

#[derive(Clone, Debug)]
pub struct CerebroTower {
    client: TowerClient,
    nextflow: NextflowConfig
}

impl CerebroTower {
    pub fn new(
        client: TowerClient,
        nextflow: NextflowConfig,
    ) -> Result<Self, TowerError> {

        Ok(Self { client, nextflow })
    }

    pub async fn watch(&self, tower_id: &str, delete_from_stage: bool) -> Result<(), TowerError> {

        let (tx, rx) = mpsc::channel(32); 

        // Limit concurrent processes launched to execute Nextflow
        let semaphore = Arc::new(
            Semaphore::new(4)
        ); 

        // Start the watchtower task
        let watchtower_handle = task::spawn(
            watchtower_task(
                self.client.clone(),
                tower_id.to_string(),
                delete_from_stage,
                tx.clone(), 
                10,
            )
        );

        // Start the process manager task
        let process_manager_handle = task::spawn(
            process_manager_task(
            self.nextflow.clone(),
            self.client.clone(),
            tower_id.to_string(),
                rx, 
                semaphore
            )
        );

        // Handle graceful shutdown
        tokio::signal::ctrl_c()
            .await
            .expect("Failed to listen for shutdown signal");

        log::info!("Shutting down...");

        // Gracefully shut down the tasks
        watchtower_handle.abort();
        process_manager_handle.abort();

        Ok(())
    }
}

#[derive(Debug)]
enum WatchtowerCommand {
    StartProcess(Vec<StagedSample>),
    PingTower,
}

async fn watchtower_task(
    client: TowerClient, 
    tower_id: String, 
    delete_from_stage: bool, 
    tx: mpsc::Sender<WatchtowerCommand>, 
    interval: u64
) {
    
    loop {
        let response = client.pull_staged_samples(
            &tower_id, 
            delete_from_stage
        ).await;

        match response {
            Ok(samples) => {

                if should_start_process(&samples) {
                    let command = WatchtowerCommand::StartProcess(samples);
                    tx.send(command).await.unwrap();

                } else if should_ping_tower(&samples) {
                    let command = WatchtowerCommand::PingTower;
                    tx.send(command).await.unwrap();
                }

            }
            Err(e) => {
                log::error!("{}", e.to_string())
            }
        }

        sleep(Duration::from_secs(interval)).await;
    }
}

async fn process_manager_task(
    nextflow: NextflowConfig,
    client: TowerClient,
    tower_id: String, 
    mut rx: mpsc::Receiver<WatchtowerCommand>,
    semaphore: Arc<Semaphore>,
) {
    while let Some(command) = rx.recv().await {
        

        match command {
            WatchtowerCommand::StartProcess(data) => {
                let permit = semaphore.clone().acquire_owned().await.unwrap();

                let mut pipeline = nextflow.clone();
                // Spawn a new task to handle the process
                task::spawn(async move {
                    if let Ok(mut child) = pipeline.launch(&data).await {

                        match child.wait().await {
                            Err(e) => {
                                log::error!("Error handling process: {}", e)
                            },
                            Ok(exit_status) => {
                                log::info!("Process completed with {exit_status}");
                                if let Err(err) = pipeline.cleanup().await {
                                    log::error!("Failed to cleanup execution directory: {}", err.to_string())
                                }
                            }
                        }

                    }
                    // Permit is automatically released when the task completes
                    drop(permit);
                });
            }
            WatchtowerCommand::PingTower => {
                if let Err(e) = client.ping_tower(&tower_id, false).await {
                    log::error!("Failed to ping tower: {}", e.to_string())
                }
            }
        }
    }
}


fn should_start_process(staged_samples: &Vec<StagedSample>) -> bool {
    !staged_samples.is_empty()
}

fn should_ping_tower(staged_samples: &Vec<StagedSample>) -> bool {
    staged_samples.is_empty()
}
