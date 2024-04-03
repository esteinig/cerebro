use std::path::PathBuf;
use chrono::Utc;
use anyhow::Result;
use reqwest::StatusCode;

use cerebro_client::client::CerebroClient;
use cerebro_model::api::files::model::WatcherConfig;
use cerebro_model::api::files::schema::RegisterFileSchema;
use crate::{error::FileSystemError, hash::fast_file_hash, weed::weed_upload};


pub struct FileSystemClient {
    api_client: CerebroClient,
    fs_url: String,
    fs_port: String
}

pub struct UploadConfig {
    data_center: Option<String>,
    max_mb: Option<i32>,
    replication: Option<String>,
    ttl: Option<String>,
}
impl Default for UploadConfig {
    fn default() -> Self {
        Self {
            data_center: None,
            max_mb: Some(16),
            replication: None,
            ttl: None
        }
    }
}

impl FileSystemClient {
     /// Creates a new instance of `FileSystemClient`.
    ///
    /// # Arguments
    ///
    /// * `api_client` - An instance of `CerebroClient`.
    /// * `weed_master` - The address of the SeaweedFS master server.
    pub fn new(api_client: &CerebroClient, fs_url: &str, fs_port: &str) -> Self {
        Self {
            api_client: api_client.clone(),
            fs_url: fs_url.to_string(),
            fs_port: fs_port.to_string()
        }
    }

    /// Pings the SeaweedFS cluster to check its health status.
    ///
    /// Makes a HEAD request to `/cluster/healthz` endpoint of the SeaweedFS master server.
    /// A successful response indicates the cluster is healthy.
    ///
    /// # Returns
    ///
    /// * `Ok(())` if the cluster is healthy.
    /// * `Err(FileSystemError)` if the cluster is unhealthy or if a network error occurs.
    pub fn ping_status(&self) -> Result<(), FileSystemError> {
        let url = format!("{}/cluster/status", self.fs_url);
        
        let response = reqwest::blocking::Client::new()
            .get(&url)
            .send()?;

        match response.status() {
            StatusCode::OK => {
                log::info!("Cerebro FS status: ok");
                Ok(())
            },
            StatusCode::SERVICE_UNAVAILABLE => Err(FileSystemError::UnhealthyCluster),
            status => Err(FileSystemError::UnexpectedResponseStatus(status)),
        }
    }
    pub fn upload_files(
        &self,
        files: &Vec<PathBuf>,
        team_name: &str,
        db_name: &str,
        upload_config: UploadConfig,
        watcher_config: WatcherConfig,
    ) -> Result<(), FileSystemError> {

        for file in files { 

            if !file.exists() {
                return Err(FileSystemError::FileNotExist(file.display().to_string()));
            }

            log::info!("Generating file hash with BLAKE3");
            let file_hash = fast_file_hash(&file)?;

            log::info!("Uploading file to SeaweedFS storage");
            let upload_response = weed_upload(
                file, 
                upload_config.data_center.clone(),
                None,
                Some(self.fs_url.clone()),
                Some(self.fs_port.clone()),
                upload_config.max_mb,
                None,
                upload_config.replication.clone(),
                upload_config.ttl.clone(),
                false
            )?;

            let file_schema = RegisterFileSchema {
                id: uuid::Uuid::new_v4().to_string(),
                date: Utc::now().to_string(),
                name: upload_response.file_name,
                hash: file_hash,
                fid: upload_response.fid,
                size: upload_response.size,
                watcher: watcher_config.clone()
            };

            log::info!("Registering file with Cerebro API");
            self.api_client.register_file(
                file_schema,
                team_name,
                db_name
            )?;
        }

        Ok(())
    }
}