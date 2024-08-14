use std::{collections::HashMap, path::PathBuf};
use cerebro_model::api::files::model::FileType;
use chrono::Utc;
use anyhow::Result;
use reqwest::StatusCode;

use cerebro_client::client::CerebroClient;
use cerebro_model::api::{pipelines::model::ProductionPipeline, watchers::model::ProductionWatcher};
use cerebro_model::api::files::schema::RegisterFileSchema;
use crate::{error::FileSystemError, hash::fast_file_hash, weed::weed_upload};


#[derive(Clone, Debug)]
pub struct FileSystemClient {
    pub api_client: CerebroClient,
    pub fs_url: String,
    pub fs_port: String
}

#[derive(Clone, Debug)]
pub struct UploadConfig {
    pub data_center: Option<String>,
    pub max_mb: Option<i32>,
    pub replication: Option<String>,
    pub ttl: Option<String>,
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
        let url = format!("{}/cluster/healthz", self.fs_url);
        
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
    // Cerebro FS delete file request
    pub fn delete_file(&self, fid: &str) -> Result<(), FileSystemError> {
        let url = format!("{}/{}", self.fs_url, fid);
        
        let response = reqwest::blocking::Client::new()
            .delete(&url)
            .send()?;

        match response.status() {
            StatusCode::OK | StatusCode::ACCEPTED => {
                log::info!("File data deleted from Cerebro FS ({fid})");
                Ok(())
            },
            StatusCode::SERVICE_UNAVAILABLE => Err(FileSystemError::UnhealthyCluster),
            status => Err(FileSystemError::UnexpectedResponseStatus(status)),
        }
    }
    // Cerebro API file entry deletion followed by Cerebro FS file deletion
    //
    // Needs improvements especially when file storage gets large!
    pub fn delete_files(
        &self,
        file_ids: &Vec<String>,
        run_id: Option<String>,
        watcher_id: Option<String>
    ) -> Result<(), FileSystemError> {

        let file_ids = match (&run_id, &watcher_id) {
            (Some(_), Some(_)) | (Some(_), None) | (None, Some(_)) => {
                self.api_client.list_files(run_id, watcher_id, 0, 1000, false)?
                    .iter()
                    .map(|file| file.id.to_owned())
                    .collect()

            },
            _ => file_ids.clone()
        };
           

        for file_id in file_ids {
            let deleted_file = self.api_client.delete_file(Some(file_id), None, None)?;
            self.delete_file(&deleted_file.fid)?;
        } 

        Ok(())
    }
    pub fn upload_files(
        &self,
        files: &Vec<PathBuf>,
        run_id: Option<String>,
        sample_id: Option<String>,
        file_type: Option<FileType>,
        upload_config: UploadConfig,
        watcher: Option<ProductionWatcher>
    ) -> Result<(), FileSystemError> {

        for file in files { 

            if !file.exists() {
                return Err(FileSystemError::FileDoesNotExist(file.to_owned()));
            }

            log::info!("Generating file hash with BLAKE3: {}", file.display());
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
                run_id: run_id.clone(),
                sample_id: sample_id.clone(),
                date: Utc::now().to_string(),
                name: upload_response.file_name,
                hash: file_hash,
                fid: upload_response.fid,
                ftype: file_type.clone(),
                size: upload_response.size,
                watcher: watcher.clone()
            };

            log::info!("Registering file with Cerebro API");
            self.api_client.register_file(
                file_schema
            )?;
        }

        Ok(())
    }
    pub fn upload_files_from_watcher(
        &self,
        files: &HashMap<String, Vec<PathBuf>>,
        run_id: String,
        file_type: Option<FileType>,
        upload_config: UploadConfig,
        watcher: ProductionWatcher,
    ) -> Result<(), FileSystemError> {

        for (sample_id, files) in files { 

            for file in files {
                if !file.exists() {
                    return Err(FileSystemError::FileDoesNotExist(file.to_owned()));
                }

                log::info!("Generating file hash with BLAKE3: {}", file.display());
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
                    run_id: Some(run_id.clone()),
                    sample_id: Some(sample_id.clone()),
                    date: Utc::now().to_string(),
                    name: upload_response.file_name,
                    hash: file_hash,
                    size: upload_response.size,
                    fid: upload_response.fid,
                    ftype: file_type.clone(),
                    watcher: Some(watcher.clone())
                };

                log::info!("{:#?}", file_schema);

                log::info!("Registering file with Cerebro API");
                self.api_client.register_file(
                    file_schema
                )?;
            }
            
        }

        Ok(())
    }
}