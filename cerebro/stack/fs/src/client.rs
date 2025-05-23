use std::{collections::HashMap, path::PathBuf};
use cerebro_model::api::{files::model::FileType, stage::model::{FileId, StagedSample}};
use chrono::Utc;
use anyhow::Result;
use reqwest::StatusCode;

use cerebro_client::client::CerebroClient;
use cerebro_model::api::watchers::model::ProductionWatcher;
use cerebro_model::api::files::schema::RegisterFileSchema;
use crate::{error::FileSystemError, hash::fast_file_hash, weed::{weed_download, weed_upload}};


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
    // Needs improvements especially when the identifiers returned becomes larger
    pub fn delete_files(
        &self,
        file_ids: &Vec<String>,
        run_id: Option<String>,
        sample_id: Option<String>,
        all: bool
    ) -> Result<(), FileSystemError> {
        
        if all {
            let confirmation = dialoguer::Confirm::new()
                .with_prompt("Do you want to delete ALL files for your team?")
                .interact()
                .unwrap();

            if !confirmation {
                return Ok(())
            } else {
                let confirmation = dialoguer::Confirm::new()
                    .with_prompt("Really?? It is the nuclear option meant for development and testing!")
                    .interact()
                    .unwrap();

                if !confirmation {
                    return Ok(())
                }
            }
        }

        if file_ids.is_empty() {
            let deleted_fids = self.api_client.delete_files(run_id, sample_id, if all { Some(all) } else { None })?;
            for fid in deleted_fids {
                self.delete_file(&fid)?
            }
        } else {
            for file_id in file_ids {
                let deleted_file = self.api_client.delete_file(&file_id)?;
                self.delete_file(&deleted_file.fid)?;
            } 
        }

        Ok(())
    }

    pub fn stage_files(
        &self,
        json: &PathBuf,
        outdir: &PathBuf,
        pipeline: Option<PathBuf>,
    ) -> Result<(), FileSystemError> {

        let staged_sample = StagedSample::from_json(&json)?;
        
        for file in &staged_sample.files {
            weed_download(
                &file.fid, 
                outdir, 
                Some(self.fs_url.clone()), 
                Some(self.fs_port.clone())
            )?
        }

        staged_sample.to_json(&outdir.join(
            format!("{}.json", staged_sample.sample_id)
        ))?;

        if let Some(file) = pipeline {
            let mut writer = csv::WriterBuilder::new()
                .has_headers(false)
                .from_path(&file)?;
            
            writer.serialize(staged_sample.pipeline)?;
            writer.flush()?;
        }

        print!("{}", staged_sample.sample_id);

        Ok(())
    }
    pub fn upload_files(
        &self,
        files: &Vec<PathBuf>,
        run_id: Option<String>,
        sample_id: Option<String>,
        pipeline_id: Option<String>,
        description: Option<String>,
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
                pipeline_id: pipeline_id.clone(),
                description: description.clone(),
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
    ) -> Result<Vec<FileId>, FileSystemError> {

        let mut file_identifiers = Vec::new();
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

                let file_id = uuid::Uuid::new_v4().to_string();
                let file_schema = RegisterFileSchema {
                    id: file_id.clone(),
                    run_id: Some(run_id.clone()),
                    sample_id: Some(sample_id.clone()),
                    pipeline_id: None,
                    description: None,
                    date: Utc::now().to_string(),
                    name: upload_response.file_name,
                    hash: file_hash,
                    size: upload_response.size,
                    fid: upload_response.fid,
                    ftype: file_type.clone(),
                    watcher: Some(watcher.clone())
                };

                log::debug!("{:#?}", file_schema);

                log::info!("Registering file with Cerebro API");
                self.api_client.register_file(
                    file_schema
                )?;

                file_identifiers.push(file_id);
            }
            
        }

        Ok(file_identifiers)
    }
}