use std::{collections::HashMap, path::PathBuf};
use cerebro_model::api::{files::model::{FileType, SeaweedFile}, stage::model::{FileId, StagedSample}};
use chrono::Utc;
use anyhow::Result;
use reqwest::StatusCode;

use cerebro_client::client::CerebroClient;
use cerebro_model::api::watchers::model::ProductionWatcher;
use cerebro_model::api::files::schema::RegisterFileSchema;
use crate::config::FsConfig;
use crate::filer::FilerClient;
use crate::{error::FileSystemError, hash::fast_file_hash, weed::{weed_download, weed_upload}};


#[derive(Clone, Debug)]
pub struct FileSystemClient {
    pub api_client: CerebroClient,
    pub config: FsConfig,
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

/// Heuristic used to route a stored identifier to the correct retrieval backend.
///
/// A SeaweedFS fid has the form `<volumeId>,<fileKey>` and never contains a `/`,
/// whereas a filer object is addressed by a path that always does. This lets
/// [`FileSystemClient::download_files`] transparently retrieve both legacy
/// fid-addressed objects (the current default) and future path-addressed filer
/// objects without the caller needing to know which is which.
///
/// # Examples
///
/// ```
/// use cerebro_fs::client::is_filer_path;
///
/// assert!(!is_filer_path("3,01637037d6"));
/// assert!(is_filer_path("/run01/sample01/reads_R1.fastq.gz"));
/// ```
pub fn is_filer_path(identifier: &str) -> bool {
    identifier.contains('/')
}

/// Recompute a file's BLAKE3 hash and compare it against the expected value.
///
/// Returns [`FileSystemError::IntegrityMismatch`] on a mismatch. This is the
/// basic verification primitive; replica-aware retry and `last_verified`
/// stamping are added in FS-6.
fn verify_file(path: &PathBuf, expected_hash: &str) -> Result<(), FileSystemError> {
    let actual = fast_file_hash(path)?;
    if actual != expected_hash {
        return Err(FileSystemError::IntegrityMismatch {
            path: path.display().to_string(),
            expected: expected_hash.to_string(),
            actual,
        });
    }
    Ok(())
}

impl FileSystemClient {
    /// Creates a new [`FileSystemClient`] using the legacy positional
    /// parameters.
    ///
    /// Retained for backwards compatibility so existing callers
    /// (`cerebro-watcher`, `cerebro-ciqa`) continue to compile unchanged. The
    /// resulting client uses [`crate::config::FsAccessMode::Weed`]. New code
    /// should prefer [`FileSystemClient::with_config`].
    ///
    /// # Arguments
    ///
    /// * `api_client` - An instance of `CerebroClient`.
    /// * `fs_url` - The SeaweedFS master base URL (e.g. `http://localhost`).
    /// * `fs_port` - The SeaweedFS master port (e.g. `9333`).
    /// * `localhost` - Whether to address the master as `host:port`.
    pub fn new(api_client: &CerebroClient, fs_url: &str, fs_port: &str, localhost: bool) -> Self {
        Self::with_config(api_client, FsConfig::weed(fs_url, fs_port, localhost))
    }

    /// Creates a new [`FileSystemClient`] from an explicit [`FsConfig`].
    pub fn with_config(api_client: &CerebroClient, config: FsConfig) -> Self {
        Self {
            api_client: api_client.clone(),
            config,
        }
    }

    /// Pings the SeaweedFS cluster to check its health status.
    ///
    /// Makes a request to the `/cluster/healthz` endpoint of the SeaweedFS
    /// master server. A successful response indicates the cluster is healthy.
    ///
    /// # Returns
    ///
    /// * `Ok(())` if the cluster is healthy.
    /// * `Err(FileSystemError)` if the cluster is unhealthy or a network error occurs.
    pub fn ping_status(&self) -> Result<(), FileSystemError> {
        let url = format!("{}/cluster/healthz", self.get_url());

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

    /// Resolve the master HTTP base URL (localhost-aware), delegating to
    /// [`FsConfig::master_health_url`].
    pub fn get_url(&self) -> String {
        self.config.master_health_url()
    }

    /// Cerebro FS delete file request (direct volume/master delete by fid).
    pub fn delete_file(&self, fid: &str) -> Result<(), FileSystemError> {
        let url = format!("{}/{}", self.get_url(), fid);

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
                Some(self.config.master_url.clone()),
                Some(self.config.master_port.clone())
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

    /// Download files from Cerebro FS to a local directory.
    ///
    /// Targets are resolved one of two ways:
    ///
    /// * **By identifier** — when `fids` is non-empty, each entry is fetched
    ///   directly. Entries are routed by [`is_filer_path`]: a value containing
    ///   `/` is treated as a filer path and retrieved via [`FilerClient`];
    ///   otherwise it is treated as a SeaweedFS fid and retrieved with
    ///   `weed download`. No Cerebro API lookup is performed, so integrity
    ///   metadata is unavailable and `verify` is skipped (with a warning).
    /// * **By run/sample** — when `fids` is empty and `run_id` is provided, the
    ///   matching [`SeaweedFile`] records are listed from the Cerebro API
    ///   (optionally filtered by `sample_id`) and each is downloaded by its
    ///   stored identifier. When `verify` is set, each file's BLAKE3 hash is
    ///   recomputed and compared against the registered hash.
    ///
    /// Paired-end Illumina read sets are returned as their two constituent
    /// `ReadPaired` files; original file names are preserved so downstream
    /// pairing by name continues to work.
    ///
    /// Returns the list of written file paths.
    pub fn download_files(
        &self,
        fids: &Vec<String>,
        run_id: Option<String>,
        sample_id: Option<String>,
        outdir: &PathBuf,
        verify: bool,
    ) -> Result<Vec<PathBuf>, FileSystemError> {

        if !outdir.exists() {
            std::fs::create_dir_all(outdir)?;
        }

        let mut written = Vec::new();

        // Direct identifier download (no API metadata available)
        if !fids.is_empty() {
            for id in fids {
                if verify {
                    log::warn!("Integrity verification skipped for direct identifier (no registered hash): {id}");
                }
                if let Some(path) = self.download_identifier(id, None, outdir)? {
                    written.push(path);
                }
            }
            return Ok(written);
        }

        // Run/sample download via the Cerebro API
        let run_id = run_id.ok_or(FileSystemError::InvalidDownloadQuery)?;

        log::info!("Listing files for run '{run_id}' from Cerebro API");
        let files = self.api_client.list_files(Some(run_id), None, 0, 100_000, false)?;

        let files: Vec<SeaweedFile> = match &sample_id {
            Some(sid) => files
                .into_iter()
                .filter(|f| f.sample_id.as_deref() == Some(sid.as_str()))
                .collect(),
            None => files,
        };

        if files.is_empty() {
            log::warn!("No files matched the requested run/sample");
        }

        for file in &files {
            log::info!("Downloading {} ({})", file.name, file.fid);
            let path = self
                .download_identifier(&file.fid, Some(file.name.as_str()), outdir)?
                .unwrap_or_else(|| outdir.join(&file.name));

            if verify {
                log::info!("Verifying BLAKE3 hash for {}", file.name);
                verify_file(&path, &file.hash)?;
                log::info!("Integrity verified: {}", file.name);
            }
            written.push(path);
        }

        Ok(written)
    }

    /// Download a single identifier (filer path or weed fid) into `outdir`.
    ///
    /// When the identifier is a filer path the resolved local file path is
    /// returned. For a weed fid the file is written by `weed download` using its
    /// stored name; the resolved path is returned only when `name` is known.
    fn download_identifier(
        &self,
        id: &str,
        name: Option<&str>,
        outdir: &PathBuf,
    ) -> Result<Option<PathBuf>, FileSystemError> {
        if is_filer_path(id) {
            let file_name = name
                .map(|s| s.to_string())
                .or_else(|| id.rsplit('/').next().map(|s| s.to_string()))
                .ok_or(FileSystemError::FileNameExtraction)?;
            let target = outdir.join(&file_name);

            let filer = FilerClient::new(
                &self.config.filer_base(),
                self.config.danger_invalid_certificate,
            )?;
            filer.download(id, &target)?;
            Ok(Some(target))
        } else {
            weed_download(
                id,
                outdir,
                Some(self.config.master_url.clone()),
                Some(self.config.master_port.clone()),
            )?;
            Ok(name.map(|n| outdir.join(n)))
        }
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
                Some(self.config.master_url.clone()),
                Some(self.config.master_port.clone()),
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
                    Some(self.config.master_url.clone()),
                    Some(self.config.master_port.clone()),
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

#[cfg(test)]
mod tests {
    use super::is_filer_path;

    #[test]
    fn weed_fid_is_not_a_filer_path() {
        assert!(!is_filer_path("3,01637037d6"));
        assert!(!is_filer_path("7,0a1b2c3d4e"));
    }

    #[test]
    fn filer_paths_are_detected() {
        assert!(is_filer_path("/run01/sample01/reads_R1.fastq.gz"));
        assert!(is_filer_path("run01/reads.fastq.gz"));
    }
}
