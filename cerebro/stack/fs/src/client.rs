use std::{collections::HashMap, path::PathBuf};
use cerebro_model::api::{files::model::{FileType, SeaweedFile}, stage::model::{FileId, StagedSample}};
use chrono::Utc;
use anyhow::Result;
use reqwest::StatusCode;

use cerebro_client::client::CerebroClient;
use cerebro_model::api::watchers::model::ProductionWatcher;
use cerebro_model::api::files::schema::RegisterFileSchema;
use cerebro_model::api::files::retention::{RestoreState, RetentionClass, RetentionPolicy, StorageTier};
use crate::config::{FsConfig, FsAccessMode};
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
    /// Storage tier to register the file under (default [`StorageTier::Hot`]).
    pub tier: StorageTier,
    /// Retention category to assign (default [`RetentionClass::Diagnostic`]).
    pub retention: RetentionClass,
    /// Whether to register the file under legal hold (exempt from expiry).
    pub legal_hold: bool,
    /// Policy used to compute `retain_until` from the retention category.
    pub retention_policy: RetentionPolicy,
}
impl Default for UploadConfig {
    fn default() -> Self {
        Self {
            data_center: None,
            max_mb: Some(16),
            replication: None,
            ttl: None,
            tier: StorageTier::default(),
            retention: RetentionClass::default(),
            legal_hold: false,
            retention_policy: RetentionPolicy::default(),
        }
    }
}

/// Result of storing a single object, abstracting over the access mode.
struct StoredObject {
    /// Weed fid (may be empty for filer-stored, path-addressed objects).
    fid: String,
    /// Filer object path, when stored via the filer.
    path: Option<String>,
    /// Stored file name.
    name: String,
    /// Object size in bytes.
    size: u64,
}

/// Build a filer object path of the form `run/sample/name`, substituting a
/// placeholder segment when the run or sample is unknown.
fn build_remote_path(run_id: Option<&str>, sample_id: Option<&str>, name: &str) -> String {
    let run = sanitize_segment(run_id.unwrap_or("_unassigned"));
    let sample = sanitize_segment(sample_id.unwrap_or("_"));
    let file = sanitize_segment(name);
    format!("{run}/{sample}/{file}")
}

/// Sanitise a single path segment: replace separators and reject empty/`.`-only
/// segments to avoid path traversal.
fn sanitize_segment(segment: &str) -> String {
    let replaced: String = segment
        .chars()
        .map(|c| if c == '/' || c == '\\' { '_' } else { c })
        .collect();
    let trimmed = replaced.trim_matches('.');
    if trimmed.is_empty() {
        "_".to_string()
    } else {
        trimmed.to_string()
    }
}

/// Outcome of a download request, separating files written to disk from those
/// that could not be retrieved because their data is archived (Glacier) and
/// requires a restore first.
#[derive(Debug, Clone)]
pub struct DownloadReport {
    /// Paths of files successfully written to the output directory.
    pub written: Vec<PathBuf>,
    /// Identifiers of files that require an archival restore before download.
    pub restore_pending: Vec<String>,
}
impl DownloadReport {
    /// Whether any requested file requires a restore before it can be retrieved.
    pub fn restore_required(&self) -> bool {
        !self.restore_pending.is_empty()
    }
}

/// Per-file result of a restore request: the object identifier, its
/// [`RestoreState`], and a human-readable message.
#[derive(Debug, Clone)]
pub struct RestoreOutcome {
    pub identifier: String,
    pub state: RestoreState,
    pub message: String,
}

/// Split resolved files into those directly retrievable and the identifiers of
/// those whose data is archived and must be restored first.
fn partition_archived(files: Vec<SeaweedFile>) -> (Vec<SeaweedFile>, Vec<String>) {
    let mut retrievable = Vec::new();
    let mut pending = Vec::new();
    for file in files {
        if file.requires_restore() {
            pending.push(file.effective_identifier().to_string());
        } else {
            retrievable.push(file);
        }
    }
    (retrievable, pending)
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
    ///   metadata is unavailable, `verify` is skipped (with a warning), and the
    ///   archival/restore state is unknown.
    /// * **By run/sample** — when `fids` is empty and `run_id` is provided, the
    ///   matching [`SeaweedFile`] records are listed from the Cerebro API
    ///   (optionally filtered by `sample_id`). Files whose data has been
    ///   archived to remote storage (Glacier, Model B) are **not** downloaded;
    ///   their identifiers are reported in [`DownloadReport::restore_pending`] so
    ///   the caller can run a restore first. The rest are downloaded, and when
    ///   `verify` is set each file's BLAKE3 hash is checked against the registered
    ///   value.
    ///
    /// Paired-end Illumina read sets are returned as their two constituent
    /// `ReadPaired` files; original file names are preserved so downstream
    /// pairing by name continues to work.
    pub fn download_files(
        &self,
        fids: &Vec<String>,
        run_id: Option<String>,
        sample_id: Option<String>,
        outdir: &PathBuf,
        verify: bool,
    ) -> Result<DownloadReport, FileSystemError> {

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
            return Ok(DownloadReport { written, restore_pending: Vec::new() });
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

        // Separate archived (Glacier) objects: these need a restore first.
        let (retrievable, restore_pending) = partition_archived(files);

        if !restore_pending.is_empty() {
            log::warn!(
                "{} object(s) are in archival storage and require a restore before download: {:?}",
                restore_pending.len(),
                restore_pending
            );
        }

        for file in &retrievable {
            let identifier = file.effective_identifier();
            log::info!("Downloading {} ({})", file.name, identifier);
            let path = self
                .download_identifier(identifier, Some(file.name.as_str()), outdir)?
                .unwrap_or_else(|| outdir.join(&file.name));

            if verify {
                log::info!("Verifying BLAKE3 hash for {}", file.name);
                verify_file(&path, &file.hash)?;
                log::info!("Integrity verified: {}", file.name);
            }
            written.push(path);
        }

        Ok(DownloadReport { written, restore_pending })
    }

    /// Report which files in a run/sample require an archival restore, and
    /// surface the restore contract for each.
    ///
    /// Resolution is by run (optionally narrowed by sample). For each archived
    /// file this returns [`RestoreState::Pending`]; for directly retrievable
    /// files it returns [`RestoreState::NotRequired`].
    ///
    /// FS-4 establishes the contract: the actual S3 Glacier `RestoreObject`
    /// execution is an operational step (or a future `s3` cargo feature). This
    /// method tells the caller precisely what must be restored and lets
    /// [`download_files`](Self::download_files) avoid blocking on Glacier
    /// objects in the meantime.
    pub fn restore_files(
        &self,
        run_id: Option<String>,
        sample_id: Option<String>,
    ) -> Result<Vec<RestoreOutcome>, FileSystemError> {
        let run_id = run_id.ok_or(FileSystemError::InvalidDownloadQuery)?;

        let files = self.api_client.list_files(Some(run_id), None, 0, 100_000, false)?;
        let files: Vec<SeaweedFile> = match &sample_id {
            Some(sid) => files
                .into_iter()
                .filter(|f| f.sample_id.as_deref() == Some(sid.as_str()))
                .collect(),
            None => files,
        };

        let mut outcomes = Vec::new();
        for file in &files {
            let identifier = file.effective_identifier().to_string();
            if file.requires_restore() {
                log::info!(
                    "Restore required for {} ({}); initiate archival (S3 Glacier) restore",
                    file.name,
                    identifier
                );
                outcomes.push(RestoreOutcome {
                    identifier,
                    state: RestoreState::Pending,
                    message: "data is archived; restore must complete before retrieval".to_string(),
                });
            } else {
                outcomes.push(RestoreOutcome {
                    identifier,
                    state: RestoreState::NotRequired,
                    message: "directly retrievable".to_string(),
                });
            }
        }
        Ok(outcomes)
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

    /// Store a single local file in SeaweedFS according to the configured
    /// access mode, returning the resulting identifier(s).
    ///
    /// * [`FsAccessMode::Weed`] (default) — upload via `weed upload`, which
    ///   chunks large files into a single manifest fid; `path` is `None`.
    /// * [`FsAccessMode::Filer`] — stream the file to the filer at a
    ///   `run/sample/name` path; `path` is set and later used for retrieval via
    ///   [`SeaweedFile::effective_identifier`](cerebro_model::api::files::model::SeaweedFile::effective_identifier).
    ///
    /// Both paths are safe for multi-gigabyte read sets: `weed` chunks on the
    /// client side and the filer streams from disk and auto-chunks server-side.
    fn store_object(
        &self,
        file: &PathBuf,
        run_id: Option<&str>,
        sample_id: Option<&str>,
        upload_config: &UploadConfig,
    ) -> Result<StoredObject, FileSystemError> {
        match self.config.access {
            FsAccessMode::Filer => {
                let name = file
                    .file_name()
                    .and_then(|s| s.to_str())
                    .ok_or(FileSystemError::FileNameExtraction)?
                    .to_string();
                let remote_path = build_remote_path(run_id, sample_id, &name);

                let filer = FilerClient::new(
                    &self.config.filer_base(),
                    self.config.danger_invalid_certificate,
                )?;
                let response = filer.upload(file, &remote_path)?;

                let size = match response.size {
                    Some(size) => size,
                    None => std::fs::metadata(file)?.len(),
                };

                Ok(StoredObject {
                    fid: response.fid.unwrap_or_default(),
                    path: Some(remote_path),
                    name,
                    size,
                })
            }
            FsAccessMode::Weed => {
                let response = weed_upload(
                    file,
                    upload_config.data_center.clone(),
                    None,
                    Some(self.config.master_url.clone()),
                    Some(self.config.master_port.clone()),
                    upload_config.max_mb,
                    None,
                    upload_config.replication.clone(),
                    upload_config.ttl.clone(),
                    false,
                )?;
                Ok(StoredObject {
                    fid: response.fid,
                    path: None,
                    name: response.file_name,
                    size: response.size,
                })
            }
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

            log::info!("Uploading file to SeaweedFS storage ({} access)", self.config.access);
            let stored = self.store_object(
                file,
                run_id.as_deref(),
                sample_id.as_deref(),
                &upload_config,
            )?;

            let now = Utc::now();
            let retain_until = upload_config
                .retention_policy
                .retain_until(upload_config.retention, now);

            let file_schema = RegisterFileSchema {
                id: uuid::Uuid::new_v4().to_string(),
                run_id: run_id.clone(),
                sample_id: sample_id.clone(),
                pipeline_id: pipeline_id.clone(),
                description: description.clone(),
                date: now.to_string(),
                name: stored.name,
                hash: file_hash,
                fid: stored.fid,
                ftype: file_type.clone(),
                size: stored.size,
                watcher: watcher.clone(),
                path: stored.path,
                tier: upload_config.tier,
                retention: upload_config.retention,
                retain_until,
                legal_hold: upload_config.legal_hold,
                replicas: None,
                archived: false,
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

                log::info!("Uploading file to SeaweedFS storage ({} access)", self.config.access);
                let stored = self.store_object(
                    file,
                    Some(run_id.as_str()),
                    Some(sample_id.as_str()),
                    &upload_config,
                )?;

                let now = Utc::now();
                let retain_until = upload_config
                    .retention_policy
                    .retain_until(upload_config.retention, now);

                let file_id = uuid::Uuid::new_v4().to_string();
                let file_schema = RegisterFileSchema {
                    id: file_id.clone(),
                    run_id: Some(run_id.clone()),
                    sample_id: Some(sample_id.clone()),
                    pipeline_id: None,
                    description: None,
                    date: now.to_string(),
                    name: stored.name,
                    hash: file_hash,
                    size: stored.size,
                    fid: stored.fid,
                    ftype: file_type.clone(),
                    watcher: Some(watcher.clone()),
                    path: stored.path,
                    tier: upload_config.tier,
                    retention: upload_config.retention,
                    retain_until,
                    legal_hold: upload_config.legal_hold,
                    replicas: None,
                    archived: false,
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
    use super::{build_remote_path, is_filer_path, sanitize_segment};

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

    #[test]
    fn remote_path_uses_run_sample_name() {
        assert_eq!(
            build_remote_path(Some("RUN01"), Some("S1"), "reads_R1.fastq.gz"),
            "RUN01/S1/reads_R1.fastq.gz"
        );
    }

    #[test]
    fn remote_path_substitutes_missing_segments() {
        assert_eq!(
            build_remote_path(None, None, "reads.fastq.gz"),
            "_unassigned/_/reads.fastq.gz"
        );
    }

    #[test]
    fn sanitize_segment_blocks_traversal() {
        assert_eq!(sanitize_segment("../etc"), "_etc");
        assert_eq!(sanitize_segment(".."), "_");
        assert_eq!(sanitize_segment("a/b"), "a_b");
    }

    fn file_with(identifier: &str, archived: bool) -> cerebro_model::api::files::model::SeaweedFile {
        use cerebro_model::api::files::retention::{RetentionClass, StorageTier};
        cerebro_model::api::files::model::SeaweedFile {
            id: "id".into(),
            date: "2025-01-01".into(),
            name: "reads.fastq.gz".into(),
            hash: "h".into(),
            size: 1,
            fid: identifier.into(),
            tags: Vec::new(),
            run_id: None,
            sample_id: None,
            ftype: None,
            watcher: None,
            path: None,
            tier: StorageTier::Cold,
            retention: RetentionClass::Diagnostic,
            retain_until: None,
            legal_hold: false,
            replicas: None,
            archived,
        }
    }

    #[test]
    fn partition_archived_separates_glacier_objects() {
        use super::partition_archived;
        let files = vec![
            file_with("3,01", false),
            file_with("3,02", true),
            file_with("3,03", false),
        ];
        let (retrievable, pending) = partition_archived(files);
        assert_eq!(retrievable.len(), 2);
        assert_eq!(pending, vec!["3,02".to_string()]);
    }
}
