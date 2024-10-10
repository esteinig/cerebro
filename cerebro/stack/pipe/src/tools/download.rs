use core::fmt;
use std::fs::{create_dir_all, remove_file, File};
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use std::time::Duration;
use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};
use tar::Archive;

use crate::error::WorkflowError;
use crate::nextflow::pathogen::{Aligner, Classifier};


/// Represents different indices available for Cerebro.
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum CerebroIndex {
    Chm13v2,
    Virus
}

impl CerebroIndex {
    /// Returns the aligner name formatted for the specified index.
    ///
    /// # Arguments
    ///
    /// * `aligner` - A reference to an `Aligner`.
    ///
    /// # Example
    ///
    /// ```
    /// let index = CerebroIndex::Virus;
    /// let aligner = Aligner::new();
    /// let name = index.aligner_name(&aligner);
    /// ```
    pub fn aligner_name(&self, aligner: &Aligner, version: Option<String>) -> String {
        match version {
            Some(version) => format!("{}-{}.{}.tar.xz", self, aligner.short_name(), version),
            None => format!("{}.{}.tar.xz", self, aligner.short_name())
        }
    }
    /// Returns the classifier name formatted for the specified index.
    ///
    /// # Arguments
    ///
    /// * `classifier` - A reference to a `Classifier`.
    ///
    /// # Example
    ///
    /// ```
    /// let index = CerebroIndex::Virus;
    /// let classifier = Classifier::new();
    /// let name = index.classifier_name(&classifier, None);
    /// ```
    pub fn classifier_name(&self, classifier: &Classifier, version: Option<String>) -> String {
        match version {
            Some(version) => format!("{}-{}.{}.tar.xz", self, classifier.short_name(), version),
            None => format!("{}.{}.tar.xz", self, classifier.short_name())
        }
        
    }
    /// Returns the reference name formatted for the specified index.
    ///
    /// # Arguments
    ///
    /// * `reference` - A reference to a `Reference`.
    ///
    /// # Example
    ///
    /// ```
    /// let index = CerebroIndex::Virus;
    /// let classifier = Classifier::new();
    /// let name = index.reference_name(&classifier, None);
    /// ```
    pub fn reference_name(&self, version: Option<String>) -> String {
        match version {
            Some(version) => format!("{}-{}.tar.xz", self, version),
            None => format!("{}.tar.xz", self)
        }
    }
}

impl fmt::Display for CerebroIndex {
    /// Formats the CerebroIndex for display.
    ///
    /// # Arguments
    ///
    /// * `f` - A mutable reference to a `fmt::Formatter`.
    ///
    /// # Example
    ///
    /// ```
    /// let index = CerebroIndex::Virus;
    /// println!("{}", index);
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CerebroIndex::Chm13v2 => write!(f, "chm13v2"),
            CerebroIndex::Virus => write!(f, "virus"),
        }
    }
}

/// Manages the download process for Cerebro indices.
pub struct CerebroDownloader {
    pub outdir: PathBuf,
    pub base_url: String,
    pub username: String,
    pub password: String,
    pub client: Client,
    pub timeout: u64,
    pub indices: Vec<CerebroIndex>,
    pub aligners: Vec<Aligner>,
    pub classifiers: Vec<Classifier>,
    pub version: Option<String>,
    pub reference: bool
}

impl CerebroDownloader {
    /// Creates a new instance of CerebroDownloader.
    ///
    /// # Arguments
    ///
    /// * `outdir` - Output directory for downloaded files.
    /// * `indices` - A list of `CerebroIndex` to download.
    ///
    /// # Errors
    ///
    /// Returns a `WorfkflowError` if the downloader cannot be created.
    ///
    /// # Example
    ///
    /// ```
    /// let outdir = PathBuf::from("/path/to/output");
    /// let indices = vec![CerebroIndex::Virus];
    /// let downloader = CerebroDownloader::new(outdir, indices);
    /// ```
    pub fn new(outdir: PathBuf, indices: Vec<CerebroIndex>) -> Result<Self, WorkflowError> {
        CerebroDownloaderBuilder::new(outdir, indices).build()
    }
    /// Creates a new instance of CerebroDownloaderBuilder.
    ///
    /// # Arguments
    ///
    /// * `outdir` - Output directory for downloaded files.
    /// * `indices` - A list of `CerebroIndex` to download.
    ///
    /// # Example
    ///
    /// ```
    /// let outdir = PathBuf::from("/path/to/output");
    /// let indices = vec![CerebroIndex::Virus];
    /// let builder = CerebroDownloader::builder(outdir, indices);
    /// ```
    pub fn builder(outdir: PathBuf, indices: Vec<CerebroIndex>) -> CerebroDownloaderBuilder {
        CerebroDownloaderBuilder::new(outdir, indices)
    }
    /// Lists the available index names for download.
    ///
    /// # Example
    ///
    /// ```
    /// let downloader = CerebroDownloader::new(outdir, indices);
    /// downloader.list();
    /// ```
    pub fn list(&self) {

        log::info!("===========================================");
        log::info!("Reference index storage for anonymous users");
        log::info!("===========================================");
        log::info!("                                           ");
        log::info!("Pre-built indices are available at:        ");
        log::info!("                                           ");
        log::info!("URL: {}                                    ", self.base_url);
        log::info!("Username: '{}' Password: '{}'              ", self.username, self.password);
        log::info!("                                           ");
        log::info!("===========================================");
        log::info!("Available index names for download (--name)");
        log::info!("===========================================");
        log::info!("                                           ");
        log::info!("Host reference => {:>16}                   ", CerebroIndex::Chm13v2);
        log::info!("Viral genomes  => {:>16}                   ", CerebroIndex::Virus);
        log::info!("                                           ");
    }
    /// Downloads the specified indices.
    ///
    /// # Errors
    ///
    /// Returns a `WorfkflowError` if any download or unpacking operation fails.
    ///
    /// # Example
    ///
    /// ```
    /// let downloader = CerebroDownloader::new(outdir, indices);
    /// downloader.download_index();
    /// ```
    pub fn download_index(&self) -> Result<(), WorkflowError> {

        if self.indices.is_empty() {
            log::warn!("No index names provided for download")
        }

        for index in &self.indices {
            for aligner in &self.aligners {
                let file_path = self.outdir.join(index.aligner_name(&aligner, self.version.clone()));
                log::info!("Downloading alignment index to file: {}", file_path.display());
                self.download(&index.aligner_name(aligner, self.version.clone()), &file_path)?;
                log::info!("Unpacking alignment index to directory: {}", self.outdir.display());
                self.unpack(&file_path, &self.outdir)?;
                log::info!("Removing download: {}", file_path.display());
                remove_file(&file_path)?;
            }
            for classifier in &self.classifiers {
                let file_path = self.outdir.join(index.classifier_name(&classifier, self.version.clone()));
                log::info!("Downloading classifier index to file: {}", file_path.display());
                self.download(&index.classifier_name(classifier, self.version.clone()), &file_path)?;
                log::info!("Unpacking classifier index to directory: {}", self.outdir.display());
                self.unpack(&file_path, &self.outdir)?;
                log::info!("Removing download: {}", file_path.display());
                remove_file(&file_path)?;
            }
            if self.reference {
                let file_path = self.outdir.join(index.reference_name(self.version.clone()));
                log::info!("Downloading reference sequences to file: {}", file_path.display());
                self.download(&index.reference_name(self.version.clone()), &file_path)?;
                log::info!("Unpacking reference sequence file to directory: {}", self.outdir.display());
                self.unpack(&file_path, &self.outdir)?;
                log::info!("Removing download: {}", file_path.display());
                remove_file(&file_path)?;
            }
        }

        Ok(())
    }
    /// Unpacks the downloaded file to the specified output directory.
    ///
    /// # Arguments
    ///
    /// * `download` - Path to the downloaded file.
    /// * `outdir` - Output directory for unpacked files.
    ///
    /// # Errors
    ///
    /// Returns a `WorfkflowError` if the unpacking operation fails.
    ///
    /// # Example
    ///
    /// ```
    /// let downloader = CerebroDownloader::new(outdir, indices);
    /// downloader.unpack(&download_path, &outdir);
    /// ```
    pub fn unpack(&self, download: &PathBuf, outdir: &PathBuf) -> Result<(), WorkflowError> {
        let file = File::open(download)?;
        let buf_reader = BufReader::new(file);
        let (reader, _compression) = niffler::get_reader(Box::new(buf_reader))?;

        let mut archive = Archive::new(reader);

        for entry in archive.entries()? {
            let mut entry = entry?;
            entry.unpack_in(outdir)?;
        }

        Ok(())
    }
    /// Downloads a file from the specified URL to the given path.
    ///
    /// # Arguments
    ///
    /// * `file_name` - The name of the file to download.
    /// * `path` - The path where the file should be saved.
    ///
    /// # Errors
    ///
    /// Returns a `WorfkflowError` if the download operation fails.
    ///
    /// # Example
    ///
    /// ```
    /// let downloader = CerebroDownloader::new(outdir, indices);
    /// downloader.download("file_name.tar.xz", &path);
    /// ```
    pub fn download(&self, file_name: &str, path: &PathBuf) -> Result<(), WorkflowError> {
        let url = format!("{}/{}", self.base_url, file_name);

        let mut response = self.client.get(&url)
            .basic_auth(&self.username, Some(&self.password))
            .timeout(Duration::from_secs(self.timeout*60))
            .send()?;

        if !response.status().is_success() {
            return Err(WorkflowError::DownloadFailedRequest(response.status()));
        }

        let mut writer = BufWriter::new(File::create(path)?);
        response.copy_to(&mut writer)?;

        Ok(())
    }
}

/// Builder for creating an instance of `CerebroDownloader`.
pub struct CerebroDownloaderBuilder {
    indices: Vec<CerebroIndex>,
    outdir: PathBuf,
    base_url: Option<String>,
    timeout: Option<u64>,
    username: Option<String>,
    password: Option<String>,
    aligners: Option<Vec<Aligner>>,
    classifiers: Option<Vec<Classifier>>,
    version: Option<String>,
    reference: bool
}

impl CerebroDownloaderBuilder {
    /// Creates a new instance of CerebroDownloaderBuilder.
    ///
    /// # Arguments
    ///
    /// * `outdir` - Output directory for downloaded files.
    /// * `indices` - A list of `CerebroIndex` to download.
    ///
    /// # Example
    ///
    /// ```
    /// let outdir = PathBuf::from("/path/to/output");
    /// let indices = vec![CerebroIndex::Virus];
    /// let builder = CerebroDownloaderBuilder::new(outdir, indices);
    /// ```
    pub fn new(outdir: PathBuf, indices: Vec<CerebroIndex>) -> Self {
        Self {
            outdir,
            indices,
            base_url: None,
            username: None,
            password: None,
            aligners: None,
            classifiers: None,
            timeout: None,
            version: None,
            reference: false
        }
    }
    /// Sets the aligners for the builder.
    ///
    /// # Arguments
    ///
    /// * `aligner` - A list of `Aligner` instances.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = CerebroDownloaderBuilder::new(outdir, indices).aligner(vec![Aligner::Bowtie2]);
    /// ```
    pub fn aligner<T: Into<Option<Vec<Aligner>>>>(mut self, aligner: T) -> Self {
        self.aligners = aligner.into();
        self
    }
    /// Sets the version for the builder.
    ///
    /// # Arguments
    ///
    /// * `version` - A version string in format {major}.{minor}.{path}
    ///
    /// # Example
    ///
    /// ```
    /// let builder = CerebroDownloaderBuilder::new(outdir, indices).version("0.1.0");
    /// ```
    pub fn version<T: Into<Option<String>>>(mut self, version: T) -> Self {
        self.version = version.into();
        self
    }
    /// Sets the classifiers for the builder.
    ///
    /// # Arguments
    ///
    /// * `classifier` - A list of `Classifier` instances.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = CerebroDownloaderBuilder::new(outdir, indices).classifier(vec![Classifier::new()]);
    /// ```
    pub fn classifier<T: Into<Option<Vec<Classifier>>>>(mut self, classifier: T) -> Self {
        self.classifiers = classifier.into();
        self
    }
    /// Sets the timeout duration for the builder.
    ///
    /// # Arguments
    ///
    /// * `timeout` - Timeout duration in seconds.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = CerebroDownloaderBuilder::new(outdir, indices).timeout(60);
    /// ```
    pub fn timeout<T: Into<Option<u64>>>(mut self, timeout: T) -> Self {
        self.timeout = timeout.into();
        self
    }
    /// Sets the base URL for the builder.
    ///
    /// # Arguments
    ///
    /// * `base_url` - The base URL for downloading files.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = CerebroDownloaderBuilder::new(outdir, indices).base_url("https://example.com");
    /// ```
    pub fn base_url<T: Into<Option<String>>>(mut self, base_url: T) -> Self {
        self.base_url = base_url.into();
        self
    }
    /// Sets the username for authentication.
    ///
    /// # Arguments
    ///
    /// * `username` - The username for basic authentication.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = CerebroDownloaderBuilder::new(outdir, indices).username("user");
    /// ```
    pub fn username<T: Into<Option<String>>>(mut self, username: T) -> Self {
        self.username = username.into();
        self
    }
    /// Sets the password for authentication.
    ///
    /// # Arguments
    ///
    /// * `password` - The password for basic authentication.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = CerebroDownloaderBuilder::new(outdir, indices).password("pass");
    /// ```
    pub fn password<T: Into<Option<String>>>(mut self, password: T) -> Self {
        self.password = password.into();
        self
    }
    /// Sets the reference to download index-associated sequences.
    ///
    /// # Arguments
    ///
    /// * `reference` - Boolean value to download sequences
    ///
    /// # Example
    ///
    /// ```
    /// let builder = CerebroDownloaderBuilder::new(outdir, indices).reference(true);
    /// ```
    pub fn reference(mut self, reference: bool) -> Self {
        self.reference = reference;
        self
    }
    /// Builds the `CerebroDownloader` instance.
    ///
    /// # Errors
    ///
    /// Returns a `WorfkflowError` if the downloader cannot be built.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = CerebroDownloaderBuilder::new(outdir, indices);
    /// let downloader = builder.timeout(60).build();
    /// ```
    pub fn build(self) -> Result<CerebroDownloader, WorkflowError> {
        
        if !self.outdir.exists() || !self.outdir.is_dir() {
            create_dir_all(&self.outdir)?;
        }
        
        let username = self.username
            .unwrap_or("u416706-sub2".to_string());
        let password = self.password
            .unwrap_or("uz2q4cyKFytNXvch".to_string());
        let base_url = self.base_url
            .unwrap_or(format!("https://{username}.your-storagebox.de/databases"));
        let aligners = self.aligners
            .unwrap_or(Vec::new());
        let classifiers = self.classifiers
            .unwrap_or(Vec::new());
        let timeout = self.timeout
            .unwrap_or(30);

        Ok(CerebroDownloader {
            outdir: self.outdir.to_owned(),
            base_url,
            username,
            password,
            client: Client::new(),
            timeout,
            indices: self.indices.clone(),
            aligners,
            classifiers,
            version: self.version,
            reference: self.reference
        })
    }
}
