use std::{path::PathBuf, process::Output};
use std::process::Command;
use reqwest::blocking::get;
use std::fs;
use std::io::Cursor;
use tar::Archive;
use anyhow::Result;

use crate::error::{WeedUploadError, WeedDownloadError};
use cerebro_model::api::files::response::WeedUploadResponse;

/// Executes the `weed` upload command for single file upload with specified options and returns a single response.
///
/// # Arguments
///
/// - `input_file`: A `PathBuf` pointing to the file to be uploaded.
/// - `data_center`: Optional data center name.
/// - `include`: File patterns to include in the upload, e.g., "*.pdf".
/// - `master`: SeaweedFS master location, defaults to "localhost:9333".
/// - `max_mb`: Split files larger than this limit, defaults to 4 MB.
/// - `options`: Additional command line options in "optionName=optionValue" format.
/// - `replication`: Replication type.
/// - `ttl`: Time to live for the uploaded file.
/// - `use_public_url`: Whether to upload to the public URL from the volume server.
///
/// # Returns
///
/// A `Result` with either a single `WeedUploadResponse` on success or an `UploadError` on failure.
///
/// # Examples
///
/// ```no_run
/// use std::path::PathBuf;
///
/// let input_file = PathBuf::from("path/to/your/file.pdf");
/// let response = weed_upload(
///     &input_file,
///     None,
///     None,
///     Some("localhost:9333"),
///     Some(4),
///     None,
///     None,
///     None,
///     false,
/// );
///
/// match response {
///     Ok(responses) => println!("{:?}", responses),
///     Err(e) => println!("An error occurred: {}", e),
/// }
/// ```
pub fn weed_upload(
    input_file: &PathBuf,
    data_center: Option<String>,
    include: Option<String>,
    master: Option<String>,
    port: Option<String>,
    max_mb: Option<i32>,
    options: Option<String>,
    replication: Option<String>,
    ttl: Option<String>,
    use_public_url: bool,
) -> Result<WeedUploadResponse, WeedUploadError> {

    let command = get_command(
        input_file,
        data_center,
        include,
        master,
        port,
        max_mb,
        options,
        replication,
        ttl,
        use_public_url
    );

    log::info!("Executing upload command: weed {command}");

    let output = run_command(&command, "weed")?;
    let output_str = String::from_utf8_lossy(&output.stdout);
    let response: Vec<WeedUploadResponse> = serde_json::from_str(&output_str)?;

    if response.len() > 1 {
        return  Err(WeedUploadError::SingleUploadResponseError);
    }

    response.first().ok_or(WeedUploadError::SingleUploadResponseError).cloned()
}




pub fn run_command(command: &str, program: &str) -> Result<Output, WeedUploadError> {
    let args: Vec<&str> = command.split_whitespace().into_iter().collect();

    let output = Command::new(program)
        .args(args)
        .output()
        .map_err(|_| WeedUploadError::ProgramExecutionFailed(program.to_string()))?;

    // Ensure command ran successfully
    if !output.status.success() {
        return Err(WeedUploadError::CommandExecutionFailed(String::from_utf8_lossy(&output.stderr).to_string()));
    }
    Ok(output)
}


pub fn get_command(
    input_file: &PathBuf,
    data_center: Option<String>,
    include: Option<String>,
    master: Option<String>,
    port: Option<String>,
    max_mb: Option<i32>,
    options: Option<String>,
    replication: Option<String>,
    ttl: Option<String>,
    use_public_url: bool,
) -> String {

    let mut command_parts = vec!["upload".to_owned()];

    if let Some(incl) = include {
        command_parts.push(format!("-include={}", incl));
    }

    if let Some(dc) = data_center {
        command_parts.push(format!("-dataCenter={}", dc));
    }

    // Strip "http://" and "https://" prefixes from the master parameter
    let master_arg = master.unwrap_or_else(|| "localhost".to_owned())
                            .replace("http://", "")
                            .replace("https://", "");

    let master_port = port.unwrap_or_else(|| "9333".to_owned());

    command_parts.push(format!("-master={}:{}", master_arg, master_port));

    if let Some(mb) = max_mb {
        command_parts.push(format!("-maxMB={}", mb));
    }

    if let Some(opts) = options {
        command_parts.push(format!("-options={}", opts));
    }

    if let Some(rep) = replication {
        command_parts.push(format!("-replication={}", rep));
    }

    if let Some(time_to_live) = ttl {
        command_parts.push(format!("-ttl={}", time_to_live));
    }

    if use_public_url {
        command_parts.push("-usePublicUrl".to_owned());
    }

    command_parts.push(input_file.display().to_string());
    
    // Combine all parts into a single command string for display
    let command_string = command_parts.join(" ");

    command_string
}

/// Downloads and installs the `weed` executable from SeaweedFS's GitHub releases into a specified path.
///
/// # Arguments
///
/// * `version` - The version of the `weed` binary to download; use "latest" for the latest release.
/// * `bin_path` - The directory where the `weed` executable will be moved after installation.
///
/// # Errors
///
/// This function returns a `Result<(), WeedDownloadError>`, which is `Ok` if the download and installation
/// process completes successfully, or an `Err` with `WeedDownloadError` detailing what went wrong.
///
/// # Examples
///
/// ```
/// let bin_path = std::path::PathBuf::from("/usr/local/bin/weed");
/// download_and_install_weed("latest", &bin_path).expect("Failed to download and install `weed`");
/// ```
pub fn download_and_install_weed(version: &str, bin_path: &PathBuf) -> Result<(), WeedDownloadError> {
    let url = format!(
        "https://github.com/seaweedfs/seaweedfs/releases/download/{version}/linux_amd64_large_disk.tar.gz"
    );

    log::info!("Downloading SeaweedFS binary from: {url}");
    let response = get(&url).map_err(|_| WeedDownloadError::DownloadError)?;

    if !response.status().is_success() {
        return Err(WeedDownloadError::DownloadError);
    }

    let bytes = response.bytes().map_err(|_| WeedDownloadError::DownloadError)?;

    log::info!("Extracting `weed` archive...");
    let cursor = Cursor::new(bytes);
    let mut archive = Archive::new(cursor);
    archive.unpack(".").map_err(|_| WeedDownloadError::ExtractError)?;

    log::info!("Setting permissions...");
    Command::new("chmod")
        .args(&["+x", "./weed"])
        .status()
        .map_err(|_| WeedDownloadError::PermissionError)?;

    log::info!("Moving binary to: {}", bin_path.display());
    fs::rename(PathBuf::from("./weed"), bin_path).map_err(|_| WeedDownloadError::MoveError)?;

    log::info!("SeaweedFS executable `weed` has been installed successfully at {}", bin_path.display());

    Ok(())
}
