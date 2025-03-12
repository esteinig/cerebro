use std::collections::HashMap;

use cerebro_pipeline::taxa::taxon::Taxon;
use futures::io::{AsyncWriteExt, AsyncReadExt};

use mongodb::gridfs::GridFsBucket;


/// Uploads the taxa HashMap to GridFS as a JSON-encoded file.
///
/// # Arguments
/// * `bucket` - A GridFsBucket instance to perform the upload.
/// * `taxa` - A reference to the taxa HashMap to be stored.
/// * `filename` - A string slice for the filename to be stored in GridFS.
///
/// # Returns
/// * On success, returns the GridFS file’s ObjectId.
///
/// **Note:** In v3 you open an upload stream, write the data to it, and then close the stream
/// to finalize the upload.
pub async fn upload_taxa_to_gridfs(
    bucket: GridFsBucket,
    taxa: &Vec<Taxon>,
    filename: &str,
) -> Result<String, Box<dyn std::error::Error>> {

    // Serialize taxa to JSON bytes.
    let taxa_bytes = serde_json::to_vec(taxa)?;

    // Open an upload stream for the given filename.
    let mut upload_stream = bucket.open_upload_stream(filename).await?;
    
    // Write all the bytes to the upload stream.
    upload_stream.write_all(&taxa_bytes).await?;
    
    // Close the stream to finalize the upload
    upload_stream.close().await?;
    
    Ok(filename.to_string())
}

/// Downloads the taxa HashMap from GridFS using the provided ObjectId.
///
/// # Arguments
/// * `bucket` - A GridFsBucket instance to perform the download.
/// * `taxa_id` - The ObjectId referencing the stored taxa data in GridFS.
///
/// # Returns
/// * On success, returns the taxa HashMap.
///
/// **Note:** open_download_stream returns a GridFsDownloadStream which implements Stream over
/// chunks of data.
pub async fn download_taxa_from_gridfs(
    bucket: GridFsBucket,
    filename: &str,
) -> Result<Vec<Taxon>, Box<dyn std::error::Error>> {
    // Open the download stream using the taxa_id.
    let mut download_stream = bucket.open_download_stream_by_name(filename).await?;
    
    let mut bytes: Vec<u8> = Vec::new();
    let _ = download_stream.read_to_end(&mut bytes).await?;
    
    // Deserialize the bytes back into a HashMap.
    let taxa: Vec<Taxon> = serde_json::from_slice(&bytes)?;
    
    Ok(taxa)
}