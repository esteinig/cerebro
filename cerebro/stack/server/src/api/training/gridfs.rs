use cerebro_model::api::cerebro::schema::PrefetchData;
use futures::io::{AsyncWriteExt, AsyncReadExt};
use futures::TryStreamExt;
use mongodb::{bson::{doc, Bson}, gridfs::GridFsBucket};


/// Uploads the PrefetchData to GridFS as a JSON-encoded file.
///
/// # Arguments
/// * `bucket` - A GridFsBucket instance to perform the upload.
/// * `data` - A reference to the PrefetchData to be stored.
/// * `filename` - A string slice for the Cerebro UUID under which the taxa are to be stored in GridFS.
///
/// # Returns
/// * On success, returns the filename
///
/// **Note:** In v3 you open an upload stream, write the data to it, and then close the stream
/// to finalize the upload.
pub async fn upload_prefetch_to_gridfs(
    bucket: &GridFsBucket,
    data: &PrefetchData,
    filename: &str,
) -> Result<Bson, Box<dyn std::error::Error>> {

    // Serialize taxa to JSON bytes.
    let taxa_bytes = serde_json::to_vec(data)?;

    // Open an upload stream for the given filename.
    let mut upload_stream = bucket.open_upload_stream(filename).await?;
    
    // Write all the bytes to the upload stream.
    upload_stream.write_all(&taxa_bytes).await?;
    
    // Close the stream to finalize the upload
    upload_stream.close().await?;
    
    Ok(upload_stream.id().clone())
}

/// Downloads the PrefetchData from GridFS using the provided ObjectId.
///
/// # Arguments
/// * `bucket` - A GridFsBucket instance to perform the download.
/// * `filename` - The Cerebro UUID referencing the stored taxa data in GridFS.
///
/// # Returns
/// * On success, returns the PrefetchData.
///
/// **Note:** open_download_stream returns a GridFsDownloadStream which implements Stream over
/// chunks of data.
pub async fn download_prefetch_from_gridfs(
    bucket: &GridFsBucket,
    filename: &str,
) -> Result<PrefetchData, Box<dyn std::error::Error>> {
    // Open the download stream using the taxa_id.
    let mut download_stream = bucket.open_download_stream_by_name(filename).await?;
    
    let mut bytes: Vec<u8> = Vec::new();
    let _ = download_stream.read_to_end(&mut bytes).await?;
    
    // Deserialize the bytes back into a HashMap.
    let data: PrefetchData = serde_json::from_slice(&bytes)?;
    
    Ok(data)
}


/// Strict: require that `filename` maps to exactly one GridFS file. Returns that file's id.
/// Errors if no match or multiple matches exist.
pub async fn find_unique_gridfs_id_by_filename(
    bucket: &GridFsBucket,
    filename: &str,
) -> Result<Bson, Box<dyn std::error::Error>> {
    let mut cursor = bucket.find(doc! { "filename": filename }).await?;
    let mut found: Option<Bson> = None;

    while let Some(file) = cursor.try_next().await? {
        if found.is_some() {
            return Err(format!("Multiple GridFS files found for filename '{}'", filename).into());
        }
        found = Some(file.id.clone());
    }

    found.ok_or_else(|| format!("No GridFS file found for filename '{}'", filename).into())
}


pub async fn delete_from_gridfs(
    bucket: &GridFsBucket,
    file_id: Bson,
) -> mongodb::error::Result<()> {
    bucket.delete(file_id).await
}