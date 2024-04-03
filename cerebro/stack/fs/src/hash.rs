use blake3::Hasher;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use crate::error::FileSystemError;

/// Computes a BLAKE3 hash of the file at the given path.
///
/// # Arguments
///
/// * `path` - The path of the file to hash.
///
/// # Returns
///
/// A `Result` which is `Ok` with a hexadecimal string representing the hash of the file,
/// or an `Err` with an `io::Error` if an error occurs during file reading or hashing.
///
/// # Examples
///
/// ```
/// use std::path::PathBuf;
///
/// let file_path = PathBuf::from("path/to/your/file.txt");
/// match fast_file_hash(&file_path) {
///     Ok(hash) => println!("File hash: {}", hash),
///     Err(e) => println!("Error computing file hash: {}", e),
/// }
/// ```
pub fn fast_file_hash(path: &PathBuf) -> Result<String, FileSystemError> {
    let file = File::open(path)?;
    let mut hasher = Hasher::new();
    let mut reader = BufReader::new(file);
    let mut buffer = [0; 8192]; // Read in chunks of 8KB.

    loop {
        let count = reader.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hasher.update(&buffer[..count]);
    }

    Ok(hasher.finalize().to_hex().to_string())
}