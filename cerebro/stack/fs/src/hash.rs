use crate::error::FileSystemError;
use blake3::Hasher;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

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
/// use cerebro_fs::hash::fast_file_hash;
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

/// Computes a BLAKE3 hash by streaming an arbitrary reader, without persisting
/// it to disk.
///
/// This is the basis for object integrity verification: an HTTP
/// response body (filer or volume GET) implements [`Read`], so the bytes are
/// hashed in flight, in constant memory, and discarded — no temp file, no
/// second pass. Returns the lowercase hex digest, matching [`fast_file_hash`].
pub fn hash_reader<R: Read>(mut reader: R) -> Result<String, FileSystemError> {
    let mut hasher = Hasher::new();
    let mut buffer = [0; 65536]; // 64KB chunks: fewer syscalls on large objects.

    loop {
        let count = reader.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hasher.update(&buffer[..count]);
    }

    Ok(hasher.finalize().to_hex().to_string())
}

/// Computes a BLAKE3 hash of an in-memory byte slice.
///
/// Used to seal the run manifest over its canonical JSON body.
pub fn hash_bytes(data: &[u8]) -> String {
    let mut hasher = Hasher::new();
    hasher.update(data);
    hasher.finalize().to_hex().to_string()
}

/// Walk `dir` recursively, computing the BLAKE3 of each file, and write a
/// b3sum-style checksum file (`<hash>  <relative-path>` per line, sorted) to
/// `output`. Produce-time integrity snapshot for the produce -> capture -> verify
/// chain (intermission-3). Returns the number of files hashed.
pub fn hash_directory(
    dir: &std::path::Path,
    output: &std::path::Path,
) -> Result<usize, FileSystemError> {
    let mut entries: Vec<PathBuf> = Vec::new();
    collect_files_recursive(dir, &mut entries)?;
    entries.sort();

    let mut out = String::new();
    let mut count = 0usize;
    for path in &entries {
        if path.as_path() == output {
            continue; // never hash the checksum file itself
        }
        let hash = fast_file_hash(path)?;
        let rel = path.strip_prefix(dir).unwrap_or(path).to_string_lossy();
        out.push_str(&format!("{}  {}\n", hash, rel));
        count += 1;
    }
    std::fs::write(output, out)?;
    Ok(count)
}

fn collect_files_recursive(
    dir: &std::path::Path,
    acc: &mut Vec<PathBuf>,
) -> Result<(), FileSystemError> {
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_dir() {
            collect_files_recursive(&path, acc)?;
        } else if path.is_file() {
            acc.push(path);
        }
    }
    Ok(())
}
