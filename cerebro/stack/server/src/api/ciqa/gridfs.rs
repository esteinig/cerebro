//! CIQA GridFS helpers. A QC dataset sample stores a `PrefetchData` exactly like training, so we
//! reuse the training GridFS helpers rather than duplicate them.
pub use crate::api::training::gridfs::{
    delete_from_gridfs, download_prefetch_from_gridfs, find_unique_gridfs_id_by_filename,
    upload_prefetch_to_gridfs,
};
