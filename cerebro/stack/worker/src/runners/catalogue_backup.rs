//! `catalogue_backup` (S4-2): back up the MongoDB control plane to the backup
//! object store, with a verifiable manifest and retention.
//!
//! Pipeline:
//! 1. **chain-verify-before** — ask the API for the audit trail and record
//!    whether the tamper-evident chain verified intact (so a restore can detect a
//!    backup taken over a broken chain). Best-effort: a failed check records
//!    `None`, it does not abort the backup.
//! 2. **dump** — shell `mongodump --archive --gzip` against a read-only backup URI.
//! 3. **checksum** — BLAKE3 the archive (matches `cerebro-fs` hashing).
//! 4. **upload** — write the archive + a `manifest.json` to the object store.
//! 5. **retain** — prune old backups (keep-last-N with an age floor).
//!
//! Heavy I/O lives here in the worker, never in an API request handler.

use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use async_trait::async_trait;
use chrono::{Duration, Utc};
use faktory::{Job, JobRunner};

use crate::backup::{
    parse_backup_refs, select_for_deletion, BackupManifest, FilesystemObjectStore, ObjectStore,
    BACKUP_ID_FORMAT, MANIFEST_FORMAT_VERSION,
};
use crate::config::BackupSettings;
use crate::context::WorkerContext;
use crate::telemetry::{JobOutcome, Metrics};

pub struct CatalogueBackup {
    ctx: Arc<WorkerContext>,
    metrics: Metrics,
    settings: Option<BackupSettings>,
}

impl CatalogueBackup {
    pub fn new(ctx: Arc<WorkerContext>, metrics: Metrics, settings: Option<BackupSettings>) -> Self {
        Self { ctx, metrics, settings }
    }

    async fn run_inner(&self, _job: Job) -> anyhow::Result<()> {
        let settings = self.settings.as_ref().ok_or_else(|| {
            anyhow::anyhow!(
                "catalogue backup is not configured — set CEREBRO_BACKUP_MONGO_URI and \
                 CEREBRO_BACKUP_STORE_PATH on the worker"
            )
        })?;

        // 1. chain-verify-before (best-effort): record the audit-chain state.
        let (audit_event_count, audit_chain_verified) = match self.ctx.api() {
            Ok(api) => match tokio::task::spawn_blocking(move || api.get_audit_trail(None, None)).await {
                Ok(Ok((events, verified))) => (events.len(), Some(verified)),
                Ok(Err(e)) => {
                    tracing::warn!("audit chain check failed; recording as unknown: {e}");
                    (0, None)
                }
                Err(e) => {
                    tracing::warn!("audit chain check task failed; recording as unknown: {e}");
                    (0, None)
                }
            },
            Err(e) => {
                tracing::warn!("API client unavailable for audit chain check: {e}");
                (0, None)
            }
        };
        if audit_chain_verified == Some(false) {
            tracing::warn!(
                "audit chain did NOT verify intact at backup time — recorded in the manifest"
            );
        }

        // 2. dump to a temporary gzip archive.
        let backup_id = Utc::now().format(BACKUP_ID_FORMAT).to_string();
        let tmp = std::env::temp_dir().join(format!("cerebro-catalogue-{backup_id}.archive.gz"));
        let args = mongodump_args(&settings.mongo_uri, &settings.mongo_db, &tmp);
        let status = tokio::process::Command::new("mongodump")
            .args(&args)
            .status()
            .await
            .map_err(|e| {
                anyhow::anyhow!(
                    "failed to run mongodump (is mongodb-database-tools installed in the worker \
                     image?): {e}"
                )
            })?;
        if !status.success() {
            anyhow::bail!("mongodump exited unsuccessfully: {status}");
        }

        // 3. checksum + size.
        let archive_blake3 = cerebro_fs::hash::fast_file_hash(&tmp)
            .map_err(|e| anyhow::anyhow!("failed to hash backup archive: {e}"))?;
        let archive_bytes = std::fs::metadata(&tmp)?.len();

        // 4. upload archive + manifest.
        let store = FilesystemObjectStore::new(&settings.store_root);
        let archive_key = format!("{}/{}/catalogue.archive.gz", settings.prefix, backup_id);
        let manifest_key = format!("{}/{}/manifest.json", settings.prefix, backup_id);

        let bytes = std::fs::read(&tmp)?;
        store.put(&archive_key, &bytes)?;

        let manifest = BackupManifest {
            format_version: MANIFEST_FORMAT_VERSION,
            backup_id: backup_id.clone(),
            created_at: Utc::now(),
            mongo_db: if settings.mongo_db.is_empty() {
                "(all databases)".to_string()
            } else {
                settings.mongo_db.clone()
            },
            archive_key: archive_key.clone(),
            archive_bytes,
            archive_blake3,
            audit_chain_verified,
            audit_event_count,
            cerebro_version: env!("CARGO_PKG_VERSION").to_string(),
        };
        store.put(&manifest_key, &serde_json::to_vec_pretty(&manifest)?)?;
        let _ = std::fs::remove_file(&tmp); // best-effort temp cleanup

        // 5. retention: keep-last-N with an age floor.
        let keys = store.list(&format!("{}/", settings.prefix))?;
        let refs = parse_backup_refs(&settings.prefix, &keys);
        let stale = select_for_deletion(
            refs,
            settings.keep_last,
            Duration::days(settings.min_age_days),
            Utc::now(),
        );
        for b in &stale {
            store.delete(&format!("{}/{}/catalogue.archive.gz", settings.prefix, b.backup_id))?;
            store.delete(&format!("{}/{}/manifest.json", settings.prefix, b.backup_id))?;
        }

        tracing::info!(
            backup_id = %backup_id,
            bytes = archive_bytes,
            chain_verified = ?audit_chain_verified,
            pruned = stale.len(),
            "catalogue backup complete"
        );
        Ok(())
    }
}

#[async_trait]
impl JobRunner for CatalogueBackup {
    type Error = io::Error;

    async fn run(&self, job: Job) -> Result<(), Self::Error> {
        let kind = job.kind().to_string();
        self.metrics.record_job(&kind, JobOutcome::Started);
        match self.run_inner(job).await {
            Ok(_) => {
                self.metrics.record_job(&kind, JobOutcome::Succeeded);
                Ok(())
            }
            Err(e) => {
                self.metrics.record_job(&kind, JobOutcome::Failed);
                tracing::error!("catalogue_backup failed: {e:#}");
                Err(io::Error::new(io::ErrorKind::Other, e.to_string()))
            }
        }
    }
}

/// Build the `mongodump` argument vector (pure, for testing).
///
/// `--archive=<file>` + `--gzip` produces a single compressed archive that
/// `mongorestore --gzip --archive=<file>` reads back. When `db` is empty the whole
/// instance is dumped (the control plane spans every team's catalogue database
/// plus users/audit/schedules); a non-empty `db` scopes the dump.
pub(crate) fn mongodump_args(uri: &str, db: &str, archive: &PathBuf) -> Vec<String> {
    let mut args = vec![format!("--uri={uri}")];
    if !db.is_empty() {
        args.push(format!("--db={db}"));
    }
    args.push(format!("--archive={}", archive.display()));
    args.push("--gzip".to_string());
    args
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mongodump_args_are_well_formed() {
        let scoped = mongodump_args(
            "mongodb://ro:pw@cerebro-database:27017",
            "cerebro",
            &PathBuf::from("/tmp/x.archive.gz"),
        );
        assert!(scoped.contains(&"--uri=mongodb://ro:pw@cerebro-database:27017".to_string()));
        assert!(scoped.contains(&"--db=cerebro".to_string()));
        assert!(scoped.contains(&"--archive=/tmp/x.archive.gz".to_string()));
        assert!(scoped.contains(&"--gzip".to_string()));

        // Empty db => whole-instance dump: no --db flag.
        let all = mongodump_args("mongodb://ro:pw@h:27017", "", &PathBuf::from("/tmp/x.gz"));
        assert!(all.iter().all(|a| !a.starts_with("--db")));
        assert!(all.contains(&"--gzip".to_string()));
    }
}
