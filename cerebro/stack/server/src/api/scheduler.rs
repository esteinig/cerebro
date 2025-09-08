// src/scheduler.rs
use std::time::Duration as StdDuration;
use cerebro_model::api::jobs::model::ScheduleJob;
use chrono::{DateTime, Duration, Utc};
use faktory::Job;
use futures::TryStreamExt;
use mongodb::{
    bson::{doc, to_bson},
    options::{FindOneAndUpdateOptions, FindOptions, ReturnDocument},
    Client as MongoClient, Collection,
};
use serde::{Deserialize, Serialize};
use uuid::Uuid;



pub fn job_with_val(kind: &str, val: serde_json::Value) -> Job {
    Job::new(kind, vec![val])  // -> args: [ { ... } ]
}

const LOCK_ID: &str = "global-scheduler";
const DEFAULT_LEASE_SECS: i64 = 30;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SchedulerLock {
    pub id: String,
    pub owner_id: String,
    pub lease_expires_at: DateTime<Utc>,
}

#[derive(Debug, Clone)]
pub struct Scheduler {
    instance_id: String,
    jobs: Collection<ScheduleJob>,
    locks: Collection<SchedulerLock>,
}

impl Scheduler {
    pub fn new(client: &MongoClient, db_name: impl Into<String>, jobs_collection_name: impl Into<String>, locks_collection_name: impl Into<String>) -> Self {
        
        let db = client.database(&db_name.into());
        let instance_id = Uuid::new_v4();

        Self {
            instance_id: instance_id.into(),
            jobs: db.collection(&jobs_collection_name.into()),
            locks: db.collection(&locks_collection_name.into()),
        }
    }

    pub fn spawn(self) {
        tokio::spawn(async move {
            let mut tick = tokio::time::interval(StdDuration::from_secs(5));
            loop {
                tick.tick().await;
                if let Err(e) = self.tick_once().await {
                    log::error!("scheduler pass failed: {e}");
                }
            }
        });
    }

    async fn tick_once(&self) -> anyhow::Result<()> {
        if !self.try_acquire_lease(DEFAULT_LEASE_SECS).await? {
            return Ok(());
        }

        let now = Utc::now();
        let filter = doc! { "enabled": true, "run_at": { "$lte": to_bson(&now)? } };
        let opts = FindOptions::builder().sort(doc!{"run_at": 1}).limit(100).build();
        
        let mut cur = self.jobs.find(filter).with_options(opts).await?;
        
        let mut fk = faktory::Client::connect().await?; // uses FAKTORY_URL

        while let Some(job) = cur.try_next().await? {
            if let Err(e) = self.enqueue_and_reschedule(&mut fk, job.clone()).await {
                log::error!("Error in job reschedule (id: {}) {}", job.id, e.to_string());

                self.jobs.update_one(
                    doc!{"id": job.id},
                    doc!{"$set": {"last_error": e.to_string(), "updated_at": to_bson(&Utc::now())?}}
                ).await?;
            }
        }

        Ok(())
    }

    async fn enqueue_and_reschedule(&self, fk: &mut faktory::Client, job: ScheduleJob) -> anyhow::Result<()> {
        
        // Enqueue object configuration args
        let mut fj = job_with_val(&job.kind, job.args).on_queue(&job.queue);

        fj.retry = Some(job.retry as isize);
        fj.reserve_for = Some(StdDuration::from_secs(job.reserve_for_seconds as u64));

        fk.enqueue(fj).await?;

        // Compute next or disable
        let now = Utc::now();

        let (enabled, next_run_at) = match job.interval_seconds {
            Some(n) if n > 0 => {
                let step = Duration::seconds(n);
                let mut next = job.run_at;
                while next <= now { next = next + step; }
                (true, next)
            }
            _ => (false, job.run_at),
        };

        self.jobs.update_one(
            doc!{"id": job.id},
            doc!{"$set": {
                "last_run_at": to_bson(&now)?,
                "run_at": to_bson(&next_run_at)?,
                "enabled": enabled,
                "last_error": mongodb::bson::Bson::Null,
                "updated_at": to_bson(&now)?,
            }}
        ).await?;

        Ok(())
    }

    // Compact lease: we can take it if we own it or it's expired.
    async fn try_acquire_lease(&self, lease_secs: i64) -> anyhow::Result<bool> {

        let now = Utc::now();
        let new_exp = now + Duration::seconds(lease_secs);

        let filter = doc! {
            "id": LOCK_ID,
            "$or": [
                { "owner_id": &self.instance_id },
                { "lease_expires_at": { "$lte": to_bson(&now)? } }
            ]
        };

        let update = doc! {
            "$set": { "owner_id": &self.instance_id, "lease_expires_at": to_bson(&new_exp)? }
        };
        let opts = FindOneAndUpdateOptions::builder()
            .upsert(true).return_document(ReturnDocument::After).build();

        let doc = self.locks
            .find_one_and_update(filter, update)
            .with_options(opts)
            .await?;

        Ok(doc.map(|d| d.owner_id == self.instance_id).unwrap_or(false))
    }
}
