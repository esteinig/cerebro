

use mongodb::bson::{doc, Document};

pub fn get_latest_files_paginated_pipeline(run_id: Option<String>, watcher_id: Option<String>, page: i64, limit: i64) -> Vec<Document> {
    
    
    if let Some(run_id) = run_id {
        return vec![
            doc! {
                "$match": {
                    "run_id": &run_id
                }  
            },
            doc! {
                "$sort": {
                    "run_id": 1
                }
            },
        ]
    }

    if let Some(watcher_id) = watcher_id {
        return vec![
            doc! {
                "$match": {
                    "watcher.id": &watcher_id
                }  
            },
            doc! {
                "$sort": {
                    "sample_id": 1
                }
            },
        ]
    }

    if let (Some(run_id), Some(watcher_id)) = (run_id, watcher_id) {
        return vec![
            doc! {
                "$match": {
                    "watcher.id": &watcher_id,
                    "run_id": &run_id
                }  
            },
            doc! {
                "$sort": {
                    "sample_id": 1
                }
            },
        ]
    }
    if limit > 0 {
        vec![
            doc! {
                "$sort": {
                    "date": -1
                }
            },
            doc! {
                // page
                "$skip": page*limit
            },
            doc! {
                "$limit": limit
            }
        ]
    } else {
        vec![
            doc! {
                "$sort": {
                    "date": -1
                }
            },
        ]
    }
    
    
}
