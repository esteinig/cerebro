

use mongodb::bson::{doc, Document};

pub fn get_latest_files_paginated_pipeline(run_id: Option<String>, page: i64, limit: i64) -> Vec<Document> {
    
    match run_id {
        Some(run_id) => {
            vec![
                doc! {
                    "$match": {
                        "run_id": &run_id
                    }  
                },
            ]
        },
        None => {
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
        }
    }
    
}
