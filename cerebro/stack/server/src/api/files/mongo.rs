

use mongodb::bson::{doc, Document};

pub fn get_latest_files_paginated_pipeline(page: i64, limit: i64) -> Vec<Document> {
    
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
