

use mongodb::bson::{doc, Document};

pub fn get_latest_logs_limit_pipeline(limit: i64, critical: bool) -> Vec<Document> {

    match critical {
        true => {
            vec![
                doc! {
                    "$match": {
                        "critical": true
                    }
                },
                doc! {
                    "$sort": {
                        "date": -1
                    }
                },
                doc! {
                    "$limit": limit
                }
            ]
        },
        false => {
            vec![
                doc! {
                    "$sort": {
                        "date": -1
                    }
                },
                doc! {
                    "$limit": limit
                }
            ]
        }
    }
}

pub fn get_latest_logs_all_pipeline(critical: bool) -> Vec<Document> {

    match critical {
        true => {
            vec![
                doc! {
                    "$match": {
                        "critical": true
                    }
                },
                doc! {
                    "$sort": {
                        "date": -1
                    }
                }
            ]
        }, 
        false => {
            vec![
                doc! {
                    "$sort": {
                        "date": -1
                    }
                }
            ]
        }
    }
    
    
}