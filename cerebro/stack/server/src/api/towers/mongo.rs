

use mongodb::bson::{doc, Document};

pub fn get_registered_towers_pipeline(id: &Option<String>) -> Vec<Document> {
    if let Some(id) = id {
        vec![
            doc! {
                "$match": {
                    "id": &id
                }
            },
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

