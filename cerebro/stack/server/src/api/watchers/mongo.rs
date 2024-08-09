

use mongodb::bson::{doc, Document};

pub fn get_registered_watchers_pipeline() -> Vec<Document> {
    
    vec![
        doc! {
            "$sort": {
                "date": -1
            }
        },
    ]
    
}
