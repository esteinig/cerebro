

use mongodb::bson::{doc, Document};

pub fn get_registered_pipelines_pipeline() -> Vec<Document> {
    
    vec![
        doc! {
            "$sort": {
                "date": -1
            }
        },
    ]
    
}
