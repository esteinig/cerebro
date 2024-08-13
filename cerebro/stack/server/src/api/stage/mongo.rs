

use cerebro_model::api::stage::schema::RegisterStagedSampleSchema;
use chrono::Utc;
use mongodb::bson::{doc, Document};
use uuid::Uuid;


pub fn get_staged_samples_pipeline(schema: &RegisterStagedSampleSchema) -> Vec<Document> {
    
    vec![
        // Match documents by the specified run_id
        doc! {
            "$match": {
                "run_id": &schema.run_id
            }
        },
        // Group the documents by sample_id
        doc! {
            "$group": {
                "_id": "$sample_id",
                "files": { "$push": "$$ROOT" },
                "run_id": { "$first": "$run_id" },
            }
        },
        // Project the grouped documents into the StagedSample structure
        doc! {
            "$project": {
                "id": Uuid::new_v4().to_string(),
                "date": Utc::now().to_string(),  
                "run_id": "$run_id",
                "sample_id": "$_id",
                "database": &schema.database,                       
                "project": &schema.project,                         
                "pipeline": format!("{}", schema.pipeline),        
                "files": "$files"
            }
        }
    ]
}



pub fn get_latest_staged_samples_pipeline(run_id: Option<String>) -> Vec<Document> {
    
    
    if let Some(run_id) = run_id {
        return vec![
            doc! {
                "$match": {
                    "run_id": &run_id
                }  
            },
            doc! {
                "$sort": {
                    "date": -1
                }
            },
        ]
    }

    vec![
        doc! {
            "$sort": {
                "date": -1
            }
        },
    ]
    
}
