

use cerebro_model::api::{stage::schema::RegisterStagedSampleSchema, towers::model::ProductionTower};
use chrono::Utc;
use mongodb::bson::{doc, Document, to_bson};
use uuid::Uuid;


pub fn create_staged_samples_pipeline(schema: &RegisterStagedSampleSchema, tower: &ProductionTower, team: &str, database: &str, project: &str) -> Vec<Document> {
    
    let mut mongo_pipeline = vec![];

    if let Some(run_id) = &schema.run_id {
        // Match documents by the specified run_id
        mongo_pipeline.push(doc! {
            "$match": {
                "run_id": run_id
            }
        })
    }

    if let Some(file_ids) = &schema.file_ids {
        // Match documents by the specified run_id
        mongo_pipeline.push(doc! {
            "$match": {
                "id": { "$in": file_ids }
            }
        })
    }

    mongo_pipeline.push(
        // Group the documents by sample_id
        doc! {
            "$group": {
                "_id": "$sample_id",
                "files": { "$push": "$$ROOT" },
                "run_id": { "$first": "$run_id" },
            }
        }
    );

    mongo_pipeline.push(
        // Project the grouped documents into the StagedSample structure
        doc! {
            "$project": {
                "id": Uuid::new_v4().to_string(),
                "date": Utc::now().to_string(),  
                "run_id": "$run_id",
                "sample_id": "$_id",
                "team": team,
                "database": database,                       
                "project": project,            
                "pipeline": format!("{}", schema.pipeline),       
                "tower": to_bson(tower).unwrap(), 
                "files": "$files"
            }
        }
    );

    mongo_pipeline
}



pub fn get_latest_staged_samples_pipeline(run_id: Option<String>, sample_id: Option<String>) -> Vec<Document> {
    
    
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

    
    if let Some(sample_id) = sample_id {
        return vec![
            doc! {
                "$match": {
                    "sample_id": &sample_id
                }  
            },
            doc! {
                "$sort": {
                    "sample_id": 1
                }
            },
        ]
    }

    if let (Some(run_id), Some(sample_id)) = (run_id, sample_id) {
        return vec![
            doc! {
                "$match": {
                    "sample_id": &sample_id,
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

    vec![
        doc! {
            "$sort": {
                "date": -1
            }
        },
    ]
    
}
