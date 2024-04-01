use actix_web::{web::Data, HttpResponse};
use mongodb::Collection;
use serde::Deserialize;
use serde::ser::Serialize;
use crate::api::server::AppState;

// Generic error response
#[derive(Debug, Deserialize)]
pub struct ErrorResponse {
    pub status: String,
    pub message: String,
}

// Get administrative database collection from admin database
pub fn get_cerebro_db_collection<T>(data: &Data<AppState>, collection: &str) -> Collection<T> {
    let db = data.db.database(&data.env.database.names.admin_database_name);
    
    let collection_name = match collection {
        "team" => &data.env.database.names.admin_database_team_collection,
        "user" => &data.env.database.names.admin_database_user_collection,
        "logs"  => &data.env.database.names.admin_database_logs_collection,
        _ => unimplemented!("Collection is not supported. Something went wrong in your code, dude!")
    };
    
    db.collection(collection_name)
}

type MongoDatabaseName = String;

// Get administrative database collection from team database
pub fn get_teams_db_collection<T>(data: &Data<AppState>, db_name: &MongoDatabaseName, collection: &str) -> Collection<T> {
    let db = data.db.database(db_name);
    
    let collection_name = match collection {
        "logs" => &data.env.database.names.team_database_logs_collection,
        "reports" => &data.env.database.names.team_database_reports_collection,
        _ => unimplemented!("Collection is not supported. Something went wrong in your code, dude!")
    };
    db.collection(collection_name)
}

pub fn as_csv_string<T>(data: Vec<T>) -> Result<String, HttpResponse> 
    where T: Serialize
{

    let mut csv_writer = csv::WriterBuilder::new()
            .delimiter(b',')
            .from_writer(Vec::new());

        for d in data {
            csv_writer.serialize(d).unwrap()
        }
        csv_writer.flush().unwrap();
        
        let csv_string = csv_writer.into_inner().unwrap().to_vec();

        let csv_result = String::from_utf8(csv_string).unwrap();

        Ok(csv_result)
}