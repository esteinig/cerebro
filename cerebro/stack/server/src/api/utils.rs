use actix_web::{web::Data, HttpResponse};
use cerebro_model::api::teams::model::{Team, TeamAdminCollection};
use cerebro_model::api::utils::AdminCollection;
use mongodb::Collection;
use serde::ser::Serialize;
use crate::api::server::AppState;

use serde::{Deserialize, Deserializer};

// Generic error response
#[derive(Debug, Deserialize)]
pub struct ErrorResponse {
    pub status: String,
    pub message: String,
}


// Get administrative database collection from admin database
pub fn get_cerebro_db_collection<T>(data: &Data<AppState>, collection: AdminCollection) -> Collection<T> {
    let db = data.db.database(&data.env.database.names.admin_database_name);
    
    let collection_name = match collection {
        AdminCollection::Teams => &data.env.database.names.admin_database_team_collection,
        AdminCollection::Users => &data.env.database.names.admin_database_user_collection,
        AdminCollection::Logs => &data.env.database.names.admin_database_logs_collection
    };
    
    db.collection(collection_name)
}

// Get administrative database collection from team admin database
pub fn get_teams_db_collection<T>(data: &Data<AppState>, team: Team, collection: TeamAdminCollection) -> Collection<T> {
    let db = data.db.database(&team.admin.database);
    db.collection(&collection.name())
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


pub fn _empty_string_as_none<'de, D>(deserializer: D) -> Result<Option<String>, D::Error>
where
    D: Deserializer<'de>,
{
    let opt = Option::<String>::deserialize(deserializer)?;
    Ok(opt.filter(|s| !s.trim().is_empty()))
}