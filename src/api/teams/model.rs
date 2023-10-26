use serde::{Deserialize, Serialize};
use crate::api::{teams::schema::{RegisterDatabaseSchema, RegisterTeamSchema, RegisterProjectSchema}, users::model::UserId};

// A type alias to better track user identifier
pub type TeamId = String;

// A team that users can belong to - teams 
// have access to specific MongoDB databases
// A user can belong to multiple teams and 
// therefore access multiple MongoDB databases
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct Team {
    pub id: TeamId,
    pub name: String,
    pub description: String,
    pub users: Vec<UserId>,
    pub databases: Vec<TeamDatabase>
}
impl Team {
    pub fn from(schema: &RegisterTeamSchema) -> Team {
        Team {
            id: uuid::Uuid::new_v4().to_string(),
            name: schema.team_name.to_owned(),
            description: schema.team_description.to_owned(),
            users: vec![schema.team_lead.to_owned()],
            databases: vec![TeamDatabase::from(&schema)],
        }
    }
}

pub type DatabaseId = String;

// TODO: ensure that `TeamDatabase.database` and 
// TODO: `ProjectCollection.collection` are unique

// Specifications of the database that users
// can access through their teams - allows for
// name validation according to MongoDB
// restrictions on database naming schemes
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct TeamDatabase {
    pub id: DatabaseId,
    // User provided alias of the database for
    // better visibility and display
    pub name: String,
    // Name with which the database is accessed
    // using the MongoDB connection handlers
    pub database: String,
    pub description: String,
    // Collections associated with this database
    // are stored as project structs -
    pub projects: Vec<ProjectCollection>,
}
impl TeamDatabase {
    pub fn from(schema: &RegisterTeamSchema) -> TeamDatabase {
        TeamDatabase {
            id: uuid::Uuid::new_v4().to_string(),
            name: schema.database_name.to_owned(),
            description: schema.database_description.to_owned(),
            database: schema.database_mongo_name.to_owned(),
            projects: vec![ProjectCollection::default(&schema)]
        }
    }
    pub fn from_database_schema(schema: &RegisterDatabaseSchema) -> TeamDatabase {
        TeamDatabase {
            id: uuid::Uuid::new_v4().to_string(),
            name: schema.database_name.to_owned(),
            description: schema.database_description.to_owned(),
            database: schema.database_mongo_name.to_owned(),
            projects: vec![ProjectCollection {
                id: uuid::Uuid::new_v4().to_string(),
                name: String::from("Data"),
                description: String::from("Default data collection"),
                collection: String::from("data")
            }]
        }
    }
}


pub type ProjectId = String;

// Specifications of the project collection that
// users can create and validate collections 
// strings faccording to MongoDB restrictions
// on collection naming schemes
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct ProjectCollection {
    pub id: ProjectId,
    pub name: String,
    pub collection: String,
    pub description: String
}
impl ProjectCollection {
    pub fn default(schema: &RegisterTeamSchema) -> ProjectCollection {
        ProjectCollection {
            id: uuid::Uuid::new_v4().to_string(),
            name: schema.project_name.to_owned(),
            description: schema.project_description.to_owned(),
            collection: schema.project_mongo_name.to_owned()
        }
    }
    pub fn from_project_schema(schema: &RegisterProjectSchema) -> ProjectCollection {
        ProjectCollection {
            id: uuid::Uuid::new_v4().to_string(),
            name: schema.project_name.to_owned(),
            description: schema.project_description.to_owned(),
            collection: schema.project_mongo_name.to_owned()
        }
    }
}