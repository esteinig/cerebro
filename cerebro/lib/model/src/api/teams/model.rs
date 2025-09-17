use serde::{Deserialize, Serialize};
use crate::api::{teams::schema::{RegisterDatabaseSchema, RegisterTeamSchema, RegisterProjectSchema}, users::model::UserId};

// A type alias to better track user identifier
pub type TeamId = String;


pub enum TeamAdminCollection {
    Logs,
    Reports,
    Files,
    Watchers,
    Towers,
    TrainingData,
    TrainingSessions
}
impl TeamAdminCollection {
    pub fn name(&self) -> String {
        match self {
            TeamAdminCollection::Logs => String::from("logs"),
            TeamAdminCollection::Reports => String::from("reports"),
            TeamAdminCollection::Files => String::from("files"),
            TeamAdminCollection::Watchers => String::from("watchers"),
            TeamAdminCollection::Towers => String::from("towers"),
            TeamAdminCollection::TrainingData => String::from("training_data"),
            TeamAdminCollection::TrainingSessions => String::from("training_sessions")
        }
    }
}

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
    pub databases: Vec<TeamDatabase>,
    // Administrative database for this team
    // stores collections of files, watchers,
    // pipelines etc.
    pub admin: TeamDatabase
}
impl Team {
    pub fn from(schema: &RegisterTeamSchema) -> Team {
        Team {
            id: uuid::Uuid::new_v4().to_string(),
            name: schema.team_name.to_owned(),
            description: schema.team_description.to_owned(),
            users: vec![schema.team_lead.to_owned()],
            databases: vec![TeamDatabase::from(&schema)],
            admin: TeamDatabase::admin()
        }
    }
    pub fn database_name_exists(&self, db_name: &str) -> bool {
        self.databases.iter().any(|db| db.name == db_name)
    }
    pub fn project_name_exists(&self, db: &str, project_name: &str) -> bool {
        self.databases.iter()
            .find(|database| database.name.as_str() == db || database.id.as_str() == db)
            .is_some_and(|database| database.projects.iter()
            .any(|project| project.name.as_str() == project_name))
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
        let uuid = uuid::Uuid::new_v4().to_string();
        TeamDatabase {
            id: uuid.clone(),
            name: schema.database_name.to_owned(),
            description: schema.database_description.to_owned(),
            database: uuid,
            projects: vec![ProjectCollection::default(&schema)]
        }
    }
    pub fn from_database_schema(schema: &RegisterDatabaseSchema) -> TeamDatabase {
        let uuid = uuid::Uuid::new_v4().to_string();
        let default_collection_uuid = uuid::Uuid::new_v4().to_string();
        TeamDatabase {
            id: uuid.clone(),
            name: schema.database_name.to_owned(),
            description: schema.database_description.to_owned(),
            database: uuid,
            projects: vec![ProjectCollection {
                id: default_collection_uuid.clone(),
                name: String::from("Data"),
                description: String::from("Default data collection"),
                collection: default_collection_uuid
            }]
        }
    }
    pub fn admin() -> TeamDatabase {
        let uuid = uuid::Uuid::new_v4().to_string();
        TeamDatabase {
            id: uuid.clone(),
            name: String::from("Admin"),
            description: String::from("Team administrative database"),
            database: uuid,
            projects: vec![
                ProjectCollection::team_files(),
                ProjectCollection::team_watchers(),
                ProjectCollection::team_pipelines(),
                ProjectCollection::team_logs(),
                ProjectCollection::team_reports(),
                ProjectCollection::team_training_data(),
                ProjectCollection::team_training_sessions(),
            ]
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
        let uuid = uuid::Uuid::new_v4().to_string();
        ProjectCollection {
            id: uuid.clone(),
            name: schema.project_name.to_owned(),
            description: schema.project_description.to_owned(),
            collection: uuid
        }
    }
    pub fn team_files() -> ProjectCollection {
        ProjectCollection {
            id: uuid::Uuid::new_v4().to_string(),
            name: String::from("Files"),
            description: String::from("File registrations for CerebroFS"),
            collection: TeamAdminCollection::Files.name()
        }
    }
    pub fn team_watchers() -> ProjectCollection {
        ProjectCollection {
            id: uuid::Uuid::new_v4().to_string(),
            name: String::from("Watchers"),
            description: String::from("Production watcher registrations for file uploads to CerebroFS"),
            collection: TeamAdminCollection::Watchers.name()
        }
    }
    pub fn team_pipelines() -> ProjectCollection {
        ProjectCollection {
            id: uuid::Uuid::new_v4().to_string(),
            name: String::from("Pipelines"),
            description: String::from("Production pipeline registrations for Cerebro"),
            collection: TeamAdminCollection::Towers.name()
        }
    }
    pub fn team_logs() -> ProjectCollection {
        ProjectCollection {
            id: uuid::Uuid::new_v4().to_string(),
            name: String::from("Logs"),
            description: String::from("Log entries for the team"),
            collection: TeamAdminCollection::Logs.name()
        }
    }
    pub fn team_reports() -> ProjectCollection {
        ProjectCollection {
            id: uuid::Uuid::new_v4().to_string(),
            name: String::from("Reports"),
            description: String::from("Report storage for team independent of database"),
            collection: TeamAdminCollection::Reports.name()
        }
    }
    pub fn team_training_data() -> ProjectCollection {
        ProjectCollection {
            id: uuid::Uuid::new_v4().to_string(),
            name: String::from("TrainingData"),
            description: String::from("Training prefetch data collection"),
            collection: TeamAdminCollection::TrainingData.name()
        }
    }
    pub fn team_training_sessions() -> ProjectCollection {
        ProjectCollection {
            id: uuid::Uuid::new_v4().to_string(),
            name: String::from("TrainingRecords"),
            description: String::from("Training session records collection"),
            collection: TeamAdminCollection::TrainingSessions.name()
        }
    }
}