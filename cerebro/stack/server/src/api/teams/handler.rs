


use cerebro_model::api::users::model::{UserId, User, Role};
use cerebro_model::api::teams::model::{TeamId, Team, TeamDatabase, ProjectCollection};
use cerebro_model::api::teams::schema::{RegisterTeamSchema, RegisterDatabaseSchema, RegisterProjectSchema, UpdateTeamSchema};
use cerebro_model::api::utils::AdminCollection;

use crate::api::auth::jwt::{self, TeamAccessQuery, TeamDatabaseAccessQuery};
use crate::api::server::AppState;
use crate::api::utils::get_cerebro_db_collection;

use actix_web::{web, get, post, HttpResponse, Responder, delete, patch};
use mongodb::{bson::doc, Collection};
use serde::Deserialize;
use futures::TryStreamExt;

// Register a new team - this is an admin operation, as access
// is currently restricted/controlled [ADMIN]. Note that Teams
// are nested within User, so we implement the CRUD operations
// in this file for user endpoints (but keep separate paths)
#[post("/teams")]
async fn register_team_handler(
    body: web::Json<RegisterTeamSchema>,
    data: web::Data<AppState>,
    auth_guard: jwt::JwtAdminMiddleware 
) -> impl Responder {

    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);  
    let team_collection: Collection<Team> = get_cerebro_db_collection(&data, AdminCollection::Teams);

    // Check if team exists already using its team name OR
    // its database field (must be unique) arguably the 
    // database field is more important as duplicates
    // of MongoDB databases would not be possible
    let team_exists = match team_collection
        .find_one(doc! { 
            "$or": [
                {"name": &body.team_name },
                {"databases.database": &body.database_mongo_name }
            ]
        }, None)
        .await
    {
        Ok(None) => false,
        Ok(Some(_)) => true,
        Err(err) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error","message": format!("Failed to register team: {}", err.to_string())}),  // not informative 
        )
    };
    if team_exists {
        return HttpResponse::Conflict().json(
            serde_json::json!({"status": "error", "message": "Team with this name and database already exists"}),
        );
    }

    // Check if the user exists (team lead)
    let user_exists = match user_collection
        .find_one(doc! { "id": &body.team_lead }, None)
        .await
    {
        Ok(None) => false,
        Ok(Some(_)) => true,
        Err(err) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error","message": format!("Failed to register team: {}", err.to_string())}),  // not informative 
        )
    };
    if !user_exists {
        return HttpResponse::Conflict().json(
            serde_json::json!({"status": "error", "message": "User (team lead) does not exist"}),
        );
    }

    // Create the new team from schema and insert, verify that the 
    let team_registration = body.into_inner();

    let new_team = Team::from(&team_registration);
    match team_collection.insert_one(&new_team, None).await {
        Ok(_) => {
            let json_response = serde_json::json!({"status": "success", "message": "Added team",  "data": serde_json::json!({"team": &new_team})});
            HttpResponse::Ok().json(json_response)
        },
        Err(err) => HttpResponse::InternalServerError().json(serde_json::json!({"status": "error", "message": format!("Failed to register team: {}", err.to_string())}))
    }

}



// Register a new database for a team
#[post("/teams/database")]
async fn register_team_database_handler(
    body: web::Json<RegisterDatabaseSchema>,
    data: web::Data<AppState>,
    auth_query: web::Query<TeamAccessQuery>,
    auth_guard: jwt::JwtDataMiddleware 
) -> impl Responder {

    let team_collection: Collection<Team> = get_cerebro_db_collection(&data, AdminCollection::Teams);

    // Get team
    let user_team = match team_collection
        .find_one(
            doc! { 
                "$and": [
                    { 
                        "$or": [
                            { "name": &auth_query.team },
                            { "id": &auth_query.team }
                        ]
                    },
                    { "users": &auth_guard.user.id },
                ]
            },
            None
        )
        .await
    {
        Ok(None) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail","message": "Failed to find user team", "data": serde_json::json!({"team_db": []})}), 
        ),
        Ok(Some(team)) => team,
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error","message": "Failed to retrieve user team", "data": serde_json::json!({"team_db": []})}),
        )
    };

    if user_team.database_name_exists(&body.database_name) {
        return HttpResponse::Conflict().json(serde_json::json!({"status": "error", "message": "Database name already exists for this team", "data": serde_json::json!({"team_db": []})}))
    }

    let new_db = TeamDatabase::from_database_schema(&body.into_inner());

    match team_collection.update_one(doc! { "id": user_team.id }, doc! { "databases": { "$push": &mongodb::bson::to_bson(&new_db).unwrap() } }, None).await {
        Ok(_) => {
            let json_response = serde_json::json!({"status": "success", "message": "Added team database",  "data": serde_json::json!({"team_db": &new_db})});
            HttpResponse::Ok().json(json_response)
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({"status": "error", "message": "Failed to register database",  "data": serde_json::json!({"team_db": &new_db})}))
    }

}


// Register a new project for a team database
#[post("/teams/project")]
async fn register_team_database_project_handler(
    body: web::Json<RegisterProjectSchema>,
    data: web::Data<AppState>,
    auth_query: web::Query<TeamDatabaseAccessQuery>,  // auth guard JWtDataMiddleware
    auth_guard: jwt::JwtDataMiddleware 
) -> impl Responder {

    // Add checks for project existence

    let team_collection: Collection<Team> = get_cerebro_db_collection(&data, AdminCollection::Teams);

    // Check if team exists already using its team name OR
    // its database name field OR using either of their 
    // unique identifiers
    let user_team = match team_collection
        .find_one(
            doc! { 
                "$and": [
                    { 
                        "$or": [
                            { "name": &auth_query.team },
                            { "id": &auth_query.team }
                        ]
                    },
                    { "users": &auth_guard.user.id },
                    { 
                        "$or": [
                            { "databases.name": &auth_query.db },
                            { "databases.id": &auth_query.db }
                        ]
                    },
                ]
            },
            None
        )
        .await
    {
        Ok(Some(team)) => team,
        Ok(None) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": format!("Failed to find team {} with database {}", &auth_query.team, &auth_query.db), "data": serde_json::json!({})}),
        ),
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error","message": "Failed to retrieve team with database", "data": serde_json::json!({})}),
        )
    };

    if user_team.project_name_exists(&auth_query.db, &body.project_name) {
        return HttpResponse::Conflict().json(serde_json::json!({"status": "error", "message": "Database name already exists for this team", "data": serde_json::json!({"team_db": []})}))
    }

    let database = match get_database_by_name(&user_team.databases, &auth_query.db) {
        Ok(db) => db,
        Err(err) => return err
    };

    let project_update_options = mongodb::options::UpdateOptions::builder().array_filters(
        vec![doc! { "db.id": &database.id }]
    ).build();

    let new_project = ProjectCollection::from_project_schema(&body.into_inner());
    let update = doc! { "$push": { "databases.$[db].projects": &mongodb::bson::to_bson(&new_project).unwrap() } };

    match team_collection.update_one(doc! { "id": user_team.id }, update, project_update_options).await {
        Ok(_) => {
            let json_response = serde_json::json!({"status": "success", "message": "Added project to team database", "data": serde_json::json!({"team_project": &new_project})});
            HttpResponse::Ok().json(json_response)
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({"status": "error", "message": "Failed to register project with team database", "data": serde_json::json!({"team_project": &new_project})}))
    }

}



#[derive(Deserialize)]
struct TeamQuery {
    team_id: Option<TeamId>,
    user_id: Option<TeamId>
}

// Get a specific team (by user or id) or all team data [ADMIN]
#[get("/teams")]
async fn get_team_handler(data: web::Data<AppState>, query: web::Query<TeamQuery>, _: jwt::JwtAdminMiddleware) -> impl Responder {
    
    let query = query.into_inner();

    let team_collection: Collection<Team> = get_cerebro_db_collection(&data, AdminCollection::Teams);

    match (&query.team_id, &query.user_id)  {
        (Some(team_id), None) => {
            match team_collection.find_one(doc! { "id": team_id }, None).await 
            {
                Ok(None) => {
                    return HttpResponse::NotFound().json(
                        serde_json::json!({"status": "fail", "message": "Team does not exist"})
                    )
                },
                Ok(Some(team)) => {
                    let json_response = serde_json::json!({
                        "status": "success",
                        "message": "Found team",
                        "data": serde_json::json!({"team": &team})
                    });
                    return HttpResponse::Ok().json(json_response)
                }
                Err(_) => return HttpResponse::InternalServerError().json(
                    serde_json::json!({"status": "error", "message": "Failed to execute team search"})
                ),
            }
        },
        (None, Some(user_id)) => {
            match team_collection.find(doc! { "users": user_id }, None).await 
            {
                Ok(cursor) => {
                    let teams = cursor.try_collect().await.unwrap_or_else(|_| vec![]);  // allow empty teams database
                    let json_response = serde_json::json!({
                        "status":  "success",
                        "message": "Listed all teams with requested user",
                        "data": serde_json::json!({"teams": teams })
                    });
                    return HttpResponse::Ok().json(json_response)
                },
                Err(err) => return HttpResponse::InternalServerError().json(
                    serde_json::json!({"status": "error", "message": format_args!("{}", err)})
                )
            }
        },
        (Some(team_id), Some(user_id)) => {
            match team_collection.find_one(doc! { "id": team_id, "users": user_id }, None).await 
            {
                Ok(None) => {
                    return HttpResponse::NotFound().json(
                        serde_json::json!({"status": "fail", "message": "Team with user does not exist"})
                    )
                },
                Ok(Some(team)) => {
                    let json_response = serde_json::json!({
                        "status": "success",
                        "message": "Found team with user",
                        "data": serde_json::json!({"team": &team})
                    });
                    return HttpResponse::Ok().json(json_response)
                }
                Err(_) => return HttpResponse::InternalServerError().json(
                    serde_json::json!({"status": "error", "message": "Failed to execute team search"})
                ),
            }
        },
        (None, None) => {
            match team_collection.find(None, None).await 
            {
                Ok(cursor) => {
                    let teams = cursor.try_collect().await.unwrap_or_else(|_| vec![]);  // allow empty teams database
                    let json_response = serde_json::json!({
                        "status":  "success",
                        "message": "Listed all teams",
                        "data": serde_json::json!({"teams": teams })
                    });
                    return HttpResponse::Ok().json(json_response)
                },
                Err(err) => return HttpResponse::InternalServerError().json(
                    serde_json::json!({"status": "error", "message": format_args!("{}", err)})
                )
            }
        }        
    }
}

// Delete a team from all members [ADMIN]
#[delete("/teams/{id}")]
async fn delete_team_handler(data: web::Data<AppState>, id: web::Path<TeamId>, _: jwt::JwtAdminMiddleware) -> impl Responder {

    let team_id = id.into_inner();
    let team_collection: Collection<Team> = get_cerebro_db_collection(&data, AdminCollection::Teams);

    // Check if the team exists at all
    let team = match team_collection
        .find_one(doc! { "id": &team_id }, None)
        .await
    {
        Ok(None) => None,
        Ok(Some(team)) => Some(team),
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error","message": "Failed to search team"}),
        )
    };

    let team = match team {
        None =>  return HttpResponse::NotFound().json(
            serde_json::json!({"status": "error", "message": "Team does not exist"}),
        ),
        Some(team) => team
    };

    // Delete team from admin collection
    match team_collection.delete_one(doc! { "id": &team_id }, None).await {
        Ok(_) => {
            // Drop all team databases
            for team_database in team.databases {
                let mongo_db = data.db.database(&team_database.database);
                if let Err(err) = mongo_db.drop(None).await {
                    return HttpResponse::InternalServerError().json(
                        serde_json::json!({"status": "error", "message": format!("Failed to drop team database: {}", err.to_string())})
                    )
                }
            }
            let json_response = serde_json::json!({"status": "success", "message": "Deleted team", "data": serde_json::json!({"id": &team_id})});
            return HttpResponse::Ok().json(json_response)
        }
        Err(err) => return HttpResponse::InternalServerError().json(serde_json::json!({"status": "error", "message": err.to_string()}))
    }
        

}


// Update a team by identifier [ADMIN]
#[patch("/teams/{id}")]
async fn update_team_handler(data: web::Data<AppState>, id: web::Path<TeamId>, body: web::Json<UpdateTeamSchema>, _: jwt::JwtAdminMiddleware) -> impl Responder {

    let team_id = id.into_inner();

    let team_collection: Collection<Team> = get_cerebro_db_collection(&data, AdminCollection::Teams);

    // Check if the team exists
    let exists = match team_collection
        .find_one(doc! { "id": &team_id }, None)
        .await
    {
        Ok(None) => false,
        Ok(Some(_)) => true,
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error","message": "Failed to search team"}),
        )
    };
    if !exists {
        return HttpResponse::NotFound().json(
            serde_json::json!({"status": "error", "message": "Team does not exist"}),
        );
    }

    // Check if the proposed team name already exists if it is provided
    if let Some(team_name) = &body.team_name {
        let name_already_exists = match team_collection
            .find_one(doc! { "name": &team_name }, None)
            .await
        {
            Ok(None) => false,
            Ok(Some(_)) => true,
            Err(_) => return HttpResponse::InternalServerError().json(
                serde_json::json!({"status": "error","message": "Failed to search team"}),
            )
        };
        if name_already_exists {
            return HttpResponse::NotFound().json(
                serde_json::json!({"status": "error", "message": "Proposed team name already exists"}),
            );
        }
    };

    let update = match &body.team_name {
        Some(team_name) => {
            doc! { 
                "$set": {
                    "name": team_name,
                    "description": &body.team_description
                }        
            }
        },
        None => {
            doc! { 
                "$set": {
                    "description": &body.team_description
                }        
            }
        }
    };

    match team_collection.update_one(
        doc! { "id": &team_id }, update, None
    ).await {
        Ok(_) => {
            let json_response = serde_json::json!(
                {"status": "success",  "message": "Updated team", "data": serde_json::json!({"id": &team_id})
            });
            HttpResponse::Ok().json(json_response)
        },
        Err(err) => HttpResponse::InternalServerError().json(serde_json::json!({"status": "error", "message": format!("Failed to register team: {}", err.to_string())}))
    }
}

#[derive(Debug, Deserialize)]
#[serde(rename_all ="lowercase")]
pub enum TeamUserOperation {
    Add,
    Remove
}

#[derive(Debug, Deserialize)]
struct TeamUserQuery {
    team_id: TeamId,
    user_id: UserId,
    update: TeamUserOperation
}

#[patch("/teams")]
async fn update_team_user_handler(data: web::Data<AppState>, query: web::Query<TeamUserQuery>,  _: jwt::JwtAdminMiddleware) -> impl Responder {

    let query = query.into_inner();

    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);  
    let team_collection: Collection<Team> = get_cerebro_db_collection(&data, AdminCollection::Teams);

    // Check if the user exists
    let user_exists = match user_collection
        .find_one(doc! { "id": &query.user_id }, None)
        .await
    {
        Ok(None) => false,
        Ok(Some(_)) => true,
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error", "message": "Failed to search user"}),
        )
    };
    if !user_exists {
        return HttpResponse::NotFound().json(
            serde_json::json!({"status": "error", "message": "User does not exist"}),
        );
    }

    // Check if the team exists
    let team = match team_collection
        .find_one(doc! { "id": &query.team_id }, None)
        .await
    {
        Ok(None) => return HttpResponse::NotFound().json(
            serde_json::json!({"status": "error", "message": "Team does not exist"}),
        ),
        Ok(Some(team)) => team,
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error","message": "Failed to search team"}),
        )
    };
    

    match query.update {
        TeamUserOperation::Add => {
            if team.users.contains(&query.user_id){
                return HttpResponse::Conflict().json(
                    serde_json::json!({"status": "error", "message": "User is already a member of the team"}),
                )
            }
            match team_collection.update_one(
                doc! { "id": &query.team_id },
                doc! { "$push": { "users": &query.user_id } }, None
            ).await {
                Ok(_) => { 
                    let json_response = serde_json::json!({
                        "status": "success", "message": "User added to team", "data": serde_json::json!({
                        "team_id": &query.team_id, "user_id": &query.user_id
                    })});
                    HttpResponse::Ok().json(json_response)
                 }, 
                Err(err) => return HttpResponse::InternalServerError().json(
                    serde_json::json!({"status": "error", "message": format!("Failed to register team: {}", err.to_string())})
                )
            }
            
        }, 
        TeamUserOperation::Remove => {
            if !team.users.contains(&query.user_id){
                return HttpResponse::Conflict().json(
                    serde_json::json!({"status": "error", "message": "User is not a member of the team"}),
                )
            }
            match team_collection.update_one(
                doc! { "id": &query.team_id },
                doc! { "$pull": { "users": &query.user_id } }, None
            ).await {
                Ok(_) => {
                    let json_response = serde_json::json!({
                        "status": "success", "message": "User removed from team", "data": serde_json::json!({
                        "team_id": &query.team_id, "user_id": &query.user_id
                    })});
                    HttpResponse::Ok().json(json_response)
                },
                Err(err) => HttpResponse::InternalServerError().json(
                    serde_json::json!({"status": "error", "message": format!("Failed to register team: {}", err.to_string())})
                )
            }
        }
    }
    
    
}


fn get_database_by_name(databases: &Vec<TeamDatabase>, db: &str) -> Result<TeamDatabase, HttpResponse>   {
    let matches: Vec<&TeamDatabase> = databases.into_iter().filter(|x| x.name == db || x.id == db ).collect();

    if matches.len() > 0 {
        Ok(matches[0].to_owned())
    } else {
        Err(HttpResponse::InternalServerError().json(serde_json::json!({"status": "fail", "message": "Failed to find database by name or identifier",  "data": serde_json::json!({"team_db": []})})))
    }
}

// Handler configuration
pub fn team_config(cfg: &mut web::ServiceConfig) {
    cfg.service(register_team_handler)
    .service(delete_team_handler)
    .service(update_team_handler)
    .service(get_team_handler)
    .service(register_team_database_handler)
    .service(register_team_database_project_handler)
    .service(update_team_user_handler);
}