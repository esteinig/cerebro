use cerebro_model::api::{teams::model::{Team, TeamId}, users::response::{FilteredUser, UserData, UserSelfResponse}, utils::AdminCollection};
use cerebro_model::api::users::model::{User, UserId, Role};
use cerebro_model::api::users::schema::{UpdateUserSchema, RegisterUserSchema};

use crate::api::utils::get_cerebro_db_collection;
use crate::api::auth::jwt;
use crate::api::server::AppState;

use argon2::{password_hash::SaltString, Argon2, PasswordHasher};
use futures::TryStreamExt;
use actix_web::{web, get, post, HttpResponse, Responder, delete, patch,};
use mongodb::{bson::doc, Collection};
use rand_core::OsRng;
use serde::Deserialize;


// Get the current user data from authentication middleware [USER]
#[get("/users/self")]
async fn get_user_self_handler(jwt_guard: jwt::JwtUserMiddleware) -> impl Responder {

    let response = UserSelfResponse {
        status: "success".to_string(),
        message: "User exists".to_string(),
        data: UserData {
            user: filter_user_record(&jwt_guard.user)
        }
    };

    HttpResponse::Ok().json(response)
}


#[derive(Deserialize)]
pub struct TeamQuery {
    id: Option<TeamId>,
    name: Option<String>
}

#[get("/users/self/teams")]
async fn get_user_self_teams_handler(
    data: web::Data<AppState>, 
    query: web::Query<TeamQuery>, 
    jwt_guard: jwt::JwtUserMiddleware,
) -> impl Responder {
    let team_collection: Collection<Team> = get_cerebro_db_collection(&data, AdminCollection::Teams);

    let mut filter = doc! { "users": &jwt_guard.user.id };

    // Add id or name to the query filter if present
    if let Some(identifier) = query.id.as_ref().map(|id| id.to_string()).or_else(|| query.name.clone()) {
        // Allow flexibility between id and name
        filter.insert("$or", vec![
            doc! { "id": &identifier },
            doc! { "name": &identifier }
        ]);
    }

    let result = if query.id.is_some() || query.name.is_some() {
        // Find one specific team
        match team_collection.find_one(filter, None).await {
            Ok(Some(team)) => HttpResponse::Ok().json(serde_json::json!({
                "status": "success",
                "message": "Obtained requested team for current user",
                "data": { "team": team }
            })),
            Ok(None) => HttpResponse::NotFound().json(serde_json::json!({
                "status": "error",
                "message": "Did not find requested team for current user"
            })),
            Err(err) => HttpResponse::InternalServerError().json(serde_json::json!({
                "status": "error",
                "message": format!("{}", err)
            })),
        }
    } else {
        // Find all teams for the current user
        match team_collection.find(filter, None).await {
            Ok(cursor) => {
                let teams: Vec<Team> = cursor.try_collect().await.unwrap_or_default(); // Allow empty teams
                HttpResponse::Ok().json(serde_json::json!({
                    "status": "success",
                    "message": "Listed all teams for current user",
                    "data": { "teams": teams }
                }))
            }
            Err(err) => HttpResponse::InternalServerError().json(serde_json::json!({
                "status": "error",
                "message": format!("{}", err)
            })),
        }
    };

    result
}

// Register a new user - this is an admin operation, as access is currently restricted/controlled [ADMIN]
#[post("/users")]
async fn register_user_handler(
    body: web::Json<RegisterUserSchema>,
    data: web::Data<AppState>,
    _: jwt::JwtAdminMiddleware 
) -> impl Responder {

    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);

    // Check if user exists by email
    let exists = match user_collection
        .find_one(doc! { "email": &body.email }, None)
        .await
    {
        Ok(None) => false,
        Ok(Some(_)) => true,
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error","message": "Failed to register user"}),  // not informative 
        )
    };
    if exists {
        return HttpResponse::Conflict().json(
            serde_json::json!({"status": "error", "message": "User with that email already exists"}),
        );
    }

    // Hash the password for storage
    let salt = SaltString::generate(&mut OsRng);
    let hashed_password = match Argon2::default().hash_password(body.password.as_bytes(), &salt) {
        Ok(password_hash) => password_hash.to_string(),
        Err(_) => {
            return HttpResponse::InternalServerError()
                .json(serde_json::json!({"status": "error", "message": "Failed to register user"}));  // not informative 
        }
    };
    

    // Create and insert the new user, return a filtered user response (no password)
    let new_user = User::from(&body.into_inner(), &hashed_password);
    match user_collection.insert_one(&new_user, None).await {
        Ok(_) => {
            let user_response = serde_json::json!({
                    "status": "success", 
                    "message": "Registers user", 
                    "data": serde_json::json!({
                        "user": filter_user_record(&new_user)
                    })
                }
            );
           HttpResponse::Ok().json(user_response)
        },
        Err(_) => HttpResponse::InternalServerError().json(serde_json::json!({"status": "error", "message": "Failed to register user"}))
    }

}

#[derive(Deserialize)]
struct UserQuery {
    id: Option<UserId>
}

// Get a specific user's data from database search [ADMIN]
#[get("/users")]
async fn get_user_handler(data: web::Data<AppState>, query: web::Query<UserQuery>, _: jwt::JwtAdminMiddleware) -> impl Responder {
    
    let query = query.into_inner();
    
    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);

    if let Some(id) = &query.id  {
        
        match user_collection.find_one(doc! { "id": id }, None).await 
        {
            Ok(None) => {
                return HttpResponse::BadRequest().json(
                    serde_json::json!({"status": "fail", "message": "Request parameters not specified correctly"})
                )
            },
            Ok(Some(user)) => {
                let json_response = serde_json::json!({
                    "status": "success",
                    "message": "Found user",
                    "data": serde_json::json!({"user": filter_user_record(&user)})
                });
                return HttpResponse::Ok().json(json_response)
            }
            Err(_) => return HttpResponse::InternalServerError().json(
                serde_json::json!({"status": "error", "message": "Failed to execute user search"})
            ),
        }
    } else {
        
        match user_collection.find(None, None).await 
        {
            Ok(cursor) => {
                let users = cursor.try_collect().await.unwrap_or_else(|_| vec![]);
                let filtered_users: Vec<FilteredUser> = users.iter().map(|user| filter_user_record(&user)).collect();
                let json_response = serde_json::json!({
                    "status":  "success",
                    "message": "Listed all users",
                    "data": serde_json::json!({"users": filtered_users})
                });
                return HttpResponse::Ok().json(json_response)
            },
            Err(err) => return HttpResponse::InternalServerError().json(
                serde_json::json!({"status": "error", "message": format_args!("{}", err)})
            )
        }
    }
}


// Delete a user from the database [ADMIN]
#[delete("/users/{id}")]
async fn delete_user_handler(data: web::Data<AppState>, id: web::Path<UserId>, _: jwt::JwtAdminMiddleware) -> impl Responder {

    let user_id =  &id.into_inner();

    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);
    let team_collection: Collection<Team> = get_cerebro_db_collection(&data, AdminCollection::Teams);

    // First remove from teams - easier fix than the other way around
    match team_collection.update_many(doc! {}, doc! { "$pull": { "users": &user_id } }, None).await 
    {
        Ok(_) => { } // continue, if the user did not exist it would have either removed them (if leftover) or not, as they would not exist in team users array
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "error", "message": "Failed to remove user from teams"})
        ),
    }

    match user_collection.find_one_and_delete(doc! { "id": &user_id }, None).await 
        {
            Ok(None) => {
                return HttpResponse::NotFound().json(
                    serde_json::json!({"status": "fail", "message": "User could not be found"})
                )
            },
            Ok(Some(user)) => {
                let json_response = serde_json::json!({
                    "status": "success", "message": "Deleted user", "data": serde_json::json!({
                        "user": filter_user_record(&user)
                    })
                });
                return HttpResponse::Ok().json(json_response)
            }
            Err(_) => return HttpResponse::InternalServerError().json(
                serde_json::json!({"status": "error", "message": "Failed to execute user deletion"})
            ),
        }
}

// Update a user (as admin access to all users is provided) [USER/ADMIN]
#[patch("/users/{id}")]
async fn patch_user_handler(data: web::Data<AppState>, id: web::Path<UserId>, body: web::Json<UpdateUserSchema>, auth_guard: jwt::JwtDataMiddleware) -> impl Responder {

    let user_id: UserId = id.into_inner();

    // Privilege safe-guard for user updates
    // if the authenticated user does not match
    // the user to be updated
    if &auth_guard.user.id != &user_id {
        // If the authenticated requester does not have
        // admin priviliges return an error
        if !auth_guard.user.roles.contains(&Role::Admin) {
            return HttpResponse::Forbidden()
                .json(serde_json::json!({"status": "error", "message": "No permissions to update user"}));  // not informative 
        }
        // Otherwise the authenticated requester has admin
        // priviliges and can update another user
    }

    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);

    // Name, email and roles must always be provided - password update is optional (admin-side, deactivated in frontend for now)
    let update = match &body.password {
        Some(pwd) => {
            
            // Hash the password for storage
            let salt = SaltString::generate(&mut OsRng);
            let hashed_password = match Argon2::default().hash_password(&pwd.as_bytes(), &salt) {
                Ok(password_hash) => password_hash.to_string(),
                Err(_) => {
                    return HttpResponse::InternalServerError()
                        .json(serde_json::json!({"status": "error", "message": "Failed to update password"})); 
                }
            };
            
            doc! {
                "$set": {
                    "name": &body.name,
                    "email": &body.email,
                    "title": &body.title,
                    "positions": &body.positions,
                    "password": hashed_password,
                    "verfied": &body.verified,
                    "roles": &mongodb::bson::to_bson(&body.roles).unwrap(),  // see if we need to handle the unwrap call, nested array
                    "updated": &mongodb::bson::to_bson(&chrono::Utc::now()).unwrap()  // see if we need to handle the unwrap call, complex type
                }
            }
        },
        None => doc! {
            "$set": {
                "name": &body.name,
                "email": &body.email,
                "title": &body.title,
                "positions": &body.positions,
                "verfied": &body.verified,
                "roles": &mongodb::bson::to_bson(&body.roles).unwrap(),  // see if we need to handle the unwrap call, required since nested array
                "updated": &mongodb::bson::to_bson(&chrono::Utc::now()).unwrap()  // see if we need to handle the unwrap call, required since complex type
            }
        }
    };


    match user_collection.find_one_and_update(doc! { "id": &user_id }, update, None).await 
        {
            Ok(None) => {
                return HttpResponse::NotFound().json(
                    serde_json::json!({"status": "fail", "message": "User does not exist"})
                )
            },
            Ok(Some(updated_user)) => {
                let json_response = serde_json::json!({
                    "status":  "success", "message": "Updated user", "data": serde_json::json!({
                        "user": filter_user_record(&updated_user)
                    })
                });
                return HttpResponse::Ok().json(json_response)
            }
            Err(_) => return HttpResponse::InternalServerError().json(
                serde_json::json!({"status": "error", "message": "Failed to update user"})
            ),
        }
}


/// Helper function to filter a user database model into
/// the response model to remove sensitive information
pub fn filter_user_record(user: &User) -> FilteredUser {
    FilteredUser {
        id: user.id.to_owned(),
        email: user.email.to_owned(),
        name: user.name.to_owned(),
        title: user.title.to_owned(),
        image: user.image.to_owned(),
        verified: user.verified.to_owned(),
        roles: user.roles.to_owned(),
        positions: user.positions.to_owned(),
        created: user.created.to_owned(),
        updated: user.updated.to_owned(),
    }
}

// Handler configuration
pub fn user_config(cfg: &mut web::ServiceConfig) {
    cfg.service(get_user_self_handler)
    .service(get_user_handler)
    .service(delete_user_handler)
    .service(patch_user_handler)
    .service(register_user_handler)
    .service(get_user_self_teams_handler);
}