//! Actix middleware for JWT authentication workflow


use core::fmt;
use std::collections::HashMap;
use std::future::{ready, Ready};
use futures::executor::block_on;

use actix_web::error::{ErrorInternalServerError, ErrorUnauthorized, ErrorForbidden, ErrorBadRequest};
use actix_web::{dev::Payload, Error as ActixWebError};
use actix_web::{http, web, FromRequest, HttpRequest};

use mongodb::Collection;
use serde::{Serialize, Deserialize};
use redis::Commands;

use crate::api::auth::token;
use crate::api::logs::model::Action;
use crate::api::logs::utils::log_admin_auth;
use crate::api::server::AppState;
use crate::api::teams::model::{Team, DatabaseId, ProjectId};
use crate::api::users::model::{User, Role};
use crate::api::utils::get_cerebro_db_collection;


#[derive(Debug, Serialize)]
struct ErrorResponse {
    status: String,
    message: String,
}
impl fmt::Display for ErrorResponse {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", serde_json::to_string(&self).unwrap_or(String::from("Failed error response serialization for display")))
    }
}

#[derive(Debug, Serialize)]
struct ErrorRefreshResponse {
    status: String,
    message: String,
    refresh: bool,
}
impl fmt::Display for ErrorRefreshResponse {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", serde_json::to_string(&self).unwrap_or(String::from("Failed error response serialization for display")))
    }
}

#[derive(Debug, Deserialize)]
pub struct DataAccessQuery {
    db: Option<DatabaseId>,
    project: Option<ProjectId>
}
impl DataAccessQuery {
    pub fn from_request(req: &HttpRequest) -> Self {
        let query_str = req.query_string();
        match web::Query::<HashMap<String, String>>::from_query(query_str) {
            Ok(query_map) => {
                match (query_map.get("db"), query_map.get("project")) {
                    (Some(db), Some(project)) => {
                        Self { db: Some(db.to_owned()), project: Some(project.to_owned())  }
                    },
                    _ => Self { db: None, project: None }
                }
            },
            Err(_) => Self { db: None, project: None }  // in case of failure do not return a valid data query struct
        }
    }
    pub fn is_valid(&self) -> bool {
        self.db.is_some() && self.project.is_some()
    }
}


// When a user attempts to access a protected route, the JWT middleware will first 
// search for the token in the Authorization header. If it is not found there, it 
// will check the Cookies object for the access_token key. If the token cannot be 
// found in either location, the middleware will send a 401 Unauthorized response 
// with the message “You are not logged in, please provide token” to the client.
//
// However, if the token is present, the middleware will call token::verify_jwt_token() 
// function to verify its authenticity. If the token is valid, the function will return 
// the token’s metadata, which will be used in the next step to query the Redis 
// database to check whether the user associated with the token has a valid session.
// 
// If the token metadata is found in the Redis database, it indicates that the user’s 
// session is still active. If the user has an active session, the middleware will use
// the user’s ID returned from the Redis query to check if the user associated
// with the token exists in the MongoDB database. If the user is found, the middleware
// will return the corresponding record obtained from the query.
//
// I'm unsure how to make variants of this middleware for role-based access - for now 
// this is implemented as the same core struct with a different role check of the user
// before the middleware returns.

fn get_token_from_auth_header(value: &http::header::HeaderValue) -> Option<String> {
    match value.to_str() {
        Ok(str) => {
            if 7 > str.len() {
                return None
            }
            Some(str.split_at(7).1.to_string())
        },
        Err(_) => return None
    }
}

pub struct JwtUserMiddleware {
    pub user: User,
    pub access_token_uuid: uuid::Uuid,
    pub team: Option<Team>
}
impl FromRequest for JwtUserMiddleware {
    type Error = ActixWebError;
    type Future = Ready<Result<Self, Self::Error>>;

    fn from_request(req: &HttpRequest, _: &mut Payload) -> Self::Future {

        // Get the application state data for access to database clients
        let data = match req.app_data::<web::Data<AppState>>() {
            Some(data) => data,
            None => {
                let json_error = ErrorResponse {
                    status: "error".to_string(),
                    message: "Failed to get application state".to_string(), 
                };
                return ready(Err(ErrorInternalServerError(json_error)));
            }
        };

        // Parse the access token from cookie or authorization header of request
        let access_token = req
            .cookie("access_token")
            .map(|c| c.value().to_string())
            .or_else(|| {

                let token = req.headers()
                    .get(http::header::AUTHORIZATION)
                    .map(|h| get_token_from_auth_header(h));  
                
                match token {
                    Some(header_token) => header_token,
                    None => None
                }

            });

        let access_token = match access_token {
            Some(token) => token,
            None => {
                let json_error = ErrorRefreshResponse {  // this may happen if client-side fetch is called after cookie expiration (and no server-side refresh initiated)
                    status: "fail".to_string(),
                    refresh: true,
                    message: "Not authorized".to_string(),
                };
                return ready(Err(ErrorUnauthorized(json_error)));
            }
        };

        // Verify the access token and provide the extracted payload as TokenDetails
        let access_token_details = match token::verify_jwt_token(
            data.env.security.token.encryption.access_public_key.to_owned(),
            &access_token,
        ) {
            Ok(token_details) => token_details,
            Err(_) => {
                let json_error = ErrorRefreshResponse {  // token may be expired (max_age)
                    status: "fail".to_string(),
                    refresh: true,
                    message: "Not authorized".to_string(),
                };
                return ready(Err(ErrorUnauthorized(json_error)));
            }
        };

          // Parse the access token Uuid string as Uuid
          let access_token_uuid = match uuid::Uuid::parse_str(&access_token_details.token_uuid.to_string()) {
            Ok(token_uuid) => token_uuid,
            Err(_) => {
                return ready(Err(ErrorInternalServerError(ErrorResponse {
                    status: "error".to_string(),
                    message: format!("Not authorized"),
                })));
            }
        };

        // Check if the access token id returns a user id result from the Redis database (session persistence)
        let user_id_redis_result = async move {
            let mut redis_client = match data.auth_session.get_connection() {
                Ok(redis_client) => redis_client,
                Err(_) => {
                    return Err(ErrorInternalServerError(ErrorResponse {
                        status: "error".to_string(),
                        message: format!("Not authorized"),
                    }));
                }
            };

            let redis_result = redis_client.get::<_, String>(access_token_uuid.clone().to_string());

            match redis_result {
                Ok(value) => Ok(value),
                Err(_) => Err(ErrorUnauthorized(ErrorRefreshResponse {   // user session may be expired (max_age)
                    status: "fail".to_string(),
                    refresh: true,
                    message: "Not authorized".to_string(),
                })),
            }
        };

        // Check that the User exists in the Mongo database and return their data
        let user_exists_result = async move {

            // Await and parse the user Uuid from the Redis database result
            let user_id = user_id_redis_result.await?;
            
            let user_id_uuid = match uuid::Uuid::parse_str(user_id.as_str()) {
                Ok(user_uuid) => user_uuid,
                Err(_) => {
                    return Err(ErrorInternalServerError(ErrorResponse {
                        status: "error".to_string(),
                        message: format!("Not authorized"),
                    }));
                }
            };

            // Database connection and user id query
            let user_collection: Collection<User> = get_cerebro_db_collection(&data, "user");
                        
            let query_result = user_collection
                .find_one(mongodb::bson::doc! { "id": format!("{}", &user_id_uuid) }, None)
                .await;
            
            match query_result {
                Ok(Some(user)) => Ok(user),
                Ok(None) => {
                    let json_error: ErrorResponse = ErrorResponse {
                        status: "fail".to_string(),
                        message: "Not authorized".to_string(),
                    };
                    Err(ErrorUnauthorized(json_error))
                }
                Err(_) => {
                    let json_error = ErrorResponse {
                        status: "error".to_string(),
                        message: "Not authorized".to_string(),
                    };
                    Err(ErrorInternalServerError(json_error))
                }
            }
        };

        // Run user exists block future to completion 
        // and return the User and access token UUID
        match block_on(user_exists_result) {
            Ok(user) => {

                // USER ROLE AUTHORIZATION
                // This is always enabled for the user-accessible endpoints

                if !user.roles.contains(&Role::User) {
                    let json_error = ErrorResponse {
                        status: "fail".to_string(),
                        message: "You do not have enough permissions".to_string(),
                    };
                    return ready(Err(ErrorForbidden(json_error)));
                };

                // If the database and project query parameters are in the header,
                // place an additional guard on the user access to team database(s)
                // and projects that contain sample data

                let data_access_query = DataAccessQuery::from_request(&req);
                let user_id = user.id.clone();

                // If we received valid data query parameters for a data access endpoint, check 
                // access of authorized user to the requested team database and project
                if data_access_query.is_valid() {

                    let team_result = async move {


                        let database_query = match &data_access_query.db {
                            Some(query) => query,
                            None => {
                                let json_error = ErrorResponse {
                                    status: "fail".to_string(),
                                    message: "No database query specified".to_string(),
                                };
                                return Err(ErrorBadRequest(json_error));
                            }
                        };

                        let project_query = match &data_access_query.project {
                            Some(query) => query,
                            None => {
                                let json_error = ErrorResponse {
                                    status: "fail".to_string(),
                                    message: "No project query specified".to_string(),
                                };
                                return Err(ErrorBadRequest(json_error));
                            }
                        };

                        // Database connection and user id query
                        let team_collection: Collection<Team> = get_cerebro_db_collection(&data, "team");

                        let query_result = team_collection
                            .find_one(mongodb::bson::doc! {
                                "$and": [
                                    { "users": &user_id },
                                    { "databases.id": &database_query }, 
                                    { "databases.projects.id": &project_query } 
                                ]
                            }, None)
                            .await;
                        
                        match query_result {
                            Ok(Some(team)) => Ok(team),
                            Ok(None) => {
                                let json_error: ErrorResponse = ErrorResponse {
                                    status: "fail".to_string(),
                                    message: "You do not have access to the requested data".to_string(),
                                };
                                Err(ErrorForbidden(json_error))
                            }
                            Err(_) => {
                                let json_error = ErrorResponse {
                                    status: "error".to_string(),
                                    message: "Failed to check user permissions for data access".to_string(),
                                };
                                Err(ErrorInternalServerError(json_error))
                            }
                        }
                    };

                    match block_on(team_result) {
                        Ok(team) => {
                            ready(Ok(JwtUserMiddleware {
                                access_token_uuid,
                                user,
                                team: Some(team)
                            }))
                        },
                        Err(error) => return ready(Err(error))
                    }
    
                } else {
                    ready(Ok(JwtUserMiddleware {
                        access_token_uuid,
                        user,
                        team: None
                    }))
                }  // if we didn't receive a team query then just go ahead - all data access endpoints should have a team query parameter

            },
            Err(error) => ready(Err(error)),
        }
    }
}

// Re-implementation with different role check at the end - 
// not the ideal solution but works for role-based guards now
pub struct JwtAdminMiddleware {
    pub user: User,
    pub access_token_uuid: uuid::Uuid,
}
impl FromRequest for JwtAdminMiddleware {
    type Error = ActixWebError;
    type Future = Ready<Result<Self, Self::Error>>;

    fn from_request(req: &HttpRequest, _: &mut Payload) -> Self::Future {

        // Get the application state data for access to database clients
        let data = match req.app_data::<web::Data<AppState>>() {
            Some(data) => data,
            None => {
                let json_error = ErrorResponse {
                    status: "error".to_string(),
                    message: "Failed to get application state data".to_string(), 
                };
                return ready(Err(ErrorInternalServerError(json_error)));
            }
        };

        // Parse the access token from cookie or authorization header of request
        let access_token = req
            .cookie("access_token")
            .map(|c| c.value().to_string())
            .or_else(|| {
                let token = req.headers()
                    .get(http::header::AUTHORIZATION)
                    .map(|h| get_token_from_auth_header(h));
                
                match token {
                    Some(header_token) => header_token,
                    None => None
                }
            });

        let access_token = match access_token {
            Some(token) => token,
            None => {
                let json_error = ErrorRefreshResponse { // this may happen if client-side fetch is called after cookie expiration (and no server-side refresh initiated)
                    status: "error".to_string(),
                    refresh: true,
                    message: "Not authorized".to_string(),
                };
                return ready(Err(ErrorUnauthorized(json_error)));
            }
        };

        // Verify the access token and provide the extracted payload as TokenDetails
        let access_token_details = match token::verify_jwt_token(
            data.env.security.token.encryption.access_public_key.to_owned(),
            &access_token,
        ) {
            Ok(token_details) => token_details,
            Err(_) => {
                let json_error = ErrorRefreshResponse {
                    status: "fail".to_string(),
                    refresh: true,
                    message: "Not authorized".to_string(),
                };
                return ready(Err(ErrorUnauthorized(json_error)));
            }
        };

        // Parse the access token Uuid string as Uuid
        let access_token_uuid = match uuid::Uuid::parse_str(&access_token_details.token_uuid.to_string()) {
            Ok(token_uuid) => token_uuid,
            Err(_) => {
                return ready(Err(ErrorInternalServerError(ErrorResponse {
                    status: "error".to_string(),
                    message: format!("Not authorized"),
                })));
            }
        };

        // Check if the access token Uuuid returns a user Uuid result from the Redis database
        let user_id_redis_result = async move {
            let mut redis_client = match data.auth_session.get_connection() {
                Ok(redis_client) => redis_client,
                Err(_) => {
                    return Err(ErrorInternalServerError(ErrorResponse {
                        status: "error".to_string(),
                        message: format!("No authorized"),
                    }));
                }
            };

            let redis_result = redis_client.get::<_, String>(access_token_uuid.clone().to_string());

            match redis_result {
                Ok(value) => Ok(value),
                Err(_) => Err(ErrorUnauthorized(ErrorRefreshResponse {
                    status: "fail".to_string(),
                    refresh: true,
                    message: "Not authorized".to_string(),
                })),
            }
        };

        // Check that the User exists in the Mongo database and return their data
        let user_exists_result = async move {

            // Await and parse the user Uuid from the Redis database result
            let user_id = user_id_redis_result.await?;

            let user_id_uuid = match uuid::Uuid::parse_str(user_id.as_str()) {
                Ok(user_uuid) => user_uuid,
                Err(_) => {
                    return Err(ErrorInternalServerError(ErrorResponse {
                        status: "error".to_string(),
                        message: format!("Failed to parse user identifier"),
                    }));
                }
            };

            // Database connection and user id query
            let user_collection: Collection<User> = get_cerebro_db_collection(&data, "user");
            
            let query_result = user_collection
                .find_one(mongodb::bson::doc! { "id": format!("{}", &user_id_uuid) }, None)
                .await;
            
            match query_result {
                Ok(Some(user)) => Ok(user),
                Ok(None) => {

                    let json_error: ErrorResponse = ErrorResponse {
                        status: "fail".to_string(),
                        message: "User belonging to this token no longer exists".to_string(),
                    };

                    match log_admin_auth(
                        &data, 
                        Some(&user_id), 
                        None,
                        Action::InvalidAdmin, 
                        "User verification failed for admin access - user may not exist", 
                        &req
                    ).await
                    {
                        Ok(_) => Err(ErrorUnauthorized(json_error)),
                        Err(_) => Err(ErrorUnauthorized(json_error))
                    }
                }
                Err(_) => {
                    let json_error = ErrorResponse {
                        status: "error".to_string(),
                        message: "Failed to check user existence".to_string(),
                    };
                    Err(ErrorInternalServerError(json_error))
                }
            }
        };

        // Run user exists block future to completion 
        // and return the User and access token UUID
        match block_on(user_exists_result) {
            Ok(user) => {

                // ADMIN GUARD
                match user.roles.contains(&Role::Admin) {
                    true => {
                        ready(Ok(JwtAdminMiddleware {
                            access_token_uuid,
                            user,
                        }))
                    },
                    false => {
                        let json_error = ErrorResponse {
                            status: "fail".to_string(),
                            message: "You do not have permissions to access this resource".to_string(),
                        };
                        ready(Err(ErrorForbidden(json_error)))
                    }
                }                
            },
            Err(error) => ready(Err(error)),
        }
    }
}