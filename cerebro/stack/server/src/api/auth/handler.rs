use crate::api::auth::{jwt, token};
use cerebro_model::api::auth::schema::{AuthLoginSchema, AuthResetPasswordSchema, AuthOneTimeCheckSchema, TokenDetails};
use cerebro_model::api::logs::model::Action;
use cerebro_model::api::users::model::{User, UserId, Role};

use cerebro_model::api::config::Config;
use cerebro_model::api::utils::AdminCollection;

use crate::api::utils::get_cerebro_db_collection;
use crate::api::email::Email;
use crate::api::server::AppState;
use crate::api::logs::utils::log_user_login;

use mongodb::{bson::doc, Collection};
use std::{fmt, str::FromStr};

use actix_web::{
    cookie::{time::Duration as ActixWebDuration, Cookie, SameSite},
    get, post, web, HttpRequest, HttpResponse, Responder, http::header::ContentType,
};
use argon2::{
    password_hash::{PasswordHash, PasswordVerifier, SaltString},
    Argon2, PasswordHasher
};

use rand_core::OsRng;
use redis::AsyncCommands;
use serde::Deserialize;
use uuid::Uuid;

// Rate limiting middleware factories

use actix_governor::{Governor, GovernorConfigBuilder, PeerIpKeyExtractor};
use governor::middleware::NoOpMiddleware;
use governor::clock::QuantaInstant;

fn login_rate_limiter_factory() -> Governor<PeerIpKeyExtractor, NoOpMiddleware<QuantaInstant>> {
    Governor::new(&GovernorConfigBuilder::default().per_second(4).burst_size(2).finish().unwrap())
}

// Make sure this is slightly more permissive as refresh tokens can be called
// frequently depending on access token expiration. This may need to be 
// configurable somehow - maybe through the service configs loaded in the
// main routine of api::server::App using actix_web::web::scopes

fn refresh_rate_limiter_factory() -> Governor<PeerIpKeyExtractor, NoOpMiddleware<QuantaInstant>> {
    Governor::new(&GovernorConfigBuilder::default().per_second(4).burst_size(2).finish().unwrap())
}

// Accoutns for scanners e.g. after emai lprovider links
fn one_time_token_rate_limiter_factory() -> Governor<PeerIpKeyExtractor, NoOpMiddleware<QuantaInstant>> {
    Governor::new(&GovernorConfigBuilder::default().per_second(8).burst_size(6).finish().unwrap())
}


#[derive(Deserialize)]
struct RoleQuery {
    role: Option<Role>
}

#[post("/auth/login", wrap="login_rate_limiter_factory()")]
async fn login_handler(
    req: HttpRequest,
    body: web::Json<AuthLoginSchema>,
    data: web::Data<AppState>,
    query: web::Query<RoleQuery>,
    // Unprotected route - critical to review policies! 
) -> impl Responder {

    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);  

    // Find the user in the database
    let user = match user_collection
        .find_one(doc! { "email": &body.email })
        .await
        .unwrap()
    {
        Some(user) => user,
        None => {
            // Log failed user email
            match log_user_login(&data, &body.email, Action::InvalidEmail, "User login attempt failed due to invalid email", &req).await {
                Ok(_) => return HttpResponse::Forbidden().json(serde_json::json!({"status": "fail", "message": "Not authorized"})),
                Err(_) => return HttpResponse::Forbidden().json(serde_json::json!({"status": "fail", "message": "Not authorized"}))
            };
        }
    };

    // Check for valid password
    let password_is_valid = PasswordHash::new(&user.password)
        .and_then(|parsed_hash| {
            Argon2::default().verify_password(body.password.as_bytes(), &parsed_hash)
        })
        .map_or(false, |_| true);

    if !password_is_valid {
        // Log failed user password
        match log_user_login(&data, &body.email, Action::InvalidPassword, "User login attempt failed due to invalid password", &req).await {
            Ok(_) => return HttpResponse::Forbidden().json(serde_json::json!({"status": "fail", "message": "Not authorized"})),
            Err(_) => return HttpResponse::Forbidden().json(serde_json::json!({"status": "fail", "message": "Not authorized"}))
        };
    }

    if !user.verified {
        // Log failed user verification status
        match log_user_login(&data, &body.email, Action::InvalidVerification, "User login attempt failed due to invalid verification", &req).await {
            Ok(_) => return HttpResponse::Forbidden().json(serde_json::json!({"status": "fail", "message": "Not authorized"})),
            Err(_) => return HttpResponse::Forbidden().json(serde_json::json!({"status": "fail", "message": "Not authorized"}))
        };
    }

    // Log successful user login attempt
    match log_user_login(&data, &body.email, Action::ValidCredentials, "User login attempt was successful", &req).await {
        Ok(_) => { }, Err(_) => return HttpResponse::InternalServerError().json(serde_json::json!({"status": "fail", "message": "Not authorized"}))
    };

    // Bot role should only be used to ge tan access token for backend pipeline automation 
    // and not continuous refresh token access. Not implemented on refresh token route!
    
    let (access_token_max_age, refresh_token_max_age) = match query.role {
        Some(Role::Bot) => {
            if user.roles.contains(&Role::Bot){
                match log_user_login(&data, &body.email, Action::ValidRole, "Bot authorization successful", &req).await {
                    Ok(_) => { }, Err(_) => return HttpResponse::InternalServerError().json(serde_json::json!({"status": "fail", "message": "Not authorized"}))
                };
                (data.env.security.token.expiration.access_max_age_bot, data.env.security.token.expiration.refresh_max_age_bot)
            } else {
                match log_user_login(&data, &body.email, Action::InvalidRole, "Bot authorization unsuccessful", &req).await {
                    Ok(_) => return HttpResponse::Forbidden().json(serde_json::json!({"status": "fail", "message": "Not authorized"})),
                    Err(_) => return HttpResponse::InternalServerError().json(serde_json::json!({"status": "fail", "message": "Not authorized"}))
                }; 
            }
        },
        _ => (data.env.security.token.expiration.access_max_age, data.env.security.token.expiration.refresh_max_age)
    };

    // Generate the access and refresh tokens using the private keys
    let access_token_details = match token::generate_jwt_token(
        uuid::Uuid::parse_str(&user.id).unwrap(),
        access_token_max_age,
        data.env.security.token.encryption.access_private_key.to_owned(),
    ) {
        Ok(token_details) => token_details,
        Err(_) => {
            return HttpResponse::BadGateway()
                .json(serde_json::json!({"status": "fail", "message": "Not authorized"})); 
        }
    };

    let refresh_token_details = match token::generate_jwt_token(
        uuid::Uuid::parse_str(&user.id).unwrap(),
        refresh_token_max_age,
        data.env.security.token.encryption.refresh_private_key.to_owned(),
    ) {
        Ok(token_details) => token_details,
        Err(_) => {
            return HttpResponse::BadGateway()
                .json(serde_json::json!({"status": "fail", "message": "Not authorized"})); 
        }
    };

    // Get a connection to the RedisDB for token storage
    let mut redis_client = match data.auth_session.get_async_connection().await {
        Ok(redis_client) => redis_client,
        Err(_) => {
            return HttpResponse::InternalServerError()
                .json(serde_json::json!({"status": "fail", "message": "Not authorized"}));
        }
    };

    // Insert the access token (key) with user identifier (value) and expiration of the key
    let access_result: redis::RedisResult<()> = redis_client
        .set_ex(
            access_token_details.token_uuid.to_string(),
            user.id.to_string(),
            (data.env.security.token.expiration.access_max_age * 60) as usize,
        )
        .await;

    if let Err(_) = access_result {
        return HttpResponse::UnprocessableEntity()
            .json(serde_json::json!({"status": "error", "message": "Not authorized"}));  // change to opaque message
    }


    // Insert the refresh token (key) with user identifier (value) and expiration of the key
    let refresh_result: redis::RedisResult<()> = redis_client
        .set_ex(
            refresh_token_details.token_uuid.to_string(),
            user.id.to_string(),
            (data.env.security.token.expiration.refresh_max_age * 60) as usize,
        )
        .await;

    if let Err(_) = refresh_result {
        return HttpResponse::UnprocessableEntity()
            .json(serde_json::json!({"status": "error", "message": "Not authorized"}));  // change to opaque message
    }

    
    // Set the access, refresh (with the encrypted token) and logged-in cookies
    let access_cookie = Cookie::build("access_token", access_token_details.token.clone().unwrap())
        .path("/")
        .max_age(ActixWebDuration::new(access_token_max_age * 60, 0))
        .domain(&data.env.security.cookies.domain) // domain must be set for cookies to be available on the application subdomain
        .http_only(true)
        .secure(data.env.security.cookies.secure)
        .same_site(SameSite::Strict)
        .finish();

    let refresh_cookie = Cookie::build("refresh_token", refresh_token_details.token.clone().unwrap())
        .path("/")
        .max_age(ActixWebDuration::new(refresh_token_max_age * 60, 0))
        .domain(&data.env.security.cookies.domain)
        .http_only(true)
        .secure(data.env.security.cookies.secure)
        .same_site(SameSite::Strict)
        .finish();

    // Return the cookies and access token for authorization
    HttpResponse::Ok()
        .cookie(access_cookie)
        .cookie(refresh_cookie)
        .json(serde_json::json!({
            "status": "success", 
            "access_token": access_token_details.token.unwrap(), 
            "refresh_token": refresh_token_details.token.unwrap()
        }))
}

#[get("/auth/refresh", wrap="refresh_rate_limiter_factory()")]
async fn refresh_access_token_handler(
    req: HttpRequest,
    data: web::Data<AppState>,
    // Unprotected route - critical to review policies! 
) -> impl Responder {

    // Get the refresh token from the cookies
    let refresh_token = match req.cookie("refresh_token") {
        Some(c) => c.value().to_string(),
        None => {
            return HttpResponse::Forbidden()
                .json(serde_json::json!({"status": "fail", "message": "You are not logged in, please provide refresh token"}));
        }
    };

    // Decode and verify the refresh token
    let refresh_token_details =
        match token::verify_jwt_token(data.env.security.token.encryption.refresh_public_key.to_owned(), &refresh_token)
        {
            Ok(token_details) => token_details,
            Err(_) => {
                return HttpResponse::Forbidden().json(
                    serde_json::json!({"status": "fail", "message": "Refresh token is invalid or session has expired"}),
                );
            }
        };
    
    // Initiate RedisDB connection and get user identfier based on refresh token UUID
    let result = data.auth_session.get_async_connection().await;
    let mut redis_client = match result {
        Ok(redis_client) => redis_client,
        Err(e) => {
            return HttpResponse::Forbidden().json(
                serde_json::json!({"status": "fail", "message": format!("Could not connect to RedisDB: {}", e)}),
            );
        }
    };
    let redis_result: redis::RedisResult<String> = redis_client
        .get(refresh_token_details.token_uuid.to_string())
        .await;

    let user_id = match redis_result {
        Ok(value) => value,
        Err(_) => {
            return HttpResponse::Forbidden()
                .json(serde_json::json!({"status": "fail", "message": "Token is invalid or session has expired"}));
        }
    };

    // Check that the user still exists in user collection
    let user_id_uuid = Uuid::parse_str(&user_id).unwrap();
    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);

    // Find the user in the database
    let user = match user_collection
        .find_one(doc! { "id": format!("{}", &user_id_uuid) })
        .await
        .unwrap()
    {
        Some(user) => user,
        None => return HttpResponse::BadRequest().json(
            serde_json::json!({"status": "fail", "message": "User belonging to this token no longer exists"}),
        )
    };

    // Generate a new access token and add it to database
    let access_token_details = match token::generate_jwt_token(
        uuid::Uuid::parse_str(&user.id).unwrap(),
        data.env.security.token.expiration.access_max_age,
        data.env.security.token.encryption.access_private_key.to_owned(),
    ) {
        Ok(token_details) => token_details,
        Err(e) => {
            return HttpResponse::BadGateway()
                .json(serde_json::json!({"status": "fail", "message": format_args!("{:?}", e)}));
        }
    };

    let redis_result: redis::RedisResult<()> = redis_client
        .set_ex(
            access_token_details.token_uuid.to_string(),
            user.id.to_string(),
            (data.env.security.token.expiration.access_max_age * 60) as usize,
        )
        .await;

    if redis_result.is_err() {
        return HttpResponse::UnprocessableEntity().json(
            serde_json::json!({"status": "error", "message": format_args!("{:?}", redis_result.unwrap_err())}),
        );
    }

    // Reset the access token and logged-in cookies to invalidate client-side
    let access_cookie = Cookie::build("access_token", access_token_details.token.clone().unwrap())
        .path("/")
        .max_age(ActixWebDuration::new(data.env.security.token.expiration.access_max_age * 60, 0))
        .domain(&data.env.security.cookies.domain)
        .http_only(true)
        .secure(data.env.security.cookies.secure)
        .same_site(SameSite::Strict)
        .finish();


    HttpResponse::Ok()
        .cookie(access_cookie)
        .json(serde_json::json!({
            "status": "success", "message": "Refresh was successful",
            "access_token": access_token_details.token.unwrap()
        }))

}


#[get("/auth/logout")]
async fn logout_handler(
    req: HttpRequest,
    data: web::Data<AppState>,
    auth_guard: jwt::JwtUserMiddleware,
) -> impl Responder {

    let message = "Token is invalid or session has expired";

    // Get the refresh token from the cookies
    let refresh_token = match req.cookie("refresh_token") {
        Some(c) => c.value().to_string(),
        None => {
            return HttpResponse::Forbidden()
                .json(serde_json::json!({"status": "fail", "message": message}));
        }
    };

    // Decode the refresh token
    let refresh_token_details =
        match token::verify_jwt_token(data.env.security.token.encryption.refresh_public_key.to_owned(), &refresh_token)
        {
            Ok(token_details) => token_details,
            Err(e) => {
                return HttpResponse::Forbidden().json(
                    serde_json::json!({"status": "fail", "message": format_args!("{:?}", e)}),
                );
            }
        };
    

    // Delete the refresh and access tokens (from guard) from the RedisDB
    let mut redis_client = data.auth_session.get_async_connection().await.unwrap();
    let redis_result: redis::RedisResult<usize> = redis_client
        .del(&[
            refresh_token_details.token_uuid.to_string(),
            auth_guard.access_token_uuid.to_string(),
        ])
        .await;

    if redis_result.is_err() {
        return HttpResponse::UnprocessableEntity().json(
            serde_json::json!({"status": "error", "message": format_args!("{:?}", redis_result.unwrap_err())}),
        );
    }

    // Consume the cookies for the client
    let access_cookie = Cookie::build("access_token", "")
        .path("/")
        .max_age(ActixWebDuration::new(-1, 0))
        .domain(&data.env.security.cookies.domain)
        .http_only(true)
        .secure(data.env.security.cookies.secure)
        .same_site(SameSite::Strict)
        .finish();

    let refresh_cookie = Cookie::build("refresh_token", "")
        .path("/")
        .max_age(ActixWebDuration::new(-1, 0))
        .domain(&data.env.security.cookies.domain)
        .http_only(true)
        .secure(data.env.security.cookies.secure)
        .same_site(SameSite::Strict)
        .finish();

    HttpResponse::Ok()
        .cookie(access_cookie)
        .cookie(refresh_cookie)
        .json(serde_json::json!({"status": "success"}))
}

#[derive(PartialEq, Clone)]
pub enum OneTimeTokenType {
    UserVerificationCheck,
    PasswordResetEmail,
    PasswordResetCheck,
    PasswordReset
}
impl fmt::Display for OneTimeTokenType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            OneTimeTokenType::UserVerificationCheck => write!(f, "user_verification_check"),
            OneTimeTokenType::PasswordResetEmail  => write!(f, "password_reset_email"),
            OneTimeTokenType::PasswordResetCheck  => write!(f, "password_reset_check"),
            OneTimeTokenType::PasswordReset  => write!(f, "password_reset")
        }
    }
}
impl FromStr for OneTimeTokenType {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {

        let parsed = match s {
            "user_verification_check" => OneTimeTokenType::UserVerificationCheck,
            "password_reset_email" => OneTimeTokenType::PasswordResetEmail,
            "password_reset_check" => OneTimeTokenType::PasswordResetCheck,
            "password_reset" => OneTimeTokenType::PasswordReset,
            _ => return Err(())
        };
        Ok(parsed)
    }
}

// Send a user verification email as admin
#[get("/auth/verification-email/{id}")]
async fn verification_email_handler(
    data: web::Data<AppState>,
    id: web::Path<UserId>,
    _: jwt::JwtAdminMiddleware,
) -> impl Responder {

    let user_id = &id.into_inner();

    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);

    // Find the user in the database
    let user = match user_collection
        .find_one(doc! { "id": &user_id })
        .await
    {   
        Ok(Some(user)) => user,
        Ok(None) => return HttpResponse::BadRequest().json(
            serde_json::json!({"status": "fail", "message": "Could not find user"}),
        ),
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": "Failed to find user"}),
        )
    };

    match create_one_time_token(
        &Uuid::from_str(&user.id).unwrap(),
        &data, 
        &OneTimeTokenType::UserVerificationCheck, 
        &data.env.security.token.expiration.onetime_max_age_email
    ).await
    {
        Ok(one_time_token_details) => {
            let email = Email::new(&user, &data.env);

            match email.send_verification(
                &format!("{}/verify?token={}", &data.env.security.cors.app_origin_public_url, &one_time_token_details.token.unwrap()), true  // dev version
            ).await {
                Ok(_) =>{
                    HttpResponse::Ok().json(
                        serde_json::json!({"status": "success", "message": "Sent verification email"}), 
                    )
                },
                Err(err) => {
                    HttpResponse::InternalServerError().json(
                        serde_json::json!({"status": "fail", "message": format!("Failed to send email: {}", err.to_string())}), 
                    )
                }
            }
        },
        Err(error_response) => error_response
    }
}


// Check validitiy of a user verification one-time token from the user verification page
#[post("/auth/verification-check", wrap="one_time_token_rate_limiter_factory()")]
async fn verification_check_handler(
    data: web::Data<AppState>,
    body: web::Json<AuthOneTimeCheckSchema>,
    // Unprotected route - critical to review policies! 
) -> impl Responder {

    // Helper to keep responses consistent and idempotent
    let ok_invalid = || HttpResponse::Ok().json(serde_json::json!({
        "status": "already_used",
        "message": "Token already used or invalid"
    }));

    let user_collection: Collection<User> =
        get_cerebro_db_collection(&data, AdminCollection::Users);

    // Consume/validate token; map any error to 200
    let (token_type, access_token_details, _user) =
        match check_one_time_token_and_user(&data, &user_collection, &body.access_token).await {
            Ok(t) => t,
            Err(_e) => return ok_invalid(),
        };

    // Only allow the expected token type
    if token_type != OneTimeTokenType::UserVerificationCheck {
        return ok_invalid();
    }

    // Issue the short-lived PasswordResetEmail token
    match create_one_time_token(
        &access_token_details.user_id,
        &data,
        &OneTimeTokenType::PasswordResetEmail,
        &data.env.security.token.expiration.onetime_max_age_reset,
    ).await {
        Ok(one_time_token_details) => HttpResponse::Ok().json(serde_json::json!({
            "status": "ok",
            "message": "Verification successful",
            "access_token": one_time_token_details.token.unwrap()
        })),
        Err(_e) => ok_invalid(),
    }
}


// Send a password reset email as admin
#[get("/auth/password-email/{id}")]
async fn password_email_handler(
    data: web::Data<AppState>,
    id: web::Path<UserId>,
    _: jwt::JwtAdminMiddleware,
) -> impl Responder {

    let user_id = &id.into_inner();
    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);

    // Find the user in the database
    let user = match user_collection
        .find_one(doc! { "id": &user_id })
        .await
    {   
        Ok(Some(user)) => user,
        Ok(None) => return HttpResponse::BadRequest().json(
            serde_json::json!({"status": "fail", "message": "Could not find user"}),
        ),
        Err(_) => return HttpResponse::InternalServerError().json(
            serde_json::json!({"status": "fail", "message": "Failed to find user"}),
        )
    };

    match create_one_time_token(
        &Uuid::from_str(&user.id).unwrap(),
        &data, 
        &OneTimeTokenType::PasswordResetCheck, 
        &data.env.security.token.expiration.onetime_max_age_email
    ).await
    {
        Ok(one_time_token_details) => {
            let email = Email::new(&user, &data.env);

            match email.send_password_reset(
                &format!("{}/password?token={}", &data.env.security.cors.app_origin_public_url, &one_time_token_details.token.unwrap())
            ).await {
                Ok(_) =>{
                    HttpResponse::Ok().json(
                        serde_json::json!({"status": "success", "message": "Sent password reset email"}), 
                    )
                },
                Err(err) => {
                    HttpResponse::InternalServerError().json(
                        serde_json::json!({"status": "fail", "message": format!("Failed to send email: {}", err.to_string())}), 
                    )
                }
            }
        },
        Err(error_response) => error_response
    }
}


// Send a password reset email using the user identified in a valid one-time token obtained during verification
#[post("/auth/password-reset-email", wrap="one_time_token_rate_limiter_factory()")]
async fn reset_password_email_handler(
    data: web::Data<AppState>,
    body: web::Json<AuthOneTimeCheckSchema>,
    // Unprotected route - critical to review policies! 
) -> impl Responder {

    let error_response = HttpResponse::Forbidden().json(
        serde_json::json!({"status": "error", "message": "Token is invalid or session has expired"})
    );

    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);

     // Check the access token, user in database from access token, delete the access token and return its details and type
     let (token_type, access_token_details, user) = match check_one_time_token_and_user(
        &data, &user_collection, &body.access_token
    ).await 
    {
        Ok((token_type, access_token_details, user)) => {
            (token_type, access_token_details, user)
        },
        Err(error_response) => { return error_response }
    };
    
    match token_type {
        OneTimeTokenType::PasswordResetEmail => {

            match create_one_time_token(
                &access_token_details.user_id,
                &data, 
                &OneTimeTokenType::PasswordResetCheck, 
                &data.env.security.token.expiration.onetime_max_age_password
            ).await
            {
                Ok(one_time_token_details) => {
                    // Send the verification email from the access token user 
                    let email = Email::new(&user, &data.env);
                    
                    match email.send_password_reset(
                        &format!("{}/password?token={}", &data.env.security.cors.app_origin_public_url, &one_time_token_details.token.unwrap())
                    ).await {
                        Ok(_) =>{
                            HttpResponse::Ok().json(
                                serde_json::json!({"status": "success", "message": "Sent password reset email"}), 
                            )
                        },
                        Err(_) => {
                            HttpResponse::InternalServerError().json(
                                serde_json::json!({"status": "fail", "message": "Failed to send password reset email"}), 
                            )
                        }
                    }
                },
                Err(error_response) => error_response
            }
            
        },
        _ => error_response
    }
}


// Check validitiy of the password reset one-time token
#[post("/auth/password-reset-check", wrap="one_time_token_rate_limiter_factory()")]
async fn reset_password_check_handler(
    data: web::Data<AppState>,
    body: web::Json<AuthOneTimeCheckSchema>,
    // Unprotected route - critical to review policies! 
) -> impl Responder {

    let ok_already = || HttpResponse::Ok().json(serde_json::json!({
        "status": "already_used",
        "message": "Token already used or invalid"
    }));

    let user_collection: Collection<User> =
        get_cerebro_db_collection(&data, AdminCollection::Users);

    let (token_type, access_token_details, _user) =
        match check_one_time_token_and_user(&data, &user_collection, &body.access_token).await {
            Ok(v) => v,
            Err(_e) => return ok_already(), // map all failures to 200
        };

    if token_type != OneTimeTokenType::PasswordResetCheck {
        return ok_already();
    }

    match create_one_time_token(
        &access_token_details.user_id,
        &data,
        &OneTimeTokenType::PasswordReset,
        &data.env.security.token.expiration.onetime_max_age_reset,
    ).await {
        Ok(one_time_token_details) => HttpResponse::Ok().json(serde_json::json!({
            "status": "ok",
            "message": "Password reset check passed",
            "access_token": one_time_token_details.token.unwrap()
        })),
        Err(_e) => ok_already(),
    }
}


// Reset a password from the password reset page
#[post("/auth/password-reset", wrap="one_time_token_rate_limiter_factory()")]
async fn reset_password_handler(
    data: web::Data<AppState>,
    body: web::Json<AuthResetPasswordSchema>,
    // Unprotected route - critical to review policies!
) -> impl Responder {

    let error_response: HttpResponse = HttpResponse::Forbidden().json(
        serde_json::json!({"status": "error", "message": "Token is invalid or session has expired"}),
    );

    // Database connections
    let user_collection: Collection<User> = get_cerebro_db_collection(&data, AdminCollection::Users);

    // Check the access token, user in database from access token, delete the access token and return its details and type
    let (token_type, _, user) = match check_one_time_token_and_user(
        &data, &user_collection, &body.access_token
    ).await 
    {
        Ok((token_type, access_token_details, user)) => {
            (token_type, access_token_details, user)
        },
        Err(error_response) => { return error_response }
    };

    match token_type {
        OneTimeTokenType::PasswordReset => {

            // Hash the password for storage
            let salt = SaltString::generate(&mut OsRng);
            
            let hashed_password = match Argon2::default().hash_password(
                &body.password.as_bytes(), 
                &salt
            ) {
                Ok(password_hash) => password_hash.to_string(),
                Err(_) => { return error_response }
            };

            // Also sets verified status, as we have confirmed or reset password by user email
            match user_collection.find_one_and_update(
                doc! { "id": &user.id },
                doc! { "$set": { "password": hashed_password, "verified": true } }
            ).await 
            {
                Ok(None) => return error_response,
                Ok(Some(_)) => {
                    return HttpResponse::Ok().json(serde_json::json!({
                        "status": "ok", "message": "Verified user and reset password",
                    }))
                }
                Err(_) => return error_response
            }

        },
        _ => return error_response
    }

}


// Helper function for checking access token and user in database
async fn check_one_time_token_and_user(
    data: &web::Data<AppState>, 
    user_collection: &Collection<User>, 
    access_token: &String,
) -> Result<(OneTimeTokenType, TokenDetails, User), HttpResponse>{

    let error_response: HttpResponse = HttpResponse::Forbidden().json(
        serde_json::json!({"status": "error", "message": "Token is invalid or session has expired"}),
    );

    // Create a local redis connection for this function
    let mut redis_client = match data.auth_onetime.get_async_connection().await {
        Ok(redis_client) => redis_client,
        Err(_) => { return Err(error_response) }
    };

    // Decode and verify the short-lived access token
    let access_token_details = match token::verify_jwt_token(
        data.env.security.token.encryption.access_public_key.to_owned(), &access_token
    )
    {
        Ok(token_details) => token_details,
        Err(_) => { return Err(error_response); }
    };
     
    // Get token type from one-time authentication database
    let redis_result: redis::RedisResult<String> = redis_client
        .get(access_token_details.token_uuid.to_string())
        .await;

    let token_type = match redis_result {
        Ok(value) => match OneTimeTokenType::from_str(&value) {
            Ok(token_type) => token_type,
            Err(_) => { return Err(error_response); }
    },
        Err(_) => { return Err(error_response); }
        
    };

    let user = match user_collection
        .find_one(doc! { "id": &access_token_details.user_id.to_string() })
        .await
    {
        Ok(None) => return Err(error_response),
        Ok(Some(user)) => user,
        Err(_) => return Err(error_response)
    };

    // Delete the one-time token from the database
    let redis_result: redis::RedisResult<usize> = redis_client
        .del(access_token_details.token_uuid.to_string())
        .await;

    if redis_result.is_err() {
        return Err(error_response)
    }

    Ok((token_type, access_token_details, user))
}


async fn create_one_time_token(user_id: &Uuid, data: &web::Data<AppState>, token_type: &OneTimeTokenType, ttl: &i64) -> Result<TokenDetails, HttpResponse> {

    let error_response: HttpResponse = HttpResponse::Forbidden().json(
        serde_json::json!({"status": "error", "message": "Token is invalid or session has expired"}),
    );

    let mut redis_client = match data.auth_onetime.get_async_connection().await {
        Ok(redis_client) => redis_client,
        Err(_) => { return Err(error_response) }
    };

    // Generate the new access token, this is a one-time token to email the 
    // password reset check link which redirects to the password reset with
    // new token
    let access_token_details = match token::generate_jwt_token(
        *user_id,
        *ttl,
        data.env.security.token.encryption.access_private_key.to_owned(),
    ) {
        Ok(token_details) => token_details,
        Err(_) => return  Err(error_response)
    };

    // Insert the access token (key) with token type (value) and expiration data
    let access_result: redis::RedisResult<()> = redis_client
        .set_ex(
            access_token_details.token_uuid.to_string(),
            token_type.to_string(),
            (ttl * 60) as usize,
        )
        .await;

    if let Err(_) = access_result {
        return Err(error_response);
    };

    Ok(access_token_details)

}


// Protected status route
#[get("/auth/status")]
async fn auth_status_handler(_: web::Data<AppState>, _: jwt::JwtUserMiddleware) -> HttpResponse {
    HttpResponse::Ok().into()
}

// Unprotected status route - critical to review policies! 
#[get("/status", wrap="login_rate_limiter_factory()")]
async fn app_status_handler(_: web::Data<AppState>) -> HttpResponse {
    HttpResponse::Ok().append_header(ContentType::plaintext()).body("ok\n").into()
}



// Handler configuration
pub fn auth_config(cfg: &mut web::ServiceConfig, config: &Config) {

    cfg.service(login_handler)
    .service(refresh_access_token_handler)
    .service(logout_handler)
    .service(auth_status_handler)
    .service(app_status_handler);

    if config.security.components.email {
        cfg.service(reset_password_email_handler)
        .service(verification_email_handler)
        .service(verification_check_handler)
        .service(reset_password_check_handler)
        .service(reset_password_handler)
        .service(password_email_handler);
    }
}