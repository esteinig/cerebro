// Utility function to generate a hashed password for 
// the admin user initialisation environment variable
// used in the `docker-compose` init script for MongoDB 

use argon2::{password_hash::SaltString, PasswordHasher, Argon2};
use rand_core::OsRng;
use crate::stack::deploy::StackConfigError;

pub fn hash_password(pwd: &str) -> Result<String, StackConfigError> {

    let salt = SaltString::generate(&mut OsRng);
    let hashed_password = Argon2::default()
        .hash_password(pwd.as_bytes(), &salt)
        .map_err(|_| StackConfigError::PasswordNotHashed)?
        .to_string();

    Ok(hashed_password)
}