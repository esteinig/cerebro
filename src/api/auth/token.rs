//! Helper functions for token signatures and verification

use base64::{engine::general_purpose, Engine as _};
use crate::api::auth::schema::{TokenDetails, TokenClaims};

/// This helper function uses the private keys to sign JSON Web Tokens using the RS256 algorithm.
/// It first converts the base64-encoded private key back to its original string format and then uses
/// it to sign the JWT. The function then returns the TokenDetails struct containing the metadata of 
/// the generated token. We will store the token metadata in the Redis database so that we can invalidate
/// the token when necessary.
pub fn generate_jwt_token(
    user_id: uuid::Uuid,
    ttl: i64,
    private_key: String,
) -> Result<TokenDetails, jsonwebtoken::errors::Error> {
    let bytes_private_key = general_purpose::STANDARD.decode(private_key).unwrap();
    let decoded_private_key = String::from_utf8(bytes_private_key).unwrap();

    let now = chrono::Utc::now();
    let mut token_details = TokenDetails {
        user_id,
        token_uuid: uuid::Uuid::new_v4(),
        expires_in: Some((now + chrono::Duration::minutes(ttl)).timestamp()),
        token: None,
    };

    let claims = TokenClaims {
        sub: token_details.user_id.to_string(),
        token_uuid: token_details.token_uuid.to_string(),
        exp: token_details.expires_in.unwrap(),
        iat: now.timestamp(),
        nbf: now.timestamp(),
    };

    let header = jsonwebtoken::Header::new(jsonwebtoken::Algorithm::RS256);
    let token = jsonwebtoken::encode(
        &header,
        &claims,
        &jsonwebtoken::EncodingKey::from_rsa_pem(decoded_private_key.as_bytes())?,
    )?;
    
    token_details.token = Some(token);

    Ok(token_details)
}


// This utility function allows for JWT authentication using the public keys. It first decodes the base64-encoded 
// public key to its original string format and then verifies the JWT using the key to extract the corresponding
// payload. The extracted payload is then used to construct the TokenDetails struct, which is returned by the function.
pub fn verify_jwt_token(
    public_key: String,
    token: &str,
) -> Result<TokenDetails, jsonwebtoken::errors::Error> {

    let bytes_public_key = general_purpose::STANDARD.decode(public_key).unwrap();
    let decoded_public_key = String::from_utf8(bytes_public_key).unwrap();

    let validation = jsonwebtoken::Validation::new(jsonwebtoken::Algorithm::RS256);

    let decoded = jsonwebtoken::decode::<TokenClaims>(
        token,
        &jsonwebtoken::DecodingKey::from_rsa_pem(decoded_public_key.as_bytes())?,
        &validation,
    )?;

    let user_id = uuid::Uuid::parse_str(decoded.claims.sub.as_str()).unwrap();
    let token_uuid = uuid::Uuid::parse_str(decoded.claims.token_uuid.as_str()).unwrap();

    Ok(TokenDetails {
        token: None,
        token_uuid,
        user_id,
        expires_in: None,
    })
}
