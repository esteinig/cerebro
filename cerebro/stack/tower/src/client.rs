use std::path::{Path, PathBuf};
use std::fs::create_dir_all;
use cerebro_model::api::auth::response::AuthLoginResponseSuccess;
use cerebro_model::api::stage::model::StagedSample;
use cerebro_model::api::stage::response::{DeleteStagedSampleResponse, ListStagedSamplesResponse};
use cerebro_model::api::utils::ErrorResponse;
use reqwest::{Client, Response, StatusCode};
use reqwest::header::AUTHORIZATION;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};

use crate::error::TowerError;



#[derive(Debug, Clone)]
pub enum Route {
    ServerStatus,
    AuthServerStatus,
    AuthLoginUser,
    AuthRefreshToken,
    DataUserSelf,
    TeamStagedSamplesRegister,
    TeamStagedSamplesList,
    TeamStagedSamplesDelete,
    TeamStagedSamplesPull,
}

impl Route {
    fn path(&self) -> &str {
        match self {
            Route::ServerStatus => "status",
            Route::AuthServerStatus => "auth/status",
            Route::AuthLoginUser => "auth/login",
            Route::AuthRefreshToken => "auth/refresh",
            Route::DataUserSelf => "users/self",
            Route::TeamStagedSamplesRegister => "stage/register",
            Route::TeamStagedSamplesList => "stage",
            Route::TeamStagedSamplesDelete => "stage",
            Route::TeamStagedSamplesPull => "stage"
        }
    }
}

#[derive(Debug, Clone)]
pub struct CerebroRoutes {
    base_url: String,
}

impl CerebroRoutes {
    pub fn new(base_url: &str) -> Self {
        Self {
            base_url: base_url.trim_end_matches('/').to_string(),
        }
    }

    pub fn url(&self, route: Route) -> String {
        format!("{}/{}", self.base_url, route.path())
    }
}

impl std::fmt::Display for Route {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.path())
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AuthLoginTokenFile {
    pub access_token: String,
    pub refresh_token: String,
}

impl From<&AuthLoginResponseSuccess> for AuthLoginTokenFile {
    fn from(response: &AuthLoginResponseSuccess) -> Self {
        Self {
            access_token: response.access_token.clone(),
            refresh_token: response.refresh_token.clone(),
        }
    }
}



#[derive(Debug, Clone)]
pub struct TowerClient {
    pub url: String,
    pub client: Client,
    pub team: Option<String>,
    pub db: Option<String>,
    pub project: Option<String>,
    pub routes: CerebroRoutes,
    pub token_data: Option<AuthLoginTokenFile>,
    pub token_file: Option<PathBuf>,
}

impl TowerClient {
    pub fn new(
        url: &str,
        token: Option<String>,
        danger_accept_invalid_certs: bool,
        token_file: Option<PathBuf>,
        team: Option<String>,
        db: Option<String>,
        project: Option<String>,
    ) -> Result<Self, TowerError> {
        let url_clean = url.trim_end_matches('/').to_string();

        let token_data = if token.is_some() {
            Some(AuthLoginTokenFile {
                access_token: token.clone().unwrap(),
                refresh_token: String::new(),
            })
        } else {
            Self::load_token(token_file.clone())?
        };

        let client = Client::builder()
            .danger_accept_invalid_certs(danger_accept_invalid_certs)
            .build()?;

        Ok(Self {
            url: url_clean.clone(),
            routes: CerebroRoutes::new(&url_clean),
            token_data,
            token_file,
            client,
            team,
            db,
            project,
        })
    }

    fn load_token(token_file: Option<PathBuf>) -> Result<Option<AuthLoginTokenFile>, TowerError> {
        if let Some(token_path) = token_file {
            if token_path.exists() {
                let token_str = std::fs::read_to_string(&token_path)?;
                let token_data: AuthLoginTokenFile = serde_json::from_str(&token_str)?;
                return Ok(Some(token_data));
            } else {
                log::error!("Token file path does not exist");
                return Err(TowerError::TokenFileNotFound);
            }
        }

        log::error!("Failed to obtain token from global argument or environment variable");
        Err(TowerError::TokenMissing)
    }

    fn get_bearer_token(&self, token: Option<String>) -> String {
        match token {
            Some(token) => format!("Bearer {}", token),
            None => {
                match &self.token_data {
                    Some(data) => format!("Bearer {}", data.access_token),
                    None => panic!("This function must only be called on a non-login route where a token is supplied!"),
                }
            }
        }
    }
    pub async fn ping_servers(&self) -> Result<(), TowerError> {
        let response = self.client
            .get(self.routes.url(Route::AuthServerStatus))
            .header(AUTHORIZATION, self.get_bearer_token(None))
            .send()
            .await?;

        if response.status().is_success() {
            log::info!("Cerebro API status: ok");
            Ok(())
        } else {
            Err(TowerError::ResponseFailure(response.status()))
        }
    }
    async fn send_request_with_team(&self, request: reqwest::RequestBuilder) -> Result<Response, TowerError> {
        let team = self.team.as_deref().ok_or(TowerError::RequireTeamNotConfigured)?;

        let response = request
            .query(&[("team", team)])
            .header(AUTHORIZATION, self.get_bearer_token(None))
            .send()
            .await?;

        Ok(response)
    }

    async fn build_request_url<T, R>(&self, route: R, params: &[(&str, T)]) -> String
    where
        T: Into<Option<String>> + Clone,
        R: Into<String> + Clone
    {
        let mut url = route.into().clone();
        let query_string: Vec<String> = params
            .iter()
            .filter_map(|(key, value)| {
                value.clone().into().map(|v| format!("{}={}", key, v))
            })
            .collect();

        if !query_string.is_empty() {
            url.push_str("?");
            url.push_str(&query_string.join("&"));
        }

        url
    }

    async fn handle_response<T: DeserializeOwned>(&self, response: Response, failure_msg: Option<&str>) -> Result<T, TowerError> {
        let status = response.status();
        if status.is_success() {
            response.json::<T>().await.map_err(|_| {
                TowerError::DataResponseFailure(status, String::from("failed to parse response data"))
            })
        } else {
            let error_response: ErrorResponse = response.json().await.map_err(|_| {
                TowerError::DataResponseFailure(status, String::from("failed to parse error response"))
            })?;
            if let Some(msg) = failure_msg {
                log::error!("{}: {}", msg, error_response.message);
            }
            Err(TowerError::ResponseFailure(status))
        }
    }

    pub async fn list_staged_samples(&self, id: &str) -> Result<Vec<StagedSample>, TowerError> {
        let url = self.routes.url(Route::TeamStagedSamplesList);
        let request = self.client.get(format!("{}/{}", url, id));

        let response = self.send_request_with_team(request).await?;
        self.handle_response::<ListStagedSamplesResponse>(response, None)
            .await?
            .data
            .ok_or_else(|| {
                TowerError::DataResponseFailure(
                    StatusCode::INTERNAL_SERVER_ERROR,
                    String::from("No staged sample data returned"),
                )
            })
    }

    pub async fn pull_staged_samples(
        &self,
        id: &str,
        delete: bool,
    ) -> Result<Vec<StagedSample>, TowerError> {

        let staged_samples = match self.list_staged_samples(id).await {
            Ok(samples) => samples,
            Err(err) => match err {
                TowerError::ResponseFailure(code) if code == StatusCode::NOT_FOUND => Vec::new(),
                _ => return Err(err),
            },
        };

        if delete {
            for staged_sample in &staged_samples {
                self.delete_staged_sample(id, Some(staged_sample.id.to_string())).await?;
            }
        }

        Ok(staged_samples)
    }

    pub async fn delete_staged_sample(
        &self,
        id: &str,
        staged_id: Option<String>,
    ) -> Result<Option<StagedSample>, TowerError> {

        let url = self.build_request_url(
            format!("{}/{id}", self.routes.url(Route::TeamStagedSamplesDelete)), 
            &[("stage_id", staged_id)]
        ).await;
        
        let request = self.client.delete(url);
        let response = self.send_request_with_team(request).await?;

        Ok(
            self.handle_response::<DeleteStagedSampleResponse>(
                response, 
                Some("Staged sample deletion failed")
            )
            .await?
            .data
        )
    }
}
