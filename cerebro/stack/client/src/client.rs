use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::create_dir_all;
use std::fs::File;
use std::io::BufWriter;
use std::time::Duration;
use std::io::Write;
use anyhow::Result;
use cerebro_model::api::cerebro::response::CerebroIdentifierResponse;
use cerebro_model::api::cerebro::response::CerebroIdentifierSummary;
use cerebro_model::api::cerebro::response::ContaminationTaxaResponse;
use cerebro_model::api::cerebro::response::FilteredTaxaResponse;
use cerebro_model::api::cerebro::response::PathogenDetectionTableResponse;
use cerebro_model::api::cerebro::response::QualityControlTableResponse;
use cerebro_model::api::cerebro::response::RetrieveModelResponse;
use cerebro_model::api::cerebro::response::TaxonHistoryResponse;
use cerebro_model::api::cerebro::schema::CerebroIdentifierSchema;
use cerebro_model::api::cerebro::schema::ContaminationSchema;
use cerebro_model::api::cerebro::schema::PathogenDetectionTableSchema;
use cerebro_model::api::cerebro::schema::QualityControlTableSchema;
use cerebro_model::api::cerebro::schema::PrevalenceContaminationConfig;
use cerebro_model::api::cerebro::schema::UpdateRunConfigSchema;
use cerebro_model::api::files::model::SeaweedFile;
use cerebro_model::api::files::model::SeaweedFileId;
use cerebro_model::api::files::response::DeleteFileResponse;
use cerebro_model::api::files::response::DeleteFilesResponse;
use cerebro_model::api::files::response::ListFilesResponse;
use cerebro_model::api::jobs::response::CreateScheduleResponse;
use cerebro_model::api::jobs::response::EnqueueJobResponse;
use cerebro_model::api::jobs::response::JobCompletionResponse;
use cerebro_model::api::jobs::response::JobsStatusResponse;
use cerebro_model::api::jobs::schema::EnqueueJobRequest;
use cerebro_model::api::jobs::schema::ScheduleJobRequest;
use cerebro_model::api::teams::schema::RegisterDatabaseSchema;
use cerebro_model::api::towers::model::ProductionTower;
use cerebro_model::api::towers::response::DeleteTowerResponse;
use cerebro_model::api::towers::response::ListTowersResponse;
use cerebro_model::api::towers::response::PingTowerResponse;
use cerebro_model::api::towers::response::RegisterTowerResponse;
use cerebro_model::api::towers::schema::RegisterTowerSchema;
use cerebro_model::api::stage::model::StagedSample;
use cerebro_model::api::stage::response::DeleteStagedSampleResponse;
use cerebro_model::api::stage::response::ListStagedSamplesResponse;
use cerebro_model::api::stage::schema::RegisterStagedSampleSchema;
use cerebro_model::api::training::response::TrainingPrefetchData;
use cerebro_model::api::training::response::TrainingResponse;
use cerebro_model::api::training::response::TrainingSessionData;
use cerebro_model::api::training::schema::CreateTrainingPrefetch;
use cerebro_model::api::watchers::model::ProductionWatcher;
use cerebro_model::api::watchers::response::DeleteWatcherResponse;
use cerebro_model::api::watchers::response::ListWatchersResponse;
use cerebro_model::api::watchers::response::PingWatcherResponse;
use cerebro_model::api::watchers::response::RegisterWatcherResponse;
use cerebro_model::api::watchers::schema::RegisterWatcherSchema;
use cerebro_pipeline::modules::pathogen::PathogenDetectionTableRecord;
use cerebro_pipeline::modules::quality::ReadQualityControl;
use cerebro_pipeline::taxa::filter::TaxonFilterConfig;
use cerebro_pipeline::taxa::taxon::Taxon;
use cerebro_pipeline::utils::read_tsv;
use cerebro_pipeline::utils::write_tsv;
use chrono::Utc;
use reqwest::blocking::RequestBuilder;
use reqwest::blocking::Response;
use reqwest::StatusCode;
use serde::de::DeserializeOwned;
use std::path::PathBuf;
use itertools::Itertools;
use reqwest::blocking::Client;
use rpassword::prompt_password;
use reqwest::header::AUTHORIZATION;
use serde::{Serialize, Deserialize};
use actix_web_httpauth::headers::authorization::Bearer;

use cerebro_model::api::utils::ErrorResponse;
use cerebro_model::api::cerebro::model::Cerebro;
use cerebro_model::api::auth::schema::AuthLoginSchema;
use cerebro_model::api::users::response::UserSelfResponse;
use cerebro_model::api::auth::response::AuthLoginResponseSuccess;
use cerebro_model::api::teams::schema::RegisterProjectSchema;
use cerebro_model::api::files::schema::RegisterFileSchema;

use crate::error::HttpClientError;
use crate::regression::RpmAnalysisResult;
use crate::regression::RpmAnalyzer;
use crate::regression::RpmConfigBuilder;
use crate::utils::DiagnosticResult;
use std::fmt;



#[derive(Debug, Clone)]
pub enum Route {
    ServerStatus,
    AuthServerStatus,
    AuthLoginUser,
    AuthRefreshToken,
    DataUserSelf,
    DataUserSelfTeams,
    DataCerebroInsertModel,
    DataCerebroRetrieveModel,
    DataCerebroUpdateModelRunConfig,
    DataCerebroIdentifiers,
    DataCerebroQualityControl,
    DataCerebroPathogenDetection,
    DataCerebroTaxaSummary,
    DataCerebroTaxaHistory,
    DataCerebroTaxaFiltered,
    DataCerebroTaxaContamination,
    TeamProjectCreate,
    TeamDatabaseCreate,
    TeamFilesRegister,
    TeamFilesList,
    TeamFilesDelete,
    TeamTowersRegister,
    TeamTowersList,
    TeamTowersDelete,
    TeamTowersPing,
    TeamWatchersRegister,
    TeamWatchersList,
    TeamWatchersDelete,
    TeamWatchersPing,
    TeamStagedSamplesRegister,
    TeamStagedSamplesList,
    TeamStagedSamplesDelete,
    TeamStagedSamplesPull,
    JobsEnqueue,
    JobsEnqueueStatus,
    JobsSchedule,
    JobsScheduleStatus,
    TrainingRegister,
    TrainingList,
    TrainingDelete,
    TrainingSessionList,
    TrainingSessionDelete
}

impl Route {
    fn path(&self) -> &str {
        match self {
            Route::ServerStatus => "status",
            Route::AuthServerStatus => "auth/status",
            Route::AuthLoginUser => "auth/login",
            Route::AuthRefreshToken => "auth/refresh",
            Route::DataUserSelf => "users/self",
            Route::DataUserSelfTeams => "users/self/teams",
            Route::DataCerebroInsertModel => "cerebro",
            Route::DataCerebroRetrieveModel => "cerebro",
            Route::DataCerebroUpdateModelRunConfig => "cerebro/run",
            Route::DataCerebroIdentifiers => "cerebro/ids",
            Route::DataCerebroQualityControl => "cerebro/table/qc",
            Route::DataCerebroPathogenDetection => "cerebro/table/pathogen",
            Route::DataCerebroTaxaSummary => "cerebro/taxa/summary",
            Route::DataCerebroTaxaHistory => "cerebro/taxa/history",
            Route::DataCerebroTaxaContamination => "cerebro/taxa/contamination",
            Route::DataCerebroTaxaFiltered => "cerebro/taxa",
            Route::TeamDatabaseCreate => "teams/database",
            Route::TeamProjectCreate => "teams/project",
            Route::TeamFilesRegister => "files/register",
            Route::TeamFilesList => "files",
            Route::TeamFilesDelete => "files",
            Route::TeamTowersRegister => "tower/register",
            Route::TeamTowersList => "tower",
            Route::TeamTowersDelete => "tower",
            Route::TeamTowersPing => "tower",
            Route::TeamWatchersRegister => "watcher/register",
            Route::TeamWatchersList => "watcher",
            Route::TeamWatchersDelete => "watcher",
            Route::TeamWatchersPing => "watcher",
            Route::TeamStagedSamplesRegister => "stage/register",
            Route::TeamStagedSamplesList => "stage",
            Route::TeamStagedSamplesDelete => "stage",
            Route::TeamStagedSamplesPull => "stage",
            Route::JobsEnqueue => "jobs/enqueue",
            Route::JobsEnqueueStatus => "jobs/enqueue",
            Route::JobsSchedule => "jobs/schedule",
            Route::JobsScheduleStatus => "jobs/schedule/status",
            Route::TrainingRegister => "training/prefetch",
            Route::TrainingList => "training/prefetch",
            Route::TrainingDelete => "training/prefetch",
            Route::TrainingSessionList => "training/sessions",
            Route::TrainingSessionDelete => "training/session",
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

// Example usage:
impl fmt::Display for Route {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
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
pub struct CerebroClient {
    pub url: String,
    pub client: Client,
    pub team: Option<String>,
    pub db: Option<String>,
    pub project: Option<String>,
    pub routes: CerebroRoutes,
    pub token_data: Option<AuthLoginTokenFile>,
    pub token_file: Option<PathBuf>,
    pub max_model_size: Option<f64>
}

impl CerebroClient {
    pub fn new(
        url: &str,
        token: Option<String>,
        login: bool,
        danger_accept_invalid_certs: bool,
        token_file: Option<PathBuf>,
        team: Option<String>,
        db: Option<String>,
        project: Option<String>,
    ) -> Result<Self, HttpClientError> {
        let url_clean = url.trim_end_matches('/').to_string();

        let token_data = if login {
            None
        } else {
            Self::load_token(token, token_file.clone())?
        };

        let client = Client::builder()
            .danger_accept_invalid_certs(danger_accept_invalid_certs)
            .timeout(Duration::from_secs(600))
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
            max_model_size: Some(100.0)
        })
    }
    pub fn log_team_warning(&self) {
        if let None = self.team {
            log::warn!("No team configured for authentication - you may need this for data access");
            log::warn!("Use team name or unique identifier with global argument '--team' or with environment variable '$CEREBRO_USER_TEAM'")
        }
    }
    pub fn log_db_warning(&self) {
        if let None = self.db {
            log::warn!("No database configured for authentication - you may need this for data access");
            log::warn!("Use database name or unique identifier with global argument '--db' or with environment variable '$CEREBRO_USER_DB'")
        }
    }
    pub fn log_project_warning(&self) {
        if let None = self.project {
            log::warn!("No database project configured for authentication - you may need this for data access");
            log::warn!("Use project name or unique identifier with global argument '--project' or with environment variable '$CEREBRO_USER_PROJECT'")
        }
    }
    pub fn set_data_auth(&mut self, team: &str, database: &str, project: &str) {
        self.team = Some(team.to_string());
        self.db = Some(database.to_string());
        self.project = Some(project.to_string());
    }
    fn load_token(
        token: Option<String>,
        token_file: Option<PathBuf>,
    ) -> Result<Option<AuthLoginTokenFile>, HttpClientError> {
        if let Some(token) = token {
            return Ok(Some(AuthLoginTokenFile {
                access_token: token,
                refresh_token: String::new(),
            }));
        }

        if let Some(token_path) = token_file {
            if token_path.exists() {
                let token_str = std::fs::read_to_string(&token_path)
                    .map_err(HttpClientError::IOFailure)?;
                let token_data: AuthLoginTokenFile = serde_json::from_str(&token_str)?;
                return Ok(Some(token_data));
            } else {
                log::error!("Token file path does not exist");
                std::process::exit(1);
            }
        }

        log::error!("Failed to obtain token from global argument (--token | --token-file) or environment variable ($CEREBRO_API_TOKEN)");
        std::process::exit(1);
    }

    // We use authorization headers for the client instead of cookies
    // as we set strict cookie policies that force same-site origin
    // requests that cannot be fulfilled by the client
    pub fn get_bearer_token(&self, token: Option<String>) -> String {
        match token {
            Some(token) => Bearer::new(token).to_string(),
            None => {
                match &self.token_data {
                    Some(data) => Bearer::new(data.access_token.clone()).to_string(),
                    None => panic!("This function must only be called on a non-login route where a token is supplied!")
                }
            }
        }
    }

    fn send_request(&self, route: Route, auth: bool) -> Result<Response, HttpClientError> {
        let mut request = self.client.get(self.routes.url(route));
        if auth {
            request = request.header(AUTHORIZATION, self.get_bearer_token(None));
        }
        request.send().map_err(|e| HttpClientError::ReqwestFailure(e))
    }

    fn check_response_status(&self, response: &Response) -> Result<(), HttpClientError> {
        if response.status().is_success() {
            Ok(())
        } else {
            Err(HttpClientError::PingServer(
                response.status(),
                format!("Server error: {}", response.status()),
            ))
        }
    }

    // Ping server status
    pub fn ping_status(&self) -> Result<(), HttpClientError> {
        let response = self.send_request(Route::ServerStatus, false)?;
        self.check_response_status(&response)?;
        log::info!("Cerebro API status: ok");
        Ok(())
    }

    // Ping server status routes - currently running a single server.
    pub fn ping_servers(&self) -> Result<(), HttpClientError> {
        let response = self.send_request(Route::AuthServerStatus, true)?;
        self.check_response_status(&response)?;
        log::info!("Cerebro API status: ok");
        Ok(())
    }
    // Login user and save tokens
    pub fn login_user(&self, email: &str, password: Option<String>) -> Result<(), HttpClientError> {
        let password = password.unwrap_or_else(|| {
            prompt_password(format!("Password [{}]: ", email)).expect("Password input failed")
        });

        let login_schema = AuthLoginSchema {
            email: email.to_owned(),
            password,
        };

        let response = self.client
            .post(self.routes.url(Route::AuthLoginUser))
            .json(&login_schema)
            .send()?;

        let status = response.status();

        if status.is_success() {
            let login_response: AuthLoginResponseSuccess = response.json()?;
            self.handle_login_success(&login_response)?;
        } else {
            let error_response: ErrorResponse = response.json()?;
            log::error!("{}", error_response.message);
            return Err(HttpClientError::ResponseFailure(status));
        }

        Ok(())
    }

    fn handle_login_success(
        &self,
        login_response: &AuthLoginResponseSuccess,
    ) -> Result<(), HttpClientError> {
        let access_token = Some(login_response.access_token.clone());
        let response = self.client
            .get(self.routes.url(Route::DataUserSelf))
            .header(AUTHORIZATION, self.get_bearer_token(access_token))
            .send()?;

        let status = response.status();

        let user_response: Result<UserSelfResponse, HttpClientError> = response.json().map_err(|_| {
            HttpClientError::ResponseFailure(status)
        });

        if status.is_success() {
            let user_response = user_response?;
            log::info!("Login successful. Welcome back, {}!", user_response.data.user.name);
        } else {
            let user_error_response = user_response.map_err(|_| {
                HttpClientError::DataResponseFailure(
                    status,
                    String::from("failed to parse error response"),
                )
            })?;
            log::error!("{}", user_error_response.message);
            return Err(HttpClientError::ResponseFailure(status));
        }

        self.save_token_to_file(login_response)?;
        Ok(())
    }

    fn save_token_to_file(&self, login_response: &AuthLoginResponseSuccess) -> Result<(), HttpClientError> {
        if let Some(path) = &self.token_file {
            log::info!("Token will be written to file: {}", path.display());
            std::fs::write(
                path,
                serde_json::to_string_pretty(&AuthLoginTokenFile::from(login_response))
                    .expect("Failed to convert token data to string"),
            )
            .map_err(HttpClientError::IOFailure)?;
        } else {
            println!("{}", login_response.access_token);
        }
        Ok(())
    }
    fn send_request_with_team(&self, request: RequestBuilder) -> Result<Response, HttpClientError> {

        self.log_team_warning();

        let team = self.team.as_deref().ok_or(HttpClientError::RequireTeamNotConfigured)?;

        let response = request
            .query(&[("team", team)])
            .header(AUTHORIZATION, self.get_bearer_token(None)).send()?;
        
        
        Ok(response)
    }
    fn send_request_with_team_db(&self, request: RequestBuilder) -> Result<Response, HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();

        let team = self.team.as_deref().ok_or(HttpClientError::RequireTeamNotConfigured)?;
        let db = self.db.as_deref().ok_or(HttpClientError::RequireDbNotConfigured)?;

        let response = request
            .query(&[("team", team)])
            .query(&[("db", db)])
            .header(AUTHORIZATION, self.get_bearer_token(None)).send()?;
        
        Ok(response)
    }
    fn send_request_with_team_db_project(&self, request: RequestBuilder) -> Result<Response, HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();
        self.log_project_warning();

        let team = self.team.as_deref().ok_or(HttpClientError::RequireTeamNotConfigured)?;
        let db = self.db.as_deref().ok_or(HttpClientError::RequireDbNotConfigured)?;
        let project = self.project.as_deref().ok_or(HttpClientError::RequireProjectNotConfigured)?;

        let response = request
            .query(&[("team", team)])
            .query(&[("db", db)])
            .query(&[("project", project)])
            .header(AUTHORIZATION, self.get_bearer_token(None)).send()?;
        
        Ok(response)
    }

    
    pub fn handle_response<S: DeserializeOwned>(
        &self,
        response: Response,
        success_msg: Option<&str>,
        failure_msg: &str,
    ) -> Result<S, HttpClientError> {
        let status = response.status();

        if status.is_success() {
            if let Some(msg) = success_msg {
                log::info!("{}", msg);
            }
            return response.json().map_err(|_| {
                HttpClientError::DataResponseFailure(status, "failed to parse response data".into())
            });
        }

        // Read body once, then try JSON -> fallback to text
        let bytes = response.bytes().unwrap_or_default();

        if let Ok(err) = serde_json::from_slice::<ErrorResponse>(&bytes) {
            if let Some(data) = &err.data {
                log::error!("{}: {} | data: {:?}", failure_msg, err.message, data);
            } else {
                log::error!("{}: {}", failure_msg, err.message);
            }
        } else {
            let text = String::from_utf8_lossy(&bytes);
            log::error!("{}: {}", failure_msg, text);
        }

        Err(HttpClientError::ResponseFailure(status))
    }

    fn build_request_url<T, R>(&self, route: R, params: &[(&str, T)]) -> String
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
    pub fn create_project(
        &self,
        name: &str,
        description: &str,
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();

        let project_schema = RegisterProjectSchema {
            project_name: name.to_owned(),
            project_description: description.to_owned()
        };

        let response = self.send_request_with_team_db(
            self.client.post(self.routes.url(Route::TeamProjectCreate)).json(&project_schema)
        )?;

        self.handle_response::<serde_json::Value>(
            response,
            Some(&format!("Project `{}` created successfully", name)),
            "Project creation failed",
        )?;
        Ok(())
    }
    // Schedule status summary (GET /jobs/schedule/status)
    // Returns a single pretty line-per-item string for "next" and "previous".
    pub fn schedule_status(&self) -> Result<String, HttpClientError> {
        let response = self.client
            .get(self.routes.url(Route::JobsScheduleStatus))
            .header(AUTHORIZATION, self.get_bearer_token(None))
            .send()?;

        // Parse standard success / error
        let status = response.status();
        if !status.is_success() {
            let err: ErrorResponse = response.json().map_err(|_| {
                HttpClientError::DataResponseFailure(
                    status,
                    String::from("failed to parse error response"),
                )
            })?;
            return Err(HttpClientError::ResponseFailureWithMessage(status, err.message));
        }

        let res: JobsStatusResponse = response.json().map_err(|_| {
            HttpClientError::DataResponseFailure(
                status,
                String::from("failed to parse schedule status response"),
            )
        })?;

        // Format summary
        if let Some(data) = res.data {
            let mut out = String::new();

            if !data.next.is_empty() {
                out.push_str("Next:\n");
                for j in &data.next {
                    // ScheduleJobSummary typically has: id, kind, queue, run_at, last_run_at (optional)
                    out.push_str(&format!(
                        "  {}  kind={}  queue={}\n",
                        j.run_at, j.kind, j.queue
                    ));
                }
            }

            if !data.previous.is_empty() {
                if !out.is_empty() { out.push('\n'); }
                out.push_str("Previous:\n");
                for j in &data.previous {
                    let last = j.last_run_at
                        .as_ref()
                        .map(|d| d.to_string())
                        .unwrap_or_else(|| "-".into());
                    out.push_str(&format!(
                        "  last_run_at={}  kind={}  queue={}\n",
                        last, j.kind, j.queue
                    ));
                }
            }

            if out.is_empty() {
                Ok(String::from("No scheduled jobs found"))
            } else {
                Ok(out.trim_end().to_string())
            }
        } else {
            Ok(String::from("No scheduled jobs found"))
        }
    }

    // Job status summary (GET /jobs/enqueue/{id})
    // Returns a compact one-line summary with completed/status and optional result/error.
    pub fn job_status(&self, id: &str) -> Result<String, HttpClientError> {
        let url = format!("{}/{}", self.routes.url(Route::JobsEnqueueStatus), id);
        let response = self.client
            .get(&url)
            .header(AUTHORIZATION, self.get_bearer_token(None))
            .send()?;

        let status = response.status();
        if !status.is_success() {
            let err: ErrorResponse = response.json().map_err(|_| {
                HttpClientError::DataResponseFailure(
                    status,
                    String::from("failed to parse error response"),
                )
            })?;
            return Err(HttpClientError::ResponseFailureWithMessage(status, err.message));
        }

        let res: JobCompletionResponse = response.json().map_err(|_| {
            HttpClientError::DataResponseFailure(
                status,
                String::from("failed to parse job completion response"),
            )
        })?;

        let mut line = format!("completed={} status={}", res.completed, res.status.unwrap_or_else(|| "unknown".into()));

        if let Some(err) = res.error {
            line.push_str(&format!(" error={}", err));
        } else if let Some(val) = res.result {
            line.push_str(&format!(" result={}", val));
        }

        Ok(line)
    }

    // Job complete? (GET /jobs/enqueue/{id})
    // Returns "yes" or "no"
    pub fn job_complete(&self, id: &str) -> Result<String, HttpClientError> {
        let url = format!("{}/{}", self.routes.url(Route::JobsEnqueueStatus), id);
        let response = self.client
            .get(&url)
            .header(AUTHORIZATION, self.get_bearer_token(None))
            .send()?;

        let status = response.status();
        if !status.is_success() {
            let err: ErrorResponse = response.json().map_err(|_| {
                HttpClientError::DataResponseFailure(
                    status,
                    String::from("failed to parse error response"),
                )
            })?;
            return Err(HttpClientError::ResponseFailureWithMessage(status, err.message));
        }

        let res: JobCompletionResponse = response.json().map_err(|_| {
            HttpClientError::DataResponseFailure(
                status,
                String::from("failed to parse job completion response"),
            )
        })?;

        Ok(if res.completed { "yes".into() } else { "no".into() })
    }
    /// Fire-and-forget enqueue to Faktory via API
    pub fn enqueue_job(
        &self,
        req: &EnqueueJobRequest,
    ) -> Result<Option<String>, HttpClientError> {
        let response = self.client
            .post(self.routes.url(Route::JobsEnqueue))
            .header(AUTHORIZATION, self.get_bearer_token(None))
            .json(req)
            .send()?;

        let res = self.handle_response::<EnqueueJobResponse>(
            response,
            Some("Job enqueued"),
            "Job enqueue failed",
        )?;

        Ok(res.data)
    }

    /// Create a schedule (one-shot or recurring)
    pub fn schedule_job(
        &self,
        req: &ScheduleJobRequest,
    ) -> Result<String, HttpClientError> {
        
        let response = self.client
            .post(self.routes.url(Route::JobsSchedule))
            .header(AUTHORIZATION, self.get_bearer_token(None))
            .json(req)
            .send()?;

        let res = self.handle_response::<CreateScheduleResponse>(
            response,
            Some("Schedule created"),
            "Schedule creation failed",
        )?;

        res.data.ok_or_else(|| {
            HttpClientError::DataResponseFailure(
                reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                String::from("Server did not return schedule id"),
            )
        })
    }
    pub fn create_database(
        &self,
        name: &str,
        description: &str,
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();

        let database_schema = RegisterDatabaseSchema {
            database_name: name.to_owned(),
            database_description: description.to_owned()
        };

        let response = self.send_request_with_team(
            self.client.post(self.routes.url(Route::TeamDatabaseCreate)).json(&database_schema)
        )?;

        self.handle_response::<serde_json::Value>(
            response,
            Some(&format!("Database `{}` created successfully", name)),
            "Database creation failed",
        )?;
        Ok(())
    }
    
    pub fn register_file(
        &self,
        register_file_schema: RegisterFileSchema,
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();

        let response = self.send_request_with_team(
            self.client
                .post(self.routes.url(Route::TeamFilesRegister))
                .json(&register_file_schema),
        )?;

        self.handle_response::<serde_json::Value>(
            response,
            Some("File registered successfully"),
            "File registration failed",
        )?;
        Ok(())
    }

    pub fn delete_file(&self, id: &str) -> Result<SeaweedFile, HttpClientError> {

        self.log_team_warning();

        let url = format!("{}/{id}", self.routes.url(Route::TeamFilesDelete));

        let response = self.send_request_with_team(
            self.client.delete(url)
        )?;

        self.handle_response::<DeleteFileResponse>(
            response,
            None,
            "File deletion failed",
        )?
        .data
        .ok_or_else(|| {
            HttpClientError::DataResponseFailure(
                reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                String::from("File data not returned after deletion"),
            )
        })
    }

    pub fn delete_files(&self, run_id: Option<String>, sample_id: Option<String>, all: Option<bool>) -> Result<Vec<SeaweedFileId>, HttpClientError> {

        self.log_team_warning();

        let url = self.build_request_url(
            self.routes.url(Route::TeamFilesDelete), &[("run_id", run_id), ("sample_id", sample_id), ("all", all.map(|b| b.to_string()))]
        );


        let response = self.send_request_with_team(
            self.client.delete(url)
        )?;

        self.handle_response::<DeleteFilesResponse>(
            response,
            None,
            "Files deletion failed",
        )?
        .data
        .ok_or_else(|| {
            HttpClientError::DataResponseFailure(
                reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                String::from("File identifiers not returned after deletion"),
            )
        })
    }

    pub fn list_files(
        &self,
        run_id: Option<String>,
        watcher_id: Option<String>,
        page: u32,
        limit: u32,
        print: bool,
    ) -> Result<Vec<SeaweedFile>, HttpClientError> {

        self.log_team_warning();

        let url = self.build_request_url(
            self.routes.url(Route::TeamFilesList),
            &[
                ("run_id", run_id),
                ("watcher_id", watcher_id),
                ("page", Some(page.to_string())),
                ("limit", Some(limit.to_string())),
            ],
        );

        let response = self.send_request_with_team(
            self.client
                .get(&url)
        )?;

        let files = self
            .handle_response::<ListFilesResponse>(
                response,
                None,
                "Failed to retrieve files",
            )?
            .data
            .ok_or_else(|| {
                HttpClientError::DataResponseFailure(
                    reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                    String::from("No file data returned"),
                )
            })?;

        if print {
            for file in &files {
                let (watcher_name, watcher_location) = file
                    .watcher
                    .as_ref()
                    .map(|w| (w.name.as_str(), w.location.as_str()))
                    .unwrap_or(("none", "none"));

                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.0} MB\t{}",
                    file.id,
                    file.date,
                    file.run_id.as_deref().unwrap_or(""),
                    file.sample_id.as_deref().unwrap_or(""),
                    watcher_name,
                    watcher_location,
                    file.fid,
                    file.size_mb(),
                    file.name
                );
            }
        }

        Ok(files)
    }

    pub fn register_tower(
        &self,
        schema: &RegisterTowerSchema,
        print: bool,
    ) -> Result<String, HttpClientError> {

        self.log_team_warning();

        let response = self.send_request_with_team(
            self.client
                .post(self.routes.url(Route::TeamTowersRegister))
                .json(schema),
        )?;

        let tower_id = self.handle_response::<RegisterTowerResponse>(
            response,
            Some("Tower registered successfully"),
            "Tower registration failed",
        )?
        .data
        .unwrap_or(schema.id.clone());

        if print {
            println!("{}", tower_id);
        }

        Ok(tower_id)
    }

    pub fn list_towers(
        &self,
        id: Option<String>,
        print: bool,
    ) -> Result<Vec<ProductionTower>, HttpClientError> {
        
        self.log_team_warning();

        let url = self.build_request_url(
            self.routes.url(Route::TeamTowersList), 
            &[("id", id)]
        );

        let response = self.send_request_with_team(
            self.client
                .get(&url)
        )?;

        let towers = self
            .handle_response::<ListTowersResponse>(
                response,
                None,
                "Failed to retrieve towers",
            )?
            .data
            .ok_or_else(|| {
                HttpClientError::DataResponseFailure(
                    reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                    String::from("No tower data returned"),
                )
            })?;

        if print {
            for towers in &towers {
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    towers.id,
                    towers.date,
                    towers.name,
                    towers.location,
                    towers.last_ping,
                    towers.stage,
                    towers.pipelines.iter().map(|p| format!("{p}")).join(", "),
                );
            }
        }

        Ok(towers)
    }

    pub fn delete_tower(
        &self, 
        id: Option<String>, 
        json: Option<PathBuf>, 
        name: Option<String>, 
        location: Option<String>
    ) -> Result<Option<ProductionTower>, HttpClientError> {

        self.log_team_warning();

        let mut url = self.routes.url(Route::TeamTowersDelete);  // deletes all

        if let Some(id) = id {
            url = format!("{url}/{id}")
        }

        if let Some(json) = json {
            let id = RegisterTowerSchema::from_json(&json)?.id;
            url = format!("{url}/{id}")
        }

        if name.is_some() | location.is_some() {
            url = self.build_request_url(
                self.routes.url(Route::TeamTowersDelete), 
                &[("name", name), ("location", location)]
            );
        }

        let response = self.send_request_with_team(
            self.client.delete(url)
        )?;

        Ok(self.handle_response::<DeleteTowerResponse>(
            response,
            None,
            "Tower deletion failed",
        )?
        .data)
    }

    pub fn ping_tower(&self, tower_id: &str, print: bool) -> Result<String, HttpClientError> {

        self.log_team_warning();

        let response = self.send_request_with_team(
            self.client
                .patch(&format!(
                    "{}/{}",
                    self.routes.url(Route::TeamTowersPing),
                    tower_id
                ))
        )?;

        let data = self
            .handle_response::<PingTowerResponse>(
                response,
                None,
                "Tower ping failed",
            )?
            .data
            .ok_or_else(|| {
                HttpClientError::DataResponseFailure(
                    reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                    String::from("No data returned after tower ping"),
                )
            })?;

        if print {
            println!("{}", data);
        }

        Ok(data)
    }

    pub fn register_watcher(
        &self,
        register_watcher_schema: &RegisterWatcherSchema,
        print: bool,
    ) -> Result<String, HttpClientError> {

        self.log_team_warning();

        let response = self.send_request_with_team(
            self.client
                .post(self.routes.url(Route::TeamWatchersRegister))
                .json(register_watcher_schema),
        )?;

        let watcher_id = self.handle_response::<RegisterWatcherResponse>(
            response,
            Some("Watcher registered successfully"),
            "Watcher registration failed",
        )?
        .data
        .unwrap_or(register_watcher_schema.id.clone());

        if print {
            println!("{}", watcher_id);
        }

        Ok(watcher_id)
    }

    pub fn list_watchers(
        &self,
        id: Option<String>,
        print: bool,
    ) -> Result<Vec<ProductionWatcher>, HttpClientError> {

        self.log_team_warning();

        let url = self.build_request_url(
            self.routes.url(Route::TeamWatchersList), 
            &[("id", id)]
        );

        let response = self.send_request_with_team(
            self.client
                .get(&url)
        )?;

        let watchers = self
            .handle_response::<ListWatchersResponse>(
                response,
                None,
                "Failed to retrieve watchers",
            )?
            .data
            .ok_or_else(|| {
                HttpClientError::DataResponseFailure(
                    reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                    String::from("No watcher data returned"),
                )
            })?;

        if print {
            for watcher in &watchers {
                println!(
                    "{}\t{} @ {}\t{}\t'{}'\t{}",
                    watcher.id,
                    watcher.name,
                    watcher.location,
                    watcher.format,
                    watcher.glob,
                    watcher.last_ping
                );
            }
        }

        Ok(watchers)
    }

    pub fn delete_watcher(
        &self, 
        id: Option<String>, 
        json: Option<PathBuf>, 
        name: Option<String>, 
        location: Option<String>
    ) -> Result<Option<ProductionWatcher>, HttpClientError> {

        self.log_team_warning();

        let mut url = self.routes.url(Route::TeamWatchersDelete); // deletes all

        if let Some(id) = id {
            url = format!("{url}/{id}")
        }

        if let Some(json) = json {
            let id = RegisterWatcherSchema::from_json(&json)?.id;
            url = format!("{url}/{id}")
        }

        if name.is_some() | location.is_some() {
            url = self.build_request_url(
                self.routes.url(Route::TeamWatchersDelete), 
                &[("name", name), ("location", location)]
            );
        }

        let response = self.send_request_with_team(
            self.client.delete(url)
        )?;

        Ok(self.handle_response::<DeleteWatcherResponse>(
            response,
            None,
            "Watcher deletion failed",
        )?
        .data)
    }

    pub fn ping_watcher(&self, id: &str, print: bool) -> Result<String, HttpClientError> {

        self.log_team_warning();

        let response = self.send_request_with_team(
            self.client
                .patch(&format!(
                    "{}/{}",
                    self.routes.url(Route::TeamWatchersPing),
                    id
                ))
        )?;

        let data = self
            .handle_response::<PingWatcherResponse>(
                response,
                None,
                "Watcher ping failed",
            )?
            .data
            .ok_or_else(|| {
                HttpClientError::DataResponseFailure(
                    reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                    String::from("No data returned after watcher ping"),
                )
            })?;

        if print {
            println!("{}", data);
        }

        Ok(data)
    }
    pub fn register_staged_samples(
        &self,
        register_staged_sample_schema: &RegisterStagedSampleSchema
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();
        self.log_project_warning();

        let response = self.send_request_with_team_db_project(
            self.client
                .post(self.routes.url(Route::TeamStagedSamplesRegister))
                .json(register_staged_sample_schema)
        )?;

        self.handle_response::<serde_json::Value>(
            response,
            Some("Samples staged successfully"),
            "Sample staging failed",
        )?;
        Ok(())
    }
    pub fn list_staged_samples(
        &self,
        id: &str,
        run_id: Option<String>,
        sample_id: Option<String>,
        print: bool,
    ) -> Result<Vec<StagedSample>, HttpClientError> {

        self.log_team_warning();

        let url = self.build_request_url(
            &format!("{}/{}", self.routes.url(Route::TeamStagedSamplesList), id), &[
                ("run_id", run_id), ("sample_id", sample_id)
            ]
        );

        let response = self.send_request_with_team(
            self.client
                .get(&url)
        )?;

        let staged_samples = self
            .handle_response::<ListStagedSamplesResponse>(
                response,
                None,
                "Failed to retrieve staged samples",
            )?
            .data
            .ok_or_else(|| {
                HttpClientError::DataResponseFailure(
                    reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                    String::from("No stage sample data returned"),
                )
            })?;

        if print {
            for staged_sample in &staged_samples {
                println!(
                    "{}\t{}\t{}\t{}\t{} @ {}\t{}\t{}",
                    staged_sample.id,
                    staged_sample.date,
                    staged_sample.sample_id,
                    staged_sample.database,
                    staged_sample.project,
                    staged_sample.tower.name,
                    staged_sample.tower.location,
                    staged_sample.tower.stage,
                );
            }
        }

        Ok(staged_samples)
    }

    pub fn delete_staged_sample(
        &self, 
        id: &str, 
        staged_id: Option<String>,
        run_id: Option<String>,
        sample_id: Option<String>,
    ) -> Result<Option<StagedSample>, HttpClientError> {

        self.log_team_warning();

        let mut url = format!("{}/{id}", self.routes.url(Route::TeamStagedSamplesDelete)); // deletes all

        if run_id.is_some() | sample_id.is_some() | staged_id.is_some() {
            url = self.build_request_url(
                format!("{}/{id}", self.routes.url(Route::TeamStagedSamplesDelete)), 
                &[("run_id", run_id), ("sample_id", sample_id), ("stage_id", staged_id)]
            );
        }

        let response = self.send_request_with_team(
            self.client
                .delete(url)
        )?;

        Ok(self.handle_response::<DeleteStagedSampleResponse>(
            response,
            None,
            "Staged sample deletion failed",
        )?
        .data)
    }

    pub fn pull_staged_samples(
        &self, 
        id: &str, 
        run_id: Option<String>,
        sample_id: Option<String>,
        outdir: &PathBuf,
        delete: bool,
    ) -> Result<(), HttpClientError> {
        
        self.log_team_warning();

        if !outdir.exists() && outdir.is_dir() {
            create_dir_all(&outdir)?;
        }

        let staged_samples = match self.list_staged_samples(id, run_id, sample_id, false){
            Ok(samples) => samples, Err(err) => {
                match err {
                    HttpClientError::ResponseFailure(code) => {
                        if code == StatusCode::NOT_FOUND {
                            Vec::new()  // no samples staged will not raise 404
                        } else {
                            return Err(err)
                        }
                    },
                    _ => return Err(err)
                }
            }
        };

        let timestamp = Utc::now().to_rfc3339_opts(chrono::SecondsFormat::Secs, true);

        for staged_sample in &staged_samples {
            staged_sample.to_json(&outdir.join(
                format!("{}.json", staged_sample.id)
            ))?;
        }

        if delete {
            for staged_sample in &staged_samples {
                self.delete_staged_sample(
                    &id, 
                    Some(staged_sample.id.to_string()), 
                    None, 
                    None
                )?;
            }
        }

        Ok(())
    }

    pub fn download_models(
        &self,
        outdir: &PathBuf
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();
        self.log_project_warning();

        let url = format!("{}", self.routes.url(Route::DataCerebroRetrieveModel));

        let response = self.send_request_with_team_db_project(
            self.client.get(url)
        )?;

        let model_response = self.handle_response::<RetrieveModelResponse>(
            response,
            Some("Models retrieved for this project"),
            "Model download failed for this project",
        )?;

        if model_response.data.is_empty() {
            log::warn!("No models found for this query")
        } else {
            for model in model_response.data {
                log::info!("Writing model to file: {}.json", model.name);
                model.write_json(
                    &outdir.join(
                        format!("{}.json", model.name)
                    )
                )?
            }
        }

        Ok(())
    }
    pub fn upload_models(
        &self,
        models: &[Cerebro]
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();
        self.log_project_warning();
        
        for model in models {
            
            if model.sample.id.is_empty() {
                return Err(HttpClientError::ModelSampleIdentifierEmpty);
            }
            
            model.check_size(
                self.max_model_size, 
                true, 
                Some(16.0)
            )?;

            let url = format!("{}", self.routes.url(Route::DataCerebroInsertModel));

            let response = self.send_request_with_team_db_project(
                self.client
                    .post(url)
                    .json(model)
            )?;

            self.handle_response::<serde_json::Value>(
                response,
                Some(&format!(
                    "Model for sample library {} uploaded successfully",
                    model.sample.id
                )),
                "Upload failed",
            )?;
        }

        Ok(())
    }

    pub fn update_models(
        &self,
        run_tsv: Option<PathBuf>
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();
        self.log_project_warning();
        

        if let Some(path) = run_tsv {

            // Run configuration update
            let update_data: Vec<UpdateRunConfigSchema> = read_tsv(&path, false, true)?;

            for update_schema in update_data {

                let url = format!("{}", self.routes.url(Route::DataCerebroUpdateModelRunConfig));

                let response = self.send_request_with_team_db_project(
                    self.client
                        .patch(url)
                        .json(&update_schema)
                )?;

                log::info!("Update schema: {:?}", &update_schema);

                match self.handle_response::<serde_json::Value>(
                    response,
                    Some(&format!(
                        "Model for sample library {} updated successfully",
                        update_schema.sample_id
                    )),
                    "Update failed",
                ) {
                    Ok(_) => {},
                    Err(err) => log::warn!("{}", err.to_string())
                };
            }
        }

        Ok(())
    }
    pub fn upload_training_prefetch(
        &self, 
        prefetch: &Vec<PathBuf>,
        collection: &str,
        description: &str
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();

        for file in prefetch {

            let schema = CreateTrainingPrefetch::from_file(file, collection, description)?;
            let url = format!("{}", self.routes.url(Route::TrainingRegister));

            log::info!("Uploading prefetch data '{}' to collection '{}'", schema.name, collection);

            let response = self.send_request_with_team(
                self.client
                    .post(url)
                    .json(&schema)
            )?;

            let _ = self.handle_response::<TrainingResponse<String>>(
                response,
                None,
                "Prefetch data upload failed",
            )?;
        }


        Ok(())
    }

    pub fn list_training_collections(
        &self, 
        collection: Option<String>
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();

        let url = if let Some(c) = collection {
            self.build_request_url(
                format!("{}", self.routes.url(Route::TrainingList)), 
                &[("collection", c)]
            )
        } else {
            self.routes.url(Route::TrainingList)
        };
        

        let response = self.send_request_with_team(
            self.client
                .get(url)
        )?;

        let training_response = self.handle_response::<TrainingResponse<Vec<TrainingPrefetchData>>>(
            response,
            None,
            "Prefetch data retrieval failed",
        )?;

        if let Some(data) = training_response.data {
            for training_prefetch in data  {
                log::info!("collection={} description='{}' name={} ", training_prefetch.collection, training_prefetch.description, training_prefetch.name);
            }
        }

        Ok(())
    }


    pub fn delete_training_data(
        &self, 
        collection: String,
        name: Option<String>
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();

        
        let url = self.build_request_url(
            format!("{}", self.routes.url(Route::TrainingList)), 
            &[("collection", collection.clone())]
        );

        let response = self.send_request_with_team(
            self.client
                .get(url)
        )?;

        let training_response = self.handle_response::<TrainingResponse<Vec<TrainingPrefetchData>>>(
            response,
            None,
            "Prefetch data retrieval failed",
        )?;


        if let Some(data) = training_response.data {
            log::info!("{} entries in training collection '{}'", data.len(), &collection);
            for training_prefetch in data  {

                if let Some(ref name) = name {
                    if name.as_str() != training_prefetch.name {
                        continue 
                    }
                }

                let url = format!("{}/{}", self.routes.url(Route::TrainingDelete), training_prefetch.id);

                log::info!("Deleting entry '{}' from collection '{collection}'", training_prefetch.name);

                let response = self.send_request_with_team(
                    self.client
                        .delete(url)
                )?;
                
                let _ = self.handle_response::<TrainingResponse<String>>(
                    response,
                    None,
                    "Deletion failed",
                )?;
            }
        } else {
            log::warn!("No data was returned for this collection")
        }

        Ok(())
    }

    pub fn list_training_sessions(
        &self, 
        completed: bool,
        user_name: Option<String>,
        results: Option<PathBuf>
    ) -> Result<Option<Vec<TrainingSessionData>>, HttpClientError> {

        self.log_team_warning();

        let url = self.routes.url(Route::TrainingSessionList);

        let response = self.send_request_with_team(
            self.client
                .get(url)
        )?;

        let training_response = self.handle_response::<TrainingResponse<Vec<TrainingSessionData>>>(
            response,
            None,
            "Training session data retrieval failed",
        )?;

        if let Some(ref data) = training_response.data {

            if data.is_empty() {
                log::warn!("No training session data available for this team")
            }

            for session in data  {
                if completed {
                    if session.completed.is_some() {
                        if let Some(ref user_name) = user_name {
                            if *user_name == session.user_name {
                                log::info!("id={} user_name='{}' user_id={} started={} completed={}", session.id, session.user_name, session.user_id, session.started, session.completed.as_deref().unwrap_or("false"));
                            }
                        } else {
                            log::info!("id={} user_name='{}' user_id={} started={} completed={}", session.id, session.user_name, session.user_id, session.started, session.completed.as_deref().unwrap_or("false"));
                        }
                        
                    }
                } else {
                    if let Some(ref user_name) = user_name {
                        if *user_name == session.user_name {
                            log::info!("id={} user_name='{}' user_id={} started={} completed={}", session.id, session.user_name, session.user_id, session.started, session.completed.as_deref().unwrap_or("false"));
                        }
                    } else {
                        log::info!("id={} user_name='{}' user_id={} started={} completed={}", session.id, session.user_name, session.user_id, session.started, session.completed.as_deref().unwrap_or("false"));
                    }
                }
                
            }

            // Output session results
            if let Some(ref results) = results {

                create_dir_all(results)?;

                for session in data {
                    if let Some(date) = &session.completed {
                        
                        let user_path = results.join(session.user_name.replace(" ", ""));
                        let session_path = user_path.join(&session.collection).join(date);
                        let session_ciqa = session_path.join("ciqa");

                        create_dir_all(&session_ciqa)?;

                        log::info!("Output complete session results to: {}", session_path.display());

                        // Output full session data:
                        session.to_json(session_path.join("session.json"))?;

                        if let Some(ref result_data) = session.result {

                            log::info!("dataset={} user_name='{}' completed={} - n={} sensitivity={:.1} specificity={:.1}", session.collection, session.user_name, date, result_data.total, result_data.sensitivity, result_data.specificity);

                            // Output evaluation results table
                            write_tsv(&result_data.data, &session_path.join("results.tsv"), true)?;
                            
                            // Output diagnostic results compatible with META-GPT and CIQA
                            for record in &result_data.data {
                                let result = DiagnosticResult::from_training_result(&record);
                                let result_filename = format!("{}.json", record.sample_name.clone().unwrap_or(record.record_id.clone()));
                                result.to_json(&session_ciqa.join(result_filename))?;
                            }
                        }

                    }
                }
            }
        }

        Ok(training_response.data)
    }


    pub fn delete_training_sessions(
        &self, 
        session_id: Option<String>,
        user_id: Option<String>,
        incomplete: bool
    ) -> Result<(), HttpClientError> {

        self.log_team_warning();

        // Delete single session by identifier
        if let Some(id) = session_id {

            let url = format!("{}/{id}", self.routes.url(Route::TrainingSessionDelete));

            let response = self.send_request_with_team(
                self.client
                    .delete(url)
            )?;

            self.handle_response::<TrainingResponse<()>>(
                response,
                None,
                "Training session data retrieval failed",
            )?;
        }

        let url = self.routes.url(Route::TrainingSessionList);

        let response = self.send_request_with_team(
            self.client
                .get(url)
        )?;

        let training_response = self.handle_response::<TrainingResponse<Vec<TrainingSessionData>>>(
            response,
            None,
            "Training session data retrieval failed",
        )?;

        if let Some(data) = training_response.data {
            for session in data  {

                if let Some(ref user_id) = user_id {
                    if *user_id != session.user_id {
                        continue;
                    }
                }

                if incomplete {
                    if session.completed.is_some() {
                        continue;
                    }
                };

                let url = format!("{}/{}", self.routes.url(Route::TrainingSessionDelete), session.id);

                let response = self.send_request_with_team(
                    self.client
                        .delete(url)
                )?;

                self.handle_response::<TrainingResponse<()>>(
                    response,
                    None,
                    "Training session data retrieval failed",
                )?;
            }
        }

        Ok(())
    }

    pub fn get_identifiers(
        &self,
        schema: &CerebroIdentifierSchema,
    ) -> Result<Option<Vec<CerebroIdentifierSummary>>, HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();
        self.log_project_warning();
        
        let url = format!("{}", self.routes.url(Route::DataCerebroIdentifiers));

        let response = self.send_request_with_team_db_project(
            self.client
                .post(url)
                .json(schema)
        )?;

        let response = self.handle_response::<CerebroIdentifierResponse>(
            response,
            None,
            "Cerebro identifiers retrieval failed",
        )?;

        Ok(response.data)
    } 
    pub fn _get_aneuploidy(&self, sample: &str) -> Result<Option<&str>, HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();
        self.log_project_warning();

        Ok(None)
    }
    pub fn get_quality_control(&self, schema: &QualityControlTableSchema, csv: Option<PathBuf>) -> Result<Vec<ReadQualityControl>, HttpClientError> {

        let url = self.build_request_url(
            format!("{}", self.routes.url(Route::DataCerebroQualityControl)), 
            &[("csv", if let Some(_) = csv { true.to_string() } else { false.to_string() })]
        );

        let response = self.send_request_with_team_db_project(
            self.client.post(url).json(schema)
        )?;

        let qc_table_response = self.handle_response::<QualityControlTableResponse>(
            response,
            None,
            "Sample summary (quality control) retrieval failed",
        )?;

        if let Some(path) = csv {
            let mut writer = BufWriter::new(File::create(path)?);
            write!(&mut writer, "{}", qc_table_response.data.csv)?;
            writer.flush()?;
        }

        Ok(qc_table_response.data.records)
    }

    pub fn get_pathogen_detection(&self, schema: &PathogenDetectionTableSchema, csv: Option<PathBuf>) -> Result<Vec<PathogenDetectionTableRecord>, HttpClientError> {

        let url = self.build_request_url(
            format!("{}", self.routes.url(Route::DataCerebroPathogenDetection)), 
            &[("csv", if let Some(_) = csv { true.to_string() } else { false.to_string() })]
        );

        let response = self.send_request_with_team_db_project(
            self.client.post(url).json(schema)
        )?;

        let pathogen_detection_table_response = self.handle_response::<PathogenDetectionTableResponse>(
            response,
            None,
            "Sample summary (pathogen detection table) retrieval failed",
        )?;

        if let Some(path) = csv {
            let mut writer = BufWriter::new(File::create(path)?);
            write!(&mut writer, "{}", pathogen_detection_table_response.data.csv)?;
            writer.flush()?;
        }

        Ok(pathogen_detection_table_response.data.records)
    }
    pub fn get_taxon_history(&self, taxon_label: String, host_label: String, regression: bool, print_regression: bool, plot: Option<PathBuf>) -> Result<Option<RpmAnalysisResult>, HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();
        self.log_project_warning();

        let url = self.build_request_url(
            format!("{}", self.routes.url(Route::DataCerebroTaxaHistory)), 
            &[("taxon_label", taxon_label), ("host_label", host_label)]
        );

        let response = self.send_request_with_team_db_project(
            self.client.get(url)
        )?;

        let response = self.handle_response::<TaxonHistoryResponse>(
            response,
            None,
            "Taxon history retrieval failed",
        )?;

        if regression & !response.data.is_empty() {

            // Build the analysis configuration.
            let config = RpmConfigBuilder::new()
            .confidence_level(0.95)
            .only_high_outliers(true)
            .build();
                
            // Create the analyzer from taxon history.
            let analyzer = RpmAnalyzer::from_taxon_history(config, response.data);

            // Run the regression and outlier detection.
            let result = match analyzer.run() {
                Ok(result) => result,
                Err(_) => return Ok(None)
            };

            if print_regression {
                log::info!("Regression result:");
                log::info!("Intercept (log₁₀ scale): {:.4}", result.intercept);
                log::info!("Slope: {:.4} ({})", result.slope, result.relationship);
                log::info!("Residual standard error: {:.4}", result.residual_standard_error);

                if !result.outliers.is_empty() {
                    log::info!("Outliers detected:");
                    for outlier in &result.outliers {
                        log::info!("Sample ID: {}, Raw host RPM: {:.2}, Raw taxon RPM: {:.2}", 
                            outlier.sample_id, outlier.raw_host, outlier.raw_taxon);
                    }
                } else {
                    log::info!("No outliers detected.");
                }
            }

            if let Some(path) = plot {
                analyzer.plot_regression(&result, &path.display().to_string())?;
            }
            

            Ok(Some(result))

        } else {
            Ok(None)
        }
    }
    fn split_taxa(&self, map: HashMap<String, (Vec<Taxon>, Vec<Taxon>)>) -> (Vec<Taxon>, Vec<Taxon>) {
        let mut sample_taxa: Vec<Taxon> = Vec::new();
        let mut contam_taxa: Vec<Taxon> = Vec::new();
        
        for (_key, (sample_vec, contam_vec)) in map.into_iter() {
            sample_taxa.extend(sample_vec);
            contam_taxa.extend(contam_vec);
        }
        
        (sample_taxa, contam_taxa)
    }
    pub fn get_prevalence_contamination(&self, 
        contam_config: &PrevalenceContaminationConfig,
        tags: Vec<String>,
    ) -> Result<HashSet<String>, HttpClientError> {

        // Obtain the collection's prevalence contamination for this tag

        let url = format!("{}", self.routes.url(Route::DataCerebroTaxaContamination));

        let contamination_schema = ContaminationSchema::from_config(
            vec![], // not limited to specific taxids present in that sample, obtain all from collection
            tags, 
            &contam_config
        );

        let response = self.send_request_with_team_db_project(
            self.client
                .post(url)
                .json(&contamination_schema)
        )?;


        let contam_response_data = self.handle_response::<ContaminationTaxaResponse>(
            response,
            None,
            "Prevalence contamination retrieval failed",
        )?;

        Ok(HashSet::from_iter(contam_response_data.data.taxid))

    }
    pub fn get_taxa(
        &self,
        schema: &CerebroIdentifierSchema,
        filter_config: &TaxonFilterConfig,
        contam_taxid_tagged: HashMap<String, HashSet<String>>, // tag: {taxid}
        contam_history: bool
    ) -> Result<(Vec<Taxon>, Vec<Taxon>), HttpClientError> {

        self.log_team_warning();
        self.log_db_warning();
        self.log_project_warning();

        // Ensure that the contam config has the same "collapse_variants" setting as the filter config,
        // since the returned 'taxid' from the filter required for the contamination detection may have
        // changed

        let identifier_response = self.get_identifiers(schema)?;

        let mut tag_data = HashMap::new();
        if let Some(identifiers) = identifier_response {

            // Initialize the DNA and RNA groups
            let mut groups: HashMap<&str, Vec<CerebroIdentifierSummary>> = HashMap::from([
                ("DNA", Vec::new()), ("RNA", Vec::new())
            ]);

            // Group each summary based on its sample tags
            for summary in identifiers {

                if summary.sample_tags.contains(&"DNA".to_string()) {
                    groups.get_mut("DNA").unwrap().push(summary.clone());
                }

                if summary.sample_tags.contains(&"RNA".to_string()) {
                    groups.get_mut("RNA").unwrap().push(summary.clone());
                }
            }

            // Log the grouped summaries
            for (tag, summaries) in groups {

                let contam_taxids = match contam_taxid_tagged.get(tag) {
                    Some(contam_taxids) => Some(contam_taxids),
                    None => {
                        log::warn!("No prevalence contamination taxonomic identifiers provided for tag: {}", tag);
                        None
                    }
                };

                if !summaries.is_empty() {

                    let ids: Vec<String> = summaries.iter().map(|s| s.cerebro_id.clone()).collect();

                    let url = self.build_request_url(
                        format!("{}", self.routes.url(Route::DataCerebroTaxaFiltered)), 
                        &[("id", ids.join(","))]
                    );

                    let response = self.send_request_with_team_db_project(
                        self.client
                            .post(url)
                            .json(filter_config)
                    )?;

                    let taxa_response_data = self.handle_response::<FilteredTaxaResponse>(
                        response,
                        None,
                        "Taxon retrieval failed",
                    )?;

                    if let Some(contam_taxids) = contam_taxids {

                        let contam_taxa: Vec<Taxon> = taxa_response_data.data.taxa
                            .iter()
                            .filter(|tax| contam_taxids.contains(&tax.taxid))
                            .cloned()
                            .collect();

                        let mut sample_control_taxa: Vec<Taxon> = taxa_response_data.data.taxa
                            .iter()
                            .filter(|tax| !contam_taxids.contains(&tax.taxid))
                            .cloned()
                            .collect();

                        let mut clear_contam_taxa = Vec::new();
                        if contam_history {
                            'contam: for contam_taxon in &contam_taxa {
                                let reg = self.get_taxon_history(
                                    format!("s__{}", contam_taxon.name), 
                                    format!("s__Homo sapiens"), // hardcoded for now
                                    true,
                                    false,
                                    None
                                )?;
                                if let Some(reg) = reg {
                                    for outlier in reg.outliers {
                                        if outlier.sample_id == schema.sample {
                                            sample_control_taxa.push(contam_taxon.to_owned());
                                            continue 'contam;
                                        }
                                    }
                                }
                                clear_contam_taxa.push(contam_taxon.clone())
                            }    
                        } else {
                            clear_contam_taxa = contam_taxa.clone();
                        }
                        tag_data.insert(tag.to_string(), (sample_control_taxa.clone(), clear_contam_taxa.clone()));
                    } else {
                        tag_data.insert(tag.to_string(), (taxa_response_data.data.taxa, Vec::new()));
                    }
                }   
            }
        }   
        
        Ok(self.split_taxa(tag_data))
    }

    
}
