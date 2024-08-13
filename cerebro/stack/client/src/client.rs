use std::fs::File;
use std::io::Write;
use anyhow::Result;
use cerebro_model::api::files::model::SeaweedFile;
use cerebro_model::api::files::response::DeleteFileResponse;
use cerebro_model::api::files::response::ListFilesResponse;
use cerebro_model::api::files::response::RegisterFileResponse;
use cerebro_model::api::pipelines::model::Pipeline;
use cerebro_model::api::pipelines::model::ProductionPipeline;
use cerebro_model::api::pipelines::response::DeletePipelineResponse;
use cerebro_model::api::pipelines::response::ListPipelinesResponse;
use cerebro_model::api::pipelines::response::PingPipelineResponse;
use cerebro_model::api::pipelines::response::RegisterPipelineResponse;
use cerebro_model::api::pipelines::schema::RegisterPipelineSchema;
use cerebro_model::api::stage::model::StagedSample;
use cerebro_model::api::stage::response::DeleteStagedSampleResponse;
use cerebro_model::api::stage::response::ListStagedSamplesResponse;
use cerebro_model::api::stage::schema::RegisterStagedSampleSchema;
use cerebro_model::api::watchers::model::ProductionWatcher;
use cerebro_model::api::watchers::response::DeleteWatcherResponse;
use cerebro_model::api::watchers::response::ListWatchersResponse;
use cerebro_model::api::watchers::response::PingWatcherResponse;
use cerebro_model::api::watchers::response::RegisterWatcherResponse;
use cerebro_model::api::watchers::schema::RegisterWatcherSchema;
use reqwest::blocking::RequestBuilder;
use reqwest::blocking::Response;
use serde::de::DeserializeOwned;
use std::path::PathBuf;
use itertools::Itertools;
use reqwest::blocking::Client;
use rpassword::prompt_password;
use reqwest::header::AUTHORIZATION;
use serde::{Serialize, Deserialize};
use actix_web_httpauth::headers::authorization::Bearer;

use cerebro_model::api::cerebro::response::TaxaSummaryDataResponse;
use cerebro_model::api::utils::ErrorResponse;
use cerebro_model::api::cerebro::model::Cerebro;
use cerebro_model::api::auth::schema::AuthLoginSchema;
use cerebro_model::api::teams::model::ProjectCollection;
use cerebro_model::api::users::response::UserSelfResponse;
use cerebro_model::api::users::response::UserSelfTeamResponse;
use cerebro_model::api::auth::response::AuthLoginResponseSuccess;
use cerebro_workflow::filters::TaxonFilterConfig;
use cerebro_model::api::cerebro::schema::{SampleSummaryQcSchema, TaxaSummarySchema};
use cerebro_model::api::teams::model::TeamDatabase;
use cerebro_model::api::teams::schema::RegisterProjectSchema;
use cerebro_model::api::files::schema::RegisterFileSchema;

use crate::error::HttpClientError;
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
    DataCerebroTaxaSummary,
    TeamProjectCreate,
    TeamFilesRegister,
    TeamFilesList,
    TeamFilesDelete,
    TeamPipelinesRegister,
    TeamPipelinesList,
    TeamPipelinesDelete,
    TeamPipelinesPing,
    TeamWatchersRegister,
    TeamWatchersList,
    TeamWatchersDelete,
    TeamWatchersPing,
    TeamStagedSamplesRegister,
    TeamStagedSamplesList,
    TeamStagedSamplesDelete,
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
            Route::DataCerebroTaxaSummary => "cerebro/taxa/summary",
            Route::TeamProjectCreate => "teams/project",
            Route::TeamFilesRegister => "files/register",
            Route::TeamFilesList => "files",
            Route::TeamFilesDelete => "files",
            Route::TeamPipelinesRegister => "pipeline/register",
            Route::TeamPipelinesList => "pipeline",
            Route::TeamPipelinesDelete => "pipeline",
            Route::TeamPipelinesPing => "pipeline",
            Route::TeamWatchersRegister => "watcher/register",
            Route::TeamWatchersList => "watcher",
            Route::TeamWatchersDelete => "watcher",
            Route::TeamWatchersPing => "watcher",
            Route::TeamStagedSamplesRegister => "stage/register",
            Route::TeamStagedSamplesList => "stage",
            Route::TeamStagedSamplesDelete => "stage"
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
    pub routes: CerebroRoutes,
    pub token_data: Option<AuthLoginTokenFile>,
    pub token_file: Option<PathBuf>
}

impl CerebroClient {
    pub fn new(
        url: &str,
        token: Option<String>,
        login: bool,
        danger_accept_invalid_certs: bool,
        token_file: Option<PathBuf>,
        team: Option<String>,
    ) -> Result<Self, HttpClientError> {
        let url_clean = url.trim_end_matches('/').to_string();

        let token_data = if login {
            None
        } else {
            Self::load_token(token, token_file.clone())?
        };

        let client = Client::builder()
            .danger_accept_invalid_certs(danger_accept_invalid_certs)
            .build()?;

        if let None = team {
            log::warn!("No team access configured - you may need this for data access");
            log::warn!("Use global argument '--team' or environment variable '$CEREBRO_USER_TEAM'")
        }

        Ok(Self {
            url: url_clean.clone(),
            routes: CerebroRoutes::new(&url_clean),
            token_data,
            token_file,
            client,
            team,
        })
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

        log::error!("Failed to obtain token from input or environmental variable");
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
        let team = self.team.as_deref().ok_or(HttpClientError::RequireTeamNotConfigured)?;

        let response = request
            .query(&[("team", team)])
            .header(AUTHORIZATION, self.get_bearer_token(None)).send()?;
        
        Ok(response)
    }

    fn handle_response<T: DeserializeOwned>(
        &self,
        response: Response,
        success_msg: Option<&str>,
        failure_msg: &str,
    ) -> Result<T, HttpClientError> {
        let status = response.status();
        if status.is_success() {
            if let Some(msg) = success_msg {
                log::info!("{}", msg);
            }
            response.json().map_err(|_| {
                HttpClientError::DataResponseFailure(
                    status,
                    String::from("failed to parse response data"),
                )
            })
        } else {
            let error_response: ErrorResponse = response.json().map_err(|_| {
                HttpClientError::DataResponseFailure(
                    status,
                    String::from("failed to parse error response"),
                )
            })?;
            log::error!("{}: {}", failure_msg, error_response.message);
            Err(HttpClientError::ResponseFailure(status))
        }
    }

    fn build_request_url<T>(&self, route: Route, params: &[(&str, T)]) -> String
    where
        T: Into<Option<String>> + Clone,
    {
        let mut url = self.routes.url(route);
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
        team_name: &str,
        db_name: &str,
        project_name: &str,
        project_description: &str,
    ) -> Result<(), HttpClientError> {
        let url = self.build_request_url(
            Route::TeamProjectCreate,
            &[("team_name", Some(team_name.to_string())), ("db_name", Some(db_name.to_string()))],
        );

        let project_schema = RegisterProjectSchema {
            project_name: project_name.to_owned(),
            project_description: project_description.to_owned(),
            project_mongo_name: project_name.split_whitespace().join("_").to_lowercase(),
        };

        let response = self
            .client
            .post(&url)
            .header(AUTHORIZATION, self.get_bearer_token(None))
            .json(&project_schema)
            .send()?;

        self.handle_response::<serde_json::Value>(
            response,
            Some(&format!("Project `{}` created successfully", project_name)),
            "Project creation failed",
        )?;
        Ok(())
    }
    
    pub fn register_file(
        &self,
        register_file_schema: RegisterFileSchema,
    ) -> Result<(), HttpClientError> {
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

    pub fn delete_file(&self, file_id: &str) -> Result<SeaweedFile, HttpClientError> {
        let response = self.send_request_with_team(
            self.client
                .delete(&format!(
                    "{}/{}",
                    self.routes.url(Route::TeamFilesDelete),
                    file_id
                ))
        )?;

        self.handle_response::<DeleteFileResponse>(
            response,
            Some(&format!("File `{}` deleted successfully", file_id)),
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

    pub fn get_files(
        &self,
        run_id: Option<String>,
        watcher_id: Option<String>,
        page: u32,
        limit: u32,
        print: bool,
    ) -> Result<Vec<SeaweedFile>, HttpClientError> {
        let url = self.build_request_url(
            Route::TeamFilesList,
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

    pub fn register_pipeline(
        &self,
        register_pipeline_schema: &RegisterPipelineSchema,
        print: bool,
    ) -> Result<String, HttpClientError> {
        let response = self.send_request_with_team(
            self.client
                .post(self.routes.url(Route::TeamPipelinesRegister))
                .json(register_pipeline_schema),
        )?;

        let pipeline_id = self.handle_response::<RegisterPipelineResponse>(
            response,
            Some("Pipeline registered successfully"),
            "Pipeline registration failed",
        )?
        .data
        .unwrap_or(register_pipeline_schema.id.clone());

        if print {
            println!("{}", pipeline_id);
        }

        Ok(pipeline_id)
    }

    pub fn get_pipelines(
        &self,
        id: Option<String>,
        print: bool,
    ) -> Result<Vec<ProductionPipeline>, HttpClientError> {
        let url = self.build_request_url(Route::TeamPipelinesList, &[("id", id)]);

        let response = self.send_request_with_team(
            self.client
                .get(&url)
        )?;

        let pipelines = self
            .handle_response::<ListPipelinesResponse>(
                response,
                None,
                "Failed to retrieve pipelines",
            )?
            .data
            .ok_or_else(|| {
                HttpClientError::DataResponseFailure(
                    reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                    String::from("No pipeline data returned"),
                )
            })?;

        if print {
            for pipeline in &pipelines {
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    pipeline.pipeline,
                    pipeline.id,
                    pipeline.date,
                    pipeline.name,
                    pipeline.location,
                    pipeline.last_ping
                );
            }
        }

        Ok(pipelines)
    }

    pub fn delete_pipeline(&self, id: &str) -> Result<ProductionPipeline, HttpClientError> {
        let response = self.send_request_with_team(
            self.client
                .delete(&format!(
                    "{}/{}",
                    self.routes.url(Route::TeamPipelinesDelete),
                    id
                ))
        )?;

        self.handle_response::<DeletePipelineResponse>(
            response,
            Some(&format!("Pipeline `{}` deleted successfully", id)),
            "Pipeline deletion failed",
        )?
        .data
        .ok_or_else(|| {
            HttpClientError::DataResponseFailure(
                reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                String::from("Pipeline data not returned after deletion"),
            )
        })
    }

    pub fn ping_pipeline(&self, pipeline_id: &str, print: bool) -> Result<String, HttpClientError> {
        let response = self.send_request_with_team(
            self.client
                .patch(&format!(
                    "{}/{}",
                    self.routes.url(Route::TeamPipelinesPing),
                    pipeline_id
                ))
        )?;

        let data = self
            .handle_response::<PingPipelineResponse>(
                response,
                None,
                "Pipeline ping failed",
            )?
            .data
            .ok_or_else(|| {
                HttpClientError::DataResponseFailure(
                    reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                    String::from("No data returned after pipeline ping"),
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

    pub fn get_watchers(
        &self,
        id: Option<String>,
        print: bool,
    ) -> Result<Vec<ProductionWatcher>, HttpClientError> {
        let url = self.build_request_url(Route::TeamWatchersList, &[("id", id)]);

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

    pub fn delete_watcher(&self, watcher_id: &str) -> Result<ProductionWatcher, HttpClientError> {
        let response = self.send_request_with_team(
            self.client
                .delete(&format!(
                    "{}/{}",
                    self.routes.url(Route::TeamWatchersDelete),
                    watcher_id
                ))
        )?;

        self.handle_response::<DeleteWatcherResponse>(
            response,
            Some(&format!("Watcher `{}` deleted successfully", watcher_id)),
            "Watcher deletion failed",
        )?
        .data
        .ok_or_else(|| {
            HttpClientError::DataResponseFailure(
                reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                String::from("Watcher data not returned after deletion"),
            )
        })
    }

    pub fn ping_watcher(&self, id: &str, print: bool) -> Result<String, HttpClientError> {
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

        let response = self.send_request_with_team(
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


    pub fn get_staged_samples(
        &self,
        id: Option<String>,
        print: bool,
    ) -> Result<Vec<StagedSample>, HttpClientError> {
        let url = self.build_request_url(Route::TeamStagedSamplesList, &[("id", id)]);

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
                    String::from("No watcher data returned"),
                )
            })?;

        if print {
            for staged_sample in &staged_samples {
                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    staged_sample.id,
                    staged_sample.sample_id,
                    staged_sample.database,
                    staged_sample.project,
                    staged_sample.pipeline,
                    staged_sample.date,
                );
            }
        }

        Ok(staged_samples)
    }

    pub fn delete_staged_sample(&self, id: &str) -> Result<StagedSample, HttpClientError> {
        let response = self.send_request_with_team(
            self.client
                .delete(&format!(
                    "{}/{}",
                    self.routes.url(Route::TeamStagedSamplesDelete),
                    id
                ))
        )?;

        self.handle_response::<DeleteStagedSampleResponse>(
            response,
            Some(&format!("Staged sample `{}` deleted successfully", id)),
            "Staged sample deletion failed",
        )?
        .data
        .ok_or_else(|| {
            HttpClientError::DataResponseFailure(
                reqwest::StatusCode::INTERNAL_SERVER_ERROR,
                String::from("Stage sample data not returned after deletion"),
            )
        })
    }

    // pub fn upload_models(
    //     &self,
    //     models: &[Cerebro],
    //     team_name: &str,
    //     project_name: &str,
    //     db_name: Option<&String>,
    // ) -> Result<(), HttpClientError> {
    //     let urls = self.get_database_and_project_queries(
    //         &self.routes.url(Route::DataCerebroInsertModel),
    //         team_name,
    //         Some(project_name),
    //         db_name,
    //     )?;

    //     for url in &urls {
    //         for model in models {
    //             if model.sample.id.is_empty() {
    //                 return Err(HttpClientError::ModelSampleIdentifierEmpty);
    //             }

    //             let response = self
    //                 .client
    //                 .post(url)
    //                 .header(AUTHORIZATION, self.get_bearer_token(None))
    //                 .json(model)
    //                 .send()?;

    //             self.handle_response::<serde_json::Value>(
    //                 response,
    //                 &format!(
    //                     "Model for sample library {} uploaded successfully",
    //                     model.sample.id
    //                 ),
    //                 "Upload failed",
    //             )?;
    //         }
    //     }
    //     Ok(())
    // }
    // pub fn taxa_summary(
    //     &self, 
    //     team_name: &str, 
    //     project_name: &str, 
    //     db_name: Option<&String>, 
    //     filter_config: Option<&PathBuf>,
    //     run_ids: Option<Vec<String>>, 
    //     sample_ids: Option<Vec<String>>, 
    //     workflow_ids: Option<Vec<String>>, 
    //     workflow_names: Option<Vec<String>>,
    //     output: &PathBuf
    // ) -> Result<(), HttpClientError> {

    //     let run_ids = run_ids.unwrap_or(Vec::new());
    //     let sample_ids = sample_ids.unwrap_or(Vec::new());
    //     let workflow_ids = workflow_ids.unwrap_or(Vec::new());
    //     let workflow_names = workflow_names.unwrap_or(Vec::new());

    //     let taxa_summary_schema = TaxaSummarySchema {
    //         run_ids: run_ids.clone(),
    //         sample_ids: sample_ids.clone(),
    //         workflow_ids: workflow_ids.clone(),
    //         workflow_names: workflow_names.clone(),
    //         filter_config: match filter_config { 
    //             Some(path) => TaxonFilterConfig::from_path(&path).map_err(|err| HttpClientError::DeserializeFilter(err))?, 
    //             None => TaxonFilterConfig::default() 
    //         }
    //     };

    //     let urls = self.get_database_and_project_queries(&self.routes.data_cerebro_taxa_summary, team_name, Some(project_name), db_name)?;

    //     if urls.len() > 1 {
    //         log::warn!("Project `{}` exists for multiple databases belonging to team `{}`", &project_name, &team_name);
    //         log::warn!("Fetching data for all projects - otherwise, specify the index of a specific database with `--db-index`:");
    //         for (i, url) in urls.iter().enumerate() {
    //             log::warn!("Index: {i} @ {url}");
    //         }
    //     }

    //     for (i, url) in urls.iter().enumerate() {
    //         log::info!(
    //             "Taxa summary query: project={} team={} run_ids={:?} sample_ids={:?} workflow_ids={:?} workflow_names={:?}", 
    //             &project_name, &team_name, &run_ids, &sample_ids, &workflow_ids, &workflow_names
    //         );

    //         let response = self.client.post(format!("{}&csv=true", url))
    //             .header(AUTHORIZATION, self.get_token_bearer(None))
    //             .json(&taxa_summary_schema)
    //             .send()?;

    //         let status = response.status();
    
    //         match status.is_success() {
    //             true => {
    //                 log::info!("Taxa summary retrieved for project {} of team {}", &project_name, &team_name);

    //                 let data_response: TaxaSummaryDataResponse = response.json().map_err(|_| {
    //                     HttpClientError::ResponseFailure(
    //                         status, 
    //                         String::from("failed to obtain data for taxa summary")
    //                     )
    //                 })?;
                    
    //                 let output_name = match i { 0 => output.clone(), _ => output.with_extension(format!("{}", &i)) };

    //                 let mut file = File::create(&output_name).unwrap();
    //                 write!(file, "{}", data_response.data.csv).unwrap();
    //             },
    //             false => {
    //                 let error_response: TaxaSummaryDataResponse = response.json().map_err(|_| {
    //                     HttpClientError::ResponseFailure(
    //                         status, 
    //                         String::from("failed to make request for taxa summary")
    //                     )
    //                 })?;
    //                 return Err(HttpClientError::ResponseFailure(
    //                     status, 
    //                     error_response.message
    //                 ))
    //             }
    //         };
    //     }
        
    //     Ok(())
    // }
    // pub fn qc_summary(
    //     &self, 
    //     team_name: &str, 
    //     project_name: &str, 
    //     db_name: Option<&String>, 
    //     cerebro_ids: Option<Vec<String>>, 
    //     sample_ids: Option<Vec<String>>, 
    //     ercc_pg: Option<f64>,
    //     output: &PathBuf
    // ) -> Result<(), HttpClientError> {

    //     let cerebro_ids = cerebro_ids.unwrap_or(Vec::new());
    //     let sample_ids = sample_ids.unwrap_or(Vec::new());

    //     let qc_summary_schema = SampleSummaryQcSchema {
    //         cerebro_ids: cerebro_ids.clone(),
    //         sample_ids: sample_ids.clone()
    //     };

    //     let urls = self.get_database_and_project_queries(&self.routes.data_cerebro_taxa_summary, team_name, Some(project_name), db_name)?;

    //     if urls.len() > 1 {
    //         log::warn!("Project `{}` exists for multiple databases belonging to team `{}`", &project_name, &team_name);
    //         log::warn!("Fetching data for all projects - otherwise, specify the index of a specific database with `--db-index`:");
    //         for (i, url) in urls.iter().enumerate() {
    //             log::warn!("Index: {i} @ {url}");
    //         }
    //     }

    //     for (i, url) in urls.iter().enumerate() {
    //         log::info!(
    //             "Quality summary query: project={} team={} cerebro_ids={:?} sample_ids={:?} ercc_pg={:?}", 
    //             &project_name, &team_name, &cerebro_ids, &sample_ids, &ercc_pg
    //         );

    //         let response = self.client.post(format!("{}&csv=true{}", url, match ercc_pg { Some(pg) => format!("&ercc={:.2}", pg), None => String::new() } ))
    //             .header(AUTHORIZATION, self.get_token_bearer(None))
    //             .json(&qc_summary_schema)
    //             .send()?;

    //         let status = response.status();
    
    //         match status.is_success() {
    //             true => {
    //                 log::info!("QC summary retrieved for project `{}` of team `{}`", &project_name, &team_name);

    //                 let data_response: TaxaSummaryDataResponse = response.json().map_err(|_| {
    //                     HttpClientError::ResponseFailure(
    //                         status, 
    //                         String::from("failed to obtain data for quality summary")
    //                     )
    //                 })?;
                    
    //                 let output_name = match i { 0 => output.clone(), _ => output.with_extension(format!("{}", &i)) };

    //                 let mut file = File::create(&output_name).unwrap();
    //                 write!(file, "{}", data_response.data.csv).unwrap();
    //             },
    //             false => {
    //                 let error_response: TaxaSummaryDataResponse = response.json().map_err(|_| {
    //                     HttpClientError::ResponseFailure(
    //                         status, 
    //                         String::from("failed to make request for quality summary")
    //                     )
    //                 })?;
    //                 return Err(HttpClientError::ResponseFailure(status, error_response.message))
    //             }
    //         };
    //     }
        
    //     Ok(())
    // }
    // pub fn get_database(&self, team_name: &str, db_name: &str) -> Result<TeamDatabase, HttpClientError> {

    //     let url = format!("{}?name={}", &self.routes.data_user_self_teams, team_name);
        
    //     // Request data on a team project for insertion of new data
    //     // log::info!("Getting user team for database ({})", &url);

    //     let response = self.client.get(url)
    //         .header(AUTHORIZATION, self.get_token_bearer(None))
    //         .send()?;
        
    //     let status = response.status();

    //     let team = match status.is_success() {
    //         true => {
    //             let team_response: UserSelfTeamResponse = response.json()?;
    //             team_response.data.team
    //         }
    //         false => {
    //             let error_response: ErrorResponse = response.json().map_err(|_| {
    //                 HttpClientError::ResponseFailure(
    //                     status, 
    //                     String::from("failed to make request")
    //                 )
    //             })?;
    //             return Err(HttpClientError::ResponseFailure(
    //                 status, 
    //                 error_response.message
    //             ))
    //         }
    //     };

    //     get_database_by_name(&team.databases, db_name)
    // }
    // pub fn get_database_and_project_queries(&self, route: &str, team_name: &str, project_name: Option<&str>, db_name: Option<&String>) -> Result<Vec<String>, HttpClientError> {

    //     let url = format!("{}?name={}", &self.routes.data_user_self_teams, team_name);
        
    //     // Request data on a team project for insertion of new data
    //     // log::info!("Getting user team for database verification ({})", &url);

    //     let response = self.client.get(url)
    //         .header(AUTHORIZATION, self.get_token_bearer(None))
    //         .send()?;
        
    //     let status = response.status();

    //     let team = match status.is_success() {
    //         true => {
    //             let team_response: UserSelfTeamResponse = response.json()?;
    //             team_response.data.team
    //         }
    //         false => {
    //             let error_response: ErrorResponse = response.json().map_err(|_| {
    //                 HttpClientError::ResponseFailure(
    //                     status, 
    //                     String::from("failed to make request")
    //                 )
    //             })?;
    //             return Err(HttpClientError::ResponseFailure(
    //                 status, 
    //                 error_response.message
    //             ))
    //         }
    //     };

    //     if team.databases.is_empty() {
    //         log::error!("No team databases exist - this is unusual, please contact system administrator");
    //         return Err(HttpClientError::TeamDatabasesNotFound)
    //     }

    //     let mut urls = Vec::new();
    //     for database in &team.databases {

    //         // If specific database name requested, check if this is it,
    //         // otherwise use all databases for this team for data insertion
    //         if let Some(name) = db_name {
    //             if &database.name != name {
    //                 log::info!("Requested database ({}) - skipping team database ({})", &name, &database.name);
    //                 continue
    //             }
    //         }
    //         match project_name {
    //             Some(name) => {
    //                 let project = get_project_by_name(&database.projects, name)?;
    //                 urls.push(format!("{}?db={}&project={}", &route, database.id, project.id))
    //             },
    //             None => {
    //                 urls.push(format!("{}?db={}", &route, database.id))
    //             }
    //         }
    //     }    

    //     Ok(urls)
    // }
}


// fn get_project_by_name(projects: &Vec<ProjectCollection>, project_name: &str) -> Result<ProjectCollection, HttpClientError>   {
//     let matches: Vec<&ProjectCollection> = projects.into_iter().filter(|x| x.name == project_name).collect();

//     if matches.len() > 0 {
//         Ok(matches[0].to_owned())
//     } else {
//         let valid_project_name_string = projects.iter()
//             .map(|project| project.name.to_owned()) // replace `field_name` with the actual field name
//             .collect::<Vec<String>>()
//             .join(", ");
    
//         Err(HttpClientError::InsertModelProjectParameter(valid_project_name_string))
//     }
// }


// fn get_database_by_name(databases: &Vec<TeamDatabase>, db_name: &str) -> Result<TeamDatabase, HttpClientError>   {
//     let matches: Vec<&TeamDatabase> = databases.into_iter().filter(|x| x.name == db_name).collect();

//     if matches.len() > 0 {
//         Ok(matches[0].to_owned())
//     } else {
//         let valid_database_name_string = databases.iter()
//             .map(|project| project.name.to_owned()) // replace `field_name` with the actual field name
//             .collect::<Vec<String>>()
//             .join(", ");
//         Err(HttpClientError::InsertModelDatabaseParameter(valid_database_name_string))
//     }
// }

