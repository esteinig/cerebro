use std::fs::File;
use std::io::Write;
use anyhow::Result;
use cerebro_model::api::files::response::ListFilesResponse;
use cerebro_model::api::files::response::RegisterFileResponse;
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

#[derive(Debug, Clone)]
pub struct CerebroRoutes {
    pub server_status: String,
    pub auth_server_status: String,
    pub auth_login_user: String,
    pub auth_refresh_token: String,

    pub data_user_self: String,
    pub data_user_self_teams: String,

    pub data_cerebro_insert_model: String,
    pub data_cerebro_taxa_summary: String,

    pub team_project_create: String,

    pub team_files_register: String,
    pub team_files_list: String
}
impl CerebroRoutes {
    pub fn new(url: &str) -> Self  {
        Self {
            server_status: format!("{}/status", url),  // unauthenticated
            auth_server_status: format!("{}/auth/status", url),
            auth_login_user: format!("{}/auth/login", url),
            auth_refresh_token: format!("{}/auth/refresh", url),
            data_user_self: format!("{}/users/self", url),
            data_user_self_teams: format!("{}/users/self/teams", url),
            data_cerebro_insert_model: format!("{}/cerebro", url),
            data_cerebro_taxa_summary: format!("{}/cerebro/taxa/summary", url),
            team_project_create: format!("{}/teams/project", url),

            team_files_register: format!("{}/files/register", url),
            team_files_list: format!("{}/files", url),
        }
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AuthLoginTokenFile {
    pub access_token: String,
    pub refresh_token: String
}
impl AuthLoginTokenFile {
    pub fn from(response_success: &AuthLoginResponseSuccess) -> Self {
        Self { access_token: response_success.access_token.to_owned(), refresh_token: response_success.refresh_token.to_owned() }
    }
}

#[derive(Debug, Clone)]
pub struct CerebroClient {
    pub url: String,
    pub client: Client,
    pub routes: CerebroRoutes,
    pub token_file: Option<PathBuf>,
    pub token_data: Option<AuthLoginTokenFile>,
}
impl CerebroClient {
    pub fn new(url: &str, token: &Option<String>, login: bool, danger_accept_invalid_certs: bool, token_file: &Option<PathBuf>) -> Result<Self, HttpClientError> {
        let url_clean = url.trim_end_matches("/");

        let token_data = match login {
            true => None,
            false => {
                match token_file {
                    Some(token_path) => {
                        let token_data = match token_path.exists() {
                            true => {
                                let token_str = std::fs::read_to_string(&token_path)
                                    .map_err(|err| HttpClientError::IOFailure(err))?;
                                let token_data: AuthLoginTokenFile = serde_json::from_str(&token_str)?;
                                Some(token_data)
                            },
                            false => {
                                log::error!("Token file path does not exist");
                                std::process::exit(1);
                            }
                        };
                        token_data
                    },
                    None => {
                        let token_data = match token {
                            Some(token) => {
                                Some(AuthLoginTokenFile {
                                    access_token: token.to_string(),
                                    refresh_token: "".to_string()
                                })
                            },
                            None => {
                                log::error!("Failed to obtain token from input (cerebro --token) or environmental variable (CEREBRO_API_TOKEN)");
                                std::process::exit(1);
                            }
                        };
                        token_data
                    }
                } 
            }
        };


        let client = reqwest::blocking::Client::builder().danger_accept_invalid_certs(danger_accept_invalid_certs).build()?;
        
        Ok(Self { 
            url: url_clean.to_owned(), 
            routes: CerebroRoutes::new(url_clean), 
            token_file: token_file.to_owned(),
            token_data,
            client
        })
    }
    // We use authorization headers for the client instead of cookies
    // as we set strict cookie policies that force same-site origin
    // requests that cannot be fulfilled by the client
    pub fn get_token_bearer(&self, token: Option<String>) -> String {
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
    // Ping server status
    pub fn ping_status(&self) -> Result<(), HttpClientError> {

        let response = self.client
            .get(&self.routes.server_status)
            .send()?;

        if response.status().is_success() {
           log::info!("Cerebro API status: ok");
        } else if response.status().is_server_error() {
            return Err(HttpClientError::PingServer(format!("{:?}", response.status()), "server error".to_string()))
        } else {
            return Err(HttpClientError::PingServer(format!("{:?}", response.status()), "failed to ping server".to_string()))
        }
        
    
        Ok(())
    }
    // Ping server status routes - currently running a single
    // server. Will separate into distinct authentication
    // and database services
    pub fn ping_servers(&self) -> Result<(), HttpClientError> {

        let response = self.client
            .get(&self.routes.auth_server_status)
            .header(AUTHORIZATION, self.get_token_bearer(None))
            .send()?;
    
        if !response.status().is_success() {
            return Err(HttpClientError::PingServer(format!("{:?}", response.status()), "failed to ping server - are you logged in?".to_string()))
        };

        log::info!("Cerebro API status: ok");
        
        Ok(())
    }
    // Login user and save tokens
    pub fn login_user(&self, email: &str, password: &Option<String>) -> Result<(), HttpClientError> {

        let password = match password {
            None => prompt_password(format!("Password [{}]: ", &email)).map_err(|_| HttpClientError::PasswordInput)?,
            Some(pwd) => pwd .to_string()
        };

        let login_schema = AuthLoginSchema {
            email: email.to_owned(),
            password
        };

        let response = self.client.post(&self.routes.auth_login_user)
            .json(&login_schema)
            .send()?;
        
        let status = response.status();

        match status.is_success() {
            true => {
                let login_response: AuthLoginResponseSuccess = response.json()?;
                let access_token = Some(login_response.access_token.clone());
                let response = self.client.get(&self.routes.data_user_self)
                    .header(AUTHORIZATION, self.get_token_bearer(access_token))
                    .send()?;

                let status = response.status();

                if status.is_success() {
                    let user_response: UserSelfResponse = response.json()?;
                    log::info!("Login successful. Welcome back, {}!", &user_response.data.user.name)
                } else {
                    log::error!("Login successful, but could not get user data");
                    let user_error_response: ErrorResponse = response.json().map_err(|_| {
                        HttpClientError::LoginUserResponseFailure(format!("{}", status), String::from("failed to make request"))
                    })?;

                    return Err(HttpClientError::LoginUserResponseFailure(
                        format!("{}", status), format!("{}", user_error_response.message)
                    ))
                }

                match &self.token_file {
                    Some(path) => {
                        log::info!("Token will be written to file: {}", path.display());

                        std::fs::write(
                            path,
                            serde_json::to_string_pretty(
                                &AuthLoginTokenFile::from(&login_response)
                            ).expect("Failed to convert token data to string"),
                        ).map_err(|err| HttpClientError::IOFailure(err))?;

                    },
                    None => {
                        print!("{}", &login_response.access_token);
                    }
                }
            },
            false => {
                let error_response: ErrorResponse = response.json()?;
                log::error!("Login failed: {} ({})", &error_response.message, &status);

                return Err(HttpClientError::LoginUserResponseFailure(
                    format!("{}", status),  error_response.message
                ))
            }
        };
        Ok(())
    }
    pub fn upload_models(&self, models: &Vec<Cerebro>, team_name: &str, project_name: &str, db_name: Option<&String>) -> Result<(), HttpClientError> {

        let urls = self.get_database_and_project_queries(&self.routes.data_cerebro_insert_model, team_name, Some(project_name), db_name)?;

        log::info!("Cerebro model upload to project `{}` for team `{}`", &project_name, &team_name);

        if urls.len() > 1 {
            log::warn!("Project `{}` exists for multiple databases belonging to team `{}`", &project_name, &team_name);
            log::warn!("Inserting data for all projects - otherwise, specify the index of a specific database with `--db-index`:");
            for (i, url) in urls.iter().enumerate() {
                log::warn!("Index: {i} @ {url}");
            }
        }

        for url in &urls {

            for model in models {
                if model.sample.id.is_empty() {
                    log::error!("Model sample identifier is an empty string - this is not allowed");
                    return Err(HttpClientError::ModelSampleIdentifierEmpty)
                }
    
                log::info!("Requesting upload of model for sample library {} (tags: {}, workflow: {})", model.sample.id, model.sample.tags.join(" "), model.workflow.id);
    
                let response = self.client.post(url)
                    .header(AUTHORIZATION, self.get_token_bearer(None))
                    .json(model)
                    .send()?;
    
                let status = response.status();
        
                match status.is_success() {
                    true => log::info!("Cerebro model uploaded successfully to project `{}` ({})", &project_name, &team_name),
                    false => {
                        let error_response: ErrorResponse = response.json().map_err(|_| {
                            HttpClientError::InsertModelResponseFailure(format!("{}", status), String::from("failed to make request"))
                        })?;

                        log::error!("Upload failed: {}", & error_response.message);
                        continue; // return Err(HttpClientError::InsertModelResponseFailure(format!("{}", status), error_response.message))
                        
                    }
                };
            }
            
        }
        Ok(())
    }
    pub fn create_project(&self, team_name: &str, db_name: &str, project_name: &str, project_description: &str) -> Result<(), HttpClientError> {

        let register_project_schema = RegisterProjectSchema {
            project_name: project_name.to_string(),
            project_description: project_description.to_string(),
            project_mongo_name: project_name.split_whitespace().join("_").to_lowercase()
        };

        let response = self.client.post(format!("{}?team_name={}&db_name={}", self.routes.team_project_create, team_name, db_name))
            .header(AUTHORIZATION, self.get_token_bearer(None))
            .json(&register_project_schema)
            .send()?;

        let status = response.status();

        match status.is_success() {
            true => log::info!("Team project created: {} - {} - {}", &team_name, &db_name, &project_name),
            false => {
                let error_response: ErrorResponse = response.json().map_err(|_| {
                    HttpClientError::InsertModelResponseFailure(format!("{}", status), String::from("failed to make request"))
                })?;
                log::error!("Project creation failed: {}", &status);
                return Err(HttpClientError::InsertModelResponseFailure(format!("{}", status), error_response.message))
            }
        };
        Ok(())
    }
    pub fn register_file(&self, register_file_schema: RegisterFileSchema, team_name: &str, db_name: &str) -> Result<(), HttpClientError> {

        let urls = self.get_database_and_project_queries(
            &self.routes.team_files_register, team_name, None, Some(&db_name.to_owned())
        )?;

        for url in &urls {
            log::info!(
                "Registering file {} in database {} for team {}", &register_file_schema.name, &db_name, &team_name
            );

            let response = self.client.post(url)
                .header(AUTHORIZATION, self.get_token_bearer(None))
                .json(&register_file_schema)
                .send()?;

            let status = response.status();
    
            match status.is_success() {
                true => {
                    log::info!("File registered successfully!");
                },
                false => {
                    let error_response: RegisterFileResponse = response.json().map_err(|_| {
                        HttpClientError::RegisterFileResponseFailure(format!("{}", status), String::from("failed to extract response data"))
                    })?;
                    log::error!("Data retrieval failed: {}", &status);
                    return Err(HttpClientError::RegisterFileResponseFailure(format!("{}", status), error_response.message))
                }
            }
        }       
        Ok(())
    }
    pub fn list_files(&self, team_name: &str, db_name: &str, page: u32, limit: u32) -> Result<(), HttpClientError> {

        let urls = self.get_database_and_project_queries(
            &self.routes.team_files_list, team_name, None, Some(&db_name.to_owned())
        )?;

        for url in &urls {
            log::info!(
                "Listing {limit} files on page {page} in database {} for team {}", &db_name, &team_name, 
            );

            let response = self.client.get(
                format!("{url}&page={page}&limit={limit}")
            )
                .header(AUTHORIZATION, self.get_token_bearer(None))
                .send()?;

            let status = response.status();
    
            match status.is_success() {
                true => {
                    let response: ListFilesResponse = response.json().map_err(|_| {
                        HttpClientError::ListFilesResponseFailure(format!("{}", status), String::from("failed to extract response data"))
                    })?;
                    for file in response.data {
                        println!("{}\t{}\t{}\t{}\t{:.0} MB\t{}", file.date, file.watcher.name, file.watcher.location, file.fid, file.size_mb(), file.name)
                    }
                },
                false => {
                    let error_response: ListFilesResponse = response.json().map_err(|_| {
                        HttpClientError::ListFilesResponseFailure(format!("{}", status), String::from("failed to extract response data"))
                    })?;
                    log::error!("Data retrieval failed: {}", &status);
                    return Err(HttpClientError::ListFilesResponseFailure(format!("{}", status), error_response.message))
                }
            }
        }       
        Ok(())
    }
    pub fn taxa_summary(
        &self, 
        team_name: &str, 
        project_name: &str, 
        db_name: Option<&String>, 
        filter_config: Option<&PathBuf>,
        run_ids: Option<Vec<String>>, 
        sample_ids: Option<Vec<String>>, 
        workflow_ids: Option<Vec<String>>, 
        workflow_names: Option<Vec<String>>,
        output: &PathBuf
    ) -> Result<(), HttpClientError> {

        let run_ids = run_ids.unwrap_or(Vec::new());
        let sample_ids = sample_ids.unwrap_or(Vec::new());
        let workflow_ids = workflow_ids.unwrap_or(Vec::new());
        let workflow_names = workflow_names.unwrap_or(Vec::new());

        let taxa_summary_schema = TaxaSummarySchema {
            run_ids: run_ids.clone(),
            sample_ids: sample_ids.clone(),
            workflow_ids: workflow_ids.clone(),
            workflow_names: workflow_names.clone(),
            filter_config: match filter_config { 
                Some(path) => TaxonFilterConfig::from_path(&path).map_err(|err| HttpClientError::DeserializeFilter(err))?, 
                None => TaxonFilterConfig::default() 
            }
        };

        let urls = self.get_database_and_project_queries(&self.routes.data_cerebro_taxa_summary, team_name, Some(project_name), db_name)?;

        if urls.len() > 1 {
            log::warn!("Project `{}` exists for multiple databases belonging to team `{}`", &project_name, &team_name);
            log::warn!("Fetching data for all projects - otherwise, specify the index of a specific database with `--db-index`:");
            for (i, url) in urls.iter().enumerate() {
                log::warn!("Index: {i} @ {url}");
            }
        }

        for (i, url) in urls.iter().enumerate() {
            log::info!(
                "Taxa summary query: project={} team={} run_ids={:?} sample_ids={:?} workflow_ids={:?} workflow_names={:?}", 
                &project_name, &team_name, &run_ids, &sample_ids, &workflow_ids, &workflow_names
            );

            let response = self.client.post(format!("{}&csv=true", url))
                .header(AUTHORIZATION, self.get_token_bearer(None))
                .json(&taxa_summary_schema)
                .send()?;

            let status = response.status();
    
            match status.is_success() {
                true => {
                    log::info!("Taxa summary retrieved for project {} of team {}", &project_name, &team_name);

                    let data_response: TaxaSummaryDataResponse = response.json().map_err(|_| {
                        HttpClientError::TaxaSummaryDataResponseFailure(format!("{}", status), String::from("failed to obtain data for taxa summary"))
                    })?;
                    
                    let output_name = match i { 0 => output.clone(), _ => output.with_extension(format!("{}", &i)) };

                    let mut file = File::create(&output_name).unwrap();
                    write!(file, "{}", data_response.data.csv).unwrap();
                },
                false => {
                    let error_response: TaxaSummaryDataResponse = response.json().map_err(|_| {
                        HttpClientError::TaxaSummaryDataResponseFailure(format!("{}", status), String::from("failed to make request for taxa summary"))
                    })?;
                    log::error!("Data retrieval failed: {}", &status);
                    return Err(HttpClientError::InsertModelResponseFailure(format!("{}", status), error_response.message))
                }
            };
        }
        
        Ok(())
    }
    pub fn qc_summary(
        &self, 
        team_name: &str, 
        project_name: &str, 
        db_name: Option<&String>, 
        cerebro_ids: Option<Vec<String>>, 
        sample_ids: Option<Vec<String>>, 
        ercc_pg: Option<f64>,
        output: &PathBuf
    ) -> Result<(), HttpClientError> {

        let cerebro_ids = cerebro_ids.unwrap_or(Vec::new());
        let sample_ids = sample_ids.unwrap_or(Vec::new());

        let qc_summary_schema = SampleSummaryQcSchema {
            cerebro_ids: cerebro_ids.clone(),
            sample_ids: sample_ids.clone()
        };

        let urls = self.get_database_and_project_queries(&self.routes.data_cerebro_taxa_summary, team_name, Some(project_name), db_name)?;

        if urls.len() > 1 {
            log::warn!("Project `{}` exists for multiple databases belonging to team `{}`", &project_name, &team_name);
            log::warn!("Fetching data for all projects - otherwise, specify the index of a specific database with `--db-index`:");
            for (i, url) in urls.iter().enumerate() {
                log::warn!("Index: {i} @ {url}");
            }
        }

        for (i, url) in urls.iter().enumerate() {
            log::info!(
                "Quality summary query: project={} team={} cerebro_ids={:?} sample_ids={:?} ercc_pg={:?}", 
                &project_name, &team_name, &cerebro_ids, &sample_ids, &ercc_pg
            );

            let response = self.client.post(format!("{}&csv=true{}", url, match ercc_pg { Some(pg) => format!("&ercc={:.2}", pg), None => String::new() } ))
                .header(AUTHORIZATION, self.get_token_bearer(None))
                .json(&qc_summary_schema)
                .send()?;

            let status = response.status();
    
            match status.is_success() {
                true => {
                    log::info!("QC summary retrieved for project `{}` of team `{}`", &project_name, &team_name);

                    let data_response: TaxaSummaryDataResponse = response.json().map_err(|_| {
                        HttpClientError::TaxaSummaryDataResponseFailure(format!("{}", status), String::from("failed to obtain data for QC summary"))
                    })?;
                    
                    let output_name = match i { 0 => output.clone(), _ => output.with_extension(format!("{}", &i)) };

                    let mut file = File::create(&output_name).unwrap();
                    write!(file, "{}", data_response.data.csv).unwrap();
                },
                false => {
                    let error_response: TaxaSummaryDataResponse = response.json().map_err(|_| {
                        HttpClientError::TaxaSummaryDataResponseFailure(format!("{}", status), String::from("failed to make request for QC summary"))
                    })?;
                    log::error!("Data retrieval failed: {}", &status);
                    return Err(HttpClientError::InsertModelResponseFailure(format!("{}", status), error_response.message))
                }
            };
        }
        
        Ok(())
    }
    pub fn get_database(&self, team_name: &str, db_name: &str) -> Result<TeamDatabase, HttpClientError> {

        let url = format!("{}?name={}", &self.routes.data_user_self_teams, team_name);
        
        // Request data on a team project for insertion of new data
        log::info!("Getting user team for database ({})", &url);

        let response = self.client.get(url)
            .header(AUTHORIZATION, self.get_token_bearer(None))
            .send()?;
        
        let status = response.status();

        let team = match status.is_success() {
            true => {
                let team_response: UserSelfTeamResponse = response.json()?;
                team_response.data.team
            }
            false => {
                let error_response: ErrorResponse = response.json().map_err(|_| {
                    HttpClientError::UserTeamsResponseFailure(format!("{}", status), String::from("failed to make request"))
                })?;
                log::error!("Failed to get team data ({})", &status);
                return Err(HttpClientError::UserTeamsResponseFailure(format!("{}", status), error_response.message))
            }
        };

        get_database_by_name(&team.databases, db_name)
    }
    pub fn get_database_and_project_queries(&self, route: &str, team_name: &str, project_name: Option<&str>, db_name: Option<&String>) -> Result<Vec<String>, HttpClientError> {

        let url = format!("{}?name={}", &self.routes.data_user_self_teams, team_name);
        
        // Request data on a team project for insertion of new data
        log::info!("Getting user team for database verification ({})", &url);

        let response = self.client.get(url)
            .header(AUTHORIZATION, self.get_token_bearer(None))
            .send()?;
        
        let status = response.status();

        let team = match status.is_success() {
            true => {
                let team_response: UserSelfTeamResponse = response.json()?;
                team_response.data.team
            }
            false => {
                let error_response: ErrorResponse = response.json().map_err(|_| {
                    HttpClientError::UserTeamsResponseFailure(format!("{}", status), String::from("failed to make request"))
                })?;
                log::error!("Failed to get team data ({})", &status);
                return Err(HttpClientError::UserTeamsResponseFailure(format!("{}", status), error_response.message))
            }
        };

        if team.databases.is_empty() {
            log::error!("No team databases exist - this is unusual, please contact system administrator");
            return Err(HttpClientError::TeamDatabasesNotFound)
        }

        let mut urls = Vec::new();
        for database in &team.databases {

            // If specific database name requested, check if this is it,
            // otherwise use all databases for this team for data insertion
            if let Some(name) = db_name {
                if &database.name != name {
                    log::info!("Requested database ({}) - skipping team database ({})", &name, &database.name);
                    continue
                }
            }
            match project_name {
                Some(name) => {
                    let project = get_project_by_name(&database.projects, name)?;
                    urls.push(format!("{}?db={}&project={}", &route, database.id, project.id))
                },
                None => {
                    urls.push(format!("{}?db={}", &route, database.id))
                }
            }
        }    

        Ok(urls)
    }
}


fn get_project_by_name(projects: &Vec<ProjectCollection>, project_name: &str) -> Result<ProjectCollection, HttpClientError>   {
    let matches: Vec<&ProjectCollection> = projects.into_iter().filter(|x| x.name == project_name).collect();

    if matches.len() > 0 {
        Ok(matches[0].to_owned())
    } else {
        let valid_project_name_string = projects.iter()
            .map(|project| project.name.to_owned()) // replace `field_name` with the actual field name
            .collect::<Vec<String>>()
            .join(", ");
    
        Err(HttpClientError::InsertModelProjectParameter(valid_project_name_string))
    }
}


fn get_database_by_name(databases: &Vec<TeamDatabase>, db_name: &str) -> Result<TeamDatabase, HttpClientError>   {
    let matches: Vec<&TeamDatabase> = databases.into_iter().filter(|x| x.name == db_name).collect();

    if matches.len() > 0 {
        Ok(matches[0].to_owned())
    } else {
        let valid_database_name_string = databases.iter()
            .map(|project| project.name.to_owned()) // replace `field_name` with the actual field name
            .collect::<Vec<String>>()
            .join(", ");
        Err(HttpClientError::InsertModelDatabaseParameter(valid_database_name_string))
    }
}
