use colored::Colorize;
use thiserror::Error;
use rust_embed::EmbeddedFile;
use crate::{utils::CRATE_VERSION, api::config::SecurityComponentsConfig};
use bcrypt::{DEFAULT_COST, hash};
use std::{path::PathBuf, io::Write};
use rcgen::generate_simple_self_signed;
use serde::{Deserialize, Serialize, Deserializer};
use base64::{engine::general_purpose, Engine as _};
use rsa::{RsaPrivateKey, RsaPublicKey, pkcs8::EncodePrivateKey, pkcs8::EncodePublicKey};
use crate::{tools::modules::password::hash_password, api::{self, config::{SmtpConfig, TokenEncryptionConfig, ConfigError}}};
use crate::stack::assets::{MongoInitTemplates, DockerTemplates, TraefikTemplates, CerebroTemplates};

#[derive(Error, Debug)]
pub enum StackConfigError {
    #[error("failed to extract project name from output directory path")]
    ProjectNameInvalid,
    #[error("failed to find embedded asset file: {0}")]
    EmbeddedFileNotFound(String),
    #[error("failed to get absolute path for output directory: {0}")]
    OutdirAbsolutePathInvalid(String),
    #[error("failed to open config file")]
    ConfigFileInputInvalid,
    #[error("failed to create config file")]
    ConfigFileOutputInvalid,
    #[error("failed to deserialize config file")]
    ConfigFileNotDeserialized(#[from] toml::de::Error),
    #[error("failed to serialize config file")]
    ConfigFileNotSerialized(#[from] toml::ser::Error),
    #[error("failed to create output path for embedded file")]
    EmbeddedFileOutputPathInvalid(#[source] std::io::Error),
    #[error("failed to write embedded file to output path")]
    EmbeddedFileOutputNotWritten(#[source] std::io::Error),
    #[error("failed to create stack tree directory")]
    StackTreeDirectoryNotCreated(#[source] std::io::Error),
    #[error("failed to create stack template directory")]
    TemplateDirectoryNotCreated(#[source] std::io::Error),
    #[error("failed to register template")]
    TemplateNotRegistered(#[from] handlebars::TemplateError),
    #[error("failed to render template")]
    TemplateNotRendered(#[from] handlebars::RenderError),
    #[error("failed to hash traefik dashboard password")]
    TraefikDashboardPasswordNotHashed(#[from] bcrypt::BcryptError),
    #[error("failed to hash cerebro database password")]
    CerebroDatabasePasswordNotHashed,
    #[error("failed to create certificate file path")]
    CertificateFilePathInvalid(#[source] std::io::Error),
    #[error("failed to write certificate file")]
    CertificateFileNotWritten(#[source] std::io::Error),
    #[error("failed to create certificate key file path")]
    CertificateKeyFilePathInvalid(#[source] std::io::Error),
    #[error("failed to write certificate key file")]
    CertificateKeyFileNotWritten(#[source] std::io::Error),
    #[error("failed to serialize certificate")]
    CertificateFileNotSerialized(#[from] rcgen::RcgenError),
    #[error("failed to create rendered template file path")]
    TemplateRenderFilePathInvalid(#[source] std::io::Error),
    #[error("failed to write rendered template file")]
    TemplateRenderFileNotWritten(#[source] std::io::Error),
    #[error("failed to generate rsa key")]
    RsaKeyNotGenerated(#[from] rsa::Error),
    #[error("failed to create pem format for private rsa key")]
    RsaPrivateKeyPemNotGenerated(#[from] rsa::pkcs8::Error),
    #[error("failed to create pem format for public rsa key")]
    RsaPublicKeyPemNotGenerated(#[from] rsa::pkcs8::spki::Error),
    #[error("failed to parse template server configuration file")]
    CerebroServerConfigNotParsed(#[from] api::config::ConfigError),
    #[cfg(feature = "libgit")]
    #[error("failed to clone remote git repository")]
    RemoteRepositoryNotCloned(#[source] git2::Error),
    #[cfg(feature = "libgit")]
    #[error("failed to checkout repository")]
    RemoteRepositoryNotCheckedOut(#[source] git2::Error)
}

fn write_embedded_file(embedded_file: &EmbeddedFile, outfile: &PathBuf) -> Result<(), StackConfigError> {
    let mut file = std::fs::File::create(outfile).map_err(|err| StackConfigError::EmbeddedFileOutputPathInvalid(err) )?;
    file.write_all(&embedded_file.data).map_err(|err| StackConfigError::EmbeddedFileOutputNotWritten(err) )
}

pub struct TemplateFiles {
    names: TemplateFileNames,
    paths: TemplateFilePaths,
}
impl TemplateFiles {
    pub fn new(outdir: &PathBuf) -> Self {
        let names = TemplateFileNames::new();
        Self {
            paths: TemplateFilePaths { 
                traefik_dynamic: outdir.join(&names.traefik_dynamic), 
                traefik_static: outdir.join(&names.traefik_static), 
                docker_compose: outdir.join(&names.docker_compose), 
                docker_compose_traefik: outdir.join(&names.docker_compose_traefik),
                docker_server: outdir.join(&names.docker_server), 
                docker_app: outdir.join(&names.docker_app), 
                mongodb_init_sh: outdir.join(&names.mongodb_init_sh), 
                mongodb_init_env: outdir.join(&names.mongodb_init_env),
                cerebro_server: outdir.join(&names.cerebro_server),
                cerebro_app: outdir.join(&names.cerebro_app)
            },
            names,
        }
    }
}

pub struct TemplateFilePaths {
    traefik_dynamic: PathBuf,
    traefik_static: PathBuf,
    docker_compose: PathBuf,
    docker_compose_traefik: PathBuf,
    docker_server: PathBuf,
    docker_app: PathBuf,
    mongodb_init_sh: PathBuf,
    mongodb_init_env: PathBuf,
    cerebro_server: PathBuf,
    cerebro_app: PathBuf,
}

pub struct TemplateFileNames {
    traefik_dynamic: String,
    traefik_static: String,
    docker_compose: String,
    docker_compose_traefik: String,
    docker_server: String,
    docker_app: String,
    mongodb_init_sh: String,
    mongodb_init_env: String,
    cerebro_server: String,
    cerebro_app: String,

}
impl TemplateFileNames {
    pub fn new() -> Self {
        Self {
            traefik_dynamic: String::from("dynamic.yml.hbs"),
            traefik_static: String::from("static.yml.hbs"),
            docker_compose: String::from("docker-compose.yml.hbs"),
            docker_compose_traefik: String::from("docker-compose.traefik.yml.hbs"),
            docker_server: String::from("Dockerfile.server.hbs"),
            docker_app: String::from("Dockerfile.app.hbs"),
            mongodb_init_sh: String::from("mongo-init.sh.hbs"),
            mongodb_init_env: String::from("database.env.hbs"),
            cerebro_server: String::from("server.toml"),
            cerebro_app: String::from("app.env"),
        }
    }
}

// Assets included in the binary to write to 
// configuration directory
pub struct StackAssets {
    mongodb: MongoInitFiles,
    docker: DockerFiles,
    traefik: TraefikFiles,
    cerebro: CerebroFiles,
    templates: TemplateFiles,
    templates_outdir: PathBuf
}
impl StackAssets {
    pub fn new(templates_outdir: &PathBuf) -> Result<StackAssets, StackConfigError> {

        Ok(Self {
            mongodb: MongoInitFiles::new()?,
            docker: DockerFiles::new()?,
            traefik: TraefikFiles::new()?,
            cerebro: CerebroFiles::new()?,
            templates: TemplateFiles::new(&templates_outdir),
            templates_outdir: templates_outdir.clone()
        })
    }
    // Writes the deployment assets for templating to the output directory
    pub fn write_asset_files(&self, traefik_deployment: TraefikDeployment) -> Result<(), StackConfigError> {

        std::fs::create_dir_all(self.templates_outdir.clone()).map_err(|err| StackConfigError::TemplateDirectoryNotCreated(err) )?;

        match traefik_deployment {
            TraefikDeployment::Web => {
                write_embedded_file(&self.traefik.web_dynamic, &self.templates.paths.traefik_dynamic)?;
                write_embedded_file(&self.traefik.web_static, &self.templates.paths.traefik_static)?;
            },
            TraefikDeployment::Localhost => {
                write_embedded_file(&self.traefik.localhost_dynamic, &self.templates.paths.traefik_dynamic)?;
                write_embedded_file(&self.traefik.localhost_static, &self.templates.paths.traefik_static)?;
            }
        }

        write_embedded_file(&self.docker.compose_traefik, &self.templates.paths.docker_compose_traefik)?;
        write_embedded_file(&self.docker.compose, &self.templates.paths.docker_compose)?;
        write_embedded_file(&self.docker.server, &self.templates.paths.docker_server)?;
        write_embedded_file(&self.docker.app, &self.templates.paths.docker_app)?;
        write_embedded_file(&self.mongodb.init_sh, &self.templates.paths.mongodb_init_sh)?;
        write_embedded_file(&self.mongodb.init_env, &self.templates.paths.mongodb_init_env)?;
        write_embedded_file(&self.cerebro.server, &self.templates.paths.cerebro_server)?;
        write_embedded_file(&self.cerebro.app, &self.templates.paths.cerebro_app)?;


        Ok(())

    }
}


pub struct MongoInitFiles {
    init_sh: EmbeddedFile,
    init_env: EmbeddedFile
}
impl MongoInitFiles {
    pub fn new() -> Result<Self, StackConfigError> {
        Ok(Self {
            init_sh: match MongoInitTemplates::get("mongodb/mongo-init.sh") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("mongodb/mongo-init.sh"))) },
            init_env: match MongoInitTemplates::get("mongodb/database.env") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("mongodb/database.env"))) }
        })
    }
}


pub struct TraefikFiles {
    web_static: EmbeddedFile,
    web_dynamic: EmbeddedFile,
    localhost_static: EmbeddedFile,
    localhost_dynamic: EmbeddedFile

}
impl TraefikFiles {
    pub fn new() -> Result<Self, StackConfigError> {
        Ok(Self {
            web_static: match TraefikTemplates::get("traefik/web/static.yml") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("traefik/web/static.yml"))) },
            web_dynamic: match TraefikTemplates::get("traefik/web/dynamic.yml") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("traefik/web/dynamic.yml"))) },
            localhost_static: match TraefikTemplates::get("traefik/localhost/static.yml") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("traefik/localhost/static.yml"))) },
            localhost_dynamic: match TraefikTemplates::get("traefik/localhost/dynamic.yml") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("traefik/localhost/dynamic.yml"))) }
        })
    }
}

pub struct DockerFiles {
    server: EmbeddedFile,
    app: EmbeddedFile,
    compose: EmbeddedFile,
    compose_traefik: EmbeddedFile
}
impl DockerFiles {
    pub fn new() -> Result<Self, StackConfigError> {
        Ok(Self {
            server: match DockerTemplates::get("docker/Dockerfile.server") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("docker/Dockerfile.server"))) },
            app: match DockerTemplates::get("docker/Dockerfile.app") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("docker/Dockerfile.app"))) },
            compose: match DockerTemplates::get("docker/docker-compose.yml") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("docker/docker-compose.yml"))) },
            compose_traefik: match DockerTemplates::get("docker/docker-compose.traefik.yml") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("docker/docker-compose.traefik.yml"))) },
        })
    }
}

pub struct CerebroFiles {
    server: EmbeddedFile,
    app: EmbeddedFile
}
impl CerebroFiles {
    pub fn new() -> Result<Self, StackConfigError> {
        
        Ok(Self {
            server: match CerebroTemplates::get("cerebro/server.toml") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("cerebro/server.toml"))) },
            app: match CerebroTemplates::get("cerebro/app.env") { Some(f) => f, None => return Err(StackConfigError::EmbeddedFileNotFound(String::from("cerebro/app.env"))) },
        })

    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MongoDbConfig {
    pub root_username: String,       
    pub root_password: String,            
    pub admin_username: String,         
    pub admin_password: String,             
    pub cerebro_admin_email: String,           
    pub cerebro_admin_name: String,         
    pub cerebro_admin_password: String, 
    #[serde(skip_deserializing)]
    pub cerebro_admin_password_hashed: String,     
}

#[derive(Debug, Serialize, Clone)]
pub enum TraefikDeployment {
    Web,
    Localhost
}
impl<'de> Deserialize<'de> for TraefikDeployment {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where D: Deserializer<'de>
    {
        let s = String::deserialize(deserializer)?;
        Ok(match s.as_str() {
            "web" => TraefikDeployment::Web,
            "localhost" => TraefikDeployment::Localhost,
            _ => unimplemented!("Traefik deployment option `{}` is not implemented", s)
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraefikConfig {
    pub deploy: TraefikDeployment,
    pub launch: bool,
    pub subdomain: TraefikSubdomain,
    pub network: TraefikNetwork,   
    pub localhost: TraefikLocalhostConfig,
    pub web: TraefikWebConfig,
    #[serde(skip_deserializing)]
    pub is_localhost: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraefikSubdomain {
    pub app: String,
    pub api: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraefikNetwork {
    pub name: String,
    pub external: bool,
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraefikLocalhostConfig {
    pub tls: bool,       
    pub domain: String,
    #[serde(skip_deserializing)]
    pub entrypoint: String,
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraefikWebConfig {
    pub email: String,
    pub domain: String,
    pub username: String,
    pub password: String,
    #[serde(skip_deserializing)]
    pub password_hashed: String

}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CerebroConfig {
    pub components: SecurityComponentsConfig,
    pub smtp: SmtpConfig,
    #[serde(skip_deserializing)]  // internal config
    pub app: AppConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AppConfig {
    pub public_cerebro_api_url: String,
    pub private_cerebro_api_access_max_age: i64,
    pub private_cerebro_api_access_cookie_secure: bool,
    pub private_cerebro_api_access_cookie_domain: String,
}
impl Default for AppConfig {
    fn default() -> Self {
        Self {
            public_cerebro_api_url: String::from("http://localhost:8080"),
            private_cerebro_api_access_max_age: 720,
            private_cerebro_api_access_cookie_secure: false,
            private_cerebro_api_access_cookie_domain: String::from("localhost"),
        }
    }
}


pub struct StackConfigTree {
    base: PathBuf,
    certs: PathBuf,
    docker: PathBuf,
    traefik: PathBuf,
    cerebro_api: PathBuf,
    mongodb: PathBuf,
    assets: PathBuf,
    keys: PathBuf
}
impl StackConfigTree {
    pub fn new(outdir: &PathBuf) -> Self {
        Self {
            base: outdir.clone(),
            certs: outdir.join("certs"),
            docker: outdir.join("docker"),
            traefik: outdir.join("traefik"),
            cerebro_api: outdir.join("cerebro_api"),
            mongodb: outdir.join("mongodb"),
            assets: outdir.join("assets"),
            keys: outdir.join("keys"),
        }
    }
    pub fn create_dir_all(&self) -> Result<(), StackConfigError> {

        for dir in [
            self.base.clone(),
            self.certs.clone(), 
            self.docker.clone(), 
            self.traefik.clone(), 
            self.cerebro_api.clone(), 
            self.mongodb.clone(), 
            self.assets.clone(),
            self.keys.clone()
        ] { 
            std::fs::create_dir_all(dir).map_err(|err| StackConfigError::StackTreeDirectoryNotCreated(err) )? 
        }

        Ok(())
    }
}


pub fn write_cert(buf: &[u8], path: &PathBuf) -> Result<(), StackConfigError> {
    let mut file = std::fs::File::create(path)
        .map_err(|err| StackConfigError::CertificateFilePathInvalid(err) )?;
    file.write_all(buf).map_err(|err| StackConfigError::CertificateFileNotWritten(err))
}

pub fn write_key(buf: &[u8], path: &PathBuf) -> Result<(), StackConfigError> {
    let mut file = std::fs::File::create(path)
        .map_err(|err| StackConfigError::CertificateKeyFilePathInvalid(err) )?;
    file.write_all(buf).map_err(|err| StackConfigError::CertificateKeyFileNotWritten(err))
}


pub fn write_rendered_template(buf: &[u8], path: &PathBuf) -> Result<(), StackConfigError> {
    let mut file = std::fs::File::create(path)
        .map_err(|err| StackConfigError::TemplateRenderFilePathInvalid(err) )?;
    file.write_all(buf).map_err(|err| StackConfigError::TemplateRenderFileNotWritten(err))
}



#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Stack {
    #[serde(skip_deserializing)]
    name: String,
    #[serde(skip_deserializing)]
    outdir: PathBuf,
    #[serde(skip_deserializing)]
    revision: String,
    #[serde(skip_deserializing)]
    dev: bool,

    cerebro: CerebroConfig,
    mongodb: MongoDbConfig,
    traefik: TraefikConfig,
    
}
impl Stack {
    pub fn from_toml(path: &PathBuf) -> Result<Self, StackConfigError> {

        log::info!("Read stack config from file: {}", path.display());

        let mut config: Self = toml::from_str(
            &std::fs::read_to_string(path).map_err(|_| StackConfigError::ConfigFileInputInvalid)?
        ).map_err(|err| StackConfigError::ConfigFileNotDeserialized(err))?;
        
        config.traefik.is_localhost = match config.traefik.deploy {
            TraefikDeployment::Web => false,
            TraefikDeployment::Localhost => true 
        };

        config.traefik.localhost.entrypoint = match config.traefik.localhost.tls {
            true => String::from("https"),
            false => String::from("http")
        };

        config.revision = CRATE_VERSION.to_string();

        log::info!("Cerebro deployment from version: {}", config.revision);

        Ok(config)

    }
    pub fn to_toml(&self, path: &PathBuf) -> Result<(), StackConfigError> {
        let config = toml::to_string_pretty(&self).map_err(|err| StackConfigError::ConfigFileNotSerialized(err))?;
        std::fs::write(path, config).map_err(|_| StackConfigError::ConfigFileOutputInvalid)
    }
    pub fn create_certs_and_keys(&self, dir_tree: &StackConfigTree) -> Result<TokenEncryptionConfig, StackConfigError> {

        if self.traefik.is_localhost && self.traefik.localhost.tls {

            log::info!("Local HTTPS deployment configured, create certificates");

            // Create a local certificate for the domain - rcgen for localhost domain
            let subject_alt_names = vec![self.traefik.localhost.domain.to_string(), "localhost".to_string()];
            let cert = generate_simple_self_signed(subject_alt_names).unwrap();

            write_cert(
                cert.serialize_pem().map_err(|err| StackConfigError::CertificateFileNotSerialized(err) )?.as_bytes(), 
                &dir_tree.certs.join("cerebro.cert")
            )?;
            write_cert(
                cert.serialize_private_key_pem().as_bytes(),
                &dir_tree.certs.join("cerebro.key")
            )?;
        };

        log::info!("Create access and refresh token keys for user authentication");

        // Create access and refresh token keys and encode to base64
        let mut rng = rand::thread_rng();
        let bits = 3072;

        let access_private_key = RsaPrivateKey::new(&mut rng, bits).map_err(|err| StackConfigError::RsaKeyNotGenerated(err) )?;
        let access_public_key = RsaPublicKey::from(&access_private_key);
        let refresh_private_key = RsaPrivateKey::new(&mut rng, bits).map_err(|err| StackConfigError::RsaKeyNotGenerated(err) )?;
        let refresh_public_key = RsaPublicKey::from(&refresh_private_key);
        
        let access_private_key_pem = access_private_key.to_pkcs8_pem(rsa::pkcs8::LineEnding::LF)
            .map_err(|err| StackConfigError::RsaPrivateKeyPemNotGenerated(err) )?;
        let access_public_key_pem = access_public_key.to_public_key_pem(rsa::pkcs8::LineEnding::LF)
            .map_err(|err| StackConfigError::RsaPublicKeyPemNotGenerated(err) )?;
        let refresh_private_key_pem = refresh_private_key.to_pkcs8_pem(rsa::pkcs8::LineEnding::LF)
            .map_err(|err| StackConfigError::RsaPrivateKeyPemNotGenerated(err) )?;
        let refresh_public_key_pem = refresh_public_key.to_public_key_pem(rsa::pkcs8::LineEnding::LF)
            .map_err(|err| StackConfigError::RsaPublicKeyPemNotGenerated(err) )?;

        write_key(
            &access_private_key_pem.as_bytes(), &dir_tree.keys.join("access_private.pem")
        )?;
        write_key(
            &access_public_key_pem.as_bytes(), &dir_tree.keys.join("access_public.pem")
        )?;
        write_key(
            &refresh_private_key_pem.as_bytes(), &dir_tree.keys.join("refresh_private.pem")
        )?;
        write_key(
            &refresh_public_key_pem.as_bytes(), &dir_tree.keys.join("refresh_public.pem")
        )?;


        Ok(TokenEncryptionConfig {
            access_private_key: general_purpose::STANDARD.encode(access_private_key_pem),
            access_public_key: general_purpose::STANDARD.encode(access_public_key_pem),
            refresh_private_key: general_purpose::STANDARD.encode(refresh_private_key_pem),
            refresh_public_key: general_purpose::STANDARD.encode(refresh_public_key_pem)
        })

    }
    pub fn configure(&mut self, outdir: &PathBuf, dev: bool, subdomain_prefix: Option<String>) -> Result<(), StackConfigError> {

        self.name = match outdir.file_name() {
            Some(os_str) => match os_str.to_str() {
                Some(dir_name) => dir_name.to_string(),
                None => return Err(StackConfigError::ProjectNameInvalid)
            },
            None => return Err(StackConfigError::ProjectNameInvalid)
        };
        log::info!("Initiate server configuration: {}", self.name);

        self.dev = dev;
        log::info!("Development deployment active: {}", dev);


        let dir_tree = StackConfigTree::new(outdir);
        dir_tree.create_dir_all()?;

        self.outdir = std::fs::canonicalize(outdir).map_err(|_| StackConfigError::OutdirAbsolutePathInvalid(outdir.display().to_string()))?;
        log::info!("Create deployment directory tree in: {}", outdir.display());

        
        log::info!("Write stack asset files to deployment directory");
        let stack_assets = StackAssets::new(&dir_tree.assets)?;
        stack_assets.write_asset_files(self.traefik.deploy.clone())?;

        if let Some(prefix) = subdomain_prefix {
            self.traefik.subdomain.api = format!("{prefix}-api");
            self.traefik.subdomain.app = format!("{prefix}-app");
            log::info!("Subdomain prefix specified, configured subdomains to `{}` and `{}`",  self.traefik.subdomain.api.blue(), self.traefik.subdomain.app.blue());
        }

        log::info!("Read default server configuration");
        // Read the default server template configuration
        let mut server_config = api::config::Config::from_toml(&stack_assets.templates.paths.cerebro_server)
            .map_err(|err| StackConfigError::CerebroServerConfigNotParsed(err))?;

        // Configure server and application
        if self.traefik.is_localhost {

            log::info!("Configure server and application deployment to localhost");
           // Localhost specific deployment configurations of Cerebro application and server
            server_config.security.cors.app_origin_public_url = format!("{}://{}.{}", self.traefik.localhost.entrypoint,  self.traefik.subdomain.app, self.traefik.localhost.domain);
            server_config.security.cookies.domain = self.traefik.localhost.domain.clone();
            if self.traefik.localhost.tls {
                server_config.security.cookies.secure = true;
            } else {
                server_config.security.cookies.secure = false;
            }
            self.cerebro.app.public_cerebro_api_url = format!("{}://{}.{}", self.traefik.localhost.entrypoint, self.traefik.subdomain.api, self.traefik.localhost.domain);

        } else {
            log::info!("Configure server and application deployment to web");
            // Web specific deployment configurations of Cerebro application and server
            server_config.security.cors.app_origin_public_url = format!("https://{}.{}", self.traefik.subdomain.app, self.traefik.web.domain);
            server_config.security.cookies.domain = self.traefik.web.domain.clone();
            server_config.security.cookies.secure = true;
            self.cerebro.app.public_cerebro_api_url =  format!("https://{}.{}", self.traefik.subdomain.api, self.traefik.web.domain);
        }

        log::info!("APP deployment configured to: {}", server_config.security.cors.app_origin_public_url.blue());
        log::info!("API deployment configured to: {}", self.cerebro.app.public_cerebro_api_url.blue());
        
        // Cookie settings for refresh access token issued by application server hook
        self.cerebro.app.private_cerebro_api_access_cookie_domain = server_config.security.cookies.domain.clone();
        self.cerebro.app.private_cerebro_api_access_max_age = server_config.security.token.expiration.access_max_age;
        self.cerebro.app.private_cerebro_api_access_cookie_secure = server_config.security.cookies.secure;
        
        
        log::info!("Configure encryption keys and certificates");
        server_config.security.token.encryption = self.create_certs_and_keys(&dir_tree)?;
        server_config.security.components = self.cerebro.components.clone();  
        

        log::info!("Hash Cerebro admin password using salted Argon2");
        // Hash the required passwords - Argon2 for Cerebro user database
        self.mongodb.cerebro_admin_password_hashed = hash_password(&self.mongodb.cerebro_admin_password).map_err(|_| StackConfigError::CerebroDatabasePasswordNotHashed)?;

        log::info!("Hash Traefik dashboard password using Bcrypt2");
        if !self.traefik.is_localhost {
            // Hash the required passwords - Bcrypt for Traefik dashboard 
            self.traefik.web.password_hashed = hash(&self.traefik.web.password, DEFAULT_COST).map_err(|err| StackConfigError::TraefikDashboardPasswordNotHashed(err))?;
        };

        log::info!("Use provided SMTP email configuration");      
        server_config.smtp = self.cerebro.smtp.clone();
 
        log::info!("Configure internal Docker network MongoDB URI"); 
        server_config.database.connections.mongodb_uri = format!(
            "mongodb://{}:{}@cerebro-database:27017/?authSource=admin", 
            self.mongodb.admin_username, 
            self.mongodb.admin_password
        );

        log::info!("Write modified server config to: {}", &dir_tree.cerebro_api.join("server.toml").display()); 
        // Write the modified server config
        server_config.to_toml(&dir_tree.cerebro_api.join("server.toml")).map_err(|_: ConfigError| StackConfigError::ConfigFileOutputInvalid)?;

        // // Write the template directories for email and report
        // self.write_templates(&stack_assets, &dir_tree)?;
        
        log::info!("Render asset templates with deployment configurations"); 
        // Render and write the asset templates
        self.render_templates(&stack_assets, &dir_tree)?;

        log::info!("Write modified stack configuration to: {}", &dir_tree.base.join("stack.toml").display()); 
        // Write the final stack configuration
        self.to_toml(&dir_tree.base.join("stack.toml"))

    }
    #[cfg(feature = "libgit")]
    pub fn clone_and_checkout_repository_libgit(&self, url: &str, branch: Option<String>, revision: Option<String>, ssh_private_key: Option<PathBuf>, passphrase: Option<String>) -> Result<(), StackConfigError> {

        let repo_subdir = self.outdir.join("cerebro");
        std::fs::create_dir_all(&repo_subdir).map_err(|err| StackConfigError::StackTreeDirectoryNotCreated(err) )?;
        log::info!("Cloning repository into deployment subdirectory: {}", &repo_subdir.display());

        // Prepare builder.
        let mut builder = git2::build::RepoBuilder::new();
        if let Some(br) = branch { log::info!("Cloning repository at branch: {}", &br); builder.branch(&br); }

        if let Some(key) = &ssh_private_key {
            // Prepare callbacks.
            let mut callbacks = RemoteCallbacks::new();
            callbacks.credentials(|_url, username_from_url, _allowed_types| {
                Cred::ssh_key(
                    username_from_url.unwrap(),
                    None,
                    key,
                    passphrase.as_deref(),
                )
            });

            // Prepare fetch options.
            let mut fo = git2::FetchOptions::new();
            fo.remote_callbacks(callbacks);
            
            builder.fetch_options(fo);

        };

        // Clone the project.
        let repo = match builder.clone(
            &url,
            &repo_subdir,
        ) {
            Ok(repo) => repo,
            Err(e) => return Err(StackConfigError::RemoteRepositoryNotCloned(e))
        };

        if let Some(rev) = revision {
            
            log::info!("Checking out revision for deployment: {}", &rev);

            // Git2 > 0.13.18 - checkout a tag (0.1.1) or a commit (8e8128)
            let (object, reference) = repo.revparse_ext(&rev).map_err(|e| StackConfigError::RemoteRepositoryNotCheckedOut(e))?;
            
            repo.checkout_tree(&object, None).map_err(|e| StackConfigError::RemoteRepositoryNotCheckedOut(e))?;

            match reference {
                // gref is an actual reference like branches or tags
                Some(gref) => repo.set_head(gref.name().unwrap()),  // TODO: unwrap call will fail on non-UTF8
                // this is a commit, not a reference
                None => repo.set_head_detached(object.id()),
            }.map_err(|e| StackConfigError::RemoteRepositoryNotCheckedOut(e))?
        }

    }
    /// Clone the repository into the deployment folder and checkout the requested revision (can be done manually after deployment)
    pub fn clone_and_checkout_repository_process(&self, url: &str, branch: Option<String>, revision: Option<String>) -> Result<(), StackConfigError> {

        let repo_dir = self.outdir.join("cerebro");
        let repo_dir_str = self.outdir.join("cerebro").display().to_string();;
        
        if repo_dir.exists() {
            log::warn!("Repository exists at: {}", repo_dir.display());
            log::warn!("Repository will be deleted and replaced with fresh clone");
            std::fs::remove_dir_all(&repo_dir).expect("Failed to remove existing deployment repository");
        }

        let git_command = "git";
        let outdir_str = self.outdir.display().to_string();

        let git_arguments = vec!["-C", &outdir_str, "clone", &url];

        log::info!("Clone repository using git system call:");
        log::info!("git {}", git_arguments.join(" "));

        let mut cmd = std::process::Command::new(git_command);

        cmd.args(&git_arguments);

        let output = cmd.output().expect("Failed to execute command");

        if output.status.success() {
            let _ = String::from_utf8_lossy(&output.stdout);
            log::info!("Cloned repository for deployment")
        } else {
            let stderr = String::from_utf8_lossy(&output.stderr);
            log::error!("Failed to clone repository for deployment: {}", stderr);
        }

        let checkout = match (branch, revision) {
            (Some(b), None) => Some(b),
            (None, Some(r)) => Some(r),
            (Some(b), Some(r)) => Some(r),  // prefer revision over branch if both supplied
            (None, None) => None
        };

        if let Some(p) = &checkout {
            let checkout_args = vec!["-C", &repo_dir_str, "checkout", p];

            log::info!("Checkout repository using git system call:");
            log::info!("git {}", checkout_args.join(" "));

            let mut cmd = std::process::Command::new(git_command);

            cmd.args(&checkout_args);

            let output = cmd.output().expect("Failed to execute command");

            if output.status.success() {
                let _ = String::from_utf8_lossy(&output.stdout);
                log::info!("Checked out repository for deployment")
            } else {
                let stderr = String::from_utf8_lossy(&output.stderr);
                log::error!("Failed to checkout repository for deployment: {}", stderr);
            }
        }
            
        Ok(())
    }
    pub fn render_templates(&self, stack_assets: &StackAssets, dir_tree: &StackConfigTree) -> Result<(), StackConfigError> {

        let mut handlebars = handlebars::Handlebars::new();
    
        handlebars.register_template_file(&stack_assets.templates.names.traefik_static, &stack_assets.templates.paths.traefik_static).map_err(|err| StackConfigError::TemplateNotRegistered(err))?;
        handlebars.register_template_file(&stack_assets.templates.names.traefik_dynamic, &stack_assets.templates.paths.traefik_dynamic).map_err(|err| StackConfigError::TemplateNotRegistered(err))?;

        let render = handlebars.render(&stack_assets.templates.names.traefik_static, &self).map_err(|err| StackConfigError::TemplateNotRendered(err))?;
        write_rendered_template(render.as_bytes(), &dir_tree.traefik.join("static.yml"))?;

        let render = handlebars.render(&stack_assets.templates.names.traefik_dynamic, &self).map_err(|err| StackConfigError::TemplateNotRendered(err))?;
        write_rendered_template(render.as_bytes(), &dir_tree.traefik.join("dynamic.yml"))?;

        handlebars.register_template_file(&stack_assets.templates.names.docker_compose, &stack_assets.templates.paths.docker_compose).map_err(|err| StackConfigError::TemplateNotRegistered(err))?;
        let render = handlebars.render(&stack_assets.templates.names.docker_compose, &self).map_err(|err| StackConfigError::TemplateNotRendered(err))?;
        write_rendered_template(render.as_bytes(), &dir_tree.base.join("docker-compose.yml"))?;

        if !self.traefik.launch {
            handlebars.register_template_file(&stack_assets.templates.names.docker_compose_traefik, &stack_assets.templates.paths.docker_compose_traefik).map_err(|err| StackConfigError::TemplateNotRegistered(err))?;
            let render = handlebars.render(&stack_assets.templates.names.docker_compose_traefik, &self).map_err(|err| StackConfigError::TemplateNotRendered(err))?;
            write_rendered_template(render.as_bytes(), &dir_tree.base.join("docker-compose.traefik.yml"))?;
        }

        handlebars.register_template_file(&stack_assets.templates.names.docker_server, &stack_assets.templates.paths.docker_server).map_err(|err| StackConfigError::TemplateNotRegistered(err))?;
        let render = handlebars.render(&stack_assets.templates.names.docker_server, &self).map_err(|err| StackConfigError::TemplateNotRendered(err))?;
        write_rendered_template(render.as_bytes(), &dir_tree.docker.join("Dockerfile.server"))?;

        handlebars.register_template_file(&stack_assets.templates.names.docker_app, &stack_assets.templates.paths.docker_app).map_err(|err| StackConfigError::TemplateNotRegistered(err))?;
        let render = handlebars.render(&stack_assets.templates.names.docker_app, &self).map_err(|err| StackConfigError::TemplateNotRendered(err))?;
        write_rendered_template(render.as_bytes(), &dir_tree.docker.join("Dockerfile.app"))?;

        handlebars.register_template_file(&stack_assets.templates.names.mongodb_init_sh, &stack_assets.templates.paths.mongodb_init_sh).map_err(|err| StackConfigError::TemplateNotRegistered(err))?;
        let render = handlebars.render(&stack_assets.templates.names.mongodb_init_sh, &self).map_err(|err| StackConfigError::TemplateNotRendered(err))?;
        write_rendered_template(render.as_bytes(), &dir_tree.mongodb.join("mongo-init.sh"))?;

        handlebars.register_template_file(&stack_assets.templates.names.mongodb_init_env, &stack_assets.templates.paths.mongodb_init_env).map_err(|err| StackConfigError::TemplateNotRegistered(err))?;
        let render = handlebars.render(&stack_assets.templates.names.mongodb_init_env, &self).map_err(|err| StackConfigError::TemplateNotRendered(err))?;
        write_rendered_template(render.as_bytes(), &dir_tree.mongodb.join("database.env"))?;

        handlebars.register_template_file(&stack_assets.templates.names.cerebro_app, &stack_assets.templates.paths.cerebro_app).map_err(|err| StackConfigError::TemplateNotRegistered(err))?;
        let render = handlebars.render(&stack_assets.templates.names.cerebro_app, &self).map_err(|err| StackConfigError::TemplateNotRendered(err))?;
        write_rendered_template(render.as_bytes(), &dir_tree.cerebro_api.join("app.env"))?;
        
        Ok(())

    }
    pub fn interactive() {

    }
}