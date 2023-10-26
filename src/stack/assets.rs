use rust_embed::RustEmbed;

#[derive(RustEmbed)]
#[folder = "templates/stack/mongodb/"]
#[prefix = "mongodb/"]
pub struct MongoInitTemplates;

#[derive(RustEmbed)]
#[folder = "templates/stack/docker/"]
#[prefix = "docker/"]
pub struct DockerTemplates;

#[derive(RustEmbed)]
#[folder = "templates/stack/traefik/"]
#[prefix = "traefik/"]
pub struct TraefikTemplates;

#[derive(RustEmbed)]
#[folder = "templates/stack/cerebro/"]
#[prefix = "cerebro/"]
pub struct CerebroTemplates;

#[derive(RustEmbed)]
#[folder = "templates/email/"]
#[prefix = "email/"]
pub struct EmailTemplates;

#[derive(RustEmbed)]
#[folder = "templates/report/"]
#[prefix = "report/"]
pub struct ReportTemplates;

