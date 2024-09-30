use rust_embed::RustEmbed;

#[derive(RustEmbed)]
#[folder = "../../templates/stack/mongodb/"]
#[prefix = "mongodb/"]
pub struct MongoInitTemplates;

#[derive(RustEmbed)]
#[folder = "../../templates/stack/docker/"]
#[prefix = "docker/"]
pub struct DockerTemplates;

#[derive(RustEmbed)]
#[folder = "../../templates/stack/traefik/"]
#[prefix = "traefik/"]
pub struct TraefikTemplates;

#[derive(RustEmbed)]
#[folder = "../../templates/stack/cerebro/"]
#[prefix = "cerebro/"]
pub struct CerebroTemplates;

