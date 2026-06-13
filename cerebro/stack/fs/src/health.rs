//! Topology health checks for Cerebro FS.
//!
//! A lightweight, dependency-free health summary of the SeaweedFS components the
//! client talks to: the **master** (cluster health) and, when a filer is in use,
//! the **filer** (path-addressed HTTP API). It is intended for readiness checks
//! and operator visibility; richer telemetry (Prometheus/Grafana) is a later
//! stage.
//!
//! Per-replica/volume status for an individual object is available via
//! [`FileSystemClient::lookup_locations`](crate::client::FileSystemClient) and is
//! surfaced per file by the integrity sweep (see [`crate::integrity`]).

use reqwest::StatusCode;

use crate::client::FileSystemClient;
use crate::config::FsAccessMode;
use crate::filer::FilerClient;

/// Health of a single storage component.
#[derive(Debug, Clone)]
pub struct ComponentHealth {
    /// Component name, e.g. `master` or `filer`.
    pub component: String,
    /// Whether the component responded successfully.
    pub reachable: bool,
    /// Human-readable detail (endpoint, status, or error).
    pub detail: String,
}

/// Aggregate health across the reachable storage components.
#[derive(Debug, Clone)]
pub struct TopologyHealth {
    pub components: Vec<ComponentHealth>,
}
impl TopologyHealth {
    /// True only when every checked component is reachable.
    pub fn healthy(&self) -> bool {
        !self.components.is_empty() && self.components.iter().all(|c| c.reachable)
    }
}

impl FileSystemClient {
    /// Check the health of the SeaweedFS topology the client can reach.
    ///
    /// Always checks the master (`/cluster/healthz`). When the access mode is
    /// [`FsAccessMode::Filer`] the filer is checked as well. Unlike
    /// [`FileSystemClient::ping_status`](crate::client::FileSystemClient::ping_status)
    /// this does not short-circuit on the first failure: it returns a status for
    /// every component so an operator can see the whole picture.
    ///
    /// The S3 archival remote (Model B) is intentionally not probed here — Glacier
    /// is not a synchronously-pingable endpoint and is the filer/volume tier-move
    /// layer's concern, not the client's.
    pub fn topology_health(&self) -> TopologyHealth {
        let mut components = Vec::new();

        // Master
        let master_url = format!("{}/cluster/healthz", self.get_url());
        let master = match reqwest::blocking::Client::new().get(&master_url).send() {
            Ok(resp) if resp.status() == StatusCode::OK => ComponentHealth {
                component: "master".to_string(),
                reachable: true,
                detail: format!("{} ok", self.get_url()),
            },
            Ok(resp) => ComponentHealth {
                component: "master".to_string(),
                reachable: false,
                detail: format!("{} returned {}", self.get_url(), resp.status()),
            },
            Err(err) => ComponentHealth {
                component: "master".to_string(),
                reachable: false,
                detail: format!("{} unreachable: {}", self.get_url(), err),
            },
        };
        components.push(master);

        // Filer (only relevant when path-addressed access is configured)
        if self.config.access == FsAccessMode::Filer {
            let base = self.config.filer_base();
            let filer = match FilerClient::new(&base, self.config.danger_invalid_certificate) {
                Ok(client) => match client.health() {
                    Ok(()) => ComponentHealth {
                        component: "filer".to_string(),
                        reachable: true,
                        detail: format!("{} ok", base),
                    },
                    Err(err) => ComponentHealth {
                        component: "filer".to_string(),
                        reachable: false,
                        detail: format!("{} unhealthy: {}", base, err),
                    },
                },
                Err(err) => ComponentHealth {
                    component: "filer".to_string(),
                    reachable: false,
                    detail: format!("{} client error: {}", base, err),
                },
            };
            components.push(filer);
        }

        TopologyHealth { components }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn component(name: &str, reachable: bool) -> ComponentHealth {
        ComponentHealth { component: name.to_string(), reachable, detail: String::new() }
    }

    #[test]
    fn healthy_requires_all_components_reachable() {
        let all_ok = TopologyHealth { components: vec![component("master", true), component("filer", true)] };
        assert!(all_ok.healthy());

        let degraded = TopologyHealth { components: vec![component("master", true), component("filer", false)] };
        assert!(!degraded.healthy());
    }

    #[test]
    fn empty_topology_is_not_healthy() {
        let empty = TopologyHealth { components: Vec::new() };
        assert!(!empty.healthy());
    }
}
