use std::collections::{HashMap};
use crate::taxa::taxon::LineageOperations;

use crate::{modules::pathogen::{AbundanceMode, ProfileTool}, taxa::{filter::TaxonFilterConfig, taxon::Taxon}};

type SampleId = String;
type Tag = String;

#[derive(Clone, Copy, Debug)]
pub enum ScoreMode {
    Hard,
    /// tau > 0. Smaller tau -> closer to hard min/max and step counts.
    Soft { tau: f64 },
}

pub trait RuleMarginScore {
    fn s_score_with(
        &self,
        cfg: &TaxonFilterConfig,
        sample_tags: &HashMap<SampleId, Vec<Tag>>,
        mode: ScoreMode,
    ) -> f64;

    #[inline]
    fn s_score(
        &self,
        cfg: &TaxonFilterConfig,
        sample_tags: &HashMap<SampleId, Vec<Tag>>,
    ) -> f64 {
        self.s_score_with(cfg, sample_tags, ScoreMode::Hard)
    }
}

impl RuleMarginScore for Taxon {
    fn s_score_with(
        &self,
        cfg: &TaxonFilterConfig,
        sample_tags: &HashMap<SampleId, Vec<Tag>>,
        mode: ScoreMode,
    ) -> f64 {
        const EPS: f64 = 1e-9;

        // ---------- helpers ----------
        #[inline] fn margin_ge(v: f64, t: f64) -> f64 { (v - t) / t.max(1e-9) }
        #[inline] fn mixed_pair_ok(tool: &ProfileTool, mode: &AbundanceMode) -> bool {
            matches!(
                (tool, mode),
                (ProfileTool::Bracken,  AbundanceMode::Profile) |
                (ProfileTool::Kraken2,  AbundanceMode::Sequence) |
                (ProfileTool::Ganon2,   AbundanceMode::Sequence) |
                (ProfileTool::Sylph,    AbundanceMode::Sequence) |
                (ProfileTool::Kmcp,     AbundanceMode::Sequence) |
                (ProfileTool::Metabuli, AbundanceMode::Sequence) |
                (ProfileTool::Blast,    AbundanceMode::Bases)    |
                (ProfileTool::Vircov,   AbundanceMode::Sequence)
            )
        }
        #[inline] fn is_kmer_tool(t: &ProfileTool) -> bool {
            !matches!(t, ProfileTool::Vircov | ProfileTool::Blast | ProfileTool::Sylph)
        }

        // Soft primitives
        let (use_soft, tau) = match mode {
            ScoreMode::Hard => (false, 0.0),
            ScoreMode::Soft { tau } => (true, tau.max(1e-6)),
        };

        let softmin = |xs: &[f64]| -> f64 {
            if !use_soft { return xs.iter().fold(f64::INFINITY, |a, b| a.min(*b)); }
            // -tau * log sum_i exp(-x_i / tau)
            let m = xs.iter().fold(f64::INFINITY, |a, &x| a.min(x));
            let z = xs.iter().map(|&x| (-(x - m) / tau).exp()).sum::<f64>();
            -(tau) * (z.ln()) + m
        };
        let softmax = |xs: &[f64]| -> f64 {
            if !use_soft { return xs.iter().fold(f64::NEG_INFINITY, |a, b| a.max(*b)); }
            let m = xs.iter().fold(f64::NEG_INFINITY, |a, &x| a.max(x));
            let z = xs.iter().map(|&x| ((x - m) / tau).exp()).sum::<f64>();
            tau * z.ln() + m
        };
        let sigmoid = |z: f64| -> f64 {
            if !use_soft { return if z >= 0.0 { 1.0 } else { 0.0 }; }
            // logistic with temperature tau
            (z / tau).tanh() * 0.5 + 0.5  // numerically stable-ish alternative to 1/(1+exp(-z/tau))
        };

        // ---------- static gates ----------
        if let Some(r) = &cfg.rank {
            let want: taxonomy::TaxRank = (*r).clone().into();
            if self.rank != want { return f64::NEG_INFINITY; }
        }
        if !cfg.domains.is_empty() {
            match self.lineage.get_domain() {
                Some(dom) if cfg.domains.contains(&dom) => {},
                _ => return f64::NEG_INFINITY,
            }
        }
        if let Some(targets) = cfg.target_set() {
            if !targets.contains(self.name.as_str()) { return f64::NEG_INFINITY; }
        }
        if let Some(bad) = &cfg.ignore_taxstr {
            if bad.iter().any(|s| self.name.contains(s)) { return f64::NEG_INFINITY; }
        }

        // ---------- NTC map ----------
        let only_ntc = sample_tags.values()
            .all(|tags| tags.contains(&"NTC".into()) || tags.contains(&"ENV".into()));
        let mut ntc_rpm: HashMap<(ProfileTool, AbundanceMode, String), f64> = HashMap::new();
        for rec in &self.evidence.profile {
            if let Some(tags) = sample_tags.get(&rec.id) {
                if tags.contains(&"NTC".into()) || tags.contains(&"ENV".into()) {
                    if let Some(na) = tags.iter().find(|t| *t == "DNA" || *t == "RNA") {
                        *ntc_rpm.entry((rec.tool.clone(), rec.mode.clone(), na.clone()))
                            .or_insert(0.0) += rec.rpm;
                    }
                }
            }
        }

        // ---------- filter profile ----------
        let mut profile = Vec::with_capacity(self.evidence.profile.len());
        'outer: for rec in self.evidence.profile.clone() {
            if let Some(tags) = sample_tags.get(&rec.id) {
                if tags.contains(&"NTC".into()) || tags.contains(&"ENV".into()) {
                    if only_ntc { profile.push(rec); }
                    continue 'outer;
                }
                if let (Some(ratio), Some(na)) = (cfg.ntc_ratio, tags.iter().find(|t| *t == "DNA" || *t == "RNA")) {
                    if let Some(ntc) = ntc_rpm.get(&(rec.tool.clone(), rec.mode.clone(), na.clone())) {
                        if (*ntc / (rec.rpm + EPS)) > ratio { continue 'outer; }
                    }
                }
            }
            let tool_ok = cfg.tools.is_empty() || cfg.tools.contains(&rec.tool);
            let mode_ok = cfg.modes.is_empty()
                || cfg.modes.contains(&rec.mode)
                || (cfg.modes.contains(&AbundanceMode::Mixed) && mixed_pair_ok(&rec.tool, &rec.mode));
            if !(tool_ok && mode_ok) { continue 'outer; }

            let pass_top = if rec.mode != AbundanceMode::Bases {
                rec.reads >= cfg.min_reads &&
                rec.rpm   >= cfg.min_rpm &&
                cfg.max_rpm.map_or(true, |mx| rec.rpm <= mx)
            } else {
                rec.bases >= cfg.min_bases &&
                rec.bpm   >= cfg.min_bpm &&
                cfg.max_bases.map_or(true, |mx| rec.bases <= mx)
            } && rec.abundance >= cfg.min_abundance;
            if !pass_top { continue 'outer; }

            profile.push(rec);
        }

        // ---------- precompute aggregates ----------
        let max_vircov_rpm = profile.iter()
            .filter(|r| matches!(r.tool, ProfileTool::Vircov))
            .map(|r| r.rpm).fold(0.0_f64, f64::max);

        let blast_count = profile.iter()
            .filter(|r| matches!(r.tool, ProfileTool::Blast))
            .count() as f64;

        let aln_tools_count = if profile.iter().any(|r| matches!(r.tool, ProfileTool::Vircov)) { 1.0 } else { 0.0 };

        let regions_margin = |min_regions: u64, cov_thresh: f64| -> f64 {
            // OR over alignment records → max margin
            if self.evidence.alignment.is_empty() {
                f64::NEG_INFINITY
            } else {
                let mut vals = Vec::with_capacity(self.evidence.alignment.len());
                for rec in &self.evidence.alignment {
                    let m = if rec.scan_coverage > cov_thresh {
                        f64::INFINITY
                    } else {
                        margin_ge(rec.scan_regions as f64, min_regions as f64)
                    };
                    vals.push(m);
                }
                softmax(&vals)
            }
        };

        let ntc_margin_min = if let Some(ratio) = cfg.ntc_ratio {
            let mut ms = Vec::new();
            for rec in &profile {
                if let Some(tags) = sample_tags.get(&rec.id) {
                    if let Some(na) = tags.iter().find(|t| *t == "DNA" || *t == "RNA") {
                        if let Some(ntc) = ntc_rpm.get(&(rec.tool.clone(), rec.mode.clone(), na.clone())) {
                            ms.push(1.0 - (ntc / (rec.rpm + EPS)) / ratio.max(EPS));
                        }
                    }
                }
            }
            if ms.is_empty() { None } else { Some(softmin(&ms)) }
        } else { None };

        // ---------- build clause scores ----------
        
        let lineage_filters = cfg.lineage.as_ref().map(|v| v.as_slice()).unwrap_or(&[]);
        let mut clause_scores = Vec::new();

        for lf in lineage_filters {
            if !lf.lineages.iter().any(|l| self.lineage.contains(l)) { continue; }

            // tag scoping
            let required_tags: Vec<String> = if !lf.tags.is_empty() { lf.tags.clone() } else { Vec::new() };
            let prof_scoped: Vec<_> = if !required_tags.is_empty() {
                profile.iter().filter(|rec| {
                    if let Some(tags) = sample_tags.get(&rec.id) {
                        tags.iter().any(|t| required_tags.contains(t))
                    } else { false }
                }).collect()
            } else { profile.iter().collect() };
            if prof_scoped.is_empty() { continue; }

            let mut atoms: Vec<f64> = Vec::new();

            // k-mer tools ≥ min_kmer_tools
            if let Some(min_tools) = lf.min_kmer_tools {
                let thr = lf.min_kmer_rpm.unwrap_or(0.0);
                if use_soft {
                    // soft count via sigmoids of per-tool margins
                    let mut ps = Vec::new();
                    for r in &prof_scoped {
                        if is_kmer_tool(&r.tool) {
                            ps.push(sigmoid(margin_ge(r.rpm, thr)));
                        }
                    }
                    let soft_c = ps.iter().sum::<f64>();
                    atoms.push(soft_c - min_tools as f64);
                } else {
                    let c = prof_scoped.iter()
                        .filter(|r| is_kmer_tool(&r.tool) && r.rpm >= thr).count() as f64;
                    atoms.push(c - min_tools as f64);
                }
            }

            // alignment tools ≥ min_alignment_tools
            if let Some(min_tools) = lf.min_alignment_tools {
                atoms.push(aln_tools_count - min_tools as f64);
            }

            // alignment rpm ≥ min_alignment_rpm (Vircov proxy)
            if let Some(min_rpm) = lf.min_alignment_rpm {
                atoms.push(margin_ge(max_vircov_rpm, min_rpm));
            }

            // alignment regions with coverage rule
            match (lf.min_alignment_regions, lf.min_alignment_regions_coverage) {
                (Some(r), Some(cov)) => atoms.push(regions_margin(r, cov)),
                (Some(r), None) => {
                    if self.evidence.alignment.is_empty() {
                        atoms.push(f64::NEG_INFINITY);
                    } else {
                        let vals: Vec<f64> = self.evidence.alignment.iter()
                            .map(|a| margin_ge(a.scan_regions as f64, r as f64)).collect();
                        atoms.push(softmax(&vals)); // OR over records
                    }
                }
                _ => {}
            }

            // assembly tools (Blast) ≥ min_assembly_tools
            if let Some(min_tools) = lf.min_assembly_tools {
                if use_soft {
                    // treat each Blast record as a 1 with tiny slope around 0 → same as hard here
                    atoms.push(blast_count - min_tools as f64);
                } else {
                    atoms.push(blast_count - min_tools as f64);
                }
            }

            if let Some(m) = ntc_margin_min { atoms.push(m); }

            if atoms.is_empty() {
                clause_scores.push(f64::INFINITY);
            } else {
                clause_scores.push(softmin(&atoms)); // AND over atoms
            }
        }

        if clause_scores.is_empty() { return f64::NEG_INFINITY; }
        softmax(&clause_scores) // OR over clauses
    }
}