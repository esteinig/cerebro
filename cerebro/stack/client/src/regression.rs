use std::collections::HashMap;
use std::collections::HashSet;

use cerebro_model::api::cerebro::response::TaxonHistoryResult;
use linregress::{FormulaRegressionBuilder, RegressionDataBuilder};
use statrs::distribution::ContinuousCDF;
use statrs::distribution::StudentsT;

use plotters::prelude::*;

use crate::error::HttpClientError;

/// Configuration for the analysis.
#[derive(Debug, Clone)]
pub struct RpmConfig {
    /// Confidence level (e.g. 0.95 for 95%).
    pub confidence_level: f64,
    /// Whether to use a prediction interval (adds extra variance) versus a confidence interval.
    pub use_prediction_interval: bool,
    /// If true, only flag as outlier those samples with higher-than-expected taxon RPM.
    pub only_high_outliers: bool,
}

impl Default for RpmConfig {
    fn default() -> Self {
        Self {
            confidence_level: 0.95,
            use_prediction_interval: true,
            only_high_outliers: false,
        }
    }
}

/// Builder for RpmConfig.
pub struct RpmConfigBuilder {
    config: RpmConfig,
}

impl RpmConfigBuilder {
    pub fn new() -> Self {
        Self {
            config: RpmConfig::default(),
        }
    }

    /// Set the confidence level.
    pub fn confidence_level(mut self, level: f64) -> Self {
        self.config.confidence_level = level;
        self
    }

    /// Specify whether to use a prediction interval.
    pub fn use_prediction_interval(mut self, use_pi: bool) -> Self {
        self.config.use_prediction_interval = use_pi;
        self
    }

    /// Specify whether to only flag outliers with higher-than-expected taxon RPM.
    pub fn only_high_outliers(mut self, only_high: bool) -> Self {
        self.config.only_high_outliers = only_high;
        self
    }

    /// Build the config.
    pub fn build(self) -> RpmConfig {
        self.config
    }
}

/// Struct to store detailed information about an outlier.
#[derive(Debug)]
pub struct Outlier {
    pub sample_id: String,
    pub raw_host: f64,
    pub raw_taxon: f64,
    pub log_host: f64,
    pub log_taxon: f64,
    pub predicted: f64,
    pub lower: f64,
    pub upper: f64,
}

/// The overall analysis result.
#[derive(Debug)]
pub struct RpmAnalysisResult {
    pub intercept: f64,
    pub slope: f64,
    pub relationship: String,
    pub residual_standard_error: f64,
    pub outliers: Vec<Outlier>,
}

impl RpmAnalysisResult {
    /// True if `sample_id` is among the detected high outliers — the shared `contam_history`
    /// rescue rule (Task B, TB-D4) used identically by the stack path (`CerebroClient::get_taxa`)
    /// and the self-contained local path (`taxon_rescued`).
    pub fn has_sample_outlier(&self, sample_id: &str) -> bool {
        self.outliers.iter().any(|o| o.sample_id == sample_id)
    }
}

/// Decide whether a prevalence-flagged contaminant taxon should be **rescued** for `sample_id`,
/// using a (run-local or hybrid) cross-sample `history` (Task B, TB-D1/TB-D4).
///
/// Mirrors the stack's per-taxon rescue exactly: filter each record's taxa to {`taxon_name`,
/// "Homo sapiens"} (the server does this via its per-taxon `s__` history query), regress
/// taxon-RPM vs host-RPM at 95% confidence detecting only high outliers, and rescue iff this
/// sample is a high outlier. `RpmAnalyzer::from_taxon_history` already drops non-positive-RPM
/// points before the log₁₀ transform, so the analyzer cannot panic on this path (TB.4); an
/// analyzer error is treated conservatively as "do not rescue".
pub fn taxon_rescued(
    history: &[TaxonHistoryResult],
    taxon_name: &str,
    sample_id: &str,
) -> bool {
    let filtered: Vec<TaxonHistoryResult> = history
        .iter()
        .map(|record| {
            let mut record = record.clone();
            record
                .taxa
                .retain(|t| t.name == taxon_name || t.name == "Homo sapiens");
            record
        })
        .collect();

    let config = RpmConfigBuilder::new()
        .confidence_level(0.95)
        .only_high_outliers(true)
        .build();

    match RpmAnalyzer::from_taxon_history(config, filtered).run() {
        Ok(result) => result.has_sample_outlier(sample_id),
        Err(_) => false,
    }
}

/// Analyzer struct that holds configuration and the data.
pub struct RpmAnalyzer {
    config: RpmConfig,
    raw_host_rpm: Vec<f64>,
    raw_taxon_rpm: Vec<f64>,
    log_host_rpm: Vec<f64>,
    log_taxon_rpm: Vec<f64>,
    sample_ids: Vec<String>,
}

impl RpmAnalyzer {
    /// Create a new analyzer with the given configuration and raw data;
    /// raw data will be log₁₀ transformed.
    /// (No sample IDs are provided in this case.)
    pub fn new(config: RpmConfig, raw_host_rpm: Vec<f64>, raw_taxon_rpm: Vec<f64>) -> Self {
        Self::new_with_sample_ids(config, raw_host_rpm, raw_taxon_rpm, Vec::new())
    }

    /// Create a new analyzer with sample IDs.
    pub fn new_with_sample_ids(
        config: RpmConfig,
        raw_host_rpm: Vec<f64>,
        raw_taxon_rpm: Vec<f64>,
        sample_ids: Vec<String>,
    ) -> Self {
        // log₁₀ is only defined for positive values. Rather than panic (these constructors are
        // public and sit on a clinical regression path — SKILLS.md: no panics on prod paths), drop
        // any point whose host OR taxon RPM is non-positive. Points are paired (host[i], taxon[i],
        // sample_ids[i]), so we filter by triple to keep the vectors index-aligned. `from_taxon_history`
        // already pre-filters to positive points, so its behaviour is unchanged; only direct callers
        // (e.g. `new`, or future callers) that would previously have crashed now skip cleanly.
        let track_ids = sample_ids.len() == raw_host_rpm.len();

        let mut kept_host_rpm = Vec::with_capacity(raw_host_rpm.len());
        let mut kept_taxon_rpm = Vec::with_capacity(raw_taxon_rpm.len());
        let mut kept_sample_ids = Vec::new();
        let mut log_host_rpm = Vec::with_capacity(raw_host_rpm.len());
        let mut log_taxon_rpm = Vec::with_capacity(raw_taxon_rpm.len());

        for (i, (&host, &taxon)) in raw_host_rpm.iter().zip(raw_taxon_rpm.iter()).enumerate() {
            if host <= 0.0 || taxon <= 0.0 {
                log::warn!(
                    "regression: skipping non-positive RPM point (host={host}, taxon={taxon}) outside the log₁₀ domain"
                );
                continue;
            }
            kept_host_rpm.push(host);
            kept_taxon_rpm.push(taxon);
            log_host_rpm.push(host.log10());
            log_taxon_rpm.push(taxon.log10());
            if track_ids {
                kept_sample_ids.push(sample_ids[i].clone());
            }
        }

        Self {
            config,
            raw_host_rpm: kept_host_rpm,
            raw_taxon_rpm: kept_taxon_rpm,
            log_host_rpm,
            log_taxon_rpm,
            sample_ids: kept_sample_ids,
        }
    }

    /// Constructs an RPMAnalyzer from a vector of TaxonHistoryResult.
    ///
    /// This method filters the history based on sample tags and computes the RPM values as:
    ///
    /// - **host_rpm:** Uses `host_reads` (or 0 if absent) and adds the average of all profile
    ///   record reads from taxa with name "Homo sapiens".
    ///
    /// - **taxon_rpm:** Uses the average of all profile record reads from taxa where the name is not "Homo sapiens".
    ///
    /// Both are normalized using `input_reads` and scaled by 1e6.
    pub fn from_taxon_history(config: RpmConfig, history: Vec<TaxonHistoryResult>) -> Self {
        let mut raw_host_rpm = Vec::new();
        let mut raw_taxon_rpm = Vec::new();
        let mut sample_ids = Vec::new();

        // Define the disqualifying tags.
        let disqualifying_tags = ["NTC", "ENV", "POS", "NEG"];

        for record in history.into_iter() {
            // Filter: must have "DNA" and not contain any disqualifying tags.
            if !record.sample_tags.contains(&"DNA".to_string()) {
                continue;
            }
            if record
                .sample_tags
                .iter()
                .any(|tag| disqualifying_tags.contains(&tag.as_str()))
            {
                continue;
            }

            let input_reads = record.input_reads as f64;
            // Skip if input_reads is zero.
            if input_reads == 0.0 {
                continue;
            }

            // Compute average profile reads for "Homo sapiens" taxa.
            let human_profile_reads: Vec<u64> = record
                .taxa
                .iter()
                .filter(|taxon| taxon.name == "Homo sapiens")
                .flat_map(|taxon| taxon.evidence.profile.iter().map(|profile| profile.reads))
                .collect();
            let avg_human_profile = if human_profile_reads.is_empty() {
                0.0
            } else {
                human_profile_reads.iter().sum::<u64>() as f64 / human_profile_reads.len() as f64
            };

            // host_reads may be missing; default to 0.
            let host_reads_val = record.host_reads.unwrap_or(0) as f64;
            // Compute host RPM.
            let host_rpm_val = ((host_reads_val + avg_human_profile) / input_reads) * 1e6;

            // Compute average profile reads for non-"Homo sapiens" taxa.
            let non_human_profile_reads: Vec<u64> = record
                .taxa
                .iter()
                .filter(|taxon| taxon.name != "Homo sapiens")
                .flat_map(|taxon| taxon.evidence.profile.iter().map(|profile| profile.reads))
                .collect();
            let avg_non_human_profile = if non_human_profile_reads.is_empty() {
                0.0
            } else {
                non_human_profile_reads.iter().sum::<u64>() as f64
                    / non_human_profile_reads.len() as f64
            };
            // Compute taxon RPM.
            let taxon_rpm_val = (avg_non_human_profile / input_reads) * 1e6;

            // Only it to the plot if both RPM are above zero for log10 conversion
            if host_rpm_val > 0. && taxon_rpm_val > 0. {
                raw_host_rpm.push(host_rpm_val);
                raw_taxon_rpm.push(taxon_rpm_val);

                // Save sample_id for later use.
                sample_ids.push(record.sample_id);
            }
        }

        // Create the analyzer from the computed raw RPM vectors and sample IDs.
        Self::new_with_sample_ids(config, raw_host_rpm, raw_taxon_rpm, sample_ids)
    }

    /// Run the regression analysis and identify outliers.
    pub fn run(&self) -> Result<RpmAnalysisResult, HttpClientError> {
        // Build the data map: keys "x" for independent variable and "y" for dependent variable.
        let mut data_map: HashMap<String, Vec<f64>> = HashMap::new();
        data_map.insert("x".to_string(), self.log_host_rpm.clone());
        data_map.insert("y".to_string(), self.log_taxon_rpm.clone());

        let data = RegressionDataBuilder::new()
            .build_from(data_map)
            .map_err(|_| HttpClientError::FailedRegressionBuild)?;

        // Fit the linear model with formula "y ~ x".
        let formula = "y ~ x";
        let model = FormulaRegressionBuilder::new()
            .data(&data)
            .formula(formula)
            .fit()
            .expect("Failed to fit regression model");

        let params = model.parameters();
        let intercept = params[0];
        let slope = params[1];
        let relationship = if slope > 0.0 {
            "Positive".to_string()
        } else {
            "Negative".to_string()
        };

        let n = self.log_host_rpm.len();
        let n_f64 = n as f64;
        let mean_x: f64 = self.log_host_rpm.iter().sum::<f64>() / n_f64;
        let sum_sq_x: f64 = self.log_host_rpm.iter().map(|x| (x - mean_x).powi(2)).sum();

        // Compute predicted values and accumulate squared residuals.
        let mut sum_sq_res = 0.0;
        let mut predictions = Vec::new();
        for (&x, &y) in self.log_host_rpm.iter().zip(self.log_taxon_rpm.iter()) {
            let pred = intercept + slope * x;
            predictions.push(pred);
            sum_sq_res += (y - pred).powi(2);
        }

        let degrees_of_freedom = n_f64 - 2.0;
        let residual_standard_error = (sum_sq_res / degrees_of_freedom).sqrt();

        // Determine the t value for the given confidence level.
        let quantile = 1.0 - (1.0 - self.config.confidence_level) / 2.0;
        let t_dist =
            StudentsT::new(0.0, 1.0, degrees_of_freedom).expect("Could not create t-distribution");
        let t_value = t_dist.inverse_cdf(quantile);

        let mut outliers = Vec::new();
        // For each observation, compute the interval and flag outliers.
        for (i, ((&x, &y), &pred)) in self
            .log_host_rpm
            .iter()
            .zip(self.log_taxon_rpm.iter())
            .zip(predictions.iter())
            .enumerate()
        {
            let interval_factor = if self.config.use_prediction_interval {
                1.0
            } else {
                0.0
            };
            let se_pred = residual_standard_error
                * (interval_factor + 1.0 / n_f64 + (x - mean_x).powi(2) / sum_sq_x).sqrt();
            let lower = pred - t_value * se_pred;
            let upper = pred + t_value * se_pred;
            // If only_high_outliers is set, only flag points with taxon RPM (y) above the upper interval.
            let is_outlier = if self.config.only_high_outliers {
                y > upper
            } else {
                y < lower || y > upper
            };
            if is_outlier {
                outliers.push(Outlier {
                    sample_id: self
                        .sample_ids
                        .get(i)
                        .cloned()
                        .unwrap_or_else(|| i.to_string()),
                    raw_host: self.raw_host_rpm[i],
                    raw_taxon: self.raw_taxon_rpm[i],
                    log_host: x,
                    log_taxon: y,
                    predicted: pred,
                    lower,
                    upper,
                });
            }
        }

        Ok(RpmAnalysisResult {
            intercept,
            slope,
            relationship,
            residual_standard_error,
            outliers,
        })
    }

    /// Creates a PNG plot of the regression.
    ///
    /// This plot shows:
    /// - The scatter points (log₁₀(host RPM) vs. log₁₀(taxon RPM)), with outliers colored red.
    /// - The regression line.
    /// - A shaded confidence (or prediction) interval.
    ///
    /// The legend has been removed.
    pub fn plot_regression(
        &self,
        result: &RpmAnalysisResult,
        filename: &str,
    ) -> Result<(), HttpClientError> {
        // Determine x-range from the log-transformed host RPM values.
        let x_min = self
            .log_host_rpm
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let x_max = self
            .log_host_rpm
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let n = self.log_host_rpm.len() as f64;
        let mean_x = self.log_host_rpm.iter().sum::<f64>() / n;
        let sum_sq_x = self
            .log_host_rpm
            .iter()
            .map(|&x| (x - mean_x).powi(2))
            .sum::<f64>();
        let dof = n - 2.0;
        let quantile = 1.0 - (1.0 - self.config.confidence_level) / 2.0;
        let t_dist = StudentsT::new(0.0, 1.0, dof)?;
        let t_value = t_dist.inverse_cdf(quantile);
        let interval_factor = if self.config.use_prediction_interval {
            1.0
        } else {
            0.0
        };

        // Generate points for regression line and confidence interval.
        let num_points = 100;
        let mut reg_line = Vec::new();
        let mut lower_bound = Vec::new();
        let mut upper_bound = Vec::new();
        for i in 0..=num_points {
            let x = x_min + (x_max - x_min) * (i as f64) / (num_points as f64);
            let pred = result.intercept + result.slope * x;
            reg_line.push((x, pred));
            let se = result.residual_standard_error
                * ((interval_factor + 1.0 / n + (x - mean_x).powi(2) / sum_sq_x).sqrt());
            lower_bound.push((x, pred - t_value * se));
            upper_bound.push((x, pred + t_value * se));
        }

        // Determine y-range based on data and computed interval.
        let y_min_data = self
            .log_taxon_rpm
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let y_max_data = self
            .log_taxon_rpm
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let y_min_interval = lower_bound
            .iter()
            .map(|&(_, y)| y)
            .fold(f64::INFINITY, f64::min);
        let y_max_interval = upper_bound
            .iter()
            .map(|&(_, y)| y)
            .fold(f64::NEG_INFINITY, f64::max);
        let y_min = y_min_data.min(y_min_interval);
        let y_max = y_max_data.max(y_max_interval);

        // Set up drawing area.
        let root_area = BitMapBackend::new(filename, (800, 600)).into_drawing_area();
        root_area.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root_area)
            .caption("RPM Regression", ("sans-serif", 30))
            .margin(20)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

        chart
            .configure_mesh()
            .x_desc("log10(host RPM)")
            .y_desc("log10(taxon RPM)")
            .draw()?;

        // Draw the confidence/prediction interval as a shaded area.
        let mut area_points: Vec<(f64, f64)> = lower_bound.clone();
        upper_bound
            .iter()
            .rev()
            .for_each(|&pt| area_points.push(pt));

        chart.draw_series(AreaSeries::new(
            area_points,
            0.0, // baseline (ignored)
            &RGBColor(200, 200, 255).mix(0.3),
        ))?;

        // Draw the regression line.
        chart.draw_series(LineSeries::new(reg_line, &BLUE))?;

        // Prepare a set of sample IDs that were flagged as outliers.
        let outlier_ids: HashSet<_> = result
            .outliers
            .iter()
            .map(|o| o.sample_id.clone())
            .collect();

        // Plot all data points; color them red if they are outliers.
        for (i, (&x, &y)) in self
            .log_host_rpm
            .iter()
            .zip(self.log_taxon_rpm.iter())
            .enumerate()
        {
            let point_color = if self
                .sample_ids
                .get(i)
                .map(|id| outlier_ids.contains(id))
                .unwrap_or(false)
            {
                &RED
            } else {
                &BLACK
            };
            chart.draw_series(std::iter::once(Circle::new(
                (x, y),
                3,
                ShapeStyle::from(point_color).filled(),
            )))?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn outlier(sample_id: &str) -> Outlier {
        Outlier {
            sample_id: sample_id.to_string(),
            raw_host: 1.0,
            raw_taxon: 1.0,
            log_host: 0.0,
            log_taxon: 0.0,
            predicted: 0.0,
            lower: 0.0,
            upper: 0.0,
        }
    }

    fn result_with(sample_ids: &[&str]) -> RpmAnalysisResult {
        RpmAnalysisResult {
            intercept: 0.0,
            slope: 1.0,
            relationship: "positive".to_string(),
            residual_standard_error: 0.0,
            outliers: sample_ids.iter().map(|s| outlier(s)).collect(),
        }
    }

    #[test]
    fn has_sample_outlier_matches_only_listed_samples() {
        let result = result_with(&["sampleA", "sampleC"]);
        assert!(result.has_sample_outlier("sampleA"));
        assert!(result.has_sample_outlier("sampleC"));
        assert!(!result.has_sample_outlier("sampleB"));
    }

    #[test]
    fn has_sample_outlier_false_when_no_outliers() {
        let result = result_with(&[]);
        assert!(!result.has_sample_outlier("sampleA"));
    }

    #[test]
    fn analyzer_drops_non_positive_points_without_panicking() {
        // Mixed positive/non-positive points: the constructor must not panic and must keep only
        // the positive points, index-aligned across host/taxon/sample_ids.
        let cfg = RpmConfig::default();
        let analyzer = RpmAnalyzer::new_with_sample_ids(
            cfg,
            vec![10.0, 0.0, 5.0, -1.0, 2.0],
            vec![100.0, 50.0, 0.0, 7.0, 20.0],
            vec![
                "s1".into(),
                "s2".into(),
                "s3".into(),
                "s4".into(),
                "s5".into(),
            ],
        );
        // Only s1 (10,100) and s5 (2,20) have BOTH coordinates positive.
        assert_eq!(analyzer.raw_host_rpm, vec![10.0, 2.0]);
        assert_eq!(analyzer.raw_taxon_rpm, vec![100.0, 20.0]);
        assert_eq!(analyzer.sample_ids, vec!["s1".to_string(), "s5".to_string()]);
        assert_eq!(analyzer.log_host_rpm.len(), 2);
        assert_eq!(analyzer.log_taxon_rpm.len(), 2);
        // log₁₀(10) == 1.0
        assert!((analyzer.log_host_rpm[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn analyzer_new_keeps_all_positive_points() {
        let analyzer =
            RpmAnalyzer::new(RpmConfig::default(), vec![10.0, 100.0], vec![1.0, 1000.0]);
        assert_eq!(analyzer.raw_host_rpm.len(), 2);
        assert!(analyzer.sample_ids.is_empty()); // new() supplies no ids
    }
}
