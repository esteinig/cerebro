use std::collections::HashMap;
use statrs::distribution::{Binomial, ChiSquared, Discrete, ContinuousCDF};
use crate::plate::{DiagnosticOutcome, DiagnosticReview};

pub trait HasSampleOutcome {
    fn sample_id(&self) -> &str;
    fn outcome(&self) -> &DiagnosticOutcome;
}

impl HasSampleOutcome for DiagnosticReview {
    fn sample_id(&self) -> &str {
        &self.sample_id
    }
    fn outcome(&self) -> &DiagnosticOutcome {
        &self.outcome
    }
}


pub fn mcnemar_from_reviews<T: HasSampleOutcome, U: HasSampleOutcome>(
    reviews_a: &[T],
    reviews_b: &[U],
) -> McNemarResult {
    let (b, c, total) = compute_bc_counts(reviews_a, reviews_b);
    mcnemar_test(b, c, total)
}

pub fn compute_bc_counts<T: HasSampleOutcome, U: HasSampleOutcome>(
    config_a: &[T],
    config_b: &[U],
) -> (u32, u32, u32) {
    let a_map: HashMap<_, _> = config_a
        .iter()
        .filter_map(|r| Some((r.sample_id().to_string(), is_correct(r.outcome())?)))
        .collect();

    let b_map: HashMap<_, _> = config_b
        .iter()
        .filter_map(|r| Some((r.sample_id().to_string(), is_correct(r.outcome())?)))
        .collect();

    let mut b = 0;
    let mut c = 0;
    let mut n = 0;

    for (sample_id, &a_correct) in &a_map {
        if let Some(&b_correct) = b_map.get(sample_id) {
            n += 1;
            match (a_correct, b_correct) {
                (true, false) => b += 1,
                (false, true) => c += 1,
                _ => {}
            }
        }
    }

    (b, c, n)
}

fn is_correct(outcome: &DiagnosticOutcome) -> Option<bool> {
    match outcome {
        DiagnosticOutcome::TruePositive | DiagnosticOutcome::TrueNegative => Some(true),
        DiagnosticOutcome::FalsePositive | DiagnosticOutcome::FalseNegative => Some(false),
        _ => None, // Indeterminate, Unknown, etc.
    }
}

#[derive(Debug)]
pub struct McNemarResult {
    pub test_type: String,
    pub p_value: f64,
    pub b: u32,
    pub c: u32,
    pub total: u32,
    pub total_discordant: u32,
}

pub fn mcnemar_test(b: u32, c: u32, total: u32) -> McNemarResult {
    let n: u32 = b + c;

    if n == 0 {
        return McNemarResult {
            test_type: "No discordant pairs".to_string(),
            p_value: 1.0,
            b,
            c,
            total,
            total_discordant: n,
        };
    }

    if n >= 25 {
        // Chi-square approximation
        let chi_sq = ((b as f64 - c as f64).powi(2)) / (b + c) as f64;
        let chi_dist = ChiSquared::new(1.0).unwrap();
        let p = 1.0 - chi_dist.cdf(chi_sq);

        McNemarResult {
            test_type: "Chi-squared approximation".to_string(),
            p_value: p,
            b,
            c,
            total,
            total_discordant: n,
        }
    } else {
        // Mid-p exact test
        let k = b.min(c);
        let binom = Binomial::new(0.5, n.into()).unwrap();

        let p_tail = (0..k).map(|i| binom.pmf(i.into())).sum::<f64>();
        let p_exact = 0.5 * binom.pmf(k.into());
        let mid_p = 2.0 * (p_tail + p_exact);

        McNemarResult {
            test_type: "Exact mid-p binomial test".to_string(),
            p_value: mid_p.min(1.0), // p-values can't exceed 1
            b,
            c,
            total,
            total_discordant: n,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-8;

    fn approx_eq(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    struct DiagnosticReviewTest {
        sample_id: String,
        outcome: DiagnosticOutcome
    }

    impl HasSampleOutcome for DiagnosticReviewTest {
        fn sample_id(&self) -> &str {
            &self.sample_id
        }
        fn outcome(&self) -> &DiagnosticOutcome {
            &self.outcome
        }
    }

    #[test]
    fn test_mcnemar_from_reviews_exact_mid_p() {
        let reviews_a = vec![
            DiagnosticReviewTest {
                sample_id: "s1".into(),
                outcome: DiagnosticOutcome::TruePositive,  // correct
            },
            DiagnosticReviewTest {
                sample_id: "s2".into(),
                outcome: DiagnosticOutcome::FalseNegative, // incorrect
            },
            DiagnosticReviewTest {
                sample_id: "s3".into(),
                outcome: DiagnosticOutcome::TrueNegative,  // correct
            },
            DiagnosticReviewTest {
                sample_id: "s4".into(),
                outcome: DiagnosticOutcome::TruePositive,  // correct
            },
        ];

        let reviews_b = vec![
            DiagnosticReviewTest {
                sample_id: "s1".into(),
                outcome: DiagnosticOutcome::TruePositive,  // correct
            },
            DiagnosticReviewTest {
                sample_id: "s2".into(),
                outcome: DiagnosticOutcome::FalsePositive, // incorrect
            },
            DiagnosticReviewTest {
                sample_id: "s3".into(),
                outcome: DiagnosticOutcome::FalsePositive, // incorrect
            },
            DiagnosticReviewTest {
                sample_id: "s4".into(),
                outcome: DiagnosticOutcome::TruePositive,  // correct
            },
        ];

        let result = mcnemar_from_reviews(&reviews_a, &reviews_b);
        assert_eq!(result.test_type, "Exact mid-p binomial test");
        assert_eq!(result.b, 1); // A correct, B wrong (s3)
        assert_eq!(result.c, 0); // B correct, A wrong (s2)
        assert_eq!(result.total_discordant, 1);
        assert!(result.p_value <= 1.0);
    }

    #[test]
    fn test_mcnemar_from_reviews_all_agree_correct() {
        let reviews_a = vec![
            DiagnosticReviewTest {
                sample_id: "x1".into(),
                outcome: DiagnosticOutcome::TruePositive,  // correct
            },
            DiagnosticReviewTest {
                sample_id: "x2".into(),
                outcome: DiagnosticOutcome::TrueNegative,  // correct
            },
        ];

        let reviews_b = vec![
            DiagnosticReviewTest {
                sample_id: "x1".into(),
                outcome: DiagnosticOutcome::TruePositive,  // correct
            },
            DiagnosticReviewTest {
                sample_id: "x2".into(),
                outcome: DiagnosticOutcome::TrueNegative,  // correct
            },
        ];

        let result = mcnemar_from_reviews(&reviews_a, &reviews_b);
        assert_eq!(result.b, 0);
        assert_eq!(result.c, 0);
        assert_eq!(result.total_discordant, 0);
        assert!(approx_eq(result.p_value, 1.0, EPSILON));
        assert_eq!(result.test_type, "No discordant pairs");
    }


    #[test]
    fn test_mcnemar_from_reviews_all_agree_incorrect() {
        let reviews_a = vec![
            DiagnosticReviewTest {
                sample_id: "x1".into(),
                outcome: DiagnosticOutcome::FalsePositive,  // incorrect
            },
            DiagnosticReviewTest {
                sample_id: "x2".into(),
                outcome: DiagnosticOutcome::FalseNegative,   // incorrect
            },
        ];

        let reviews_b = vec![
            DiagnosticReviewTest {
                sample_id: "x1".into(),
                outcome: DiagnosticOutcome::FalseNegative,  // incorrect
            },
            DiagnosticReviewTest {
                sample_id: "x2".into(),
                outcome: DiagnosticOutcome::FalsePositive,   // incorrect
            },
        ];

        let result = mcnemar_from_reviews(&reviews_a, &reviews_b);
        assert_eq!(result.b, 0);
        assert_eq!(result.c, 0);
        assert_eq!(result.total_discordant, 0);
        assert!(approx_eq(result.p_value, 1.0, EPSILON));
        assert_eq!(result.test_type, "No discordant pairs");
    }

    #[test]
    fn test_mcnemar_realistic_small_model_improvement() {
        let mut reviews_a = vec![];
        let mut reviews_b = vec![];

        // 20 samples where both configs agree (true or false)
        for i in 0..10 {
            let id = format!("same_correct_{}", i);
            reviews_a.push(DiagnosticReviewTest {
                sample_id: id.clone(),
                outcome: DiagnosticOutcome::TruePositive,
            });
            reviews_b.push(DiagnosticReviewTest {
                sample_id: id,
                outcome: DiagnosticOutcome::TruePositive,
            });
        }
        for i in 0..10 {
            let id = format!("same_incorrect_{}", i);
            reviews_a.push(DiagnosticReviewTest {
                sample_id: id.clone(),
                outcome: DiagnosticOutcome::FalseNegative,
            });
            reviews_b.push(DiagnosticReviewTest {
                sample_id: id,
                outcome: DiagnosticOutcome::FalseNegative,
            });
        }

        // 6 samples where A is correct, B is incorrect → b += 6
        for i in 0..6 {
            let id = format!("a_correct_b_wrong_{}", i);
            reviews_a.push(DiagnosticReviewTest {
                sample_id: id.clone(),
                outcome: DiagnosticOutcome::TrueNegative,
            });
            reviews_b.push(DiagnosticReviewTest {
                sample_id: id,
                outcome: DiagnosticOutcome::FalsePositive,
            });
        }

        // 3 samples where B is correct, A is incorrect → c += 3
        for i in 0..3 {
            let id = format!("b_correct_a_wrong_{}", i);
            reviews_a.push(DiagnosticReviewTest {
                sample_id: id.clone(),
                outcome: DiagnosticOutcome::FalsePositive,
            });
            reviews_b.push(DiagnosticReviewTest {
                sample_id: id,
                outcome: DiagnosticOutcome::TrueNegative,
            });
        }

        let result = mcnemar_from_reviews(&reviews_a, &reviews_b);

        println!("{result:#?}");

        assert_eq!(result.test_type, "Exact mid-p binomial test");
        assert_eq!(result.b, 6);
        assert_eq!(result.c, 3);
        assert_eq!(result.total_discordant, 9);
        assert!(result.p_value >= 0.0 && result.p_value <= 1.0);
    }
    
    #[test]
    fn test_mcnemar_realistic_chi_squared_case() {
        let mut reviews_a = vec![];
        let mut reviews_b = vec![];

        // 10 samples where both are correct
        for i in 0..10 {
            let id = format!("same_correct_{}", i);
            reviews_a.push(DiagnosticReviewTest {
                sample_id: id.clone(),
                outcome: DiagnosticOutcome::TrueNegative,
            });
            reviews_b.push(DiagnosticReviewTest {
                sample_id: id,
                outcome: DiagnosticOutcome::TrueNegative,
            });
        }

        // 10 samples where both are incorrect
        for i in 0..10 {
            let id = format!("same_incorrect_{}", i);
            reviews_a.push(DiagnosticReviewTest {
                sample_id: id.clone(),
                outcome: DiagnosticOutcome::FalsePositive,
            });
            reviews_b.push(DiagnosticReviewTest {
                sample_id: id,
                outcome: DiagnosticOutcome::FalsePositive,
            });
        }

        // 30 samples where A correct, B wrong → b += 30
        for i in 0..30 {
            let id = format!("a_correct_b_wrong_{}", i);
            reviews_a.push(DiagnosticReviewTest {
                sample_id: id.clone(),
                outcome: DiagnosticOutcome::TruePositive,
            });
            reviews_b.push(DiagnosticReviewTest {
                sample_id: id,
                outcome: DiagnosticOutcome::FalseNegative,
            });
        }

        // 5 samples where B correct, A wrong → c += 5
        for i in 0..5 {
            let id = format!("b_correct_a_wrong_{}", i);
            reviews_a.push(DiagnosticReviewTest {
                sample_id: id.clone(),
                outcome: DiagnosticOutcome::FalsePositive,
            });
            reviews_b.push(DiagnosticReviewTest {
                sample_id: id,
                outcome: DiagnosticOutcome::TrueNegative,
            });
        }

        let result = mcnemar_from_reviews(&reviews_a, &reviews_b);

        assert_eq!(result.test_type, "Chi-squared approximation");
        assert_eq!(result.b, 30);
        assert_eq!(result.c, 5);
        assert_eq!(result.total_discordant, 35);
        assert!(result.p_value >= 0.0 && result.p_value <= 1.0);
    }

    #[test]
    fn test_chi_squared_case() {
        let result = mcnemar_test(20, 10, 0);
        assert_eq!(result.test_type, "Chi-squared approximation");
        assert!(approx_eq(result.total_discordant as f64, 30.0, EPSILON));
        assert!(result.p_value >= 0.0 && result.p_value <= 1.0);
    }

    #[test]
    fn test_mid_p_case_symmetric() {
        // Perfect symmetry: b == c → mid-p ≈ 1.0
        let result = mcnemar_test(4, 4, 0);
        assert_eq!(result.test_type, "Exact mid-p binomial test");
        assert!(approx_eq(result.p_value, 1.0, EPSILON));
        assert_eq!(result.total_discordant, 8);
    }

    #[test]
    fn test_mid_p_case_asymmetric() {
        // Highly asymmetric → lower p-value
        let result = mcnemar_test(5, 0, 0);
        assert_eq!(result.test_type, "Exact mid-p binomial test");
        assert!(result.p_value < 0.05);
        assert_eq!(result.total_discordant, 5);
    }

    #[test]
    fn test_no_discordant_pairs() {
        let result = mcnemar_test(0, 0, 0);
        assert_eq!(result.test_type, "No discordant pairs");
        assert!(approx_eq(result.p_value, 1.0, EPSILON));
        assert_eq!(result.total_discordant, 0);
    }

    #[test]
    fn test_upper_bound() {
        let result = mcnemar_test(12, 12, 0);
        assert!(result.p_value <= 1.0);
        assert!(result.p_value >= 0.0);
    }
}