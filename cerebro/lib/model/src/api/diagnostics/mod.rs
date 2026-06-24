//! Shared diagnostic evaluation, META-GPT output reading, and regression-against-baseline.
//!
//! The single home for sensitivity/specificity scoring (`eval`), for reading META-GPT output
//! without depending on the `meta-gpt` crate (`metagpt`), and for regression of a run's
//! statistics against a versioned baseline (`regression`). Training and CIQA delegate their
//! counting/statistics here so they report identical numbers from one implementation.

pub mod eval;
pub mod metagpt;
pub mod regression;
