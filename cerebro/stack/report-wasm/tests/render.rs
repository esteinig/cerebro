//! Stage 11 — Typst render validation. Renders the regression template (PASS / REGRESSED) and
//! the gated META-GPT appendix to a PDF and asserts content. This exercises the part Stage 08
//! could not validate headlessly: that the templates actually *compile* and that the verdict /
//! colour-coding branches and the appendix gate behave.
//!
//! Requires the typst toolchain + fonts at runtime, so it is a CI integration test (it does not
//! run in a headless sandbox without fonts). Build feature: `lib` (default).

use std::path::Path;

use cerebro_report::{
    LibraryReportCompiler, MetaGptAssessment, PathogenDetectionReport, ReportType, TemplateFormat,
};

/// Render a config file to the (tera-substituted) Typst source, then compile it to PDF bytes.
fn render_pdf(report_type: ReportType, fixture: &str) -> (String, Vec<u8>) {
    let mut c = LibraryReportCompiler::new(".".to_string(), report_type)
        .expect("compiler init");
    let src = c
        .report_from_path(Path::new(fixture), TemplateFormat::Json)
        .expect("template render (tera)");
    let pdf = c.pdf_data(&src, "report.typ".to_string()).expect("typst compile");
    (src, pdf)
}

fn assert_pdf(bytes: &[u8]) {
    assert!(!bytes.is_empty(), "empty PDF");
    assert_eq!(&bytes[..5], b"%PDF-", "not a PDF");
}

#[test]
fn regression_pass_renders_and_reads_pass() {
    let (src, pdf) = render_pdf(ReportType::Regression, "tests/fixtures/regression_pass.json");
    assert_pdf(&pdf);
    // The rendered Typst source carries the substituted verdict + stats.
    assert!(src.contains("RUN-PASS"));
    assert!(src.contains("96.2%")); // current sensitivity, formatted
    assert!(src.contains("+0.4 pp")); // sensitivity delta
    // PASS fixture must not present as regressed.
    assert!(!src.to_lowercase().contains("regressed = true"));
}

#[test]
fn regression_regressed_renders_and_reads_regressed() {
    let (src, pdf) = render_pdf(
        ReportType::Regression,
        "tests/fixtures/regression_regressed.json",
    );
    assert_pdf(&pdf);
    assert!(src.contains("RUN-REGRESSED"));
    assert!(src.contains("-8.2 pp")); // sensitivity delta below threshold
    assert!(src.contains("McNemar significant") || src.contains("0.0120"));
}

#[test]
fn pathogen_appendix_present_when_enabled_absent_when_not() {
    // Build configs in-test (no fixture file): default report has no meta_gpt block.
    let mut off = PathogenDetectionReport::default();
    off.meta_gpt = None;

    let mut on = PathogenDetectionReport::default();
    on.meta_gpt = Some(MetaGptAssessment {
        enabled: true,
        ..Default::default()
    });

    let mut c = LibraryReportCompiler::new(".".to_string(), ReportType::PathogenDetection)
        .expect("compiler init");

    let src_off = c.report(&off).expect("render off");
    let src_on = c.report(&on).expect("render on");

    // Appendix C is gated: present only when enabled (the back-compat guarantee, rendered).
    assert!(!src_off.contains("Appendix C"), "appendix leaked when disabled");
    assert!(src_on.contains("Appendix C"), "appendix missing when enabled");

    // And both still compile to PDF.
    assert_pdf(&c.pdf_data(&src_off, "off.typ".to_string()).expect("compile off"));
    assert_pdf(&c.pdf_data(&src_on, "on.typ".to_string()).expect("compile on"));
}
