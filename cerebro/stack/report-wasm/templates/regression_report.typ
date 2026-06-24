// regression_report.typ — Cerebro CIQA regression report (QC dataset / baseline)
//
// Rendered by `cerebro-report` (the WASM/Typst module) via tera, then compiled to PDF by the
// `typst` crate — the same pipeline as templates/pathogen_detection.typ. The idioms below
// (page frame, `outlinebox`, 4-column header tables, the colour-coded indicator, and the
// `table.header(...)` multi-page table) mirror that template deliberately.
//
// STATUS: correct-by-inspection against templates/pathogen_detection.typ; not compiled here.
// Validate by compiling with the typst crate (native + wasm32) over a PASS and a REGRESSED
// fixture and inspecting the PDF, per 08-STAGE-REPORTING.md §S8.5.

#set text(font: "IBM Plex Sans")

// ── VARIABLES (filled by build_context / tera) ──────────────────────────────

#let report_header_logo_enabled = {{ report_header_logo_enabled }};

#let regression_header_dataset        = "{{ regression_header_dataset }}";
#let regression_header_dataset_version= "{{ regression_header_dataset_version }}";
#let regression_header_baseline_id    = "{{ regression_header_baseline_id }}";
#let regression_header_model          = "{{ regression_header_model }}";
#let regression_header_quantization   = "{{ regression_header_quantization }}";
#let regression_header_run_id         = "{{ regression_header_run_id }}";
#let regression_header_evaluated      = "{{ regression_header_evaluated }}";
#let regression_header_config_hash    = "{{ regression_header_config_hash }}";

#let regression_verdict_regressed = {{ regression_verdict_regressed }};   // bool
#let regression_verdict_passed    = {{ regression_verdict_passed }};      // bool (threshold floor)
#let regression_verdict_reasons   = "{{ regression_verdict_reasons }}";

#let stats_current_sensitivity  = "{{ stats_current_sensitivity }}";
#let stats_current_specificity  = "{{ stats_current_specificity }}";
#let stats_current_ppv          = "{{ stats_current_ppv }}";
#let stats_current_npv          = "{{ stats_current_npv }}";
#let stats_baseline_sensitivity = "{{ stats_baseline_sensitivity }}";
#let stats_baseline_specificity = "{{ stats_baseline_specificity }}";
#let stats_baseline_ppv         = "{{ stats_baseline_ppv }}";
#let stats_baseline_npv         = "{{ stats_baseline_npv }}";
#let stats_delta_sensitivity    = "{{ stats_delta_sensitivity }}";
#let stats_delta_specificity    = "{{ stats_delta_specificity }}";
#let stats_threshold_sensitivity= "{{ stats_threshold_sensitivity }}";
#let stats_threshold_specificity= "{{ stats_threshold_specificity }}";

#let mcnemar_p_value     = "{{ mcnemar_p_value }}";
#let mcnemar_discordant  = "{{ mcnemar_discordant }}";
#let mcnemar_adjustment  = "{{ mcnemar_adjustment }}";
#let mcnemar_panel_size  = "{{ mcnemar_panel_size }}";

// per-sample decisions: (sample_id, reference_result, predicted_result, decision)
#let regression_decisions = (
  {% for row in regression_decisions %}("{{ row.sample_id }}","{{ row.reference }}","{{ row.predicted }}","{{ row.decision }}"),{% endfor %}
);

// ── PAGE FRAME (same structure as pathogen_detection.typ) ───────────────────

#let header_title = [
    #set text(11pt, weight: "bold");
    #smallcaps[METAGENOMICS];
    #v(-5pt);
    #smallcaps[CIQA REGRESSION REPORT];
    #v(-5pt);
    #set text(8pt, weight: "medium")
    #overline(offset: -1.2em)[#smallcaps[NOT FOR DIAGNOSTIC USE]]
]

#let logo_title = [ #set text(11pt, weight: "bold"); #smallcaps[REPORT] ]
#let logo = [
  #if report_header_logo_enabled { image("logo.png", width: {{ report_header_logo_width }}, fit: "stretch") } else { logo_title }
]

#set par(justify: true)
#set page(
  margin: (top: 96pt, bottom: 120pt),
  header: grid(columns: (1fr, 1fr), align(left)[#header_title], align(right)[#logo]),
  footer: context [
    #set text(8pt);
    #line(length: 100%);
    Page #counter(page).display("1 / 1", both: true)
    #h(1fr) Dataset: #regression_header_dataset (#regression_header_dataset_version)
    #h(1fr) Baseline: #regression_header_baseline_id
    #h(1fr) Evaluated: #regression_header_evaluated
  ],
  numbering: "1 of 1"
)

// ── BOX FUNCTION (identical to pathogen_detection.typ) ──────────────────────

#let outlinebox(title: "title", body) = {
  return block(stroke: 0.25pt + luma(40), radius: 4pt, width: 100%, above: 26pt)[
    #place(top + start, dy: -12pt, dx: 10pt)[
      #block(stroke: 0.25pt + luma(40), fill: rgb(255,255,255), inset: 8pt, radius: 4pt)[
        #text(fill: rgb(0,0,0), weight: "bold")[#title]
      ]
    ]
    #block(width: 100%, inset: (top: 30pt, x: 10pt, bottom: 25pt))[ #body ]
  ]
}

// ── HEADER TABLE (4-col, same idiom as the appendix headers) ────────────────

#let regression_header = {
  show table.cell.where(x: 0): set text(weight: 500)
  show table.cell.where(x: 2): set text(weight: 500)
  table(
    columns: 4,
    [Dataset], [#regression_header_dataset], [Run ID], [#regression_header_run_id],
    [Dataset Version], [#regression_header_dataset_version], [Model], [#regression_header_model],
    [Baseline], [#regression_header_baseline_id], [Quantization], [#regression_header_quantization],
    [Evaluated], [#regression_header_evaluated], [Config Hash], [#regression_header_config_hash],
  )
}

// ── VERDICT (colour-coded, same styling family as the result block) ─────────

#let regression_verdict = [
  #if regression_verdict_regressed {
    [
      #set text(size: 18pt, weight: "extrabold", fill: rgb(150, 0, 0));
      #smallcaps[REGRESSED]
    ]
  } else {
    [
      #set text(size: 18pt, weight: "extrabold", fill: rgb(0, 90, 0));
      #smallcaps[PASS]
    ]
  }
  #v(4pt)
  #text(size: 10pt)[#regression_verdict_reasons]
  #v(2pt)
  #text(size: 8pt, weight: "thin", stroke: 0.5pt + rgb(170,170,170))[
    A regression is recorded when the absolute sensitivity/specificity floor is breached or a
    significant adverse paired shift is detected. This report does not modify the baseline.
  ]
]

// ── STATISTICS (current vs baseline + deltas + thresholds) ──────────────────

#let regression_statistics = {
  show table.cell.where(y: 0): set text(weight: "bold")
  table(
    columns: 5,
    [Metric], [Current], [Baseline], [Δ], [Threshold],
    [Sensitivity], [#stats_current_sensitivity], [#stats_baseline_sensitivity], [#stats_delta_sensitivity], [#stats_threshold_sensitivity],
    [Specificity], [#stats_current_specificity], [#stats_baseline_specificity], [#stats_delta_specificity], [#stats_threshold_specificity],
    [PPV], [#stats_current_ppv], [#stats_baseline_ppv], [—], [—],
    [NPV], [#stats_current_npv], [#stats_baseline_npv], [—], [—],
  )
}

// ── SIGNIFICANCE (McNemar + power caveat) ───────────────────────────────────

#let regression_significance = [
  McNemar p-value: #text(weight: "bold")[#mcnemar_p_value]
  #h(2fr) Discordant pairs: #mcnemar_discordant
  #h(2fr) Adjustment: #mcnemar_adjustment
  #h(2fr) Panel size: #mcnemar_panel_size

  #v(4pt)
  #text(size: 8pt, weight: "thin", stroke: 0.5pt + rgb(100,100,100))[
    A non-significant result is not evidence of no change on a small or unbalanced panel; the
    absolute threshold floor is the primary gate.
  ]
]

// ── PER-SAMPLE DECISIONS (multi-page; header repeats) ───────────────────────

#let regression_decision_table = [
  #set text(size: 9pt)
  #table(
    columns: (1.6fr, 1fr, 1fr, auto),
    table.header([Sample], [Reference], [Predicted], [Decision]),
    ..regression_decisions.map(row => (
      row.at(0), row.at(1), row.at(2),
      // emphasise discordant outcomes
      if row.at(3) == "FP" or row.at(3) == "FN" { text(weight: "bold", fill: rgb(150,0,0))[#row.at(3)] } else { [#row.at(3)] },
    )).flatten()
  )
  #v(4pt)
  #text(size: 8pt, weight: "thin", stroke: 0.5pt + rgb(170,170,170))[
    TP/TN/FP/FN/Excluded per sample, from the shared diagnostic evaluator. This table may continue
    across multiple pages.
  ]
]

// ── COMPOSITION ─────────────────────────────────────────────────────────────

#regression_header
#v(16pt)
#outlinebox(title: "Verdict")[#regression_verdict]
#outlinebox(title: "Statistics — Current vs Baseline")[#regression_statistics]
#outlinebox(title: "Paired Significance")[#regression_significance]
#pagebreak()
= Per-sample Decisions
#regression_decision_table
