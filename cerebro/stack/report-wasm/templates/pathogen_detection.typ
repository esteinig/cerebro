#set text(font: "IBM Plex Sans")

// VARIABLES

#let report_header_logo_enabled = {{ report_header_logo_enabled }};

#let report_footer_patient_id = "{{ report_footer_patient_id }}";
#let report_footer_date_collected = "{{ report_footer_date_collected }}";
#let report_footer_date_reported = "{{ report_footer_date_reported }}";
#let report_footer_reporting_location = "{{ report_footer_reporting_location }}";

#let patient_header_patient_name = "{{ patient_header_patient_name }}";
#let patient_header_patient_dob = "{{ patient_header_patient_dob }}";
#let patient_header_patient_urn = "{{ patient_header_patient_urn }}";
#let patient_header_requested_doctor = "{{ patient_header_requested_doctor }}";
#let patient_header_hospital_site = "{{ patient_header_hospital_site }}";
#let patient_header_laboratory_number = "{{ patient_header_laboratory_number }}";
#let patient_header_specimen_id = "{{ patient_header_specimen_id }}";
#let patient_header_date_collected = "{{ patient_header_date_collected }}";
#let patient_header_date_received = "{{ patient_header_date_received }}";
#let patient_header_specimen_type = "{{ patient_header_specimen_type }}";
#let patient_header_reporting_laboratory = "{{ patient_header_reporting_laboratory }}";
#let patient_header_reporting_date = "{{ patient_header_reporting_date }}";

#let patient_result_pathogen_detected = {{ patient_result_pathogen_detected }};
#let patient_result_pathogen_reported = "{{ patient_result_pathogen_reported }}";
#let patient_result_review_date = "{{ patient_result_review_date }}";
#let patient_result_comments = "{{ patient_result_comments }}";
#let patient_result_orthogonal_tests = "{{ patient_result_orthogonal_tests }}";
#let patient_result_clinical_notes = "{{ patient_result_clinical_notes}}";
#let patient_result_actions = "{{ patient_result_actions }}";
#let patient_result_contact_name = "{{ patient_result_contact_name }}";
#let patient_result_contact_email = "{{ patient_result_contact_email }}";

#let report_legal_disclaimer = "{{ report_legal_disclaimer }}";
#let report_legal_disclosure = "{{ report_legal_disclosure }}";
#let report_legal_liability = "{{ report_legal_liability }}";

#let report_authorisation_signatures = (
  {% for signature in report_authorisation_signatures %} ("{{ signature.name }}", "{{ signature.position }}", "{{ signature.institution }}"),{% endfor %}
);

#let report_authorisation_laboratory = "{{ report_authorisation_laboratory }}";
#let report_authorisation_identifier = "{{ report_authorisation_identifier }}";


#let appendix_laboratory_enabled = {{ appendix_laboratory_enabled }};
#let appendix_laboratory_description = "{{ appendix_laboratory_description }}";
#let appendix_laboratory_comments = "{{ appendix_laboratory_comments }}";

#let appendix_laboratory_header_protocol = "{{ appendix_laboratory_header_protocol }}";
#let appendix_laboratory_header_sample_id = "{{ appendix_laboratory_header_sample_id }}";
#let appendix_laboratory_header_version = "{{ appendix_laboratory_header_version }}";
#let appendix_laboratory_header_run_id = "{{ appendix_laboratory_header_run_id }}";
#let appendix_laboratory_header_extraction = "{{ appendix_laboratory_header_extraction }}";
#let appendix_laboratory_header_extraction_control = "{{ appendix_laboratory_header_extraction_control }}";
#let appendix_laboratory_header_rna_depletion = "{{ appendix_laboratory_header_rna_depletion }}";
#let appendix_laboratory_header_library_control = "{{ appendix_laboratory_header_library_control }}";
#let appendix_laboratory_header_adapter = "{{ appendix_laboratory_header_adapter }}";
#let appendix_laboratory_header_sequencing_control = "{{ appendix_laboratory_header_sequencing_control }}";
#let appendix_laboratory_header_library = "{{ appendix_laboratory_header_library }}";
#let appendix_laboratory_header_negative_control = "{{ appendix_laboratory_header_negative_control }}";
#let appendix_laboratory_header_sequencer = "{{ appendix_laboratory_header_sequencer }}";
#let appendix_laboratory_header_positive_control = "{{ appendix_laboratory_header_positive_control }}";


#let appendix_bioinformatics_enabled = {{ appendix_bioinformatics_enabled }};
#let appendix_bioinformatics_description = "{{ appendix_bioinformatics_description }}";
#let appendix_bioinformatics_comments = "{{ appendix_bioinformatics_comments }}";
#let appendix_bioinformatics_libraries = "{{ appendix_bioinformatics_libraries }}";
#let appendix_bioinformatics_evidence = "{{ appendix_bioinformatics_evidence }}";

#let appendix_bioinformatics_header_pipeline = "{{ appendix_bioinformatics_header_pipeline }}";
#let appendix_bioinformatics_header_version = "{{ appendix_bioinformatics_header_version }}";
#let appendix_bioinformatics_header_pipeline_id = "{{ appendix_bioinformatics_header_pipeline_id }}";
#let appendix_bioinformatics_header_configuration = "{{ appendix_bioinformatics_header_configuration }}";
#let appendix_bioinformatics_header_started = "{{ appendix_bioinformatics_header_started }}";
#let appendix_bioinformatics_header_completed = "{{ appendix_bioinformatics_header_completed }}";
#let appendix_bioinformatics_header_sample_id = "{{ appendix_bioinformatics_header_sample_id }}";
#let appendix_bioinformatics_header_libraries = "{{ appendix_bioinformatics_header_libraries }}";
#let appendix_bioinformatics_header_databases = "{{ appendix_bioinformatics_header_databases }}";
#let appendix_bioinformatics_header_taxonomy = "{{ appendix_bioinformatics_header_taxonomy }}";

// PAGE HEADER

#let header_title = [
    #set text(11pt, weight: "bold");
    #smallcaps[METAGENOMICS];
    #v(-5pt);
    #smallcaps[CLINICAL REPORT];
    #v(-5pt);
    #set text(8pt, weight: "medium")
    #overline(offset: -1.2em)[#smallcaps[NOT FOR DIAGNOSTIC USE]]
]

#let logo_title = [
    #set text(11pt, weight: "bold");
    #smallcaps[REPORT]
]

#let logo = [
  #if report_header_logo_enabled {
    image("logo.png", width: 11%, fit: "stretch")
  } else {
    logo_title
  }
]

// PAGE FOOTER

#let page_disclaimer = [
  *Disclaimer*: #report_legal_disclaimer
]

// REPORT LAYOUT

#set par(justify: true)
#set page(
  margin: (top: 96pt, bottom: 192pt),
  header: grid(
    columns: (1fr, 1fr),
    align(left)[#header_title],
    align(right)[#logo]
  ),
  footer:context [
    #set text(10pt);
    #page_disclaimer;
    #line(length: 100%);
    #set text(8pt);
    Page
    #counter(page).display(
      "1 / 1",
      both: true,
    )
    #h(1fr) Patient ID: #report_footer_patient_id
    #h(1fr) Collected: #report_footer_date_collected
    #h(1fr) Reported: #report_footer_date_reported
    #h(1fr) Location: #report_footer_reporting_location
      
    
  ],
  numbering: "1 of 1"
)


// BOX FUNCTION

#let outlinebox(title: "title", body) = {
  let strokeColor = luma(40)
  return block(
      stroke: 0.25pt + luma(40),
      radius: 4pt,
      width: 100%,
      above: 26pt
    )[
      #place(top + start, dy: -12pt, dx:10pt)[
        #block(
          stroke: 0.25pt + luma(40),
          fill: rgb(255, 255, 255),
          inset: 8pt,
          radius: 4pt,
        )[
          #text(fill: rgb(0, 0, 0), weight: "bold")[#title]
        ]
      ]
      #block(
        width: 100%,
        inset: ( top: 30pt, x: 10pt, bottom: 25pt)
      )[
        #body
      ]
    ]
}

// FIRST PAGE

#let first_page_header = {
  show table.cell.where(x: 0): set text(weight: 500)
  show table.cell.where(x: 2): set text(weight: 500)
  table(
    columns: 4,
    [Patient Name], [#patient_header_patient_name], [Specimen ID], [#patient_header_specimen_id],
    [Date of Birth], [#patient_header_patient_dob], [Specimen Type], [#patient_header_specimen_type],
    [URN Number], [#patient_header_patient_urn], [Date Collected], [#patient_header_date_collected],
    [Requested Doctor], [#patient_header_requested_doctor], [Date Received], [#patient_header_date_received],
    [Hospital Site], [#patient_header_hospital_site], [Reporting Laboratory], [#patient_header_reporting_laboratory],
    [Laboratory Number], [#patient_header_laboratory_number], [Reporting Date], [#patient_header_reporting_date],
  )
}

#let first_page_referral = [
  
    #set text(11pt);
    #set par(justify: true);
    
    #let referral = [
      Thank you for referring your #text(11pt)[PATIENT] #text(weight: "bold")[#patient_header_patient_name] #text(11pt)[DOB] #text(11pt, weight: "bold")[#patient_header_patient_dob] to the metagenomic next‐generation sequencing for improving infectious disease diagnostics and pathogen discovery clinical research project (#text(weight: "medium")[MetaGP ‐ HREC/74254/2021.110]).
    ]
    #let results = [
      Metagenomic sequencing and analysis is complete. Your patient's case and data were reviewed and discussed during a multi‐disciplinary research team meeting on #text(weight: "medium")[#patient_result_review_date].
    ]
    #let note = [
      Please refer to the #text(10pt, weight: "medium")[APPENDIX] for details regarding the assay and sequencing results.
    ]

    #referral;
    
    #results;
    #if (appendix_bioinformatics_enabled or appendix_laboratory_enabled) {
      note;
    }
]

#let first_page_results = [ 
  
  #if patient_result_pathogen_detected {
    [
      #set text(size: 16pt, weight: "extrabold", style: "italic");
      #smallcaps[#patient_result_pathogen_reported]
    ]
  } else {
    [
      #set text(size: 16pt, weight: "extrabold");
      #smallcaps[No pathogen detected];
      
      #text(size: 8pt, weight: "light", stroke: 0.5pt + rgb(100, 100, 100))[No pathogen detected] #text(size: 8pt, weight: "light", stroke: 0.5pt + rgb(170, 170, 170))[does not exclude the possibility of a pathogen being present.]
    ]
  }
  
  #text(size: 8pt, weight: "thin", stroke: 0.5pt + rgb(170, 170, 170))[
    A pathogen is reported when a detection at the taxonomic rank of species has been reviewed.
  ]
  
]

#let first_page_comments = [
  #patient_result_comments
]

// SECOND PAGE COMPONENTS

#let second_page_actions = [
  #patient_result_actions
]

#let second_page_orthogonal = [
  #patient_result_orthogonal_tests

  #text(size: 8pt, weight: "thin", stroke: 0.5pt + rgb(100, 100, 100))[
    Orthogonal test results are supplementary information not covered by this diagnostic test.
  ]
]

#let second_page_clinical = [
  #patient_result_clinical_notes


  #text(size: 8pt, weight: "thin", stroke: 0.5pt + rgb(100, 100, 100))[
    Patient clinical notes are supplementary information not covered by this diagnostic test.
  ]
]

#let second_page_contact = [
  For any further questions regarding these results, please contact:
  
  #text(size: 10pt, weight: "extrabold")[#patient_result_contact_name (#patient_result_contact_email)] 
]


#let second_page_disclosure = [
  #report_legal_disclosure;
]

#let second_page_liability = [
  #report_legal_liability;
]

// THIRD PAGE COMPONENTS

#let third_page_authorisation = [
  Your patient's case and data were reviewed and discussed at a multi-disciplinary meeting on #text(weight: "bold")[#patient_result_review_date]. Results are based on data meeting specific reporting thresholds and clinical relevance. 
]

#let third_page_signatures = [
  #for signature in report_authorisation_signatures {
    [ 
      #v(64pt);
      #grid(columns: (3fr, 1fr), stroke: (_, y) => if y == 0 { (bottom: 1pt) }, gutter: 16pt)[]
      #grid(columns: (3fr, 1fr), stroke: none, gutter: 16pt)[#signature.at(0)][Date]
      
      #v(-5pt)
      #signature.at(1), #signature.at(2)
      
    ]
  }
]

#let third_page_issuing_laboratory = [
  Issued by: #report_authorisation_laboratory

  Unique report identifier: #report_authorisation_identifier
]

// PAGE COMPOSITION

#let first_page = [
  #first_page_header
  #v(16pt)
  #first_page_referral
  #v(20pt)
  #outlinebox(title: "Test Result")[#first_page_results]
  
  #outlinebox(title: "Comments")[#first_page_comments]
  
  #pagebreak()
]

#let second_page = [  
  
  #outlinebox(title: "Orthogonal Tests")[#second_page_orthogonal]

  #outlinebox(title: "Clinical Notes")[#second_page_clinical]

  #outlinebox(title: "Actions")[#second_page_actions]
  
  #outlinebox(title: "Contact")[#second_page_contact]
  
  #pagebreak()
]


#let legal_page = [

  #outlinebox(title: "Disclosure Statement")[#second_page_disclosure]
  
  #outlinebox(title: "Liability Statement")[#second_page_liability]

  #pagebreak()
]

#let authorisation_page = [
  #outlinebox(title: "Issuing Laboratory")[#third_page_issuing_laboratory]

  #outlinebox(title: "Authorisation")[
    #third_page_authorisation
    #third_page_signatures
  ]

  #if (appendix_laboratory_enabled or appendix_bioinformatics_enabled) {
    pagebreak()
  }
]

// APPENDICES

#let appendix_laboratory_header = {
  show table.cell.where(x: 0): set text(weight: 500)
  show table.cell.where(x: 2): set text(weight: 500)
  table(
    columns: 4,
    [Protocol], [#appendix_laboratory_header_protocol], [Sample ID], [#appendix_laboratory_header_sample_id],
    [Version], [#appendix_laboratory_header_version], [Run ID], [#appendix_laboratory_header_run_id],
    [Extraction], [#appendix_laboratory_header_extraction], [Extraction Control], [#appendix_laboratory_header_extraction_control],
    [RNA Depletion], [#appendix_laboratory_header_rna_depletion], [Library Control], [#appendix_laboratory_header_library_control],
    [Adapters], [#appendix_laboratory_header_adapter], [Sequencing Control], [#appendix_laboratory_header_sequencing_control],
    [Library], [#appendix_laboratory_header_library], [Negative Control], [#appendix_laboratory_header_negative_control],
    [Sequencer], [#appendix_laboratory_header_sequencer], [Positive Control], [#appendix_laboratory_header_positive_control],
  )
}

#let appendix_laboratory = [
  = Appendix A: Laboratory Assay

  #appendix_laboratory_header
  #v(16pt)
  #outlinebox(title: "Description")[#appendix_laboratory_description]

  #outlinebox(title: "Comments")[#appendix_laboratory_comments]

  #if (appendix_bioinformatics_enabled) {
    pagebreak()
  }
]

#let appendix_bioinformatics_header = {
  show table.cell.where(x: 0): set text(weight: 500)
  show table.cell.where(x: 2): set text(weight: 500)
  table(
    columns: 4,
    [Pipeline], [#appendix_bioinformatics_header_pipeline], [Sample ID], [#appendix_bioinformatics_header_sample_id],
    [Version], [#appendix_bioinformatics_header_version], [Libraries], [#appendix_bioinformatics_header_libraries],
    [Pipeline ID], [#appendix_bioinformatics_header_pipeline_id], [Configuration], [#appendix_bioinformatics_header_configuration],
    [Started], [#appendix_bioinformatics_header_started], [Databases], [#appendix_bioinformatics_header_databases],
    [Completed], [#appendix_bioinformatics_header_completed], [Taxonomy], [#appendix_bioinformatics_header_taxonomy],
  )
}

#let appendix_bioinformatics = [
  = Appendix B: Bioinformatics Analysis

  #appendix_bioinformatics_header
  #v(16pt)
  #outlinebox(title: "Description")[#appendix_bioinformatics_description]

  #outlinebox(title: "Comments")[#appendix_bioinformatics_comments]

  #pagebreak()

  #outlinebox(title: "Libraries")[#appendix_bioinformatics_libraries]

  #outlinebox(title: "Evidence")[#appendix_bioinformatics_evidence]

];

// PAGES

#first_page

#second_page

#legal_page

#authorisation_page

#if appendix_laboratory_enabled {
  [#appendix_laboratory]
}
#if appendix_bioinformatics_enabled {
  [#appendix_bioinformatics]
}

