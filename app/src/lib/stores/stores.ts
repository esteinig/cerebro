import { browser } from "$app/environment";
import { CNS_BACTERIA, getCurrentDate, PATHOGENS, POSITIVE_CONTROLS, VERTEBRATE_VIRUSES } from "$lib/utils/helpers";

import { 
    ProfileTool, 
    PathogenDetectionRank, 
} from "$lib/utils/types";

import type {
    TaxonOverviewRecord, 
    Cerebro, 
    TaxonFilterConfig, 
    ClientFilterConfig, 
    ClientFilterMinimum, 
    ClientFilterModules, 
    HighlightConfig, 
    QualityControlSummary, 
    TaxonHighlightConfig, 
    WorkflowConfig, 
    ClientFilterContam, 
    PrevalenceContaminationConfig, 
    PathogenDetectionReport,
    ReportHeader,
    PatientHeader,
    PatientResult,
    ReportFooter,
    ReportAuthorisation,
    ReportLegal,
    AppendixBioinformatics,
    AppendixLaboratory,
    AppendixBioinformaticsHeader,
    AppendixLaboratoryHeader,
    ClientFilterEvidence, 
} from "$lib/utils/types";

import { writable, type Writable } from "svelte/store";

// Session based theme store. Grabs the current theme from the current body.
export const storeTheme = writable(browser ? document.body.getAttribute('data-theme') ?? '' : 'dali');


interface SessionSettings {
    acceptedCommentRisk: boolean
}

// Session based settings store. Keeps user settings like risk acceptance checks.
export const storeSettings: Writable<SessionSettings> = writable({
    acceptedCommentRisk: false
} satisfies SessionSettings);


// Stores for '/data/samples/[sample]'

export const selectedWorkflowConfiguration = writable<WorkflowConfig | null>(null);
export const selectedModels = writable<Cerebro[]>([]);
export const selectedIdentifiers = writable<string[]>([]);
export const selectedWorkflowIdentifier = writable<string>("0"); // 0 = latest in loading function

export const navigationLoading = writable<boolean | null>(null);


const defaultClientFilterConfig: ClientFilterConfig = {
    domains: [],
    genera: [],
    species: [],
    tools: [ProfileTool.Vircov, ProfileTool.Bracken, ProfileTool.Metabuli, ProfileTool.Ganon2, ProfileTool.Blast],
    contam: {
        display: true,
        opacity: 20
    } satisfies ClientFilterContam,
    evidence: {
        display: false,
        opacity: 20
    } satisfies ClientFilterEvidence,
    modules: {
        alignment: false,
        profile: false,
        assembly: false
    } satisfies ClientFilterModules,
    minimum: {
        rpm: 0,
        rpm_kmer: 0,
        rpm_alignment: 0,
        contigs: 0,
        bases: 0
    } satisfies ClientFilterMinimum
}

const defaultServerFilterConfig: TaxonFilterConfig = {
    rank: PathogenDetectionRank.Species,
    domains: [],            // Filter by domain names
    tools: [],              // Filter by specific detection tools
    modes: [],              // Filter by detection modes (Sequence/Profile)
    min_bases: 0,           // Minimum contig length for inclusion
    max_bases: null,
    min_bpm: 0.0,           // Minimum contig length per million for inclusion
    min_reads: 0,           // Minimum read count for inclusion
    min_rpm: 0.0,           // Minimum RPM for inclusion
    max_rpm: null,          // Minimum RPM for inclusion
    min_abundance: 0.0,     // Minimum abundance for inclusion
    ntc_ratio: 1,
    lineage: null,
    targets: null,
    collapse_variants: false
}

export const primaryThresholdServerFilterConfig: TaxonFilterConfig = {
    rank: PathogenDetectionRank.Species,
    domains: [],
    tools: ["Vircov", "Bracken", "Metabuli", "Ganon2", "Blast"], 
    modes: ["Mixed"], 
    min_bases: 0,   // Check this threshold TODO V75
    max_bases: null,
    min_bpm: 0.0,      
    min_reads: 0,       
    min_rpm: 0.0,      
    max_rpm: null,
    min_abundance: 0,   
    ntc_ratio: 1,
    lineage: [
      {
        lineages: ["d__Bacteria", "d__Archaea"],
        tags: ["DNA"],
        min_alignment_tools: null,
        min_alignment_rpm: null,
        min_alignment_regions: null,
        min_alignment_regions_coverage: null,
        min_kmer_tools: 3,
        min_kmer_rpm: 10.0,
        min_assembly_tools: null,
      },
      {
        lineages: ["d__Viruses"],
        tags: ["DNA", "RNA"],
        min_alignment_tools: 1,
        min_alignment_rpm: 10.0,
        min_alignment_regions: null,
        min_alignment_regions_coverage: null,
        min_kmer_tools: 3,
        min_kmer_rpm: 10.0,
        min_assembly_tools: null,
      },
      {
        lineages: ["d__Eukaryota"],
        tags: ["DNA"],
        min_alignment_tools: null,
        min_alignment_rpm: null,
        min_alignment_regions: null,
        min_alignment_regions_coverage: null,
        min_kmer_tools: 3,
        min_kmer_rpm: 10.0,
        min_assembly_tools: 1,
      },
    ],
    targets: null,
    collapse_variants: false
  };
  
export const secondaryThresholdServerFilterConfig: TaxonFilterConfig = {
    rank: PathogenDetectionRank.Species,
    domains: [],
    tools: ["Vircov", "Bracken", "Metabuli", "Ganon2", "Blast"], 
    modes: ["Mixed"], 
    min_bases: 0, 
    max_bases: null,
    min_bpm: 0.0, 
    min_reads: 0, 
    min_rpm: 0.0, 
    max_rpm: null,
    min_abundance: 0.0,
    ntc_ratio: 1,
    lineage: [
        {
            lineages: ["d__Bacteria", "d__Archaea"],
            tags: ["DNA"],
            min_alignment_tools: null,
            min_alignment_rpm: null,
            min_alignment_regions: null,
            min_alignment_regions_coverage: null,
            min_kmer_tools: 2,
            min_kmer_rpm: 3.0,
            min_assembly_tools: null,
        },
        {
            lineages: ["d__Viruses"],
            tags: ["DNA", "RNA"],
            min_alignment_tools: 1,
            min_alignment_rpm: 0.0,
            min_alignment_regions: null,
            min_alignment_regions_coverage: null,
            min_kmer_tools: 2,
            min_kmer_rpm: 3.0,
            min_assembly_tools: null,
        },
        {
            lineages: ["d__Eukaryota"],
            tags: ["DNA"],
            min_alignment_tools: null,
            min_alignment_rpm: null,
            min_alignment_regions: null,
            min_alignment_regions_coverage: null,
            min_kmer_tools: 3,
            min_kmer_rpm: 5.0,
            min_assembly_tools: 1,
        },
    ],
    targets: null,
    collapse_variants: false
};

export const validationThresholdServerFilterConfig: TaxonFilterConfig = {
    rank: PathogenDetectionRank.Species,
    domains: [],
    tools: ["Vircov", "Bracken", "Metabuli", "Ganon2", "Blast"], 
    modes: ["Mixed"], 
    min_bases: 0, 
    max_bases: null,
    min_bpm: 0.0, 
    min_reads: 0, 
    min_rpm: 0.0, 
    max_rpm: null,
    min_abundance: 0.0,
    ntc_ratio: 1,
    lineage: [
        {
            lineages: ["d__Viruses"],
            tags: ["DNA", "RNA"],
            min_alignment_tools: 1,
            min_alignment_rpm: 0.0,
            min_alignment_regions: 2,
            min_alignment_regions_coverage: 20.0,
            min_kmer_tools: 1,
            min_kmer_rpm: 0.0,
            min_assembly_tools: null,
        },
    ],
    targets: VERTEBRATE_VIRUSES,
    collapse_variants: false
};

const defaultPrevalenceContamConfig: PrevalenceContaminationConfig = {
    min_rpm: 0.0,
    threshold: 0.6,
    sample_type: null,
    collapse_variants: false
}

// Contamination highlights

const contamHighlightConfig: HighlightConfig = {
    species: [],
    taxid: [],
    color: "secondary"
}

const syndromeHighlightConfig: HighlightConfig = {
    species: PATHOGENS,
    taxid: [],
    color: "tertiary"
}

const validationHighlightConfig: HighlightConfig = {
    species: CNS_BACTERIA,
    taxid: [],
    color: "secondary"
}

const positiveControlsHighlightConfig: HighlightConfig = {
    species: POSITIVE_CONTROLS,
    taxid: [],
    color: "secondary"
}

const defaultTaxonHighlightConfig: TaxonHighlightConfig = {
    contamination: contamHighlightConfig,
    syndrome: syndromeHighlightConfig,
    validation: validationHighlightConfig,
    positive_controls: positiveControlsHighlightConfig
}

const defaultReport: PathogenDetectionReport = {
    header: {
        logo_enabled: true,
        logo: null
    } satisfies ReportHeader,
    footer: {
        reporting_location: "Melbourne"
    } satisfies ReportFooter,
    authorisation: {
        laboratory: "Victorian Infectious Diseases Reference Laboratory (VIDRL), 792 Elizabeth St, Melbourne VIC 3004",
        identifier: crypto.randomUUID(),
        signatures: [
            {"name": "Dr. Chuan Lim", "position": "Acting Director", "institution": "VIDRL"},
            {"name": "Dr. Prashanth Ramachandran", "position": "Clinical Metagenomics Research Fellow", "institution": "University of Melbourne"}
        ]
    } satisfies ReportAuthorisation,
    legal: {
        disclosure: "This report overrides all previous reports. The results relate solely to the sample/s as received and are limited to the specific tests undertaken as listed on the report. The results of this report are confidential and are not to be used or disclosed to any other person or used for any other purpose, whether directly or indirectly, unless that use is disclosedor the purpose is expressly authorised in writing by the Victorian Infectious Diseases Reference Laboratory (VIDRL) and the named recipient on this report.",
        liability: "To the fullest extent permitted by law, VIDRL will not be liable for any loss or claim (including legal costs calculated on an indemnity basis) which arise because of (a) problems related to the merchantability, fitness or quality of the sample/s, or (b) any negligent or unlawful act or omissions by VIDRL that is connected with any activities or services provided by VIDRL under this agreement (including the timing and/or method under which the sample/s were taken, stored or transported).",
        disclaimer: "This metagenomics assay is a research‐only assay, under ongoing development, and should only be used in addition to usual diagnostic testing algorithms. The test has been performed using an unapproved IVD medical device that has not been evaluated for quality, safety or performance by the Therapeutic Goods Administration. Results are not sufficient to rule in or out causative pathogens, and should not be used as the basis for clinical management decisions. Methods used for this analysis have not yet been validated according to the current NPAAC standards and results should be interpreted accordingly. NATA/RCPA accreditation does not cover the performance of this service."
    } satisfies ReportLegal,
    patient_header: {
        patient_name: "",
        patient_urn: "",
        patient_dob: "",
        requested_doctor: "",
        hospital_site: "",
        laboratory_number: "",
        specimen_id: "",
        date_collected: "",
        date_received: "",
        specimen_type: "",
        reporting_laboratory: "VIDRL",
        reporting_date: getCurrentDate(),
    } satisfies PatientHeader,
    patient_result: {
        pathogen_detected: false,
        pathogen_reported: "No pathogen detected",
        review_date: getCurrentDate(),
        comments: "No further comments.",
        actions: "No further actions taken.",
        orthogonal_tests: "No orthogonal test results provided in this report.",
        clinical_notes: "No clinical notes provided in this report.",
        contact_name: "Dr. Chuan Lim",
        contact_email: "chuan.lim@mh.org.au"
    } satisfies PatientResult,
    appendix_laboratory: {
        enabled: false,
        description: "",
        comments: "",
        header: {
            protocol: "",
            sample_id: "",
            version: "",
            run_id: "",
            extraction: "",
            extraction_control: "",
            rna_depletion: "",
            library_control: "",
            adapter: "",
            sequencing_control: "",
            library: "",
            negative_control: "",
            sequencer: "",
            positive_control: ""
        } satisfies AppendixLaboratoryHeader
    } satisfies AppendixLaboratory,
    appendix_bioinformatics: {
        enabled: false,
        description: "",
        comments: "",
        header: {
            pipeline: "",
            version: "",
            pipeline_id: "",
            configuration: "",
            started: "",
            completed: "",
            sample_id: "",
            libraries: "",
            databases: "",
            taxonomy: ""
        } satisfies AppendixBioinformaticsHeader,
        libraries: "",
        evidence: ""
    } satisfies AppendixBioinformatics
}

export const selectedReportSchema = writable<PathogenDetectionReport>(defaultReport)
export const selectedClientFilterConfig = writable<ClientFilterConfig>(defaultClientFilterConfig)
export const selectedServerFilterConfig = writable<TaxonFilterConfig>(defaultServerFilterConfig)
export const selectedPrevalenceContamConfig = writable<PrevalenceContaminationConfig>(defaultPrevalenceContamConfig)
export const selectedTaxonHighlightConfig = writable<TaxonHighlightConfig>(defaultTaxonHighlightConfig)
export const selectedTaxa = writable<TaxonOverviewRecord[]>([])

