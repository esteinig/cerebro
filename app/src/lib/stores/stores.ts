import { browser } from "$app/environment";
import { PathogenDetectionRank, type TaxonOverviewRecord, type Cerebro, type CerebroFilterConfig, type ClientFilterConfig, type ClientFilterMinimum, type ClientFilterModules, type HighlightConfig, type QualityControlSummary, type TaxonHighlightConfig, type WorkflowConfig } from "$lib/utils/types";
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


const defaultClientFilterConfig: ClientFilterConfig = {
    domains: [],
    genera: [],
    species: [],
    modules: {
        alignment: false,
        kmer: false,
        assembly: false
    } satisfies ClientFilterModules,
    minimum: {
        rpm: 10,
        rpm_kmer: 0,
        rpm_alignment: 0,
        contigs: 0,
        bases: 0
    } satisfies ClientFilterMinimum
}

const defaultServerFilterConfig: CerebroFilterConfig = {
    rank: PathogenDetectionRank.Species,
    domains: [],
    tags: [],
    // K-mer data
    kmer_min_reads: 3,
    kmer_databases: [],
    // Scanning part of alignment pipelines
    alignment_min_reads: 3,
    alignment_min_bases: 0,
    alignment_min_regions: 0,
    alignment_min_coverage: 0,
    // General alignment section of the pipeline
    alignment_min_ref_length: 2000,
    // LCA BLAST/Diamond on assembled contigs
    assembly_min_contig_length: 200,
    assembly_min_contig_identity: 60.0,
    assembly_min_contig_coverage: 60.0
}


// Contamination highlights

const contamHighlightConfig: HighlightConfig = {
    species: [],
    taxid: [],
    color: "secondary"
}

const syndromeHighlightConfig: HighlightConfig = {
    species: [
        "Cryptococcus",
        "Neisseria",
        "Human betaherpesvirus",
        "Streptococcus",
        "Haemophilus",
        "Lymphocryptovirus"

    ],
    taxid: [],
    color: "tertiary"
}

const defaultTaxonHighlightConfig: TaxonHighlightConfig = {
    contamination: contamHighlightConfig,
    syndrome: syndromeHighlightConfig
}

export const selectedClientFilterConfig = writable<ClientFilterConfig>(defaultClientFilterConfig)
export const selectedServerFilterConfig = writable<CerebroFilterConfig>(defaultServerFilterConfig)
export const selectedTaxonHighlightConfig = writable<TaxonHighlightConfig>(defaultTaxonHighlightConfig)
export const selectedTaxa = writable<TaxonOverviewRecord[]>([])

