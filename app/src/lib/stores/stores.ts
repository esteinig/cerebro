import { browser } from "$app/environment";
import { PATHOGENS } from "$lib/utils/helpers";
import { PathogenDetectionRank, type TaxonOverviewRecord, type Cerebro, type TaxonFilterConfig, type ClientFilterConfig, type ClientFilterMinimum, type ClientFilterModules, type HighlightConfig, type QualityControlSummary, type TaxonHighlightConfig, type WorkflowConfig, type ClientFilterContam, type PrevalenceContaminationConfig, ProfileTool } from "$lib/utils/types";
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
    min_bases: 200,         // Minimum contig length for inclusion
    min_bpm: 0.0,           // Minimum contig length per million for inclusion
    min_reads: 0,           // Minimum read count for inclusion
    min_rpm: 1.0,           // Minimum RPM for inclusion
    min_abundance: 0,       // Minimum abundance for inclusion
    ntc_ratio: 3
}

const defaultPrevalenceContamConfig: PrevalenceContaminationConfig = {
    min_rpm: 0.0,
    threshold: 0.5
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

const defaultTaxonHighlightConfig: TaxonHighlightConfig = {
    contamination: contamHighlightConfig,
    syndrome: syndromeHighlightConfig
}

export const selectedClientFilterConfig = writable<ClientFilterConfig>(defaultClientFilterConfig)
export const selectedServerFilterConfig = writable<TaxonFilterConfig>(defaultServerFilterConfig)
export const selectedPrevalenceContamConfig = writable<PrevalenceContaminationConfig>(defaultPrevalenceContamConfig)
export const selectedTaxonHighlightConfig = writable<TaxonHighlightConfig>(defaultTaxonHighlightConfig)
export const selectedTaxa = writable<TaxonOverviewRecord[]>([])

