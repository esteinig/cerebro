<script lang="ts">
	import { page } from "$app/stores";
	import DataViewSelection from "$lib/components/data/sample/selections/DataViewSelection.svelte";
	import WorkflowSelection from "$lib/components/data/sample/selections/WorkflowSelection.svelte";
	import SampleSelection from "$lib/components/data/sample/selections/SampleSelection.svelte";
    import QualityControl from "$lib/components/data/sample/QualityControl.svelte";
	import type { Cerebro, CerebroFilterConfig, ClientFilterConfig, ClientFilterMinimum, ClientFilterModules, HighlightConfig, QualityControlSummary, TaxonHighlightConfig, WorkflowConfig } from "$lib/utils/types";
	import CommentBox from "$lib/components/data/sample/CommentBox.svelte";
	import { getToastStore, SlideToggle } from "@skeletonlabs/skeleton";
    import type { ToastSettings } from "@skeletonlabs/skeleton";
	import Classification from "$lib/components/data/sample/Classification.svelte";
    import ServerFilterConfiguration from "$lib/components/data/sample/taxa/filters/ServerFilterConfiguration.svelte";
    import ClientFilterConfiguration from "$lib/components/data/sample/taxa/filters/ClientFilterConfiguration.svelte";
	import Candidates from "$lib/components/data/sample/Candidates.svelte";
	import Reports from "$lib/components/data/sample/Reports.svelte";
	import ContamHighlights from "$lib/components/data/sample/taxa/filters/TaxonomyHighlights.svelte";
	import TaxonomyHighlights from "$lib/components/data/sample/taxa/filters/TaxonomyHighlights.svelte";

    let toastStore = getToastStore();

    let selectedView: string = "qc";

    let selectedWorkflowConfiguration: WorkflowConfig;
    let selectedWorkflowIdentifier: string = $page.data.requestedWorkflow;

    let selectedIdentifiers: string[] = $page.data.sampleCerebro?.length ? [$page.data.sampleCerebro[0].id] : [];
    let deselectedIdentifiers: string[] = [];

    let selectedModels: Cerebro[] = $page.data.sampleCerebro?.length ? [$page.data.sampleCerebro[0]] : [];
    let selectedQualityControlSummaries: QualityControlSummary[] = $page.data.qualityControlSummaries?.length ? $page.data.qualityControlSummaries[0] : [];

    
    const changeWorkflow = async() => {

    }

    const reloadTable = () => {
        // To reload classficiation table we re-assign the current 
        // selected identfiers to trigger the reactive load in child 
        // component
        selectedIdentifiers = selectedIdentifiers;
    }


    // Change selected model on change of library/control selection
    $: {
        selectedWorkflowConfiguration = $page.data.sampleWorkflows.filter(
            (workflow: WorkflowConfig) => selectedWorkflowIdentifier === workflow.id
        )[0]
        selectedModels = [...$page.data.sampleCerebro, ...$page.data.controlCerebro].filter(
            cerebro => selectedIdentifiers.includes(cerebro.id)
        ).filter((value, index, self) => 
            self.findIndex(v => v.id === value.id) === index
        ); // filter unique in case we select controls as main samples
        
        let summaries: QualityControlSummary[] = $page.data.qualityControlSummaries.filter(
            (summary: QualityControlSummary) => summary.model_id ? selectedIdentifiers.includes(summary.model_id) : false
        );

        let sortOrder = selectedModels.map(cerebro => cerebro.id);
        summaries.sort(function (a, b) {
            if (a.model_id === null || b.model_id === null){
                toastStore.trigger({
                    message: "Failed to sort quality control data for some reason. Please report this bug to the issues section of Cerebro",
                    background: "variant-filled-tertiary"
                } satisfies ToastSettings)
            }
            return sortOrder.indexOf(a.model_id ?? "") - sortOrder.indexOf(b.model_id ?? "");
        })
        selectedQualityControlSummaries = summaries;

    }

    $: {
        // Deselect workflow libraries and controls to make clear
        // that comments refer to the overall sample rther than a
        // specific variant
        if (selectedView === "comments") {
            selectedWorkflowIdentifier = "";
            deselectedIdentifiers = selectedIdentifiers;
            selectedIdentifiers = [];
        } else {
            // Restore the previous selections when changing view again
            selectedWorkflowIdentifier = $page.data.requestedWorkflow;
            // Guard on page load to initiate selection
            if (!selectedIdentifiers.length){
                selectedIdentifiers = deselectedIdentifiers;
            }
        }
    }

    let clientFilterConfig: ClientFilterConfig = {
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

    // Default filter configuration in settings
    let serverFilterConfig: CerebroFilterConfig = {
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

    let contamHighlightConfig: HighlightConfig = {
        species: [
            "Herbaspirillum huttiense",
            "Stutzerimonas stutzeri",
            "Cutibacterium acnes",
            "Delftia acidovorans",
            "Micrococcus luteus",
            "Ralstonia insidiosa",
            "Tepidimonas taiwanensis",
            "Eimeria maxima",
            "Babesia bigemina",
            "Nannochloropsis gaditana",
            "Pseudomonas brenneri",
            "Methylorubrum populi",
            "Pseudomonas alcaligenes",
            "Corynebacterium kefirresidentii",
            "Malassezia restricta",
            "Cutibacterium",
            "Corynebacterium",
            "Geobacillus",
            "Cloacibacterium",
            "Bifidobacterium",
            "Delftia",
            "Finegoldia",
            "Cupriavidus"


        ],
        taxid: [],
        color: "secondary"
    }

    let syndromeHighlightConfig: HighlightConfig = {
        species: [
            "Cryptococcus",
            "Neisseria",
            "Human betaherpesvirus",
            "Streptococcus",
            "Streptococcus",
            "Haemophilus",

        ],
        taxid: [],
        color: "tertiary"
    }

    let taxonHighlightConfig: TaxonHighlightConfig = {
        contamination: contamHighlightConfig,
        syndrome: syndromeHighlightConfig
    }

   $: serverFiltersActive = serverFilterConfig.domains.length > 0 ||
        serverFilterConfig.tags.length > 0 ||
        serverFilterConfig.kmer_min_reads > 0 ||
        serverFilterConfig.kmer_databases.length > 0 ||
        serverFilterConfig.alignment_min_reads > 0 ||
        serverFilterConfig.alignment_min_bases > 0 ||
        serverFilterConfig.alignment_min_regions > 0 ||
        serverFilterConfig.alignment_min_coverage > 0 ||
        serverFilterConfig.alignment_min_ref_length > 0 ||
        serverFilterConfig.assembly_min_contig_length > 0 ||
        serverFilterConfig.assembly_min_contig_identity > 0 ||
        serverFilterConfig.assembly_min_contig_coverage > 0

    $: clientFiltersActive = clientFilterConfig.domains.length > 0 ||
        clientFilterConfig.genera.length > 0 ||
        clientFilterConfig.species.length > 0 ||
        clientFilterConfig.modules.alignment ||
        clientFilterConfig.modules.kmer ||
        clientFilterConfig.modules.assembly ||
        clientFilterConfig.minimum.rpm > 0 ||
        clientFilterConfig.minimum.rpm_kmer > 0 ||
        clientFilterConfig.minimum.rpm_alignment > 0 ||
        clientFilterConfig.minimum.contigs > 0 ||
        clientFilterConfig.minimum.bases > 0


    let showServerSideFilters: boolean = false;

</script>

<div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-3 gap-16 pt-10">
    <div>
        <div class="w-full">
            <div class="mb-1">
                <ol class="breadcrumb justify-start">
                    <li class="crumb opacity-60">{$page.params.db_name}</li>
                    <li class="crumb-separator" aria-hidden>&rsaquo;</li>
                    <li class="crumb opacity-60">{$page.params.project_name}</li>
                    <li class="crumb-separator" aria-hidden>&rsaquo;</li>
                    <li class="crumb">{$page.params.sample}</li>
                </ol>
            </div>
            <div class="mb-4 p-4 border border-primary-500 rounded-md">
                <DataViewSelection bind:selectedView={selectedView}/>
            </div>
        </div>
        <p class="mb-1"><span class="opacity-60">Data selection</span></p>
        <div class="w-full border border-primary-500 rounded-md p-4">
            <p class="mb-1"><span class="opacity-40">Workflows</span></p>
            <div class="mb-4 p-4 ">
                <WorkflowSelection workflows={$page.data.sampleWorkflows} bind:selectedWorkflowIdentifier={selectedWorkflowIdentifier} />
            </div>
            <p class="mb-1"><span class="opacity-40">Libraries</span></p>
            <div class="mb-4 p-4 text-sm">
                <SampleSelection models={$page.data.sampleCerebro} bind:selectedIdentifiers={selectedIdentifiers}  variant="sample" variantColor="primary"/>
            </div>
            {#if selectedView !== "reports"}
            <p class="mb-1"><span class="opacity-40">Controls</span></p>
            <div class="mb-4 p-4 text-sm">
                <SampleSelection models={$page.data.controlCerebro} bind:selectedIdentifiers={selectedIdentifiers} variant="control" variantColor="primary"/>
            </div>
            {/if}
        </div>
        {#if selectedView === "classification"}
            <p class="mb-1 mt-4">
                <span class="opacity-60">Taxonomy filters</span>
            </p>
            <div class="w-full border border-primary-500 rounded-md p-4" >
                
                <div class="grid grid-cols-3 sm:grid-cols-3 md:grid-cols-3 gap-x-4 align-center">
                    <div class="col-span-2">
                        <p class="opacity-60 text-xs pb-8 pr-4">
                            {#if showServerSideFilters}
                                Apply classification table filters by switching to client-side filters
                            {:else}
                                Apply fine-grained evidence filters by switching to server-side filters
                            {/if}
                            
                        </p>
                    </div>
                    <div class="col-span-1 text-xs">
                        <SlideToggle name="server-side-filters" bind:checked={showServerSideFilters} active="variant-filled-tertiary dark:variant-filled-tertiary" size="sm">Server-side</SlideToggle>
                    </div>
                </div>
                {#if showServerSideFilters}
                    <ServerFilterConfiguration bind:serverFilterConfig={serverFilterConfig}></ServerFilterConfiguration>
                    <div class="text-center py-8">
                        <button type="button" class="btn variant-outline-primary w-3/4" on:click={() => reloadTable()}>
                            <svg aria-hidden="true" fill="none" stroke="currentColor" class="w-5 h-5 mr-3" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                <path d="M16.023 9.348h4.992v-.001M2.985 19.644v-4.992m0 0h4.992m-4.993 0l3.181 3.183a8.25 8.25 0 0013.803-3.7M4.031 9.865a8.25 8.25 0 0113.803-3.7l3.181 3.182m0-4.991v4.99" stroke-linecap="round" stroke-linejoin="round"></path>
                            </svg>
                            Reload Table
                        </button>
                    </div>
                {:else}
                    <ClientFilterConfiguration bind:clientFilterConfig={clientFilterConfig}></ClientFilterConfiguration>
                {/if}
                
            </div>

            <p class="mb-1 mt-4">
                <span class="opacity-60">Taxonomy highlights</span>
            </p>
            <div class="w-full border border-primary-500 rounded-md p-4" >
                
                <TaxonomyHighlights title="Common contamination" bind:highlightConfig={contamHighlightConfig}></TaxonomyHighlights>
                <TaxonomyHighlights title="Syndromic pathogens" bind:highlightConfig={syndromeHighlightConfig}></TaxonomyHighlights>

            </div>

        {/if}
    </div>
    {#if selectedView === "qc"}
        <div class="col-span-2">
            <div class="w-full mb-5">
                <p class="mb-1"><span class="opacity-60">Quality Control</span>
                <span class="opacity-60 text-xs ml-5"> Select rows in the table to view detailed results and configurations</span>
                </p>
                <div class="mb-4 border border-primary-500 p-8 rounded-md">
                    <QualityControl selectedWorkflowConfiguration={selectedWorkflowConfiguration} selectedModels={selectedModels} selectedQualityControlSummaries={selectedQualityControlSummaries} />
                </div>
            </div>
        </div>
    {:else if selectedView === "classification"}
        <div class="col-span-2">
            <div>
                <div class="w-full mb-5">
                    <div class="">
                        <span class="opacity-60">Classification overview</span>
                        {#if clientFiltersActive}
                            <span class="text-warning-500 opacity-60 text-xs ml-5">
                                Client-side table filters are active
                            </span>
                         {/if}
                        {#if serverFiltersActive}
                            <span class="text-warning-500 opacity-60 text-xs ml-5">
                                Server-side evidence filters are active
                            </span>
                        {/if}
                    </div>
                    <div class="mb-4 border border-primary-500 p-4 rounded-md">
                        <Classification selectedIdentifiers={selectedIdentifiers} selectedModels={selectedModels} clientFilterConfig={clientFilterConfig} serverFilterConfig={serverFilterConfig} taxonHighlightConfig={taxonHighlightConfig}></Classification>
                    </div>
                </div>
            </div>
        </div>
    {:else if selectedView === "candidates"}
        <div class="col-span-2">
            <div class="w-full mb-5">
                <p class="mb-1"><span class="opacity-60">Candidates</span></p>
                <div class="mb-4 border border-primary-500 p-4 rounded-md">
                    <Candidates selectedModels={selectedModels}></Candidates>
                </div>
            </div>
        </div>
    {:else if selectedView === "comments"}
        <div class="col-span-2">
            <div class="w-full mb-5">
                <p class="mb-1">
                    <span class="opacity-60">Sample comments</span>
                    <span class="opacity-60 text-xs ml-5">Comments are submitted for all libraries of <span class="font-semibold">{$page.params.sample}</span>. Other team members are able to view and submit comments.
                </p>
                <div class="mb-4 border border-primary-500 p-4 rounded-md">
                    <CommentBox selectedDatabaseId={$page.params.db} selectedProjectId={$page.params.project}></CommentBox>
                </div>
            </div>
        </div>
    {:else if selectedView === "reports"}
        <div class="col-span-2">
            <div class="w-full mb-5">
                <p class="mb-1"><span class="opacity-60">Sample reports</span>
                    <span class="opacity-60 text-xs ml-5"></p>
                <div class="mb-4 border border-primary-500 p-4 rounded-md">
                    <Reports selectedWorkflowConfiguration={selectedWorkflowConfiguration} selectedIdentifiers={selectedIdentifiers} selectedModels={selectedModels}></Reports>
                </div>
            </div>
        </div>
    {/if}
</div>
