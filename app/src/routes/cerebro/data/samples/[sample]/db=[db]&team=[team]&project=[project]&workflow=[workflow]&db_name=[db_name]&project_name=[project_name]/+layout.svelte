<script lang="ts">


    import {
        selectedModels,
        selectedIdentifiers,
        selectedServerFilterConfig,
        selectedClientFilterConfig,
        selectedTaxonHighlightConfig,

		selectedPrevalenceContamConfig

    } from '$lib/stores/stores';

	import { page } from "$app/stores";
	import { ListBox, ListBoxItem, SlideToggle } from "@skeletonlabs/skeleton";
	import { FileTag, type Cerebro, type WorkflowConfig } from "$lib/utils/types";

	import WorkflowSelection from "$lib/components/data/sample/selections/WorkflowSelection.svelte";
	import SampleSelection from "$lib/components/data/sample/selections/SampleSelection.svelte";
    import ServerFilterConfiguration from "$lib/components/data/sample/taxa/filters/ServerFilterConfiguration.svelte";
    import ClientFilterConfiguration from "$lib/components/data/sample/taxa/filters/ClientFilterConfiguration.svelte";
	import TaxonomyHighlights from "$lib/components/data/sample/taxa/filters/TaxonomyHighlights.svelte";
    
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import { goto } from "$app/navigation";

    enum DataView {
        QualityControl = "qc",
        TaxonomicProfile = "profile",
        Candidates = "candidates",
        Comments = "comments",
        Reports = "reports"
    }

    function getDnaLibraryIdentifiers(): string[] {
        const dnaFilter: Cerebro[] = $page.data.sampleCerebro.filter((model: Cerebro) => model.sample.tags.includes(FileTag.DNA));
        return dnaFilter.map(model => model.id)
    }
    function getDnaNegativeTemplateControlIdentifiers(): string[] {
        const dnaFilter: Cerebro[] = $page.data.controlCerebro.filter((model: Cerebro) => model.sample.tags.includes(FileTag.DNA));
        return dnaFilter.map(model => model.id)
    }

    function getLibraryModels(selectedLibraryIdentifiers: string[]): Cerebro[] {
        return $page.data.sampleCerebro.filter((model: Cerebro) => selectedLibraryIdentifiers.includes(model.id))
    }

    function getNegativeTemplatecontrolIdentifiers(selectedNegativeTemplateControlIdentifiers: string[]): Cerebro[] {
        return $page.data.controlCerebro.filter((model: Cerebro) => selectedNegativeTemplateControlIdentifiers.includes(model.id))
    }

    
    let showServerSideFilters: boolean = false;
    let selectedView: string = DataView.QualityControl;
    
    let workflowIdentifier: string = $page.data.requestedWorkflow;

    let selectedLibraryIdentifiers: string[] = getDnaLibraryIdentifiers();
    let selectedNegativeTemplateControlIdentifiers: string[] = getDnaNegativeTemplateControlIdentifiers();

    $: {
        const selectedLibraryModels = getLibraryModels(selectedLibraryIdentifiers);
        const selectedNegativeTemplateControlModels = getNegativeTemplatecontrolIdentifiers(selectedNegativeTemplateControlIdentifiers);

        selectedIdentifiers.set([...selectedLibraryIdentifiers, ...selectedNegativeTemplateControlIdentifiers]);
        selectedModels.set([...selectedLibraryModels, ...selectedNegativeTemplateControlModels]);
    }


    function reloadTable() {
        selectedIdentifiers.set($selectedIdentifiers);
    }

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
                <ListBox>
                    <ListBoxItem bind:group={selectedView} name="medium" value={DataView.QualityControl} active='variant-soft' rounded='rounded-token' on:click={() => goto(`./${DataView.QualityControl}`)}>Quality Control</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value={DataView.TaxonomicProfile}  active='variant-soft' rounded='rounded-token' on:click={() => goto(`./${DataView.TaxonomicProfile}`)}>Taxonomic Profile</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value={DataView.Candidates}  active='variant-soft' rounded='rounded-token' on:click={() => goto(`./${DataView.Candidates}`)}>Pathogen Candidates</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value={DataView.Comments}  active='variant-soft' rounded='rounded-token' on:click={() => goto(`./${DataView.Comments}`)}>Team Comments</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value={DataView.Reports}  active='variant-soft' rounded='rounded-token' on:click={() => goto(`./${DataView.Reports}`)}>Clinical Reports</ListBoxItem>
                </ListBox>
            </div>
        </div>
        <p class="mb-1"><span class="opacity-60">Data selection</span></p>
        <div class="w-full border border-primary-500 rounded-md p-4">
            <p class="mb-1"><span class="opacity-40">Workflows</span></p>
            <div class="mb-4 p-4 ">
                <WorkflowSelection workflows={$page.data.sampleWorkflows} bind:selectedWorkflowIdentifier={workflowIdentifier} />
            </div>
            <p class="mb-1"><span class="opacity-40">Libraries</span></p>
            <div class="mb-4 p-4 text-sm">
                <SampleSelection models={$page.data.sampleCerebro} bind:identifiers={selectedLibraryIdentifiers} variant="sample" variantColor="primary"/>
            </div>
            <p class="mb-1"><span class="opacity-40">Negative template controls</span></p>
            <div class="mb-4 p-4 text-sm">
                <SampleSelection models={$page.data.controlCerebro} bind:identifiers={selectedNegativeTemplateControlIdentifiers}  variant="control" variantColor="secondary"/>
            </div>
        </div>
        {#if selectedView === DataView.TaxonomicProfile}
            <p class="mb-1 mt-4">
                <span class="opacity-60">Client-side filters</span>
            </p>
            <div class="w-full border border-primary-500 rounded-md p-4">
                <div class="grid grid-cols-3 sm:grid-cols-3 md:grid-cols-3 gap-x-4 align-center">
                    <div class="col-span-2">
                        <p class="opacity-60 text-xs pb-8 pr-4">
                            Apply fine-grained evidence filters and negative template control comparison thresholds using server-side filters below
                        </p>
                    </div>
                </div>
                <ClientFilterConfiguration bind:clientFilterConfig={$selectedClientFilterConfig}></ClientFilterConfiguration>
            </div>
            <p class="mb-1 mt-4">
                <span class="opacity-60">Server-side filters</span>
            </p>
            <div class="w-full border border-primary-500 rounded-md p-4" >
                <ServerFilterConfiguration bind:prevalenceContamConfig={$selectedPrevalenceContamConfig} bind:serverFilterConfig={$selectedServerFilterConfig}></ServerFilterConfiguration>
                <div class="text-center py-8">
                    <button type="button" class="btn variant-outline-primary w-3/4" on:click={() => reloadTable()}>
                        <svg aria-hidden="true" fill="none" stroke="currentColor" class="w-5 h-5 mr-3" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M16.023 9.348h4.992v-.001M2.985 19.644v-4.992m0 0h4.992m-4.993 0l3.181 3.183a8.25 8.25 0 0013.803-3.7M4.031 9.865a8.25 8.25 0 0113.803-3.7l3.181 3.182m0-4.991v4.99" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg>
                        Reload Table
                    </button>
                </div>
            </div>
            <p class="mb-1 mt-4">
                <span class="opacity-60">Taxonomy highlights</span>
            </p>
            <div class="w-full border border-primary-500 rounded-md p-4" >
                <TaxonomyHighlights title="Common contamination" bind:highlightConfig={$selectedTaxonHighlightConfig.contamination}></TaxonomyHighlights>
                <TaxonomyHighlights title="Syndromic pathogens" bind:highlightConfig={$selectedTaxonHighlightConfig.syndrome}></TaxonomyHighlights>
            </div>
        {/if}
    </div>

    <div class="col-span-2">
        {#if $selectedIdentifiers.length}
            <slot />
        {:else}
            <div class="flex justify-center py-16">
                <ErrorAnimation></ErrorAnimation>
            </div>
            <p class="flex justify-center text-lg pb-4">No libraries or controls selected</p>
        {/if}
    </div>
        
</div>
