<script lang="ts">
	import { page } from "$app/stores";
	import DataViewSelection from "$lib/components/data/sample/selections/DataViewSelection.svelte";
	import WorkflowSelection from "$lib/components/data/sample/selections/WorkflowSelection.svelte";
	import SampleSelection from "$lib/components/data/sample/selections/SampleSelection.svelte";
    import QualityControl from "$lib/components/data/sample/QualityControl.svelte";
	import type { Cerebro, CerebroFilterConfig, ClientFilterConfig, ClientFilterMinimum, ClientFilterModules, HighlightConfig, QualityControlSummary, TaxonHighlightConfig, WorkflowConfig } from "$lib/utils/types";
	import CommentBox from "$lib/components/data/sample/CommentBox.svelte";
	import { getToastStore, ListBox, ListBoxItem, SlideToggle } from "@skeletonlabs/skeleton";
    import type { ToastSettings } from "@skeletonlabs/skeleton";
	import Classification from "$lib/components/data/sample/TaxonomicProfile.svelte";
    import ServerFilterConfiguration from "$lib/components/data/sample/taxa/filters/ServerFilterConfiguration.svelte";
    import ClientFilterConfiguration from "$lib/components/data/sample/taxa/filters/ClientFilterConfiguration.svelte";
	import Candidates from "$lib/components/data/sample/Candidates.svelte";
	import Reports from "$lib/components/data/sample/Reports.svelte";
	import ContamHighlights from "$lib/components/data/sample/taxa/filters/TaxonomyHighlights.svelte";
	import TaxonomyHighlights from "$lib/components/data/sample/taxa/filters/TaxonomyHighlights.svelte";
    
    import {
        selectedWorkflowConfiguration,
        selectedModels,
        selectedIdentifiers,
        selectedServerFilterConfig,
        selectedClientFilterConfig,
        selectedTaxonHighlightConfig
    } from '$lib/stores/stores';

	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import { goto } from "$app/navigation";

    let selectedView: string = "Quality Control";

    let identifiers = $page.data.sampleCerebro?.length ? [$page.data.sampleCerebro[0].id] : [];

    let selectedWorkflowIdentifier: string = $page.data.requestedWorkflow;

    $: {
        // Update the workflow configuration based on selected workflow
        selectedWorkflowConfiguration.set(
            $page.data.sampleWorkflows.find(
                (workflow: WorkflowConfig) => workflow.id === selectedWorkflowIdentifier
            )
        );
        // Filter selected models based on selected sample identifiers
        const models = [...$page.data.sampleCerebro, ...$page.data.controlCerebro].filter((cerebro) =>
            $selectedIdentifiers.includes(cerebro.id)
        );
        selectedModels.set(models);
    }

    let showServerSideFilters: boolean = false;


    const reloadTable = () => {
        // To reload classficiation table we re-assign the current 
        // selected identfiers to trigger the reactive load in child 
        // component
        selectedIdentifiers.set(
            $selectedIdentifiers
        );
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
                    <ListBoxItem bind:group={selectedView} name="medium" value="Quality Control" active='variant-soft' rounded='rounded-token' on:click={() => goto('./qc')}>Quality Control</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="Taxonomic Profile", active='variant-soft' rounded='rounded-token' on:click={() => goto('./profile')}>Taxonomic Profile</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="Pathogen Candidates", active='variant-soft' rounded='rounded-token' on:click={() => goto('./candidates')}>Pathogen Candidates</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="Team Comments", active='variant-soft' rounded='rounded-token' on:click={() => goto('./comments')}>Team Comments</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="Clinical Reports", active='variant-soft' rounded='rounded-token' on:click={() => goto('./reports')}>Clinical Reports</ListBoxItem>
                </ListBox>
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
                <SampleSelection models={$page.data.sampleCerebro} identifiers={identifiers} variant="sample" variantColor="primary"/>
            </div>
            {#if selectedView !== "Clinical Reports"}
                <p class="mb-1"><span class="opacity-40">Controls</span></p>
                <div class="mb-4 p-4 text-sm">
                    <SampleSelection models={$page.data.controlCerebro} identifiers={identifiers} variant="control" variantColor="primary"/>
                </div>
            {/if}
        </div>
        {#if selectedView === "Taxonomic Profile"}
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
                    <ServerFilterConfiguration bind:serverFilterConfig={$selectedServerFilterConfig}></ServerFilterConfiguration>
                    <div class="text-center py-8">
                        <button type="button" class="btn variant-outline-primary w-3/4" on:click={() => reloadTable()}>
                            <svg aria-hidden="true" fill="none" stroke="currentColor" class="w-5 h-5 mr-3" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                <path d="M16.023 9.348h4.992v-.001M2.985 19.644v-4.992m0 0h4.992m-4.993 0l3.181 3.183a8.25 8.25 0 0013.803-3.7M4.031 9.865a8.25 8.25 0 0113.803-3.7l3.181 3.182m0-4.991v4.99" stroke-linecap="round" stroke-linejoin="round"></path>
                            </svg>
                            Reload Table
                        </button>
                    </div>
                {:else}
                    <ClientFilterConfiguration bind:clientFilterConfig={$selectedClientFilterConfig}></ClientFilterConfiguration>
                {/if}
                
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
