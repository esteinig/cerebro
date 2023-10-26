<script lang="ts">
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
    import type { Cerebro, CerebroFilterConfig, PriorityTaxon, TaxonOverview } from "$lib/utils/types";
	import { SlideToggle } from "@skeletonlabs/skeleton";
	import SpeciesOverviewTable from "./taxa/SpeciesOverviewTable.svelte";
	import CandidateTaxaTable from "./taxa/candidates/CandidateTaxaTable.svelte";
	import CandidateTaxonCard from "./taxa/candidates/CandidateTaxonCard.svelte";
	import CandidateTaxonComments from "./taxa/candidates/CandidateTaxonComments.svelte";

    export let selectedModels: Cerebro[];

    // Unique priority taxa for this sample
    let priorityTaxa: PriorityTaxon[] = []; 

    let selectedIdentifiers: string[];
    let taxaOverview: TaxonOverview[];
    let modelNameTags: Map<string, string[]>;
    let serverFilterConfigs: CerebroFilterConfig[];

    let showEvidenceTable: boolean = false;

    let selectedPriorityTaxon: PriorityTaxon | null = null;

    $: {
        priorityTaxa = selectedModels.flatMap(cerebro => cerebro.sample.priority).filter((value, index, self) => 
            self.findIndex(v => v.id === value.id) === index
        )

        selectedIdentifiers = selectedModels.map(cerebro => cerebro.id);
        taxaOverview = priorityTaxa.map(priorityTaxon => priorityTaxon.taxon_overview);
        modelNameTags = new Map(selectedModels.map(cerebro => [cerebro.name, cerebro.sample.tags]));
        serverFilterConfigs = priorityTaxa.map(priorityTaxon => priorityTaxon.filter_config);

    }


</script>

<div class="p-4">
    {#if !priorityTaxa.length}

        <div class="flex justify-center py-16 ">
            <ErrorAnimation></ErrorAnimation>
        </div>
        <p class="flex justify-center text-lg pb-4">No candidate taxa available</p>
        
    {:else}
        {#if showEvidenceTable}
            <div class="pb-10">
                <p class="pb-3 opacity-60 text-xl">Evidence</p>
                <SpeciesOverviewTable selectedIdentifiers={selectedIdentifiers} taxonOverview={taxaOverview} modelNameTags={modelNameTags} clientFilterConfig={null} serverFilterConfig={serverFilterConfigs} candidateButton={false} pagination={false}></SpeciesOverviewTable>
            </div>
        {/if}
        <div>
            <div class="flex items-center align-center justify-between pb-2">
                <p class="opacity-60 text-xl">Submissions</p>
                <div class="text-xs mr-1">
                    <SlideToggle name="taxa-table" bind:checked={showEvidenceTable} active="variant-filled-tertiary dark:variant-filled-tertiary" size="sm">Evidence</SlideToggle>
                </div>
            </div>
            <div class="pb-10">
                <CandidateTaxaTable priorityTaxa={priorityTaxa} bind:selectedCandidate={selectedPriorityTaxon}></CandidateTaxaTable>
            </div>
            <div>
                {#if selectedPriorityTaxon}
                    <p class="opacity-60 text-xl pb-2">Decisions</p>

                    <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-4 w-full">
                        <CandidateTaxonCard priorityTaxon={selectedPriorityTaxon}></CandidateTaxonCard>
                        <CandidateTaxonComments priorityTaxon={selectedPriorityTaxon}></CandidateTaxonComments>
                    </div>

                {/if}
            </div>
        </div>
    {/if}

</div>