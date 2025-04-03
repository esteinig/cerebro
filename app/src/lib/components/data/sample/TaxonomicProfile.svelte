<script lang="ts">
	import { AbundanceMode, DisplayData, DisplayVisualisation } from "$lib/utils/types";
	import SpeciesOverviewTable from "./taxa/SpeciesOverviewTable.svelte";
	import VisualisationSelection from "./taxa/visualisations/VisualisationSelection.svelte";
	import { primaryThresholdServerFilterConfig, secondaryThresholdServerFilterConfig, selectedModels, selectedTaxa, validationThresholdServerFilterConfig } from "$lib/stores/stores";
	import TableDataSelection from "./taxa/visualisations/TableDataSelection.svelte";
	import TableModeSelection from "./taxa/visualisations/TableModeSelection.svelte";
	import TaxonHistory from "$lib/components/visualisations/taxa/history/TaxonHistory.svelte";
    
    let selectedVisualisation: DisplayVisualisation = DisplayVisualisation.Threshold;

    let selectedTableData: DisplayData = DisplayData.Rpm;
    let selectedTableMode: AbundanceMode = AbundanceMode.Mixed;

    $: numSelectedModels = $selectedModels.length;

    $: {
        selectedVisualisation // on any change in selected visualisation reset selected taxa to empty
        $selectedTaxa = [];
    }

    $: {
        $selectedModels // on any change in selected models reset selected taxa to empty
        $selectedTaxa = [];
    }
</script>


<div class="flex justify-center pt-4 pb-8 gap-8">
    <VisualisationSelection bind:selectedVisualisation={selectedVisualisation} numberModels={numSelectedModels}></VisualisationSelection>
    {#if selectedVisualisation === "Table"}
        <TableDataSelection bind:selectedTableData={selectedTableData}></TableDataSelection>
        {#if selectedVisualisation === "Table"}
            <TableModeSelection bind:selectedTableMode={selectedTableMode}></TableModeSelection>
        {/if}
    {/if}
</div>
<div>
    {#if selectedVisualisation === DisplayVisualisation.Run}
        <p class="flex justify-center py-12 opacity-60">Select taxa from the table for run comparison</p>
    {:else if selectedVisualisation === DisplayVisualisation.History}
        {#if $selectedTaxa.length > 0}
            {#each $selectedTaxa as selectedTaxon}
                <TaxonHistory selectedTaxon={selectedTaxon}></TaxonHistory>
            {/each}
        {:else}
            <p class="flex justify-center py-12 opacity-60">Select taxa from the table for historical comparison</p>
        {/if}
    {:else if selectedVisualisation === DisplayVisualisation.Threshold}
    {#if $selectedTaxa.length > 0}
            {#each $selectedTaxa as selectedTaxon}
                <TaxonHistory selectedTaxon={selectedTaxon}></TaxonHistory>
            {/each}
        {:else}
            <p class="flex justify-center py-12 opacity-60">Select taxa from the table for historical comparison</p>
        {/if}
        <p class="text-lg opacity-60 py-12">Threshold taxonomic profiling</p>
        <SpeciesOverviewTable displayData={selectedTableData} displayMode={selectedTableMode} selectedVisualisation={selectedVisualisation} serverFilterConfig={primaryThresholdServerFilterConfig}></SpeciesOverviewTable>
        <p class="text-lg opacity-60 py-12">Sub-threshold taxonomic profiling</p>
        <SpeciesOverviewTable displayData={selectedTableData} displayMode={selectedTableMode} selectedVisualisation={selectedVisualisation} serverFilterConfig={secondaryThresholdServerFilterConfig}></SpeciesOverviewTable>
        <p class="text-lg opacity-60 py-12">Validation plate</p>
        <SpeciesOverviewTable displayData={DisplayData.Rpm} displayMode={selectedTableMode} selectedVisualisation={selectedVisualisation} serverFilterConfig={validationThresholdServerFilterConfig} disablePrevalenceContamination={true}></SpeciesOverviewTable>
    {:else}

        <SpeciesOverviewTable displayData={selectedTableData} displayMode={selectedTableMode} selectedVisualisation={selectedVisualisation}></SpeciesOverviewTable>
    {/if}
</div>


