<script lang="ts">
	import { AbundanceMode, DisplayData } from "$lib/utils/types";
	import SpeciesOverviewTable from "./taxa/SpeciesOverviewTable.svelte";
	import TaxonHeatmap from "$lib/components/visualisations/taxa/heatmap/TaxonHeatmap.svelte";
	import VisualisationSelection from "./taxa/visualisations/VisualisationSelection.svelte";
	import { selectedModels, selectedIdentifiers, selectedTaxa } from "$lib/stores/stores";
	import TableDataSelection from "./taxa/visualisations/TableDataSelection.svelte";
	import TableModeSelection from "./taxa/visualisations/TableModeSelection.svelte";
	import TaxonHistory from "$lib/components/visualisations/taxa/history/TaxonHistory.svelte";
    
    let selectedVisualisation: string = "Table";

    let selectedTableData: DisplayData = DisplayData.Rpm;
    let selectedTableMode: AbundanceMode = AbundanceMode.Mixed;

    $: numSelectedModels = $selectedModels.length;

    $: {
        selectedVisualisation
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
    {#if selectedVisualisation === "Heatmap"}
        <TaxonHeatmap selectedModels={$selectedModels} selectedIdentifiers={$selectedIdentifiers} displayData={DisplayData.Rpm}></TaxonHeatmap>
    {:else if selectedVisualisation === "History"}
        {#each $selectedTaxa as selectedTaxon}
            <TaxonHistory selectedTaxon={selectedTaxon}></TaxonHistory>
        {/each}
    {:else}
        <p>Not found</p>
    {/if}
    <SpeciesOverviewTable displayData={selectedTableData} displayMode={selectedTableMode} disableDrawer={selectedVisualisation !== "Table"}></SpeciesOverviewTable>
</div>


