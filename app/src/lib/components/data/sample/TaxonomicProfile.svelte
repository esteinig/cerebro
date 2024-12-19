<script lang="ts">
	import { type Taxon, DisplayData } from "$lib/utils/types";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import SpeciesOverviewTable from "./taxa/SpeciesOverviewTable.svelte";
	import TaxonHeatmap from "$lib/components/visualisations/taxa/heatmap/TaxonHeatmap.svelte";
	import VisualisationSelection from "./taxa/visualisations/VisualisationSelection.svelte";
	import { selectedModels, selectedIdentifiers, selectedServerFilterConfig } from "$lib/stores/stores";
    
    let selectedVisualisation: string = "Table";


</script>

<div class="flex justify-center pt-4 pb-8">
    <VisualisationSelection bind:selectedVisualisation={selectedVisualisation}></VisualisationSelection>
</div>
<div>
    
    {#if selectedVisualisation === "Heatmap"}
        <TaxonHeatmap selectedModels={$selectedModels} selectedIdentifiers={$selectedIdentifiers} displayData={DisplayData.Rpm}></TaxonHeatmap>
    {/if}
    
    <SpeciesOverviewTable></SpeciesOverviewTable>
    
</div>


