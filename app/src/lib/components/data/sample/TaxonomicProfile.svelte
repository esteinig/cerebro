<script lang="ts">
    import { page } from "$app/stores";
    import CerebroApi, { ApiResponse } from "$lib/utils/api";
    import { getToastStore, ProgressRadial } from "@skeletonlabs/skeleton";
	import { type Taxon, DisplayData } from "$lib/utils/types";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import SpeciesOverviewTable from "./taxa/SpeciesOverviewTable.svelte";
	import TaxonHeatmap from "$lib/components/visualisations/taxa/heatmap/TaxonHeatmap.svelte";
	import VisualisationSelection from "./taxa/visualisations/VisualisationSelection.svelte";
	import { selectedModels, selectedIdentifiers, selectedServerFilterConfig } from "$lib/stores/stores";
    
    let publicApi = new CerebroApi();
    let toastStore = getToastStore();

    let loading: boolean = false;
    let taxa: Taxon[] = [];
    let selectedVisualisation: string = "Table";

    const getAggregatedTaxaOverview = async(selectedIdentifiers: string[]) => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.taxa}?team=${$page.params.team}&db=${$page.params.db}&project=${$page.params.project}&id=${selectedIdentifiers.join(",")}&overview=true`,
            { 
                method: 'POST',  
                mode: 'cors',
                credentials: 'include', 
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify($selectedServerFilterConfig) 
            } as RequestInit,
            $page.data.refreshToken, toastStore, "Taxa loaded"
        )
        
        loading = false;
        
        if (response.ok){
            taxa = response.json.data.taxa;
        }
    }    

    $: if ($selectedIdentifiers.length > 0) {
        getAggregatedTaxaOverview($selectedIdentifiers);
    }

</script>

<div class="flex justify-center pt-4 pb-8">
    <VisualisationSelection bind:selectedVisualisation={selectedVisualisation}></VisualisationSelection>
</div>
<div>
    {#if loading}
        <div class="flex justify-center py-24">
            <ProgressRadial width="sm:w-12 md:w-24" stroke={20} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
        </div>
    {:else}
        {#if !taxa.length}
            <div class="flex justify-center py-16 "><ErrorAnimation></ErrorAnimation></div>
            <p class="flex justify-center text-lg pb-4">No taxa available</p>
        {:else}
        
        {#if selectedVisualisation === "Heatmap"}
            <TaxonHeatmap selectedModels={$selectedModels} selectedIdentifiers={$selectedIdentifiers} displayData={DisplayData.Rpm}></TaxonHeatmap>
        {/if}
        <SpeciesOverviewTable taxa={taxa}></SpeciesOverviewTable>
        {/if}
    {/if}

</div>


