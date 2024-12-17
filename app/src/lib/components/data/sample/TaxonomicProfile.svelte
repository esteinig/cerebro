<script lang="ts">
    import { page } from "$app/stores";
    import CerebroApi, { ApiResponse } from "$lib/utils/api";
    import { getToastStore, ProgressRadial } from "@skeletonlabs/skeleton";
	import { type Cerebro, type TaxonFilterConfig, type ClientFilterConfig, type Taxon, DisplayData } from "$lib/utils/types";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import SpeciesOverviewTable from "./taxa/SpeciesOverviewTable.svelte";
	import TaxonHeatmap from "$lib/components/visualisations/taxa/heatmap/TaxonHeatmap.svelte";
	import VisualisationSelection from "./taxa/visualisations/VisualisationSelection.svelte";


    export let selectedIdentifiers: string[] = [];
    export let selectedModels: Cerebro[] = [];

    
    export let serverFilterConfig: TaxonFilterConfig;

    // export let clientFilterConfig: ClientFilterConfig;
    // export let taxonHighlightConfig: TaxonHighlightConfig;

    let publicApi = new CerebroApi();
    let toastStore = getToastStore();

    let loading: boolean = false;
    let taxa: Taxon[] = [];
    let modelNameTags: Map<string, string[]> = new Map();
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
                body: JSON.stringify(serverFilterConfig) 
            } as RequestInit,
            $page.data.refreshToken, toastStore, "Taxa loaded"
        )
        
        loading = false;
        
        if (response.ok){
            taxa = response.json.data.taxa;
        }
    }    

    $: if (selectedIdentifiers.length > 0) {
        // Clear the selected taxa store whenever we select new libraries
        // selectedTaxa.set([]);
        getAggregatedTaxaOverview(selectedIdentifiers);
    }

    $: {
        modelNameTags = new Map(selectedModels.map(cerebro => [cerebro.name, cerebro.sample.tags]));
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
            <TaxonHeatmap selectedModels={selectedModels} selectedIdentifiers={selectedIdentifiers} displayData={DisplayData.Rpm}></TaxonHeatmap>
        {/if}
        <SpeciesOverviewTable taxa={taxa}></SpeciesOverviewTable>
     
            <!-- <DatatableTaxonOverview data={taxaOverview}></DatatableTaxonOverview> -->
        {/if}
    {/if}

</div>


