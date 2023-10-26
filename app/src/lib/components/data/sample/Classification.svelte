<script lang="ts">
    import { page } from "$app/stores";
    import CerebroApi, { ApiResponse } from "$lib/utils/api";
    import { getToastStore, ProgressRadial, type ToastSettings } from "@skeletonlabs/skeleton";
	import type { Cerebro, CerebroFilterConfig, ClientFilterConfig, GenusOverview, TaxonOverview } from "$lib/utils/types";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import SpeciesOverviewTable from "./taxa/SpeciesOverviewTable.svelte";
	import GenusOverviewTable from "./taxa/GenusOverviewTable.svelte";
    
    export let selectedIdentifiers: string[] = [];
    export let selectedModels: Cerebro[] = [];
    export let clientFilterConfig: ClientFilterConfig;
    export let serverFilterConfig: CerebroFilterConfig;

    let publicApi = new CerebroApi();
    let toastStore = getToastStore();

    let loading: boolean = false;
    let taxaOverview: TaxonOverview[] = [];
    let genusOverview: GenusOverview[] = [];
    let modelNameTags: Map<string, string[]> = new Map();

    const getAggregatedTaxaOverview = async(selectedIdentifiers: string[]) => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.taxa}?db=${$page.params.db}&project=${$page.params.project}&id=${selectedIdentifiers.join(",")}&overview=true`,
            { 
                method: 'POST',  
                mode: 'cors',
                credentials: 'include', 
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify(serverFilterConfig) 
            } as RequestInit,
            $page.data.refreshToken, toastStore, "Taxonomic overview loaded"
        )
        
        loading = false;
        
        if (response.ok){
            // Get taxon overview and sort by total RPM
            taxaOverview = response.json.data.taxa;
            taxaOverview.sort((a, b) => a.rpm < b.rpm ? 1 : -1)

            // Get genus specific view groupings
            let genusGroups = taxaOverview.reduce((group: {[key: string]: GenusOverview}, item) => {
                
                // TODO: Make sure we are not excluding unusual ones!
                if (item.genus === null) return group;

                if (!group[item.genus]) {
                    group[item.genus] = {
                        taxid: "",
                        genus: item.genus,
                        domain: item.domain,
                        species: []
                    } satisfies GenusOverview;
                }

                group[item.genus].species.push(item);
                return group;
                }, {}
            );
            genusOverview = Object.entries(genusGroups).map(([_, overview]) => overview);
        }
    }


    

    $: {
        getAggregatedTaxaOverview(selectedIdentifiers);
        modelNameTags = new Map(selectedModels.map(cerebro => [cerebro.name, cerebro.sample.tags]));
    }
    let showGenusOverview: boolean = false;

</script>

<div>
    {#if loading}
        <div class="flex justify-center py-24">
            <ProgressRadial width="sm:w-12 md:w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
        </div>
    {:else}
        {#if !taxaOverview.length}

            <div class="flex justify-center py-16 ">
                <ErrorAnimation></ErrorAnimation>
            </div>
            <p class="flex justify-center text-lg pb-4">No taxa available</p>
            
        {:else}
            {#if showGenusOverview}
                <GenusOverviewTable genusOverview={genusOverview}></GenusOverviewTable>
            {:else}
                <SpeciesOverviewTable selectedIdentifiers={selectedIdentifiers} taxonOverview={taxaOverview} modelNameTags={modelNameTags} clientFilterConfig={clientFilterConfig} serverFilterConfig={serverFilterConfig}></SpeciesOverviewTable>
            {/if}
        {/if}
    {/if}

</div>


