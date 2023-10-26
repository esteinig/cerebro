<script lang="ts">
	import type { SampleOverviewData } from "$lib/utils/types";
    import { PriorityTaxonType } from "$lib/utils/types";
	import { Paginator, type PaginationSettings } from "@skeletonlabs/skeleton";
    import { baseTags, getDateTimeString } from "$lib/utils/helpers";
	import CandidateIcon from "$lib/general/icons/CandidateIcon.svelte";

    export let sampleOverviewData: Array<SampleOverviewData>;
    export let selectedSampleOverview: SampleOverviewData;
    export let selectedSamples: string[] = [];

    let tableData: Array<SampleOverviewData>;
    
    let checkAll: boolean = false;
    
    let paginationSettings: PaginationSettings = {
        page: 0,
        limit: 50,
        size: sampleOverviewData.length,
        amounts: [5, 10, 50, 100, 500],
    };

    const changeSelectedSample = (selectedSampleId: string) => {

        selectedSamples = selectedSamples.includes(selectedSampleId) ? 
            selectedSamples.filter(sampleId => sampleId !== selectedSampleId) : 
            [...selectedSamples, selectedSampleId]

    }

    const changeSelectAll = () => {
        selectedSamples = checkAll ? [] : sampleOverviewData.map(sample => sample.id);
        checkAll = selectedSamples.length > 0;
    }

	$: {
        paginationSettings.size = sampleOverviewData.length;
		tableData = sampleOverviewData.slice(
			paginationSettings.page * paginationSettings.limit,
			paginationSettings.page * paginationSettings.limit + paginationSettings.limit
		);
        selectedSampleOverview = tableData[0];
    }

    const getUniquePriorityTaxonTypes = (sample: SampleOverviewData): string[] => {
        let taxonTypes = sample.priority.map(x => x.taxon_type);
        return taxonTypes.filter((v,i,a)=>a.indexOf(v)==i)
    }

</script>

<div class="table-container">
    <table class="table table-hover table-compact">
        <thead>
			<tr>
                <th>Sample</th>
                <th>Libraries</th>
                <th>Specimen</th>
                <th>Group</th>
                <th>Candidates</th>
                <th>Completed</th>
                <th>
                    <div class="flex justify-end mr-2">
                        <input class="checkbox" type="checkbox" on:click={changeSelectAll}/>
                    </div>
                </th>
            </tr>
        </thead>
        <tbody>
            {#each tableData as sample}
                <tr class="hover:cursor-pointer" on:click={() => {selectedSampleOverview = sample}}>
                    <td><span class="ml-1">{sample.id}</span></td>
                    <td class="truncate"><span class="ml-1">{baseTags(sample.tags, true).join(", ")}</span></td>
                    <td class="truncate"><span class="ml-1">{sample.types.join(", ")}</span></td>
                    <td class="truncate"><span class="ml-1">{sample.groups.join(", ")}</span></td>
                    <td>
                        <div class="grid grid-cols-3 sm:grid-cols-3 md:grid-cols-3 gap-1 items-center w-2/3 ml-1">
                        
                            {#if getUniquePriorityTaxonTypes(sample).includes(PriorityTaxonType.Pathogen)}
                                <CandidateIcon priorityTaxonType={PriorityTaxonType.Pathogen}></CandidateIcon>
                            {:else}
                                <div></div>
                            {/if}

                            {#if getUniquePriorityTaxonTypes(sample).includes(PriorityTaxonType.Unknown)}
                                <CandidateIcon priorityTaxonType={PriorityTaxonType.Unknown}></CandidateIcon>
                            {:else}
                                <div></div>
                            {/if}

                            {#if getUniquePriorityTaxonTypes(sample).includes(PriorityTaxonType.Contaminant)}
                                <CandidateIcon priorityTaxonType={PriorityTaxonType.Contaminant}></CandidateIcon>
                            {:else}
                                <div></div>
                            {/if}
                        </div>
                    </td>
                    <td><span class="ml-1">{getDateTimeString(sample.latest_workflow)}</span></td>
                    <td>
                        <div class="flex justify-end mr-3">
                            <input class="checkbox" type="checkbox" checked={checkAll} on:click={() => changeSelectedSample(sample.id)}/>
                        </div>
                    </td>
                </tr>
            {/each}
        </tbody>
    </table>
</div>

<Paginator bind:settings={paginationSettings} showFirstLastButtons={false} showPreviousNextButtons={true} class="mt-2" select="paginator-select select text-xs" regionControl="opacity-30 dark:variant-filled-surface dark:opacity-60 text-xs"/>

<style lang="postcss">
    .pathogen-icon {
        fill: turquoise
    }
</style>