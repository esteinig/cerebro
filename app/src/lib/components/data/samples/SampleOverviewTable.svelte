<script lang="ts">
	import type { ProjectCollection, RunConfig, SampleOverviewData, Team, TeamDatabase } from "$lib/utils/types";
    import { PriorityTaxonType } from "$lib/utils/types";
	import { Paginator, type PaginationSettings } from "@skeletonlabs/skeleton";
    import { baseTags, getDateTimeString } from "$lib/utils/helpers";
	import CandidateIcon from "$lib/general/icons/CandidateIcon.svelte";
	import FileTagChip from "$lib/general/icons/FileTagChip.svelte";

    export let sampleOverviewData: Array<SampleOverviewData>;
    export let selectedSampleOverview: SampleOverviewData;
    export let selectedSamples: string[] = [];
    
    export let selectedTeam: Team;
    export let selectedDatabase: TeamDatabase;
    export let selectedProject: ProjectCollection;

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

    const getLatestRunId = (sample: SampleOverviewData): string | undefined => {
        return sample.runs
            .sort((a, b) => b.date.localeCompare(a.date)) // Sort by date descending
            [0]?.id; // Get the id of the first (latest) run
    };

</script>

<div class="table-container">
    <table class="table table-hover table-compact">
        <thead>
			<tr>
                <th>Sample</th>
                <th>Libraries</th>
                <th>Tags</th>
                <th>Specimen</th>
                <th>Group</th>
                <th>Candidates</th>
                <th>RunID</th>
                <th>Completed</th>
                <th></th>
                <th>
                    <div class="flex justify-end mr-2">
                        <input class="checkbox" type="checkbox" on:click={changeSelectAll}/>
                    </div>
                </th>
            </tr>
        </thead>
        <tbody>
            {#each tableData as sample}
                <tr class="hover:cursor-pointer items-center align-center" on:click={() => {selectedSampleOverview = sample}}>

                        <td class="pt-1"><span class="ml-1 text-base">{sample.id}</span></td>
                        <td class="truncate"><span class="ml-1 text-base"><FileTagChip tags={baseTags(sample.tags, true)} join={false}></FileTagChip></td>
                        <td class="truncate"><span class="ml-1 text-base"><FileTagChip tags={baseTags(sample.tags, true, ["ENV", "NTC", "PS", "S", "NS"])}></FileTagChip></td>
                        <td class="truncate"><span class="ml-1 text-base">{sample.types.join(", ")}</span></td>
                        <td class="truncate"><span class="ml-1 text-base">{sample.groups.join(", ")}</span></td>
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
                        <td><span class="ml-1 text-base">{getLatestRunId(sample)}</span></td>
                        <td><span class="ml-1 text-base">{getDateTimeString(sample.latest_workflow)}</span></td>
                        <td>
                            <div class="text-base -mt-1">
                                <!-- <button class="btn btn-sm variant-outline-primary mr-1">
                                    <div class="w-4 h-4 mr-2 ">
                                        <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="size-6">
                                            <path stroke-linecap="round" stroke-linejoin="round" d="m16.862 4.487 1.687-1.688a1.875 1.875 0 1 1 2.652 2.652L6.832 19.82a4.5 4.5 0 0 1-1.897 1.13l-2.685.8.8-2.685a4.5 4.5 0 0 1 1.13-1.897L16.863 4.487Zm0 0L19.5 7.125" />
                                        </svg>  
                                    </div>
                                    Edit
                                </button> -->
                                <a href={`${selectedSampleOverview?.id}/workflow=0`} class="btn btn-sm variant-outline-secondary mr-1 opacity-40">
                                    <div class="w-4 h-4 mr-2 ">
                                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                            <path d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" stroke-linecap="round" stroke-linejoin="round"></path>
                                        </svg>
                                    </div>
                                    Report
                                </a>
                                <a href={`team=${selectedTeam.name}&db=${selectedDatabase.name}&project=${selectedProject.name}/${selectedSampleOverview?.id}/workflow=0/qc`} class="btn btn-sm variant-outline-primary">
                                    <div class="w-4 h-4 mr-2">
                                        <svg aria-hidden="true" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg">
                                            <path clip-rule="evenodd" d="M3 10a.75.75 0 01.75-.75h10.638L10.23 5.29a.75.75 0 111.04-1.08l5.5 5.25a.75.75 0 010 1.08l-5.5 5.25a.75.75 0 11-1.04-1.08l4.158-3.96H3.75A.75.75 0 013 10z" fill-rule="evenodd"></path>
                                        </svg>
                                    </div>
                                    Results
                                </a>
                            </div>
                        </td>
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