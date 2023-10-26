<script lang="ts">
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import type { PriorityTaxon, ProjectCollection, SampleDescriptionSchema, SampleOverviewData, TeamDatabase } from "$lib/utils/types";
	import { ProgressRadial, getToastStore } from "@skeletonlabs/skeleton";
	import { baseTags, getDateTimeString, getDateTimeStringUtc } from "$lib/utils/helpers";
	import { page } from "$app/stores";
	import CandidateIcon from "$lib/general/icons/CandidateIcon.svelte";

    export let selectedSampleOverview: SampleOverviewData;

    export let selectedDatabaseId: string;
    export let selectedProjectId: string;

    export let selectedDatabase: TeamDatabase;
    export let selectedProject: ProjectCollection;

    let edit: boolean = false;
    let loading: boolean = false;

    const publicApi: CerebroApi = new CerebroApi();
    const toastStore = getToastStore();

    let uniquePriorityTaxa: PriorityTaxon[] = [];
    let uniqueTags: string[] = [];
    let descriptions: string[] = [];

    let description: string = "";
    let sampleGroup: string = "";
    let sampleType: string = "";
    
    const updateSample = async() => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.updateSampleDescription}?db=${selectedDatabaseId}&project=${selectedProjectId}&id=${selectedSampleOverview.id}`, {
                method: 'PATCH',
                mode: 'cors',
                credentials: 'include',
                body: JSON.stringify({
                    description: description,
                    sample_group: sampleGroup,
                    sample_type: sampleType,
                } satisfies SampleDescriptionSchema),
                headers:  { 'Content-Type': 'application/json' }
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, "Description updated"
        );
        
        if (response.ok) {
            selectedSampleOverview.description = [description];
            edit = false;
        }

        loading = false;
    }  

    $: {
        uniquePriorityTaxa = selectedSampleOverview?.priority.filter((x, i, a) => a.indexOf(x) == i);
        descriptions = selectedSampleOverview?.description.filter(descr => descr !== null);

        // Display unique base tags (DNA, RNA, NTC) in title
        uniqueTags = baseTags(selectedSampleOverview?.tags, true);
    }

</script>

{#if selectedSampleOverview}
    <div class="card max-w-xl2">
        <header class="card-header">
            <div class="flex p-4">
                <h3 class="h3">{selectedSampleOverview.id}</h3>
                <div>
                    {#each uniqueTags as tag}
                        <span class="code text-xs ml-2">{tag}</span>
                    {/each}
                </div>
                <div class="ml-auto">
                    {#if edit}
                        <button type="button" class="btn btn-md variant-outline-tertiary m-2" on:click={() => { edit = false; } }>
                            <div class="w-5 h-5 mr-2">
                                <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                    <path d="M9.75 9.75l4.5 4.5m0-4.5l-4.5 4.5M21 12a9 9 0 11-18 0 9 9 0 0118 0z" stroke-linecap="round" stroke-linejoin="round"></path>
                                </svg>
                            </div>
                            Cancel
                        </button>
                        <button type="button" class="btn btn-md variant-outline-secondary m-2" on:click={updateSample}>
                            <div class="w-5 h-5 mr-2">
                                <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                    <path d="M9 8.25H7.5a2.25 2.25 0 00-2.25 2.25v9a2.25 2.25 0 002.25 2.25h9a2.25 2.25 0 002.25-2.25v-9a2.25 2.25 0 00-2.25-2.25H15m0-3l-3-3m0 0l-3 3m3-3V15" stroke-linecap="round" stroke-linejoin="round"></path>
                                </svg>
                            </div>
                            Submit
                        </button>
                    {:else}
                        <button type="button" class="btn btn-md variant-outline-secondary m-2" on:click={() => edit = true }>
                            <div class="w-5 h-5 mr-2">
                                <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                    <path d="M16.862 4.487l1.687-1.688a1.875 1.875 0 112.652 2.652L10.582 16.07a4.5 4.5 0 01-1.897 1.13L6 18l.8-2.685a4.5 4.5 0 011.13-1.897l8.932-8.931zm0 0L19.5 7.125M18 14v4.75A2.25 2.25 0 0115.75 21H5.25A2.25 2.25 0 013 18.75V8.25A2.25 2.25 0 015.25 6H10" stroke-linecap="round" stroke-linejoin="round"></path>
                                </svg>
                            </div>
                            Edit
                        </button>
                        <a href={`samples/${selectedSampleOverview.id}/db=${selectedDatabaseId}&project=${selectedProjectId}&workflow=0&db_name=${selectedDatabase.name}&project_name=${selectedProject.name}`} class="btn btn-md variant-outline-primary m-2" on:click={() => {} }>
                            <div class="w-5 h-5 mr-2">
                                <svg aria-hidden="true" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg">
                                    <path clip-rule="evenodd" d="M3 10a.75.75 0 01.75-.75h10.638L10.23 5.29a.75.75 0 111.04-1.08l5.5 5.25a.75.75 0 010 1.08l-5.5 5.25a.75.75 0 11-1.04-1.08l4.158-3.96H3.75A.75.75 0 013 10z" fill-rule="evenodd"></path>
                                </svg>
                            </div>
                            Results
                        </a>
                    {/if}
                </div>
            </div>
        </header>
        <section class="p-4">
            <div class="p-4">
                <p class="text-sm opacity-60 pb-5">
                    Description
                </p>
                {#if loading}
                    <div class="flex justify-center">
                        <ProgressRadial width="sm:w-12 md:w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
                    </div>
                {:else}
                    {#if edit}
                        <textarea class="textarea" rows="2" placeholder="Enter a short sample description..." bind:value={description}/>
                    {:else}
                        {#if descriptions?.length}
                            <p class="text-sm opacity-30">
                                {descriptions.join(";")}
                            </p>
                        {:else}
                            <p class="text-sm opacity-30">
                                No description available for this sample
                            </p>
                        {/if}
                    {/if}
                {/if}
            </div>
            <div class="p-4">
                <p class="text-sm opacity-60">
                    Attributes
                </p>
                    <div class="pt-5">
                        {#if edit}
                            <input class="input w-1/2" type="text" placeholder="Enter a sample type..." bind:value={sampleType}/>
                        {:else}
                            {#each selectedSampleOverview.types as specimen}
                                <span class="code text-sm ml ml-2 mr-1 text-gray-300 bg-gray-600 dark:bg-gray-500 dark:text-gray-300">Specimen</span><span class="code text-sm">{specimen}</span>
                            {/each}
                        {/if}
                        {#if edit}
                            <input class="input w-1/2 mt-5" type="text" placeholder="Enter a sample group..." bind:value={sampleGroup}/>
                        {:else}
                            {#each selectedSampleOverview.groups as group}
                                <span class="code text-sm ml-2 mr-1 text-gray-300 bg-gray-600 dark:bg-gray-500 dark:text-gray-300">Group</span><span class="code text-sm">{group}</span>
                            {/each}
                        {/if}
                    </div>
            </div>
            <div class="p-4">
                <hr class="mt-4">
            </div>
            <div class="p-4">
                <p class="text-sm opacity-60">
                    Candidates
                </p>
                <div class="table-container pt-5 text-xs">
                    {#if uniquePriorityTaxa.length} 
                    <table class="table table-hover table-compact">
                        <thead>
                            <tr>
                                <th>Candidate</th>
                                <th>Domain</th>
                                <th>Scientific Name</th>
                                <th>Submitted By</th>
                                <th>Submitted On</th>
                            </tr>
                        </thead>
                        <tbody>
                            {#each uniquePriorityTaxa as priorityTaxon}
                                <tr class="">
                                    <td>
                                        <div class="flex items-center ml-1"> 
                                            <CandidateIcon priorityTaxonType={priorityTaxon.taxon_type}></CandidateIcon> 
                                            <span class="ml-2 text-sm">{priorityTaxon.taxon_type} </span>
                                        </div>
                                    </td>
                                    <td><span class="ml-1">{priorityTaxon.taxon_overview.domain}</span></td>
                                    <td><span class="ml-1 italic">{priorityTaxon.taxon_overview.name}</span></td>
                                    <td><span class="ml-1">{priorityTaxon.user_name}</span></td>
                                    <td><span class="ml-1">{getDateTimeStringUtc(priorityTaxon.date)}</span></td>
                                </tr>
                            {/each}
                        </tbody>
                    </table>
                    {:else}
                    <p class="flex ml-2 text-sm opacity-30">
                        No candidate taxa have been selected
                    </p>
                    {/if}
                </div>
            </div>

            <div class="p-4">
                <p class="text-sm opacity-60">
                    Reports
                </p>
                <div class="table-container pt-5 text-sm">
                    <table class="table table-hover table-compact">
                        <thead>
                            <tr>
                                <th class='text-center'>Date</th>
                                <th class='text-center'>ID</th>
                                <th class='text-center'>Negative</th>
                                <th class='text-center'>Reviewed</th>
                            </tr>
                        </thead>
                        <tbody>
                            {#each selectedSampleOverview.reports as report}
                                <tr class="">
                                    <td class='text-center'><span class="ml-1">{getDateTimeStringUtc(report.date)}</span></td>
                                    <td class='text-center'><span class="ml-1">{report.id.substring(0, 8)}</span></td>
                                    <td class='text-center'><span class="ml-1">{report.negative}</span></td>
                                    <td class='text-center'><span class="ml-1">{report.review_date}</span></td>
                                </tr>
                            {/each}
                        </tbody>
                    </table>
                </div>
            </div>

            <div class="p-4">
                <hr class="mt-4">
            </div>
            <div class="p-4">
                <p class="text-sm opacity-60">
                    Runs
                </p>
                <div class="table-container pt-5 text-sm">
                    <table class="table table-hover table-compact">
                        <thead>
                            <tr>
                                <th>Run ID</th>
                                <th>Date</th>
                                <th>Description</th>
                                <th>Protocol</th>
                                <th>Version</th>
                            </tr>
                        </thead>
                        <tbody>
                            {#each selectedSampleOverview.runs as run}
                                <tr>
                                    <td><span class="ml-1">{run.id}</span></td>
                                    <td><span class="ml-1">{run.date}</span></td>
                                    <td><span class="ml-1">Phage spike-in experiment</span></td>
                                    <td><span class="ml-1">CNS UMI</span></td>
                                    <td><span class="ml-1">v0.9.0</span></td>
                                </tr>
                            {/each}
                        </tbody>
                    </table>
                </div>
            </div>
            <div class="p-4">
                <p class="text-sm opacity-60">
                    Libraries
                </p>
                <div class="table-container pt-5 text-sm">
                    <table class="table table-hover table-compact">
                        <thead>
                            <tr>
                                <th>Tags</th>
                                <th>Sample</th>
                                <th>Specimen</th>
                                <th>Group</th>
                            </tr>
                        </thead>
                        <tbody>
                            {#each selectedSampleOverview.samples as sample}
                                <tr>
                                    <td><span class="ml-1">
                                    <span class="code text-xs?">{sample.tags.join("-")}</span>
                                    </span></td>
                                    <td><span class="ml-1">{sample.id}</span></td>
                                    <td><span class="ml-1">{sample.sample_type}</span></td>
                                    <td><span class="ml-1">{sample.sample_group}</span></td>
                                    
                                </tr>
                            {/each}
                        </tbody>
                    </table>
                </div>
            </div>
            <div class="p-4">
                <p class="text-sm opacity-60">
                    Workflows
                </p>
                <div class="table-container pt-5 text-sm">
                    <table class="table table-hover table-compact">
                        <thead>
                            <tr>
                                <th>UUID</th>
                                <th>Mnemonic</th>
                                <th>Started</th>
                                <th>Completed</th>
                                <th>Version</th>
                            </tr>
                        </thead>
                        <tbody>
                            {#each selectedSampleOverview.workflows as workflow}
                                <tr class="">
                                    <td><span class="ml-1">{workflow.id.substring(0, 8)}</span></td>
                                    <td><span class="ml-1">{workflow.name}</span></td>
                                    <td><span class="ml-1">{getDateTimeString(workflow.started)}</span></td>
                                    <td><span class="ml-1">{getDateTimeString(workflow.completed)}</span></td>
                                    <td><span class="ml-1">v{workflow.version}</span></td>
                                </tr>
                            {/each}
                        </tbody>
                    </table>
                </div>
            </div>
           
        </section>
        <footer class="card-footer">
            
        </footer>
    </div>
{/if}