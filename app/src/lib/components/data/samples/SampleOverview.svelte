<script lang="ts">
	import { page } from "$app/stores";
	import { InputChip, getModalStore,type ModalSettings } from "@skeletonlabs/skeleton";
	import SampleOverviewTable from "./SampleOverviewTable.svelte";
    import { Role } from "$lib/utils/types";
	import type { SampleDeleteSchema, SampleOverviewData, SampleSummarySchema } from "$lib/utils/types";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
    import { ProgressRadial, getToastStore } from '@skeletonlabs/skeleton';
	import { getDateTimeString } from "$lib/utils/helpers";
	import SampleViewSelection from "./SampleViewSelection.svelte";
	import { goto } from "$app/navigation";

    const modalStore = getModalStore();
    const toastStore = getToastStore();

    const publicApi = new CerebroApi();

    export let sampleOverviewData: SampleOverviewData[] = [];

    // All need to be reactive since navigating to a new page with changes
    // in query parameters invalidates server side data, but does not 
    // update client-side data\

    $: selectedTeamName = $page.data.selectedTeam.name;
    $: selectedDatabaseName = $page.data.selectedDatabase.name;
    $: selectedProjectName = $page.data.selectedProject.name;

    $: selectedTeam = $page.data.selectedTeam;
    $: selectedDatabase = $page.data.selectedDatabase;
    $: selectedProject = $page.data.selectedProject;

    let searchTerms: string[] = [];
    let selectedSamples: string[] = [];

    let loading: boolean = false;

    function clearData() {
        selectedSamples = []
        searchTerms = []
    }

    async function changeTeam() {
        clearData()
        await goto(`/cerebro/data/samples/team=${selectedTeamName}&db=0&project=0`, { invalidateAll: true })
    }
    async function changeDatabase() {
        clearData()
        await goto(`/cerebro/data/samples/team=${selectedTeamName}&db=${selectedDatabaseName}&project=0`, { invalidateAll: true })
    }
    async function changeProject() {
        clearData()
        await goto(`/cerebro/data/samples/team=${selectedTeamName}&db=${selectedDatabaseName}&project=${selectedProjectName}`, { invalidateAll: true })
    }

    async function deleteSelectedSamples() {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.deleteSamples}?team=${selectedTeam.id}&db=${selectedDatabase.id}&project=${selectedProject.id}`, {
                method: 'DELETE',
                mode: 'cors',
                credentials: 'include',
                body: JSON.stringify({sample_id: selectedSamples} satisfies SampleDeleteSchema),
                headers:  { 'Content-Type': 'application/json' }
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, `${selectedSamples.length > 1 ? 'Samples' : 'Sample'} deleted`
        );
        
        if (response.ok) {
            sampleOverviewData = sampleOverviewData.filter((overview: SampleOverviewData) => !selectedSamples.includes(overview.id))
            selectedSamples = []
        }
        
        loading = false;

    }

    async function getSelectedSampleCsv() {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.getSampleSummary}?team=${selectedTeam.id}&db=${selectedDatabase.id}&project=${selectedProject.id}&csv=true`, {
                method: 'POST',
                mode: 'cors',
                credentials: 'include',
                body: JSON.stringify({sample_ids: selectedSamples, cerebro_ids: []} satisfies SampleSummarySchema),
                headers:  { 'Content-Type': 'application/json' }
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, `Sample summary retrieved`
        );

        if (response.ok) {
            selectedSamples = []
            if (response.json.data?.csv) {
                exportCsv(response.json.data.csv)
            }
        }

        loading = false;

    }

    const exportCsv = (csv: string) => {
        var blob = new Blob([csv], {type: 'text/csv'});
        var link = document.createElement("a");
        if (link.download !== undefined) { // feature detection
            // Browsers that support HTML5 download attribute
            let url = URL.createObjectURL(blob);
            link.setAttribute("href", url);
            let now = getDateTimeString(new Date().toISOString(), false);
            link.setAttribute("download", `${now}_QualityControlSummary_${selectedTeam.name}_${selectedDatabase.name}_${selectedProject.name}.csv`);
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }
    }

    const openSampleDeleteVerification = async() => {

        const modal: ModalSettings = {
            type: 'confirm',
            title: 'Confirm Deletion',
            body: `<p class="opacity-60 mt-2">Delete results of selected samples in <span>${selectedDatabase.name}</span> <span>${selectedProject.name}</span></p>`,
            response: async(confirmed: boolean) => {
                if (confirmed) {
                    await deleteSelectedSamples()
                }
            }
        };
        modalStore.trigger(modal);
    }

    async function createTeamDatabase() {

        // loading = true;

        // let response: ApiResponse = await publicApi.fetchWithRefresh(
        //     `${publicApi.routes.teams.createProject}?team=${selectedTeamId}`, {
        //         method: 'POST',
        //         mode: 'cors',
        //         credentials: 'include',
        //         body: JSON.stringify({

        //         } satisfies ),
        //         headers:  { 'Content-Type': 'application/json' }
        //     } satisfies RequestInit,
        //     $page.data.refreshToken, toastStore, `Sample summary retrieved`
        // );

        // if (response.ok) {
        //     selectedSamples = []
        //     if (response.json.data?.csv) {
        //         exportCsv(response.json.data.csv)
        //     }
        // }
    }

    enum TeamDataCrud {
        CreateProject = 'Project',
        CreateDatabase = 'Database'
    }

    let newProjectName: string = "";
    let newProjectDesription: string = "";
    let newDatabaseName: string = "";
    let newDatabaseDescription: string = "";

    let create: TeamDataCrud = TeamDataCrud.CreateProject;

    let selectedView: string = "samples";

</script>

    
<div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-4 gap-16 pt-10">
    <div class="gap-y-4 md:mt-5 col-span-1 text-lg">
        <p class="opacity-80">Database collections</p>
        <div class="mb-4 p-4 border border-primary-500">
            <SampleViewSelection bind:selectedView={selectedView}></SampleViewSelection>
        </div>
        <div class="mb-4 p-4 pl-6 border border-primary-500">
            <div class="mb-4">
                <ol class="breadcrumb justify-lstart">
                    <li class="crumb opacity-70">Team {selectedTeam.name}</li>
                    <li class="crumb-separator" aria-hidden>&rsaquo;</li>
                    <li class="crumb opacity-80">Database {selectedDatabase.name}</li>
                </ol>
            </div>
            <div class="mb-8">
                <div class="flex justify-start">
                    <div class="w-full">
                        <InputChip name="sampleSearch" placeholder=" Search and subset table ... " bind:value={searchTerms} allowUpperCase/> 
                    </div>
                </div>
            </div>
            {#if create == "Project"}
                <div class="mt-10 mb-4">
                    <p class="opacity-80 mb-4">Create a new project</p>
                    <p class="opacity-40 mb-4 text-xs">New projects are the primary way to store workflow outputs in a collection.</p>

                    <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 gap-4 items-center w-full mb-4">
                        <div class="">
                            <label class="label text-sm">
                                <span class="opacity-60">Project name</span>
                                <input class="input w-full" type="text" placeholder="..." bind:value={newProjectName}/>
                            </label>
                        </div>
                        <div class="">
                            <label class="label text-sm">
                                <span class="opacity-60">Project description</span>
                                <textarea class="textarea" rows="3" placeholder="..." bind:value={newProjectDesription} />
                            </label>
                        </div>
                    </div>
                    <div class="flex items-center justify-between">
                        <div class="">
                            {#each [TeamDataCrud.CreateProject, TeamDataCrud.CreateDatabase] as c}
                                <button
                                    class="chip {create === c ? 'variant-filled-primary' : 'variant-soft'} ml-2"
                                    on:click={() => { create = c }}
                                    on:keypress
                                >
                                    <span>{c}</span>
                                </button>
                            {/each} 
                        </div>
                        <button type="button" class="btn btn-md variant-outline-primary" on:click={createTeamDatabase}>
                            <div class="w-5 h-5 mr-2">
                                <svg data-slot="icon" aria-hidden="true" fill="none" stroke-width="1.5" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                    <path d="M2.25 12.75V12A2.25 2.25 0 0 1 4.5 9.75h15A2.25 2.25 0 0 1 21.75 12v.75m-8.69-6.44-2.12-2.12a1.5 1.5 0 0 0-1.061-.44H4.5A2.25 2.25 0 0 0 2.25 6v12a2.25 2.25 0 0 0 2.25 2.25h15A2.25 2.25 0 0 0 21.75 18V9a2.25 2.25 0 0 0-2.25-2.25h-5.379a1.5 1.5 0 0 1-1.06-.44Z" stroke-linecap="round" stroke-linejoin="round"></path>
                                </svg>
                            </div>            
                            Create   
                        </button>
                    </div>
                </div>
            {:else}
                <div class="mt-10 mb-4">
                    <p class="opacity-80 mb-4">Create a new database</p>
                    <p class="opacity-40 mb-4 text-xs">New databases are the primary way to organise project environments e.g. collecting data in production, from validation experiments or clinical research cohorts.</p>

                    <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 gap-4 items-center w-full mb-4">
                        <div class="">
                            <label class="label text-sm">
                                <span class="opacity-60">Database name</span>
                                <input class="input w-full" type="text" placeholder="..." bind:value={newDatabaseName}/>
                            </label>
                        </div>
                        <div class="">
                            <label class="label text-sm">
                                <span class="opacity-60">Database description</span>
                                <textarea class="textarea" rows="3" placeholder="..." bind:value={newDatabaseDescription} />
                            </label>
                        </div>
                    </div>
                    <div class="flex items-center justify-between">
                        <div class="">
                            {#each [TeamDataCrud.CreateProject, TeamDataCrud.CreateDatabase] as c}
                                <button
                                    class="chip {create === c ? 'variant-filled-primary' : 'variant-soft'} ml-2"
                                    on:click={() => { create = c }}
                                    on:keypress
                                >
                                    <span>{c}</span>
                                </button>
                            {/each} 
                        </div>
                        <button type="button" class="btn btn-md variant-outline-primary" on:click={getSelectedSampleCsv}>
                            <div class="w-5 h-5 mr-2">
                                <svg data-slot="icon" aria-hidden="true" fill="none" stroke-width="1.5" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                    <path d="M20.25 6.375c0 2.278-3.694 4.125-8.25 4.125S3.75 8.653 3.75 6.375m16.5 0c0-2.278-3.694-4.125-8.25-4.125S3.75 4.097 3.75 6.375m16.5 0v11.25c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125V6.375m16.5 0v3.75m-16.5-3.75v3.75m16.5 0v3.75C20.25 16.153 16.556 18 12 18s-8.25-1.847-8.25-4.125v-3.75m16.5 0c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125" stroke-linecap="round" stroke-linejoin="round"></path>
                                </svg>
                            </div>            
                            Create   
                        </button>
                    </div>
                </div>
            {/if}
        </div>
    </div>
    <div class="gap-y-4 md:mt-1 col-span-3">
        
        {#if loading}
            <div class="flex justify-center pt-36">
                <ProgressRadial width="sm:w-12 md:w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
            </div>
        {:else}

            <div class="flex gap-x-4 align-center items-end">

                <div class="w-1/4">
                    <p class="mb-1"><span class="opacity-60">Team</span></p>
                    <select id="teamSelect" class="select" bind:value={selectedTeamName} on:change={changeTeam}>
                        {#each $page.data.userTeams as team}
                            <option value={team.name}>{team.name}</option>
                        {/each}
                    </select>
                </div>

                <div class="w-1/4">
                    <p class="mb-1"><span class="opacity-60">Database</span></p>
                    <select id="projectSelect"  class="select" bind:value={selectedDatabaseName} on:change={changeDatabase}>
                        {#each $page.data.selectedTeam.databases as db}
                            <option value={db.name}>{db.name}</option>
                        {/each}
                    </select>
                </div>
                <div class="w-1/4">
                    <p class="mb-1"><span class="opacity-60">Project</span></p>
                    <select id="projectSelect"  class="select" bind:value={selectedProjectName} on:change={changeProject}>
                        {#each $page.data.selectedDatabase.projects as project}
                            <option value={project.name}>{project.name}</option>
                        {/each}
                    </select>
                </div>
                {#if selectedSamples.length > 0}
                    <div class="flex items-center gap-x-4 ml-auto">
                        {#if $page.data.userData.roles.includes(Role.Data)}
                            <div class="mt-7">
                                <button type="button" class="btn btn-large variant-outline-tertiary" on:click={openSampleDeleteVerification}>
                                    <div class="w-6 h-6 mr-2">
                                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                            <path d="M14.74 9l-.346 9m-4.788 0L9.26 9m9.968-3.21c.342.052.682.107 1.022.166m-1.022-.165L18.16 19.673a2.25 2.25 0 01-2.244 2.077H8.084a2.25 2.25 0 01-2.244-2.077L4.772 5.79m14.456 0a48.108 48.108 0 00-3.478-.397m-12 .562c.34-.059.68-.114 1.022-.165m0 0a48.11 48.11 0 013.478-.397m7.5 0v-.916c0-1.18-.91-2.164-2.09-2.201a51.964 51.964 0 00-3.32 0c-1.18.037-2.09 1.022-2.09 2.201v.916m7.5 0a48.667 48.667 0 00-7.5 0" stroke-linecap="round" stroke-linejoin="round"></path>
                                        </svg>
                                    </div>    
                                    Delete
                                </button>
                            </div>
                        {/if}
                        <div class="mt-7">
                            <button type="button" class="btn btn-md variant-outline-secondary" on:click={getSelectedSampleCsv}>

                                <div class="w-5 h-5 mr-2">
                                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="size-6">
                                        <path stroke-linecap="round" stroke-linejoin="round" d="m16.862 4.487 1.687-1.688a1.875 1.875 0 1 1 2.652 2.652L6.832 19.82a4.5 4.5 0 0 1-1.897 1.13l-2.685.8.8-2.685a4.5 4.5 0 0 1 1.13-1.897L16.863 4.487Zm0 0L19.5 7.125" />
                                      </svg>                                      
                                </div>            
                                Edit    
                            </button>
                        </div>
                        <div class="mt-7">
                            <button type="button" class="btn btn-md variant-outline-primary" on:click={getSelectedSampleCsv}>

                                <div class="w-5 h-5 mr-2">
                                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" stroke-linecap="round" stroke-linejoin="round"></path>
                                    </svg>
                                </div>            
                                Data     
                            </button>
                        </div>
                    </div>
                {/if}
            </div>
            <div class="mt-5">
                {#if selectedView === "samples"}
                <SampleOverviewTable 
                    sampleOverviewData={sampleOverviewData} 
                    bind:selectedSamples 
                    selectedTeam={selectedTeam}
                    selectedDatabase={selectedDatabase}
                    selectedProject={selectedProject}
                ></SampleOverviewTable>
                {:else if selectedView === "runs"}
                    
                {:else}

                {/if}
            </div>
        {/if}
    </div>
</div>

