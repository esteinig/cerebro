<script lang="ts">
	import { page } from "$app/stores";
	import { InputChip, getModalStore, type ToastSettings, type ModalSettings } from "@skeletonlabs/skeleton";
	import SampleOverviewTable from "./SampleOverviewTable.svelte";
    import { Role } from "$lib/utils/types";
	import type { ProjectCollection, SampleDeleteSchema, SampleOverviewData, SampleSummarySchema, Team, TeamDatabase } from "$lib/utils/types";
	import SampleOverviewCard from "./SampleOverviewCard.svelte";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
    import { ProgressRadial, getToastStore } from '@skeletonlabs/skeleton';
	import { getDateTimeString } from "$lib/utils/helpers";
	import SampleViewSelection from "./SampleViewSelection.svelte";

    const modalStore = getModalStore();
    const toastStore = getToastStore();

    const publicApi = new CerebroApi();

    let sampleOverviewData: Array<SampleOverviewData> = $page.data.sampleOverview;

    let selectedTeamId: string = $page.data.defaultTeam.id;
    let selectedDatabaseId: string = $page.data.defaultDatabase.id;
    let selectedProjectId: string = $page.data.defaultProject.id;

    let selectedTeam: Team = $page.data.defaultTeam;
    let selectedDatabase: TeamDatabase = $page.data.defaultDatabase;
    let selectedProject: ProjectCollection = $page.data.defaultProject;
    
    let databases: Array<TeamDatabase> = $page.data.defaultTeam.databases;
    let projects: Array<ProjectCollection> = $page.data.defaultDatabase.projects;

    let searchTerms: string[] = [];

    let selectedSampleOverview: SampleOverviewData;
    let selectedSamples: string[] = [];

    let loading: boolean = false;


    // Helpers
    const getSelectedTeam = (): Team | null => {
        let matchedTeams = $page.data.userTeams.filter((team: Team) => team.id === selectedTeamId);
        if (!matchedTeams.length) {
            toastStore.trigger(
                { message: "Team could not be selected", background: "variant-filled-error" } satisfies ToastSettings
            )
            return null
        }
        return matchedTeams[0]
    }
    const getSelectedDatabase = (team: Team, first: boolean = false): TeamDatabase | null => {
        
        let fail: boolean;

        if (first) {
            fail = !team.databases.length
        } else {
            let matchedDatabases = team.databases.filter((db: TeamDatabase) => db.id === selectedDatabaseId);
            fail = !matchedDatabases.length
        }

        if (fail) {
            toastStore.trigger(
                { message: "No databases available for selected team", background: "variant-filled-error" } satisfies ToastSettings
            )
            return null
        }
        return team.databases[0]
    }
    const getSelectedProject = (database: TeamDatabase, first: boolean = false): ProjectCollection | null => {

        let fail: boolean;

        if (first) {
            fail = !database.projects.length
        } else {
            let matchedProjects = database.projects.filter((project: ProjectCollection) => project.id === selectedProjectId);
            fail = !matchedProjects.length
        }

        if (fail) {
            toastStore.trigger(
                { message: "No projects available for selected database", background: "variant-filled-error" } satisfies ToastSettings
            )
            return null
        }
        return database.projects[0]
    }


    // Selection helpers
    const changeTeam = async() => {

        let team: Team | null = getSelectedTeam();

        if (team !== null) {
            let teamDatabase = getSelectedDatabase(team, true);

            if (teamDatabase !== null) {
                let teamProject = getSelectedProject(teamDatabase, true);
            
                if (teamProject !== null) {
                    selectedDatabaseId = teamDatabase.id;
                    selectedDatabase = teamDatabase;
                    databases = team.databases;

                    selectedProjectId = teamProject.id;
                    selectedProject = teamProject;
                    projects = teamDatabase.projects;

                    await loadData();
                }
            }
        }
    }

    const changeDatabase = async() => {

        let team: Team | null = getSelectedTeam();

        if (team !== null) {
            let teamDatabase = getSelectedDatabase(team);

            if (teamDatabase !== null) {
                let teamProject = getSelectedProject(teamDatabase, true);
            
                if (teamProject !== null) {
                    selectedProjectId = teamProject.id;
                    selectedProject = teamProject;
                    projects = teamDatabase.projects;

                    await loadData();
                }
            }
        }
    }

    const changeProject = async() => {

        let team: Team | null = getSelectedTeam();

        if (team !== null) {
            let teamDatabase = getSelectedDatabase(team);

            if (teamDatabase !== null) {
                let teamProject = getSelectedProject(teamDatabase);
            
                if (teamProject !== null) {
                    selectedProject = teamProject;

                    await loadData();
                }
            }
        }
    }

    // Data request functions

    const loadData = async() => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.sampleOverview}?db=${selectedDatabaseId}&project=${selectedProjectId}&page=0&limit=${$page.data.defaultPageLimit}&notag=${$page.data.defaultNegativeTemplateControl}`, {
                method: 'GET',
                mode: 'cors',
                credentials: 'include'
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, null
        );
        
        if (response.ok) {
            sampleOverviewData = response.json.data?.sample_overview;

            // We check for an ok but empty response - the request
            // was successful, but there are no samples in the 
            // requested database and project
            if (sampleOverviewData.length == 0){
                toastStore.trigger(
                    { message: "No samples available", background: "variant-filled-tertiary" } satisfies ToastSettings
                )
            }
        }

        loading = false;

    }

    const deleteSelectedSamples = async() => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.deleteSamples}?db=${selectedDatabaseId}&project=${selectedProjectId}`, {
                method: 'DELETE',
                mode: 'cors',
                credentials: 'include',
                body: JSON.stringify({sample_id: selectedSamples} satisfies SampleDeleteSchema),
                headers:  { 'Content-Type': 'application/json' }
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, `${selectedSamples.length > 1 ? 'Samples' : 'Sample'} deleted`
        );
        
        if (response.ok) {
            sampleOverviewData = sampleOverviewData.filter(overview => !selectedSamples.includes(overview.id))
            selectedSamples = []
        }
        
        loading = false;

    }

    const getSelectedSampleCsv = async() => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.getSampleSummary}?db=${selectedDatabaseId}&project=${selectedProjectId}&csv=true`, {
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

    let newProjectName: string = "";
    let newProjectDesription: string = "";
    let newDatabaseName: string = "";
    let newDatabaseDescription: string = "";
    let create: string = 'Project';

</script>

    
<div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-4 gap-16 pt-10">
    <div class="gap-y-4 md:mt-5 col-span-1 text-lg">
        <div class="mb-4">
            <SampleViewSelection selectedView={"samples"}></SampleViewSelection>
        </div>
        <div class="mb-4 p-4 pl-6 border-token border-primary-500">
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
                            {#each ['Project', 'Database'] as c}
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
                            {#each ['Project', 'Database'] as c}
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
                    <select id="teamSelect" class="select" bind:value={selectedTeamId} on:change={changeTeam}>
                        {#each $page.data.userTeams as team}
                            <option value={team.id}>{team.name}</option>
                        {/each}
                    </select>
                </div>

                <div class="w-1/4">
                    <p class="mb-1"><span class="opacity-60">Database</span></p>
                    <select id="projectSelect"  class="select" bind:value={selectedDatabaseId} on:change={changeDatabase}>
                        {#each databases as db}
                            <option value={db.id}>{db.name}</option>
                        {/each}
                    </select>
                </div>
                <div class="w-1/4">
                    <p class="mb-1"><span class="opacity-60">Project</span></p>
                    <select id="projectSelect"  class="select" bind:value={selectedProjectId} on:change={changeProject}>
                        {#each projects as project}
                            <option value={project.id}>{project.name}</option>
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
                                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" stroke-linecap="round" stroke-linejoin="round"></path>
                                    </svg>
                                </div>            
                                Taxa    
                            </button>
                        </div>
                        <div class="mt-7">
                            <button type="button" class="btn btn-md variant-outline-secondary" on:click={getSelectedSampleCsv}>

                                <div class="w-5 h-5 mr-2">
                                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" stroke-linecap="round" stroke-linejoin="round"></path>
                                    </svg>
                                </div>            
                                Quality     
                            </button>
                        </div>
                    </div>
                {/if}
            </div>
            <div class="mt-5">
                <SampleOverviewTable 
                sampleOverviewData={sampleOverviewData} 
                bind:selectedSampleOverview 
                bind:selectedSamples 
                selectedTeamId={selectedTeamId} 
                selectedTeamName={selectedTeam.name} 
                selectedDatabaseId={selectedDatabaseId} 
                selectedProjectId={selectedProjectId} 
                selectedDatabaseName={selectedDatabase.name} 
                selectedProjectName={selectedProject.name}
                ></SampleOverviewTable>
            </div>
        {/if}
    </div>
</div>

