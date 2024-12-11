<script lang="ts">
	import { goto } from "$app/navigation";
    import { page } from "$app/stores";
    import { onMount } from 'svelte';
    import { invalidate } from '$app/navigation';
    import { getToastStore } from "@skeletonlabs/skeleton";

	import { 
        type RegisterStagedSampleSchema, 
        type ProductionTower, 
        type ProductionWatcher, 
        type ProjectCollection, 
        type TeamDatabase, 
		Pipeline

    } from "$lib/utils/types";

	import { getDateTimeStringUtc, isWithinTimeLimit } from "$lib/utils/helpers";
	import ActiveIndicator from "$lib/general/icons/ActiveIndicator.svelte";
	import FileSelection from "$lib/components/production/watchers/FileSelection.svelte";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import CerebroApi from "$lib/utils/api.js";

    export let data;
    
    let toastStore = getToastStore();
    const publicApi = new CerebroApi()


    // Invalidates page data every three minute to update 
    // the watcher/pipeline activity status indicators
    onMount(() => {
        const interval = setInterval(() => {
            invalidate('watchers:data');
        }, 3 * 60 * 1000);

        return () => clearInterval(interval);
    });
    
    async function changeTeam() {
        await goto(`/cerebro/production/watchers/team=${selectedTeamId}&watcher=0`)
    }

    async function changeWatcher() {
        await goto(`/cerebro/production/watchers/team=${selectedTeamId}&watcher=${selectedWatcher?.id}`)
    }

    let selectedFilesIds: string[] = [];
    let selectedTeamId: string = $page.params.team;

    let selectedWatcher: ProductionWatcher | undefined = data.defaultWatcher;

    $: watcherSelection = data.registeredWatchers;
    $: selectedWatcher = watcherSelection.find(watcher => watcher.id === selectedWatcher?.id) || watcherSelection[0];

    let selectedTower: ProductionTower | undefined = data.registeredTowers[0];

    $: towerSelection = data.registeredTowers;
    $: selectedTower = towerSelection.find(tower => tower.id === selectedTower?.id) || towerSelection[0];

    let selectedPipeline: Pipeline | undefined = data.registeredTowers[0]?.pipelines[0];

    $: pipelineSelection = selectedTower?.pipelines;
    $: selectedPipeline= pipelineSelection?.find(pipeline => pipeline === selectedPipeline);

    let selectedDatabase: TeamDatabase = data.defaultDatabase;
    
    $: databaseSelection = data.selectedTeam.databases;
    $: selectedDatabase = databaseSelection.find(db => db.id === db.id) || databaseSelection[0];

    let selectedProject: ProjectCollection = data.defaultProject;

    $: projectSelection = selectedDatabase.projects;
    $: selectedProject = projectSelection.find(project => project.id === selectedProject.id) || projectSelection[0];

    $: towerIsActive = isWithinTimeLimit(
        data.registeredTowers.find((tower) => tower.id === selectedTower?.id)?.last_ping, 5
    );
    $: watcherIsActive = isWithinTimeLimit(
        data.registeredWatchers.find((watcher) => watcher.id === selectedWatcher?.id)?.last_ping, 5
    );    


    const stageLibraries = async() => {

        if (!selectedTower || !selectedPipeline || !selectedFilesIds) {
            return
        }

        let registerStagedSampleSchema: RegisterStagedSampleSchema = {
            tower_id: selectedTower.id,
            pipeline: selectedPipeline,
            file_ids: selectedFilesIds,
            run_id: null
        }
        await publicApi.fetchWithRefresh(
            `${publicApi.routes.stage.registerSamples}?team=${$page.params.team}&db=${selectedDatabase.id}&project=${selectedProject.id}`,
            { 
                method: 'POST',  
                mode: 'cors',
                credentials: 'include', 
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify(registerStagedSampleSchema) 
            } as RequestInit,
            $page.data.refreshToken, toastStore, "Staged libraries for pipeline execution"
        )
    }
</script>


<div class="grid sm:grid-cols-1 md:grid-cols-4 gap-16 pt-10">
    <div class="col-span-1">
        <div class="grid grid-rows-2 md:grid-rows-2 gap-y-6">
            <div class="row-span-1">
                <p class="opacity-60 mb-2">Data selection</p>
                <div class="card border border-primary-500 !bg-transparent p-6 gap-4">
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Team</span></p>
                            <select id="teamSelect" class="select" bind:value={selectedTeamId} on:change={changeTeam}>
                                {#each $page.data.userTeams as team}
                                    <option value={team.id}>{team.name}</option>
                                {/each}
                            </select>
                        </div>
                    </div>
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Watcher</span></p>
                            <select id="watcherSelect" class="select" bind:value={selectedWatcher} on:change={changeWatcher} placeholder="No watchers registered">
                                {#each watcherSelection as watcher}
                                    <option value={watcher}>{watcher.name} @ {watcher.location}</option>
                                {/each}
                            </select>
                            <div class="my-2 ml-1 inline-flex items-center">
                                <ActiveIndicator size=0.6 divClass="mt-0.5" active={watcherIsActive}></ActiveIndicator>
                                <span class="text-sm opacity-60 ml-1.5">{selectedWatcher ? getDateTimeStringUtc(selectedWatcher.last_ping, true, true) : "No watcher selected"}</span>
                            </div>
                        </div>
                    </div><div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Filters</span></p>
                            
                        </div>
                    </div>
                </div>
            </div>
            <div class="row-span-1">
                <p class="opacity-60 mb-2">Pipeline selection</p>
                <div class="card border border-primary-500 !bg-transparent p-6 gap-4">
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Tower</span></p>
                            <select id="towerSelect" class="select" bind:value={selectedTower} placeholder="No towers registered">
                                {#each towerSelection as tower}
                                    <option value={tower}>{tower.name} @ {tower.location}</option>
                                {/each}
                            </select>
                            <div class="my-2 ml-1 inline-flex items-center">
                                <ActiveIndicator size=0.6 divClass="mt-0.5" active={towerIsActive}></ActiveIndicator>
                                <span class="text-sm opacity-60 ml-1.5">{selectedTower ? getDateTimeStringUtc(selectedTower.last_ping, true, true) : "No towers selected"}</span>
                            </div>
                        </div>
                    </div>
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Pipeline</span></p>
                            <select id="pipelineSelect" class="select" bind:value={selectedPipeline}>
                                {#each pipelineSelection ?? [] as pipeline}
                                    <option value={pipeline}>{pipeline}</option>
                                {/each}
                            </select>
                        </div>
                    </div>
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Database</span></p>
                            <select id="databaseSelect" class="select" bind:value={selectedDatabase}>
                                {#each databaseSelection as db}
                                    <option value={db}>{db.name}</option>
                                {/each}
                            </select>
                        </div>
                    </div>
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Project</span></p>
                            <select id="locationSelect" class="select" bind:value={selectedProject}>
                                {#each projectSelection as project}
                                    <option value={project}>{project.name}</option>
                                {/each}
                            </select>
                        </div>
                    </div>
                </div>
                <div class="flex gap-4 justify-center mt-4">
                    <button type="button" class="btn btn-md variant-outline-primary align-center w-3/4" disabled={!towerIsActive} on:click={stageLibraries}>
                        <div class="w-4 h-4 mr-4 -mt-1.5">
                            <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-brain-circuit"><path d="M12 5a3 3 0 1 0-5.997.125 4 4 0 0 0-2.526 5.77 4 4 0 0 0 .556 6.588A4 4 0 1 0 12 18Z"/><path d="M9 13a4.5 4.5 0 0 0 3-4"/><path d="M6.003 5.125A3 3 0 0 0 6.401 6.5"/><path d="M3.477 10.896a4 4 0 0 1 .585-.396"/><path d="M6 18a4 4 0 0 1-1.967-.516"/><path d="M12 13h4"/><path d="M12 18h6a2 2 0 0 1 2 2v1"/><path d="M12 8h8"/><path d="M16 8V5a2 2 0 0 1 2-2"/><circle cx="16" cy="13" r=".5"/><circle cx="18" cy="3" r=".5"/><circle cx="20" cy="21" r=".5"/><circle cx="20" cy="8" r=".5"/></svg>
                        </div>            
                        Stage libraries  
                        <div class="w-4 h-4 mr-4 -mt-1.5">
                            <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-brain-circuit transform scale-x-[-1]"><path d="M12 5a3 3 0 1 0-5.997.125 4 4 0 0 0-2.526 5.77 4 4 0 0 0 .556 6.588A4 4 0 1 0 12 18Z"/><path d="M9 13a4.5 4.5 0 0 0 3-4"/><path d="M6.003 5.125A3 3 0 0 0 6.401 6.5"/><path d="M3.477 10.896a4 4 0 0 1 .585-.396"/><path d="M6 18a4 4 0 0 1-1.967-.516"/><path d="M12 13h4"/><path d="M12 18h6a2 2 0 0 1 2 2v1"/><path d="M12 8h8"/><path d="M16 8V5a2 2 0 0 1 2-2"/><circle cx="16" cy="13" r=".5"/><circle cx="18" cy="3" r=".5"/><circle cx="20" cy="21" r=".5"/><circle cx="20" cy="8" r=".5"/></svg>
                        </div>   
                    </button>
                </div>
            </div>
        </div>
    </div>
    <div class="col-span-1 md:col-span-3">
        {#if selectedWatcher}
            <FileSelection title={`${selectedWatcher?.name} @ ${selectedWatcher?.location}`} files={data.files} bind:selected={selectedFilesIds}></FileSelection>
        {:else}
            <div class="mt-[10%] flex justify-center items-center">
                <div class="text-center space-y-4">
                    <div class="flex justify-center">
                        <ErrorAnimation />
                    </div>
                    <h3 class="h5">
                        No watchers configured for team 
                    </h3>   
                    <p>{#if $page.error} {$page.error.message} {/if}</p>
                </div>
            </div>
        {/if}
        
    </div>
</div>
