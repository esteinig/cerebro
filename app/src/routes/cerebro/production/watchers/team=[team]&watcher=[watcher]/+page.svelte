<script lang="ts">
	import { goto } from "$app/navigation";
    import { page } from "$app/stores";
    import { onMount } from 'svelte';
    import { invalidate } from '$app/navigation';

	import { type ProductionPipeline, type ProductionWatcher, type ProjectCollection, type TeamDatabase } from "$lib/utils/types";
	import { isWithinTimeLimit } from "$lib/utils/helpers";
	import ActiveIndicator from "$lib/general/icons/ActiveIndicator.svelte";
	import FileSelection from "$lib/components/production/watchers/FileSelection.svelte";

    export let data;

    // Invalidates page data every five minutes to update 
    // the watcher/pipeline activity status indicators
    onMount(() => {
        const interval = setInterval(() => {
            invalidate('watchers:data');
        }, 5 * 60 * 1000);

        return () => clearInterval(interval);
    });
    
    async function changeTeam() {
        await goto(`/cerebro/production/watchers/team=${selectedTeamId}&watcher=0`)
    }

    async function changeWatcher() {
        await goto(`/cerebro/production/watchers/team=${selectedTeamId}&watcher=${selectedWatcher?.id}`)
    }


    let selectedWatcher: ProductionWatcher | undefined = data.defaultWatcher;

    $: watcherSelection = data.registeredWatchers;
    $: selectedWatcher = watcherSelection.find(watcher => watcher.id === selectedWatcher?.id) || watcherSelection[0];

    let selectedPipeline: ProductionPipeline | undefined = data.registeredPipelines[0];

    $: pipelineSelection = data.registeredPipelines;
    $: selectedPipeline = pipelineSelection.find(pipeline => pipeline.id === selectedPipeline?.id) || pipelineSelection[0];

    $: pipelineIsActive = isWithinTimeLimit(
        data.registeredPipelines.find((pipeline) => pipeline.id === selectedPipeline?.id)?.last_ping, 5
    );
    $: watcherIsActive = isWithinTimeLimit(
        data.registeredWatchers.find((watcher) => watcher.id === selectedWatcher?.id)?.last_ping, 5
    );    

    let selectedTeamId: string = $page.params.team;
    let selectedDatabaseId: string = $page.params.db;
    let selectedProjectId: string = data.defaultProject.id;

    let databases: Array<TeamDatabase> = data.defaultTeam.databases;
    let projects: Array<ProjectCollection> = data.defaultDatabase.projects;



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
                            <select id="watcherSelect" class="select" bind:value={selectedWatcher} on:change={changeWatcher}>
                                {#each watcherSelection as watcher}
                                    <option value={watcher}>{watcher.name} @ {watcher.location}</option>
                                {/each}
                            </select>
                            <div class="my-2 ml-1 inline-flex items-center">
                                <ActiveIndicator size=0.6 divClass="mt-0.5" active={watcherIsActive}></ActiveIndicator>
                                <span class="text-sm opacity-60 ml-1.5">{selectedWatcher ? selectedWatcher.format : "No watcher selected"}</span>
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
                            <p class="mb-1"><span class="opacity-60">Pipeline</span></p>
                            <select id="pipelineSelect" class="select" bind:value={selectedPipeline}>
                                {#each pipelineSelection as pipeline}
                                    <option value={pipeline}>{pipeline.name} @ {pipeline.location}</option>
                                {/each}
                            </select>
                            <div class="my-2 ml-1 inline-flex items-center">
                                <ActiveIndicator size=0.6 divClass="mt-0.5" active={pipelineIsActive}></ActiveIndicator>
                                <span class="text-sm opacity-60 ml-1.5">{selectedPipeline ? selectedPipeline.pipeline : "No pipeline selected"}</span>
                            </div>
                        </div>
                    </div>
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Database</span></p>
                            <select id="databaseSelect" class="select" bind:value={selectedDatabaseId}>
                                {#each databases as db}
                                    <option value={db.id}>{db.name}</option>
                                {/each}
                            </select>
                        </div>
                    </div>
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Project</span></p>
                            <select id="locationSelect" class="select" bind:value={selectedProjectId}>
                                {#each projects as project}
                                    <option value={project.id}>{project.name}</option>
                                {/each}
                            </select>
                        </div>
                    </div>
                </div>
                <div class="flex gap-4 justify-center mt-4">
                    <button type="button" class="btn btn-md variant-outline-primary align-center w-1/2" disabled={!pipelineIsActive}>
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
        <FileSelection title={`${selectedWatcher?.name} @ ${selectedWatcher?.location}`} files={data.files}></FileSelection>
    </div>
</div>
