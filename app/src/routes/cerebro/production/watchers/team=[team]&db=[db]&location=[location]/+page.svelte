<script lang="ts">
	import { goto } from "$app/navigation";
    import { page } from "$app/stores";
    import { onMount } from 'svelte';
    import { invalidate } from '$app/navigation';

	import type { ProductionPipeline, ProjectCollection, SeaweedFile, TeamDatabase } from "$lib/utils/types";
	import { getDateTimeStringUtc, isWithinTimeLimit } from "$lib/utils/helpers";
    import { TreeView, TreeViewItem, type TreeViewNode } from '@skeletonlabs/skeleton';
	import PipelineIndicator from "$lib/general/icons/PipelineIndicator.svelte";

    export let data;

    // Invalidates page data every five minutes to update 
    // the pipeline activity status and indicator
    onMount(() => {
        const interval = setInterval(() => {
            invalidate('watchers:data');
        }, 5 * 60 * 1000);

        return () => clearInterval(interval);
    });
    
    /**
     * Group SeaweedFile objects by a specified key (watcher_location, run_id, or date).
     * Adjusted to handle nullable run_id by excluding those files from run_id grouping
     * or placing them under a specific key.
     *
     * @param files Array of SeaweedFile objects to be grouped.
     * @param groupByKey The key by which the files should be grouped ('watcher_location', 'run_id', or 'date').
     * @param dateRange (Optional) Specifies the date range for grouping if 'date' is chosen as groupByKey. 
     *                  Format: { start: 'YYYY-MM-DD', end: 'YYYY-MM-DD' }
     * @returns An object with keys as the group by value and values as arrays of SeaweedFiles belonging to that group.
     *
     * @example
     * // Example usage for grouping by watcher location
     * const files: SeaweedFile[] = [...];
     * const groupedByLocation = groupSeaweedFiles(files, 'watcher_location');
     *
     * // Example usage for grouping by run_id, excluding null run_ids
     * const groupedByRunId = groupSeaweedFiles(files, 'run_id');
     *
     * // Example usage for grouping by date with a specified range
     * const groupedByDateRange = groupSeaweedFiles(files, 'date', { start: '2022-01-01', end: '2022-01-31' });
     */
    function groupSeaweedFiles(files: SeaweedFile[] | undefined, groupByKey: 'watcher_location' | 'watcher_name' | 'run_id' | 'sample_id' | 'date', dateRange?: { start: string, end: string }): Record<string, SeaweedFile[]> {
        
        const grouped: Record<string, SeaweedFile[]> = {};

        if (files === undefined) {
            return grouped 
        }

        files.forEach(file => {
            let key: string | null = null;

            switch (groupByKey) {
                case 'watcher_location':
                    key = file.watcher.location;
                    break;
                case 'watcher_name':
                    key = file.watcher.name;
                    break;
                case 'run_id':
                    // Handle nullable run_id by setting a default key or excluding
                    key = file.run_id ? file.run_id : 'no_run_id'; // Or `null` to exclude
                    break;
                case 'sample_id':
                    // Handle nullable run_id by setting a default key or excluding
                    key = file.sample_id ? file.sample_id : 'no_sample_id'; // Or `null` to exclude
                    break;
                case 'date':
                    if (dateRange && file.date >= dateRange.start && file.date <= dateRange.end) {
                        key = file.date;
                    } else if (!dateRange) {
                        key = file.date;
                    } // Files outside the specified range are excluded
                    break;
                default:
                    throw new Error(`Unsupported groupByKey: ${groupByKey}`);
            }

            if (key) {
                grouped[key] = grouped[key] || [];
                grouped[key].push(file);
            }
        });

        return grouped;
    }

    let selectedWatcherLocation = data.defaultWatcherLocation;
    let selectedWatcherLocationFiles: Record<string, SeaweedFile[]> = groupSeaweedFiles(data.files[selectedWatcherLocation], "run_id");  // files can be undefined

    async function changeWatcherSurface() {
        await goto(`/cerebro/production/watchers/team=${selectedTeamId}&db=${selectedDatabaseId}&location=${selectedLocation}`)
    }

    let selectedTeamId: string = $page.params.team;
    let selectedDatabaseId: string = $page.params.db;
    let selectedProjectId: string = data.defaultProject.id;
    let selectedLocation: string = data.defaultWatcherLocation;

    let databases: Array<TeamDatabase> = data.defaultTeam.databases;
    let projects: Array<ProjectCollection> = data.defaultDatabase.projects;

    let runTree = ['runs'];
    let sampleTree = ['sample'];

    let controlTag: string = "";
    let controlTags: string[] = ["TPL", "ENV", "NEG", "POS"];
    
    let nucleicAcidTag: string = "";
    let nucleicAcidTags: string[] = ["DNA", "RNA"]

    let specimenTag: string = "";
    let specimenTags: string[] = ["CSF", "BLD", "RES"]
    
    let pipelineSelection: ProductionPipeline[] = data.registeredPipelines;
    let selectedPipeline: ProductionPipeline = pipelineSelection[0];
    
    $: pipelineIsActive = isWithinTimeLimit(selectedPipeline.last_ping, 30);

</script>


<div class="grid sm:grid-cols-1 md:grid-cols-5 gap-16 pt-10">
    <div class="col-span-1">
        <div class="grid grid-rows-2 md:grid-rows-2 gap-y-6">
            <div class="row-span-1">
                <p class="opacity-60 mb-2">Sentinel selection</p>
                <div class="border p-4 border-primary-500 gap-4 rounded-2xl">
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Team</span></p>
                            <select id="teamSelect" class="select" bind:value={selectedTeamId} on:change={changeWatcherSurface}>
                                {#each $page.data.userTeams as team}
                                    <option value={team.id}>{team.name}</option>
                                {/each}
                            </select>
                        </div>
                    </div>
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Database</span></p>
                            <select id="databaseSelect" class="select" bind:value={selectedDatabaseId} on:change={changeWatcherSurface}>
                                {#each databases as db}
                                    <option value={db.id}>{db.name}</option>
                                {/each}
                            </select>
                        </div>
                    </div>
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Sentinel</span></p>
                            <select id="watcherSelect" class="select" bind:value={selectedLocation} on:change={changeWatcherSurface}>
                                {#each Object.keys(data.files) as location }
                                    <option value={location}>{location}</option>
                                {/each}
                            </select>
                        </div>
                    </div>
                </div>
            </div>
            <div class="row-span-1">
                <p class="opacity-60 mb-2">Cerebro configuration</p>
                <div class="border p-4 border-primary-500 gap-4 rounded-2xl">
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Pipeline</span></p>
                            <select id="pipelineSelect" class="select" bind:value={selectedPipeline}>
                                {#each pipelineSelection as pipeline}
                                    <option value={pipeline}>{pipeline.name} @ {pipeline.location}</option>
                                {/each}
                            </select>
                            <div class="my-2 ml-1 inline-flex items-center">
                                <PipelineIndicator size=0.6 divClass="mt-0.5" active={pipelineIsActive}></PipelineIndicator>
                                <span class="text-sm opacity-60 ml-1.5">{selectedPipeline.pipeline} Workflow</span>
                            </div>
                            
                        </div>
                    </div>
                    <div class="flex gap-4 mb-4">
                        <div class="w-full">
                            <p class="mb-1"><span class="opacity-60">Database Project</span></p>
                            <select id="locationSelect" class="select" bind:value={selectedProjectId}>
                                {#each projects as project}
                                    <option value={project.id}>{project.name}</option>
                                {/each}
                            </select>
                        </div>
                    </div>
                </div>
                <div class="flex gap-4 justify-center mt-4">
                    <button type="button" class="btn btn-md variant-outline-primary align-center w-3/4" disabled={!pipelineIsActive}>
                        <div class="w-4 h-4 mr-4 -mt-1.5">
                            <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-brain-circuit"><path d="M12 5a3 3 0 1 0-5.997.125 4 4 0 0 0-2.526 5.77 4 4 0 0 0 .556 6.588A4 4 0 1 0 12 18Z"/><path d="M9 13a4.5 4.5 0 0 0 3-4"/><path d="M6.003 5.125A3 3 0 0 0 6.401 6.5"/><path d="M3.477 10.896a4 4 0 0 1 .585-.396"/><path d="M6 18a4 4 0 0 1-1.967-.516"/><path d="M12 13h4"/><path d="M12 18h6a2 2 0 0 1 2 2v1"/><path d="M12 8h8"/><path d="M16 8V5a2 2 0 0 1 2-2"/><circle cx="16" cy="13" r=".5"/><circle cx="18" cy="3" r=".5"/><circle cx="20" cy="21" r=".5"/><circle cx="20" cy="8" r=".5"/></svg>
                        </div>            
                        Stage files for Cerebro  
                        <div class="w-4 h-4 mr-4 -mt-1.5">
                            <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-brain-circuit transform scale-x-[-1]"><path d="M12 5a3 3 0 1 0-5.997.125 4 4 0 0 0-2.526 5.77 4 4 0 0 0 .556 6.588A4 4 0 1 0 12 18Z"/><path d="M9 13a4.5 4.5 0 0 0 3-4"/><path d="M6.003 5.125A3 3 0 0 0 6.401 6.5"/><path d="M3.477 10.896a4 4 0 0 1 .585-.396"/><path d="M6 18a4 4 0 0 1-1.967-.516"/><path d="M12 13h4"/><path d="M12 18h6a2 2 0 0 1 2 2v1"/><path d="M12 8h8"/><path d="M16 8V5a2 2 0 0 1 2-2"/><circle cx="16" cy="13" r=".5"/><circle cx="18" cy="3" r=".5"/><circle cx="20" cy="21" r=".5"/><circle cx="20" cy="8" r=".5"/></svg>
                        </div>   
                    </button>
                </div>
            </div>
        </div>
    </div>
    <div class="col-span-1 md:col-span-3">

        <p class="opacity-60 mb-2">Sequence runs @ {selectedLocation}</p>
        <TreeView class="border p-4 border-primary-500 rounded-2xl" selection multiple>
            {#each Object.entries(selectedWatcherLocationFiles) as [run_id, files] }
                <TreeViewItem bind:group={runTree} name="runs" value="runs">
                        <div class="grid grid-cols-7 sm:grid-cols-7 md:grid-cols-7 gap-4">
                            <div class="opacity-60 col-span-2">{run_id}</div> 
                            <div class="opacity-60 col-span-1">{files.length} files</div>
                         </div>
                        <svelte:fragment slot="children">
                            {#each Object.entries(groupSeaweedFiles(files, "sample_id"))  as [sample_id, files] }
                            <TreeViewItem bind:group={sampleTree} name="samples" value="sample">

                                <div class="grid grid-cols-7 sm:grid-cols-7 md:grid-cols-7 gap-4 items-center">
                                    <div class="p-4">
                                        <p class="opacity-60">{sample_id}</p>
                                    </div>
                                    <div class="col-span-1 flex justify-center">
                                        {#each nucleicAcidTags as tag}
                                            <button
                                                class="chip {nucleicAcidTag === tag ? 'variant-filled-primary' : 'variant-soft'} ml-2"
                                                on:click={() => { nucleicAcidTag = tag }}
                                                on:keypress
                                            >
                                                <span>{tag}</span>
                                            </button>
                                        {/each}
                                    </div>
                                    <div class="col-span-1 flex justify-center">
                                        {#each specimenTags as tag}
                                            <button
                                                class="chip {specimenTag === tag ? 'variant-filled-primary' : 'variant-soft'} ml-2"
                                                on:click={() => { specimenTag = tag }}
                                                on:keypress
                                            >
                                                <span>{tag}</span>
                                            </button>
                                        {/each}
                                    </div>
                                    <div class="col-span-2 flex justify-center">
                                        {#each controlTags as tag}
                                            <button
                                                class="chip {controlTag === tag ? 'variant-filled-primary' : 'variant-soft'} ml-2"
                                                on:click={() => { controlTag = tag }}
                                                on:keypress
                                            >
                                                <span>{tag}</span>
                                            </button>
                                        {/each}
                                    </div>
                                </div>
                                <svelte:fragment slot="children">
                                {#each files as file}
                                        <TreeViewItem>
                                            <div class="grid grid-cols-7 sm:grid-cols-7 md:grid-cols-7 gap-4 items-center">

                                                <div class="col-span-1 text-center">
                                                    {file.watcher.name}
                                                </div>
                                                <div class="col-span-1 text-center">
                                                    {getDateTimeStringUtc(file.date, false)}
                                                </div>
                                                <div class="col-span-1">
                                                    {file.name}
                                                </div>
                                                <div class="col-span-2 text-center">
                                                    {file.id}
                                                </div>
                                                <div class="col-span-1 text-center">
                                                    {file.fid}
                                                </div>
                                                
                                            </div>
                                        </TreeViewItem>
                                {/each}
                            </svelte:fragment>
                            </TreeViewItem>
                        {/each}
                        </svelte:fragment>
                </TreeViewItem>
            {/each}
        </TreeView>
    </div>
</div>
