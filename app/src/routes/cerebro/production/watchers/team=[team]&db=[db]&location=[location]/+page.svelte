<script lang="ts">
	import { goto } from "$app/navigation";
    import { page } from "$app/stores";
	import type { SeaweedFile, Team, TeamDatabase } from "$lib/utils/types";
	import { getDateTimeStringUtc } from "$lib/utils/helpers";
    import { TreeView, TreeViewItem, type TreeViewNode } from '@skeletonlabs/skeleton';
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    export let data;

    // let selectedView: string = "team";

    // let groupedFiles: Record<string, SeaweedFile[]> = groupSeaweedFiles(data.files, "watcher_location");
    
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
    let watcherLocationFiles: Record<string, SeaweedFile[]> = groupSeaweedFiles(data.files[selectedWatcherLocation], "run_id");

    async function changeWatcherSurface() {
        await goto(`/cerebro/production/watchers/team=${selectedTeamId}&db=${selectedDatabaseId}&location=${selectedLocation}`)
    }


    let selectedTeamId: string = $page.params.team;
    let selectedDatabaseId: string = $page.params.db;
    let selectedLocation: string = $page.params.location;

    
    let databases: Array<TeamDatabase> = $page.data.defaultTeam.databases;


    let mediums = ['movies'];
    let books: string[] = [];
    let booksChildren: TreeViewItem[] = [];
</script>


<div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-4 gap-16 pt-10 ">
    <div class="col-span-1">
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
            <div class="flex gap-4 mb-6">
                <div class="w-full">
                    <p class="mb-1"><span class="opacity-60">Database</span></p>
                    <select id="projectSelect"  class="select" bind:value={selectedDatabaseId} on:change={changeWatcherSurface}>
                        {#each databases as db}
                            <option value={db.id}>{db.name}</option>
                        {/each}
                    </select>
                </div>
            </div>
            <div class="flex gap-4 mb-4">
                <div class="w-full">
                    <p class="mb-1"><span class="opacity-60">Watcher Location</span></p>
                    <select id="projectSelect"  class="select" bind:value={selectedLocation} on:change={changeWatcherSurface}>
                        {#each Object.keys(data.files) as location }
                            <option value={location}>{location}</option>
                        {/each}
                    </select>
                </div>
            </div>
            <div class="flex gap-4">
                <div class="">
                    <p class="mb-1"><span class="opacity-60">Date Range</span></p>
                    
                </div>
            </div>
        </div>
    </div>
    <div class="col-span-1 md:col-span-3">



        <TreeView class="border p-4 border-primary-500 rounded-2xl" selection multiple>
            {#each Object.entries(watcherLocationFiles) as [run_id, files] }
                <TreeViewItem bind:group={mediums} name="medium" value="books">
                        <div class="grid grid-cols-4 sm:grid-cols-4 md:grid-cols-4 gap-4">
                            <div class="opacity-60 col-span-2">{run_id}</div> 
                            <div class="opacity-60 col-span-1">{files.length} files</div>
                         </div>
                        <svelte:fragment slot="children">
                            {#each Object.entries(groupSeaweedFiles(files, "sample_id"))  as [sample_id, files] }
                            <TreeViewItem bind:this={booksChildren[0]} bind:group={books} name="books" value="Clean Code">
                                <div class="p-4">
                                <p class="opacity-60">{sample_id}</p>
                                </div>
                                <svelte:fragment slot="children">
                                {#each files as file}
                                        <TreeViewItem>
                                            <div class="grid grid-cols-7 sm:grid-cols-7 md:grid-cols-7 gap-4">
                                                <div class="col-span-1">
                                                    {getDateTimeStringUtc(file.date)}
                                                </div>
                                                <div class="col-span-1 text-center">
                                                    {file.id.substring(0, 8)}
                                                </div>
                                                <div class="col-span-1">
                                                    {file.watcher.name}
                                                </div>
                                                <div class="col-span-1">
                                                    {file.name}
                                                </div>
        
                                                <div class="col-span-1">
                                                    <div class="flex justify-end mr-2">
                                                        <input class="checkbox" type="checkbox" />
                                                    </div>
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
<!--         
        <ListBox class="border p-4 border-primary-500">

            {#each Object.entries(watcherLocationFiles) as [run_id, files] }
                <div class="flex justify-between align-middle items-center">
                    <p class="opacity-60">{run_id}</p>
                    <div>
                        <div class="flex justify-end mr-10">
                            <input class="checkbox" type="checkbox" />
                        </div>
                    </div>
                </div>
                {#each Object.entries(groupSeaweedFiles(files, "sample_id"))  as [sample_id, files] }
                    <div class="p-4">
                        <p class="opacity-60">{sample_id}</p>
                        {#each files as file}
                                <ListBoxItem bind:group={selectedView} name="medium" value="{file.id}" active='variant-soft' rounded='rounded-token'>

                                    
                                </ListBoxItem>
                        {/each}
                    </div>
                {/each}
            {/each}

            
        </ListBox> -->
    </div>
</div>
