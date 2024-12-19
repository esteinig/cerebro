<script lang="ts">
	import {  RecursiveTreeView, type TreeViewNode } from "@skeletonlabs/skeleton";
	import {  type SeaweedFile, SeaweedFileType, WatcherFileGrouping } from "$lib/utils/types";
	import FileSelectionRun from "./FileSelectionRun.svelte";
	import FileSelectionSample from "./FileSelectionSample.svelte";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import { page } from "$app/stores";

    export let title: string = "Sequence Data";
    export let files: SeaweedFile[];
    export let selected: string[] = [];

    export let firstLevel: WatcherFileGrouping = WatcherFileGrouping.RunId;
    export let secondLevel: WatcherFileGrouping = WatcherFileGrouping.SampleId;

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
     function groupSeaweedFiles(files: SeaweedFile[] | undefined, group: WatcherFileGrouping, dateRange?: { start: string, end: string }): Record<string, SeaweedFile[]> {
        
        const grouped: Record<string, SeaweedFile[]> = {};

        if (files === undefined) {
            return grouped 
        }

        files.forEach(file => {
            let key: string | null = null;

            switch (group) {
                case WatcherFileGrouping.WatcherLocation:
                    key = file.watcher? file.watcher.location : null;
                    break;
                case WatcherFileGrouping.WatcherName:
                    key = file.watcher? file.watcher.name : null;
                    break;
                case WatcherFileGrouping.RunId:
                    // Handle nullable run_id by setting a default key or excluding
                    key = file.run_id ? file.run_id : 'Missing Run ID'; // Or `null` to exclude
                    break;
                case WatcherFileGrouping.SampleId:
                    // Handle nullable run_id by setting a default key or excluding
                    key = file.sample_id ? file.sample_id : 'Missing Sample ID'; // Or `null` to exclude
                    break;
                case WatcherFileGrouping.Date:
                    if (dateRange && file.date >= dateRange.start && file.date <= dateRange.end) {
                        key = file.date;
                    } else if (!dateRange) {
                        key = file.date;
                    } // Files outside the specified range are excluded
                    break;
                default:
                    throw new Error(`Unsupported file grouping: ${group}`);
            }

            if (key) {
                grouped[key] = grouped[key] || [];
                grouped[key].push(file);
            }
        });

        return grouped;
    }

    let groupedFilesFirstLevel: Record<string, SeaweedFile[]> = {};

    $: groupedFilesFirstLevel = groupSeaweedFiles(files, firstLevel); 

    let numberRuns =  new Set(files.map(file => file.run_id)).size;

    let showInvalidSamples: boolean = false;
    
    $: treeView = Object.entries(groupedFilesFirstLevel).map(([runId, groupedFiles1]) => {
        const children = Object.entries(groupSeaweedFiles(groupedFiles1, secondLevel))
            .filter(([_, groupedFiles2]) => {

                if (showInvalidSamples) {
                    return true
                } else {
                    // Filter out unpaired nodes
                    const unpaired = groupedFiles2.filter(file => file.ftype === SeaweedFileType.ReadsPaired).length % 2 !== 0;
                    return !unpaired;
                }
            })
            .map(([sampleId, groupedFiles2]) => {
                return {
                    id: sampleId,
                    content: FileSelectionSample,
                    contentProps: {
                        id: sampleId,
                        files: groupedFiles2
                    }
                } as TreeViewNode;
            });

        return {
            id: runId,
            content: FileSelectionRun,
            contentProps: {
                id: runId,
                files: groupedFiles1
            },
            children: children
        } as TreeViewNode;
    });

    let checkedNodes: string[] = [];
    let indeterminateNodes: string[] = [];

    $: selected = checkedNodes.reduce((acc: string[], node) => {
        const matchingFiles = files.filter(file => file.sample_id === node);
        
        const unpaired = matchingFiles.filter(file => file.ftype === SeaweedFileType.ReadsPaired).length % 2 !== 0;
        if (!unpaired){
            matchingFiles.forEach(file => acc.push(file.id));
        }
        return acc;
    }, []);


    $: numberSelectedLibraries = checkedNodes.filter(node => !Object.keys(groupedFilesFirstLevel).includes(node)).length
        

</script>

<p class="opacity-60 mb-2">{title} 
    <span class="ml-6 text-xs opacity-40">{numberRuns} {numberRuns == 1 ? "run": "runs"} detected</span>
    
    <span class="ml-6 text-xs opacity-40">{numberSelectedLibraries} {numberSelectedLibraries == 1 ? "library": "libraries"} selected</span>
</p>
<div class="card border border-primary-500 !bg-transparent p-6">
    {#if numberRuns == 0}
        <div class="mt-[3%] flex justify-center items-center">
            <div class="text-center space-y-4">
                <div class="flex justify-center">
                    <ErrorAnimation />
                </div>
                <h3 class="h5">
                    No sequence runs detected for this watcher
                </h3>   
                <p>{#if $page.error} {$page.error.message} {/if}</p>
            </div>
        </div>
    {:else}

        <RecursiveTreeView 
            nodes={treeView} 
            selection 
            multiple
            relational
            bind:checkedNodes={checkedNodes} 
            bind:indeterminateNodes={indeterminateNodes}
        ></RecursiveTreeView>
    {/if}
</div>