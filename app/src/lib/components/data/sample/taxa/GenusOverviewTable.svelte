<script lang="ts">
	import type { GenusOverview } from "$lib/utils/types";
	import { ListBox, ListBoxItem, Paginator, type PaginationSettings } from "@skeletonlabs/skeleton";
  
    export let genusOverview: GenusOverview[];
    export let modelNameTags: Map<string, string[]> = new Map();

    let tableData: GenusOverview[] = [];

    const isTagData = (item: string[] | undefined): item is string[] => { return !!item }
    
    const getSelectedTags = (names: string[]) => {
        return names.map(name => modelNameTags.get(name)).filter(isTagData);
    }

    let paginationSettings: PaginationSettings = {
        page: 0,
        limit: 50,
        size: genusOverview.length,
        amounts: [5, 10, 20, 50, 100, 500],
    };

    $: {
        paginationSettings.size = genusOverview.length;
		tableData = genusOverview.slice(
			paginationSettings.page * paginationSettings.limit,
			paginationSettings.page * paginationSettings.limit + paginationSettings.limit
		);
    }

    let selectedTaxid: string;

</script>

<div>
    <div>
        <ListBox>
            <ListBoxItem group="header" name="header" value="qc" active='variant-soft' hover='hover:cursor-default' rounded='rounded-token'>
                <div class="grid grid-cols-9 sm:grid-cols-9 md:grid-cols-9 gap-x-1 gap-y-4 w-full text-sm opacity-60">
                    <div class="col-span-1">Domain</div>
                    <div class="col-span-2">Genus</div>
                    <div class="col-span-2">Tags</div>
                    <div class="text-right">t-RPM</div>
                    <div class="text-right">Contigs</div>
                    <div class="text-right">Modules</div>
                    <div class="flex justify-end"></div>
                </div>
            </ListBoxItem>
            {#each tableData as overview}
                <ListBoxItem bind:group={selectedTaxid} name={overview.genus} value={overview.taxid} active='variant-soft' rounded='rounded-token'>
                    
                    <div class="grid grid-cols-9 sm:grid-cols-9 md:grid-cols-9 gap-x-1 gap-y-4 w-full text-sm">
                        <div class="col-span-1 opacity-70">{overview.domain}</div>
                        <div class="col-span-2 truncate italic">{overview.genus} spp.</div>
                        <div class="col-span-2 truncate">
                            <!-- {#each getSelectedTags(overview.names) as tags}
                            <span class="code bg-primary-500/30 text-primary-700 dark:bg-primary-500/20 dark:text-primary-400 mr-1">{tags.join("-")}</span>
                            {/each} -->
                        </div>
                        <!-- <div class="text-right">{overview.rpm.toFixed(1)}</div>
                        <div class="text-right">{overview.contigs.toFixed(0)}</div>
                        <div class="flex justify-end gap-x-1 items-center align-center pt-1">

                            <div class="grid grid-cols-3 sm:grid-cols-3 md:grid-cols-3 gap-x-1 w-1/3 text-sm">
                                {#if overview.alignment}
                                    <div class="rounded-full bg-primary-500 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                {#if overview.kmer}
                                    <div class="rounded-full bg-secondary-500 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                {#if overview.assembly}
                                    <div class="rounded-full bg-tertiary-500 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                            </div>
                        </div>
                        <div class="flex justify-end gap-x-1 items-center align-center">
                            <div class="grid grid-cols-3 sm:grid-cols-3 md:grid-cols-3 gap-x-1 w-1/2 text-sm">
                               Test
                            </div>
                        </div> -->
                    </div>
                </ListBoxItem>
            {/each}
        </ListBox>
    </div>
    <div class="mt-8">
        <Paginator bind:settings={paginationSettings} showFirstLastButtons={false} showPreviousNextButtons={true} class="mt-2" select="paginator-select select text-xs" regionControl="opacity-30 dark:variant-filled-surface dark:opacity-60 text-xs"/>
    </div>
</div>