<script lang="ts">
	import { formatAsThousands } from "$lib/utils/helpers";
	import type { QualityControlSummary } from "$lib/utils/types";
	import { Tab, TabGroup } from "@skeletonlabs/skeleton";

    export let selectedQualityControlSummary: QualityControlSummary;
    export let selectedStage: string = "Input";

    let tabSet: number = 0;

</script>

<div class="">
    <TabGroup>
        <Tab bind:group={tabSet} name="overview" value={0}>
            <span class="text-sm">Overview</span>
        </Tab>
        <Tab bind:group={tabSet} name="description" value={1}>
            <span class="text-sm">Description</span>
        </Tab>
        <Tab bind:group={tabSet} name="configuration" value={2}>
            <span class="text-sm">Configuration</span>
        </Tab>
        <svelte:fragment slot="panel">
            {#if tabSet === 0}
                
                <div class="card max-w-xl2">
                    <header class="card-header">
                        <div class="p-1">
                            <div class="flex">
                                <div class="flex items-center">
                                    <div class="flex-none">
                                        <div class="md:h4 lg:h3">{selectedStage}</div>
                                        <div class="flex items-center">
                                            <p class="flex opacity-60 text-xs pr-1">
                                                    Unprocessed paired-end reads
                                            </p>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </header>
                    <section class="p-4">
                        <div class="p-1">
                            
                                <div class="grid grid-cols-2 sm:grid-cols-2 md:grid-cols-2 gap-x-8 gap-y-4 w-3/4">
                                    <div><span class="opacity-60 mr-4">Total reads</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.total_reads)}</div>
                                    <div><span class="opacity-60 mr-4">Total biomass</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.total_biomass)} pg</div>
                                </div>
                        </div>
                    </section>
                    <footer class="card-footer">
                
                        
                    </footer>
                </div>
            {/if}
        </svelte:fragment>
    </TabGroup>
</div>
