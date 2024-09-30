<script lang="ts">
	import { formatAsThousands } from "$lib/utils/helpers";
	import type { QualityControlSummary, WorkflowParamsQcHostDepletion } from "$lib/utils/types";
	import { Tab, TabGroup } from "@skeletonlabs/skeleton";

    export let selectedQualityControlSummary: QualityControlSummary;
    export let selectedStageParams: WorkflowParamsQcHostDepletion | undefined;
    export let selectedStage: string = "Host Depletion";

    let tabSet: number = 0;

</script>

<div class="">
    <TabGroup>
        <Tab bind:group={tabSet} name="overview" value={0}>
            <span class="text-sm">Results</span>
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
                                                Processed paired-end reads
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
                                    <div><span class="opacity-60 mr-4">Host reads</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.host_reads)}</div>
                                    <div><span class="opacity-60 mr-4">Host biomass</span></div>
                                    <div>{selectedQualityControlSummary.host_biomass} pg</div>
                                    <div><span class="opacity-60 mr-4">Other reads</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.other_reads)}</div>
                                    <div><span class="opacity-60 mr-4">Other biomass</span></div>
                                    <div>{selectedQualityControlSummary.other_biomass} pg</div>
                                </div>
                        </div>
                    </section>
                    <footer class="card-footer">
                
                        
                    </footer>
                </div>
            {:else if tabSet === 1}
                <div class="card max-w-xl2">
                    <header class="card-header">
                        <div class="p-1">
                            <div class="flex">
                                <div class="flex items-center">
                                    <div class="flex-none">
                                        <div class="md:h4 lg:h3">{selectedStage}</div>
                                        <div class="flex items-center">
                                            <p class="flex opacity-60 text-xs pr-1">
                                                Background information
                                            </p>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </header>
                    <section class="p-4">
                    </section>
                </div>
                {:else if tabSet === 2}
                    <div class="card max-w-xl2">
                        <header class="card-header">
                            <div class="p-1">
                                <div class="flex">
                                    <div class="flex items-center">
                                        <div class="flex-none">
                                            <div class="md:h4 lg:h3">{selectedStage}</div>
                                            <div class="flex items-center">
                                                <p class="flex opacity-60 text-xs pr-1">
                                                    Parameters and commands
                                                </p>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </header>
                        <section class="p-4">
                            
                            <div class="p-1 pb-8">
                                <p class="text-sm opacity-80 pb-5">
                                    Workflow configuration
                                </p>
                                <div class="grid grid-cols-2 sm:grid-cols-2 md:grid-cols-2 gap-x-8 gap-y-4 w-3/4">
                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">databases</span></div>
                                    <div>{selectedStageParams?.databases.split(",").join(", ")}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">references</span></div>
                                    <div>{selectedStageParams?.references.split(",").join(", ")}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">taxa</span></div>
                                    <div>{selectedStageParams?.taxa.split(" ").join(", ")}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">direct</span></div>
                                    <div>{selectedStageParams?.direct.split(" ").join(", ")}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">min_cov</span></div>
                                    <div>{selectedStageParams?.min_cov}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">min_len</span></div>
                                    <div class="text-sm">{selectedStageParams?.min_len}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">min_mapq</span></div>
                                    <div class="text-sm">{selectedStageParams?.min_mapq}</div>
                                    
                                </div>
                            </div>
                            <div class="p-1">
                                <p class="text-sm opacity-80 pb-5">
                                    Workflow commands
                                </p>
                            </div>
                        </section>
                    </div>
            {/if}
        </svelte:fragment>
    </TabGroup>
</div>
