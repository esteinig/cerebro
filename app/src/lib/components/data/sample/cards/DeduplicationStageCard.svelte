<script lang="ts">
	import { formatAsThousands } from "$lib/utils/helpers";
	import type { QualityControlSummary, WorkflowParamsQcDeduplication } from "$lib/utils/types";
	import { Tab, TabGroup } from "@skeletonlabs/skeleton";

    export let selectedQualityControlSummary: QualityControlSummary;
    export let selectedStageParams: WorkflowParamsQcDeduplication;
    export let selectedStage: string = "Deduplication";

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
                                                Removal of duplicate read (PCR amplicons)
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
                                    <div><span class="opacity-60 mr-4">Duplicate reads</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.deduplicated_reads)}</div>
                                </div>
                        </div>
                    </section>
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
                        <div class="p-1 pb-8">
                            <p class="text-sm opacity-80 pb-5">
                                Description
                            </p>
                        
                            <p class="text-sm font-bold opacity-80 my-2">
                                UMI-based assays and deduplication methods must be implemented for
                                accurate biomass estimation when using ERCC/EDCC.
                            </p>
                        
                            <p class="text-sm opacity-60 my-2">
                                Deduplication is conducted as the first step in the quality control module. This is because
                                we want to avoid trimming reads, as nearly all implemented options for deduplication use
                                read similarity assessments. 
                            </p>
                              
                            <p class="text-sm opacity-60 my-2">
                                If sufficient memory is available we recommend the naive clustering method, which also does not use
                                hashing algorithms with high collision probability as when using Fastp deduplication modes. 
                            </p>
                        
                        
                            <p class="text-sm opacity-60 my-2">
                                UMI-tools is not recommended due to apparent non-deterministic outcomes despite seed value parameter, and because
                                only host and synthetic control reads can be deduplicated using alignment/position-based methods in 
                                metagenomic data.
                            </p>
                        </div>
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
                                <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">method</span></div>
                                <div>{selectedStageParams.method}</div>
                                <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">umi_tools.seed</span></div>
                                <div>{selectedStageParams.umi_tools.seed}</div>
                                <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">umi_tools.reference</span></div>
                                <div>{selectedStageParams.umi_tools.reference}</div>
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
