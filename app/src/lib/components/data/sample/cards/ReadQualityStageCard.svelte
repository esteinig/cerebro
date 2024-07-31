<script lang="ts">
	import { calculateTotalQuality, formatAsPercentage, formatAsThousands } from "$lib/utils/helpers";
	import type { QualityControlSummary, WorkflowParamsQcFastp } from "$lib/utils/types";
	import { Tab, TabGroup } from "@skeletonlabs/skeleton";
    
    export let selectedQualityControlSummary: QualityControlSummary;
    export let selectedStageParams: WorkflowParamsQcFastp | undefined;
    export let selectedStage: string = "Read Quality";

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
                                                Removal of low-quality reads
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

                                    <div><span class="opacity-60 mr-4">Total removed</span></div>
                                    <div>{formatAsThousands(calculateTotalQuality(
                                        selectedQualityControlSummary.low_complexity_reads,
                                        selectedQualityControlSummary.qc_reads
                                    ))}</div>

                                    <div><span class="opacity-60 mr-4 ml-4">Low complexity</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.low_complexity_reads)}</div>

                                    <div><span class="opacity-60 mr-4 ml-4">Low quality</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.qc_reads)}</div>

                                    <div><span class="opacity-60 mr-4 ml-8">Minimum length</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.qc_min_length_reads)}</div>

                                    <div><span class="opacity-60 mr-4 ml-8">Low read quality</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.qc_low_quality_reads)}</div>

                                    <div><span class="opacity-60 mr-4 ml-8">Missing bases</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.qc_missing_bases_reads)}</div>

                                    <div><span class="opacity-60 mr-4">{">"} Q20</span></div>
                                    <div>{formatAsPercentage(selectedQualityControlSummary.q20_percent)}</div>

                                    <div><span class="opacity-60 mr-4">{">"} Q30</span></div>
                                    <div>{formatAsPercentage(selectedQualityControlSummary.q30_percent)}</div>

                                    <div><span class="opacity-60 mr-4">Mean length</span></div>
                                    <div>{selectedQualityControlSummary.mean_length_r1} / {selectedQualityControlSummary.mean_length_r2}</div>


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
                            <p class="text-sm opacity-60">
                                Read quality control is conducted with 
                                <span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90">Fastp v.</span> .
                            </p>
                            <p class="text-sm opacity-60 my-2">
                                This step includes various optional configurations such as length and adapter trimming, as well as estimation of 
                                remaining duplication levels.
                            </p>
                            <p class="text-sm opacity-60 my-2">
                                <span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90">Fastp v.</span> 
                                is never used for deduplication due to a high probability of hash collisions and resulting non-deterministic outcomes.
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
                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">min_read_length</span></div>
                                    <div>{selectedStageParams?.min_read_length}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">cut_tail_quality</span></div>
                                    <div>{selectedStageParams?.cut_tail_quality}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">complexity_threshold</span></div>
                                    <div>{selectedStageParams?.complexity_threshold}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">adapter_auto_detect</span></div>
                                    <div>{selectedStageParams?.adapter_auto_detect}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">adapter_file</span></div>
                                    <div>{selectedStageParams?.adapter_file}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">adapter_seq_1</span></div>
                                    <div class="text-sm">{selectedStageParams?.adapter_seq_1}</div>

                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">adapter_seq_2</span></div>
                                    <div class="text-sm">{selectedStageParams?.adapter_seq_2}</div>
                                    
                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">trim_poly_g</span></div>
                                    <div class="text-base">{selectedStageParams?.trim_poly_g}</div>
                                </div>
                            </div>
                            <div class="p-1">
                                <p class="text-sm opacity-80 pb-5">
                                    Workflow command
                                </p>
                                
                            </div>
                        </section>
                    </div>
                {/if}

        </svelte:fragment>
    </TabGroup>
</div>
