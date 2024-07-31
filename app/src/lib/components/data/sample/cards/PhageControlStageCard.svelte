<script lang="ts">
	import { formatAsPercentage, formatAsThousands } from "$lib/utils/helpers";
	import type { QualityControlSummary, WorkflowParamsQcPhage } from "$lib/utils/types";
	import { Tab, TabGroup } from "@skeletonlabs/skeleton";

    export let selectedStageParams: WorkflowParamsQcPhage | undefined;
    export let selectedQualityControlSummary: QualityControlSummary;
    export let selectedStage: string = "Phage Control";

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
                                                <span class="mr-1">Extraction and sequencing controls</span>
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

                                    <div><span class="opacity-60 mr-4">DNA phage</span></div>
                                    <div>{selectedQualityControlSummary.dna_phage_id}</div>
                                    <div><span class="opacity-60 mr-4">DNA phage reads</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.dna_phage_reads)}</div>
                                    <div><span class="opacity-60 mr-4">DNA phage coverage</span></div>
                                    <div>{formatAsPercentage(selectedQualityControlSummary.dna_phage_coverage_percent)}</div>
                                    <div><span class="opacity-60 mr-4">DNA phage biomass</span></div>
                                    <div>{selectedQualityControlSummary.dna_phage_biomass} pg</div>
                                </div>
                                <div class="grid grid-cols-2 sm:grid-cols-2 md:grid-cols-2 gap-x-8 gap-y-4 w-3/4 mt-12">
                                    <div><span class="opacity-60 mr-4">RNA phage</span></div>
                                    <div>{selectedQualityControlSummary.rna_phage_id}</div>
                                    <div><span class="opacity-60 mr-4">RNA phage reads</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.rna_phage_reads)}</div>
                                    <div><span class="opacity-60 mr-4">RNA phage coverage</span></div>
                                    <div>{formatAsPercentage(selectedQualityControlSummary.rna_phage_coverage_percent)}</div>
                                    <div><span class="opacity-60 mr-4">RNA phage biomass</span></div>
                                    <div>{selectedQualityControlSummary.rna_phage_biomass} pg</div>
                                </div>

                                <div class="grid grid-cols-2 sm:grid-cols-2 md:grid-cols-2 gap-x-8 gap-y-4 w-3/4 mt-12">
                                    <div><span class="opacity-60 mr-4">SEQ phage</span></div>
                                    <div>{selectedQualityControlSummary.seq_phage_id}</div>
                                    <div><span class="opacity-60 mr-4">SEQ phage reads</span></div>
                                    <div>{formatAsThousands(selectedQualityControlSummary.seq_phage_reads)}</div>
                                    <div><span class="opacity-60 mr-4">SEQ phage coverage</span></div>
                                    <div>{formatAsPercentage(selectedQualityControlSummary.seq_phage_coverage_percent)}</div>
                                    <div><span class="opacity-60 mr-4">SEQ phage biomass</span></div>
                                    <div>{selectedQualityControlSummary.seq_phage_biomass} pg</div>
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
                            <p class="text-sm opacity-60 my-2">
                                T4 phage spike-ins are used as internal extraction controls. Coverage should be 100% for 
                                high confidence in DNA extractions. More information can be found in the assay validation
                                and documentation.
                            </p>
                            <p class="text-sm opacity-60">
                                Phage controls are processed in a sub-workflow using 
                                <span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90">minimap2 v.</span> 
                                alignments against the <span class="italic">Escherichia</span> T4 phage genome. Aligned reads are depleted with 
                                <span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90">Scrubby v.</span> 
                                 and a summary of alignments and coverage is produced with 
                                <span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90">Vircov v.</span> 
                            </p>
                            <p class="text-sm opacity-60 my-2">
                                This step is conducted after host genome depletion to avoid costly alignment of host reads and to keep the phage alignment
                                and evaluation distinct from other steps in the quality control module. 
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
                                    <div><span class="code text-xs bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 mr-4">fasta</span></div>
                                    <div>{selectedStageParams?.fasta}</div>
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
