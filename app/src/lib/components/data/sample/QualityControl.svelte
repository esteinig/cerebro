<script lang="ts">
	import type { Cerebro, QualityControlSummary, WorkflowConfig } from "$lib/utils/types";
   	import SyntheticControlStageCard from "./cards/SyntheticControlStageCard.svelte";
	import PhageControlStageCard from "./cards/PhageControlStageCard.svelte";
	import ReadQualityStageCard from "./cards/ReadQualityStageCard.svelte";

    export let selectedWorkflowConfiguration: WorkflowConfig;
    export let selectedModels: Cerebro[] = [];
    export let selectedQualityControlSummaries: QualityControlSummary[] = [];

</script>

<div class="w-full">
    <div class="flex pb-10 opacity-90 text-sm">
        <ol class="breadcrumb justify-center">
            <li class="crumb opacity-70">Input</li>
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            {#if selectedWorkflowConfiguration?.params?.qc.controls.ercc.enabled}
                <li class="crumb">Synthetic Controls</li>
            {:else}
                <li class="crumb opacity-10">Synthetic Controls</li>
            {/if}
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            {#if selectedWorkflowConfiguration?.params?.qc.host.depletion.enabled}
                <li class="crumb">Internal Controls</li>
            {:else}
                <li class="crumb opacity-10">Phage Controls</li>
            {/if}
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            {#if selectedWorkflowConfiguration?.params?.qc.deduplication.enabled}
                <li class="crumb">Deduplication</li>
            {:else}
                <li class="crumb opacity-10">Deduplication</li>
            {/if}
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            {#if selectedWorkflowConfiguration?.params?.qc.reads.fastp.enabled}
                <li class="crumb">Read Quality</li>
            {:else}
                <li class="crumb opacity-10">Read Quality</li>
            {/if}
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            {#if selectedWorkflowConfiguration?.params?.qc.host.depletion.enabled}
                <li class="crumb">Host Depletion</li>
            {:else}
                <li class="crumb opacity-10">Host Depletion</li>
            {/if}
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            {#if selectedWorkflowConfiguration?.params?.qc.host.depletion.enabled}
                <li class="crumb">Background Depletion</li>
            {:else}
                <li class="crumb opacity-10">Background Depletion</li>
            {/if}
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            <li class="crumb opacity-70">Output</li>
        </ol>
    </div>
    {#each selectedModels as model, i}  
        <div class="grid grid-cols-1 md:grid-cols-2 grid-rows-2 gap-12 py-8">
            <!-- ReadQualityCard: First column spanning two rows -->
            <div class="p-4 row-span-2 rounded-xl bg-surface-500/5">
                <p class="flex align-center"><span class="code bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 text-sm">{model.sample.id}</span> <span class="code ml-2 text-sm">{model.sample.tags.join("-")}</span></p>
                    
                <ReadQualityStageCard selectedQualityControlSummary={selectedQualityControlSummaries[i]}></ReadQualityStageCard>
            </div>
    
            <!-- SyntheticControlCard: Second column, second row -->
            <div class="p-8 rounded-xl bg-surface-500/5">
                <SyntheticControlStageCard selectedModel={model} selectedQualityControlSummary={selectedQualityControlSummaries[i]}></SyntheticControlStageCard>
            </div>
            
            <!-- PhageControlCard: Second column, first row -->
            <div class="p-8 rounded-xl bg-surface-500/5">
                <PhageControlStageCard selectedModel={model}></PhageControlStageCard>
            </div>
    
        </div>
    {/each}
</div>
