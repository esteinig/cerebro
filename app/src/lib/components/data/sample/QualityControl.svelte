<script lang="ts">
	import type { Cerebro, WorkflowConfig } from "$lib/utils/types";
   	import SyntheticControlStageCard from "./cards/SyntheticControls.svelte";
	import PhageControlStageCard from "./cards/InternalControls.svelte";
	import ReadQualityStageCard from "./cards/ReadQuality.svelte";

    export let selectedWorkflowConfiguration: WorkflowConfig | null;
    export let selectedModels: Cerebro[] = [];

</script>

<div class="w-full">
    <div class="flex pb-10 opacity-90 text-sm">
        <ol class="breadcrumb justify-center">
            <li class="crumb opacity-70">Input</li>
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            {#if selectedWorkflowConfiguration?.params?.qc.host.depletion.enabled}
                <li class="crumb">Internal Controls</li>
            {:else}
                <li class="crumb opacity-10">Internal Controls</li>
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
                <ReadQualityStageCard selectedModel={model}></ReadQualityStageCard>
            </div>
    
            <!-- SyntheticControlCard: Second column, second row -->
            <div class="p-8 rounded-xl bg-surface-500/5">
                <SyntheticControlStageCard selectedModel={model}></SyntheticControlStageCard>
            </div>
            
            <!-- PhageControlCard: Second column, first row -->
            <div class="p-8 rounded-xl bg-surface-500/5">
                <PhageControlStageCard selectedModel={model}></PhageControlStageCard>
            </div>
    
        </div>
    {/each}
</div>
