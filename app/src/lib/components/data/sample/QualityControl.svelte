<script lang="ts">
	import type { Cerebro, QualityControlSummary, WorkflowConfig } from "$lib/utils/types";
    import { calculateOutput, calculateTotalQuality, formatAsPercentage, formatAsThousands } from "$lib/utils/helpers";
	import InputStageCard from "./cards/InputStageCard.svelte";
	import DeduplicationStageCard from "./cards/DeduplicationStageCard.svelte";
	import SyntheticControlStageCard from "./cards/SyntheticControlStageCard.svelte";
	import ReadQualityStageCard from "./cards/ReadQualityStageCard.svelte";
	import HostDepletionStageCard from "./cards/HostDepletionStageCard.svelte";
	import PhageControlStageCard from "./cards/PhageControlStageCard.svelte";

    export let selectedWorkflowConfiguration: WorkflowConfig;
    export let selectedModels: Cerebro[] = [];
    export let selectedQualityControlSummaries: QualityControlSummary[] = [];
    
    let selectedStages: string[] = [];


    const changeStage = (stage: string, i: number) => {
        selectedStages[i] = stage;
    }

    $: {
        selectedStages = selectedModels.map(_ => "Input");
    }

    const getTotalPhageReads = (selectedQualityControlSummaries: QualityControlSummary[], i: number): number => {
        
        let dnaReads: number = selectedQualityControlSummaries[i].dna_phage_reads || 0;
        let rnaReads: number = selectedQualityControlSummaries[i].rna_phage_reads || 0;
        let seqReads: number = selectedQualityControlSummaries[i].seq_phage_reads || 0;

        return dnaReads+rnaReads+seqReads
    }
    const getTotalPhagePercent = (selectedQualityControlSummaries: QualityControlSummary[], i: number): number  => {
        let qc = selectedQualityControlSummaries[i];

        let dnaPercent: number = qc.dna_phage_percent ? +qc.dna_phage_percent : 0.;
        let rnaPercent: number = qc.rna_phage_percent ? +qc.rna_phage_percent : 0.;
        let seqPercent: number = qc.seq_phage_percent ? +qc.seq_phage_percent : 0.;

        return  dnaPercent+rnaPercent+seqPercent
    }

</script>

<div class="w-full">
    <div class="flex pb-10 opacity-90 text-sm">
        <ol class="breadcrumb justify-center">
            <li class="crumb opacity-70">Input</li>
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            {#if selectedWorkflowConfiguration?.params?.qc.deduplication.enabled}
                <li class="crumb">Deduplication</li>
            {:else}
                <li class="crumb opacity-10">Deduplication</li>
            {/if}
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            {#if selectedWorkflowConfiguration?.params?.qc.controls.ercc.enabled}
                <li class="crumb">Synthetic Controls</li>
            {:else}
                <li class="crumb opacity-10">Synthetic Controls</li>
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
                <li class="crumb">Phage Controls</li>
            {:else}
                <li class="crumb opacity-10">Phage Controls</li>
            {/if}
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            <li class="crumb opacity-70">Output</li>
        </ol>
    </div>
    <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-12">
        
        {#each selectedModels as model, i}
            <div>
                <p class="flex align-center"><span class="code bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90 text-sm">{model.sample.id}</span> <span class="code ml-2 text-sm">{model.sample.tags.join("-")}</span></p>
                <div class="table-container py-3 text-sm">
                    <table class="table table-hover table-compact">
                        <thead>
                            <tr>
                                <th class="text-center">Stage</th>
                                <th class="text-end">Reads</th>
                                <th class="text-end">Percent</th>
                                <th class="text-end">QC Pass</th>
                        </thead>
                        <tbody> 
                                <tr class="opacity-60 hover:cursor-pointer"  on:click={() => changeStage("Input", i)}>
                                    <td class="text-center"><span class="ml-1">Input</span></td>
                                    <td class="text-end"><span class="ml-1">{formatAsThousands(selectedQualityControlSummaries[i].total_reads)}</span></td>
                                    <td class="flex justify-end">{formatAsPercentage(100)}</td>
                                    <td class="text-end"><span class="ml-1">

                                    </span></td>
                                </tr>
                                <tr class="hover:cursor-pointer" on:click={() => changeStage("Deduplication", i)}>
                                    <td class="text-center"><span class="ml-1">Deduplication</span></td>
                                    <td class="text-end"><span class="ml-1">{formatAsThousands(selectedQualityControlSummaries[i].deduplicated_reads)}</span></td>
                                    <td class="text-end">{formatAsPercentage(selectedQualityControlSummaries[i].deduplicated_percent)}</td>
                                </tr>
                                <tr class="hover:cursor-pointer" on:click={() => changeStage("Synthetic Control", i)}>
                                    <td class="text-center"><span class="ml-1">Synthetic Controls</span></td>
                                    <td class="text-end"><span class="ml-1">{formatAsThousands(selectedQualityControlSummaries[i].ercc_reads)}</span></td>
                                    <td class="flex justify-end">{formatAsPercentage(selectedQualityControlSummaries[i].ercc_percent)}</td>
                                </tr>
                                <tr class="hover:cursor-pointer" on:click={() => changeStage("Read Quality", i)}>
                                    <td class="text-center"><span class="ml-1">Read Quality</span></td>
                                    <td class="text-end"><span class="ml-1">{formatAsThousands(
                                        calculateTotalQuality(
                                            selectedQualityControlSummaries[i].low_complexity_reads, 
                                            selectedQualityControlSummaries[i].qc_reads
                                        )
                                        )}</span></td>
                                    <td class="flex justify-end">{formatAsPercentage(
                                        calculateTotalQuality(
                                            selectedQualityControlSummaries[i].low_complexity_percent,
                                            selectedQualityControlSummaries[i].qc_percent
                                        )
                                        )}</td>
                                </tr>
                                <tr class="hover:cursor-pointer" on:click={() => changeStage("Host Depletion", i)}>
                                    <td class="text-center"><span class="ml-1">Host Depletion</span></td>
                                    <td class="text-end"><span class="ml-1">{formatAsThousands(selectedQualityControlSummaries[i].host_reads)}</span></td>
                                    <td class="flex justify-end">{formatAsPercentage(selectedQualityControlSummaries[i].host_percent)}</td>
                                </tr>
                                <tr class="hover:cursor-pointer" on:click={() => changeStage("Phage Control", i)}>
                                    <td class="text-center"><span class="ml-1">Phage Control</span></td>
                                    <td class="text-end"><span class="ml-1">{formatAsThousands(getTotalPhageReads(selectedQualityControlSummaries, i))}</span></td>
                                    <td class="flex justify-end">{formatAsPercentage(getTotalPhagePercent(selectedQualityControlSummaries, i))}</td>
                                </tr>
                                <tr class="opacity-60 hover:cursor-pointer" on:click={() => changeStage("Output", i)}>
                                    <td class="text-center"> <span class="ml-1">Output</span></td>
                                    <td class="text-end"><span class="ml-1">{formatAsThousands(selectedQualityControlSummaries[i].output_reads)}</span></td>
                                    <td class="text-end">{formatAsPercentage(selectedQualityControlSummaries[i].output_percent)}</td>
                                </tr>
                        </tbody>
                    </table>
                </div>
            </div>
            <div>
                {#if selectedStages[i] === "Input"}
                    <InputStageCard selectedQualityControlSummary={selectedQualityControlSummaries[i]}></InputStageCard>
                {:else if selectedStages[i] === "Deduplication"}
                    <DeduplicationStageCard selectedQualityControlSummary={selectedQualityControlSummaries[i]} selectedStageParams={selectedWorkflowConfiguration?.params?.qc.deduplication}></DeduplicationStageCard>
                {:else if selectedStages[i] === "Synthetic Control"}
                    <SyntheticControlStageCard selectedModel={selectedModels[i]} selectedQualityControlSummary={selectedQualityControlSummaries[i]} selectedStageParams={selectedWorkflowConfiguration?.params?.qc.controls.ercc}></SyntheticControlStageCard>
                {:else if selectedStages[i] === "Read Quality"}
                    <ReadQualityStageCard selectedQualityControlSummary={selectedQualityControlSummaries[i]} selectedStageParams={selectedWorkflowConfiguration?.params?.qc.reads.fastp}></ReadQualityStageCard>
                {:else if selectedStages[i] === "Host Depletion"}
                    <HostDepletionStageCard selectedQualityControlSummary={selectedQualityControlSummaries[i]} selectedStageParams={selectedWorkflowConfiguration?.params?.qc.host.depletion}></HostDepletionStageCard>
                {:else if selectedStages[i] === "Phage Control"}
                    <PhageControlStageCard selectedQualityControlSummary={selectedQualityControlSummaries[i]} selectedStageParams={selectedWorkflowConfiguration?.params?.qc.controls.phage}></PhageControlStageCard>
                {/if}
            </div>
        {/each}
    </div>
</div>
