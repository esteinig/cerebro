<script lang="ts">
	import { type QualityControlSummary } from "$lib/utils/types";
    import { MeterChart, ChartTheme } from '@carbon/charts-svelte'
    import { Statuses } from '@carbon/charts';

    export let selectedQualityControlSummary: QualityControlSummary;

    let data = [
        {
            group: 'Depth',
            value: 10
        },
        {
            group: 'Deduplication',
            value: 10
        },
        {
            group: 'Controls',
            value: 10
        },
        {
            group: 'Read QC',
            value: 10
        }
    ]
    
	let groupColors: Record<string, string> = {};

    if (typeof window !== "undefined") {
        const style = getComputedStyle(document.body);
        groupColors = {
            'Depth': `rgb(${style.getPropertyValue("--color-primary-100").split(" ").join(",")})`,
            'Deduplication': `rgb(${style.getPropertyValue("--color-primary-200").split(" ").join(",")})`,
            'Controls': `rgb(${style.getPropertyValue("--color-primary-300").split(" ").join(",")})`,
            'Read QC': `rgb(${style.getPropertyValue("--color-primary-400").split(" ").join(",")})`,
        };
    }
    
    $: options = {
		theme: ChartTheme.G100,
        toolbar: { enabled: false }, 
        title: '',
        color: { scale: groupColors },
        height: '200px',
        meter: {
            peak: 60,
            proportional: {
                total: 100,
            },
            status: {
                ranges: [
                    {
                        range: [0, 50] as [number, number],
                        status: 'danger' as Statuses
                    },
                    {
                        range: [50, 60] as [number, number],
                        status: 'warning' as Statuses
                    },
                    {
                        range: [60, 100] as [number, number],
                        status: 'success' as Statuses
                    }
                ]
            }
        },
    };


</script>

<div class="p-6 rounded-lg space-y-2">
    <!-- Highlight Input Reads -->
    <div class="text-center space-y-1">
        <p class="text-lg font-semibold text-primary-600 dark:text-primary-300 opacity-40">Input Reads</p>
        <p class="text-2xl font-bold text-gray-900 dark:text-white">
            <span class="text-primary-800 dark:text-primary-400">{selectedQualityControlSummary.input_reads.toLocaleString()}</span>
            <span class="ml-1 text-sm text-gray-500 dark:text-gray-300">(100%)</span>
        </p>
    </div>

    <div class="grid grid-cols-2 px-32">
    <!-- Reads Section -->
    <div>
        <div class="grid grid-cols-2 sm:grid-cols-1 gap-4 text-medium text-center pt-8">
            <div class="space-y-1">
                <p class="font-medium opacity-40">Synthetic Controls</p>
                <p class="text-gray-800 dark:text-gray-100">
                    <span>{selectedQualityControlSummary.ercc_reads?.toLocaleString() ?? 'N/A'}</span>
                    <span class="ml-1 text-sm text-primary-600 dark:text-primary-400">
                        ({selectedQualityControlSummary.ercc_reads_percent?.toFixed(2) ?? 'N/A '}%)
                    </span>
                </p>
            </div>
            <div class="space-y-1">
                <p class="font-medium opacity-40">Internal Controls</p>
                <p class="text-gray-800 dark:text-gray-100">
                    <span>{selectedQualityControlSummary.control_reads?.toLocaleString() ?? 'N/A'}</span>
                    <span class="ml-1 text-sm text-primary-600 dark:text-primary-400">
                        ({selectedQualityControlSummary.control_reads_percent?.toFixed(2) ?? 'N/A '}%)
                    </span>
                </p>
            </div>
            <div class="space-y-1">
                <p class="font-medium opacity-40">Deduplicated</p>
                <p class="text-gray-800 dark:text-gray-100">
                    <span>{selectedQualityControlSummary.deduplicated_reads?.toLocaleString() ?? 'N/A'}</span>
                    <span class="ml-1 text-sm text-primary-600 dark:text-primary-400">
                        ({selectedQualityControlSummary.deduplicated_percent?.toFixed(2) ?? 'N/A '}%)
                    </span>
                </p>
            </div>
            <div class="space-y-1">
                <p class="font-medium opacity-40">Quality Control</p>
                <p class="text-gray-800 dark:text-gray-100">
                    <span>{selectedQualityControlSummary.qc_reads?.toLocaleString() ?? 'N/A'}</span>
                    <span class="ml-1 text-sm text-primary-600 dark:text-primary-400">
                        ({selectedQualityControlSummary.qc_reads_percent?.toFixed(2) ?? 'N/A '}%)
                    </span>
                </p>
            </div>
            <div class="space-y-1">
                <p class="font-medium opacity-40">Host Depletion</p>
                <p class="text-gray-800 dark:text-gray-100">
                    <span>{selectedQualityControlSummary.host_reads?.toLocaleString() ?? 'N/A'}</span>
                    <span class="ml-1 text-sm text-primary-600 dark:text-primary-400">
                        ({selectedQualityControlSummary.host_reads_percent?.toFixed(2) ?? 'N/A '}%)
                    </span>
                </p>
            </div>
            <div class="space-y-1">
                <p class="font-medium opacity-40">Background Depletion</p>
                <p class="text-gray-800 dark:text-gray-100">
                    <span>{selectedQualityControlSummary.background_reads?.toLocaleString() ?? 'N/A'}</span>
                    <span class="ml-1 text-sm text-primary-600 dark:text-primary-400">
                        ({selectedQualityControlSummary.background_reads_percent?.toFixed(2) ?? 'N/A '}%)
                    </span>
                </p>
            </div>
        </div>
    </div>

    <!-- Length and Quality Section -->
    <div>
        <div class="grid grid-cols-1 sm:grid-cols-1 gap-4 text-medium text-center pt-8">
            <div class="space-y-1">
                <p class="font-medium opacity-40">Low Quality</p>
                <p class="text-gray-800 dark:text-gray-100">
                    <span>{selectedQualityControlSummary.low_quality_reads?.toLocaleString() ?? 'N/A'}</span>
                    <span class="ml-1 text-sm text-primary-600 dark:text-primary-400">
                        ({selectedQualityControlSummary.low_quality_percent?.toFixed(2) ?? 'N/A '}%)
                    </span>
                </p>
            </div>
            <div class="space-y-1">
                <p class="font-medium opacity-40">Low Complexity</p>
                <p class="text-gray-800 dark:text-gray-100">
                    <span>{selectedQualityControlSummary.low_complexity_reads?.toLocaleString() ?? 'N/A'}</span>
                    <span class="ml-1 text-sm text-primary-600 dark:text-primary-400">
                        ({selectedQualityControlSummary.low_complexity_percent?.toFixed(2) ?? 'N/A '}%)
                    </span>
                </p>
            </div>
            <div class="space-y-1">
                <p class="font-medium opacity-40">Adapter Trimmed</p>
                <p class="text-gray-800 dark:text-gray-100">
                    <span>{selectedQualityControlSummary.adapter_trimmed_reads?.toLocaleString() ?? 'N/A'}</span>
                    <span class="ml-1 text-sm text-primary-600 dark:text-primary-400">
                        ({selectedQualityControlSummary.adapter_trimmed_percent?.toFixed(2) ?? 'N/A '}%)
                    </span>
                </p>
            </div>
            <div class="space-y-1">
                <p class="font-medium opacity-40">Min Length</p>
                <p class="text-gray-800 dark:text-gray-100">
                    <span>{selectedQualityControlSummary.min_length_reads?.toLocaleString() ?? 'N/A'}</span>
                    <span class="ml-1 text-sm text-primary-600 dark:text-primary-400">
                        ({selectedQualityControlSummary.min_length_percent?.toFixed(2) ?? 'N/A '}%)
                    </span>
                </p>
            </div>
            <div class="space-y-1">
                <p class="font-medium opacity-40">Mean Length</p>
                <p class="text-gray-800 dark:text-gray-100">
                    {selectedQualityControlSummary.mean_read_length_r1?.toLocaleString() ?? 'N/A'} bp
                    {#if selectedQualityControlSummary.mean_read_length_r2}
                    | {selectedQualityControlSummary.mean_read_length_r2?.toLocaleString() ?? 'N/A'} bp
                    {/if}
                </p>
                
            </div>
            <div class="space-y-1">
                <p class="font-medium opacity-40">Q20 | Q30</p>
                <p class="text-gray-800 dark:text-gray-100">
                    {selectedQualityControlSummary.q20_percent?.toFixed(2) ?? 'N/A '}%
                    | {selectedQualityControlSummary.q30_percent?.toFixed(2) ?? 'N/A '}%
                </p>
            </div>
        </div>
    </div>
    </div>

    <!-- Highlight Output Reads -->
    <div class="text-center space-y-1 pt-4">
        <p class="text-lg font-semibold text-primary-600 dark:text-primary-300 opacity-40">Output Reads</p>
        <p class="text-2xl font-bold text-gray-900 dark:text-white">
            <span class="text-primary-800 dark:text-primary-400">{selectedQualityControlSummary.output_reads.toLocaleString()}</span>
            <span class="ml-1 text-sm text-gray-500 dark:text-gray-300">
                ({selectedQualityControlSummary.output_reads_percent.toFixed(2)}%)
            </span>
        </p>
    </div>

    <div class="pt-16 px-32">
        <MeterChart {data} {options} />
    </div>

</div>

<style lang="postcss">
    :root {
         --cds-grid-bg: rgb(0, 0, 0, 0);
    }
    /* Scoped globally but applied only within this component */
    :global(.proportional-meter-total) {
        display: none;
    }
    /* Scoped globally but applied only within this component */
    :global(.proportional-meter-title) {
        display: none;
    }
    /* Scoped globally but applied only within this component */
    :global(.status-indicator) {
        display: none;
    }
</style>
 
<link rel="stylesheet" href="/carbon.css">