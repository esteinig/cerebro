<script lang="ts">
	import { formatAsThousands, ERCC_CONCENTRATIONS, ERCC_GROUPS } from "$lib/utils/helpers";
	import type { Cerebro, Ercc, QualityControlSummary, WorkflowParamsQcErcc } from "$lib/utils/types";
	import { Tab, TabGroup } from "@skeletonlabs/skeleton";
    import { sampleCorrelation } from "simple-statistics";
    
	import { ScatterChart, ChartTheme, ScaleTypes } from '@carbon/charts-svelte';

    export let selectedModel: Cerebro;
    export let selectedQualityControlSummary: QualityControlSummary;
    export let selectedStageParams: WorkflowParamsQcErcc | undefined;
    export let selectedStage: string = "Synthetic controls";
    
    interface CorrelationData {
        group: string,
        alignments: number,
        concentration: number
    }

    const isCorrelationData = (item: CorrelationData | undefined): item is CorrelationData => { return !!item }
    
    const getErccCorrelationData = (ercc: Ercc | null): CorrelationData[] => {
        if (ercc === null) {
            return []
        } 
        return ercc.records.map(record => {
            let concentration = ERCC_CONCENTRATIONS.get(record.reference);
            let group = ERCC_GROUPS.get(record.reference);
            if (concentration !== undefined && record.alignments > 0){  // LOG scale plot guard
                return {
                    group: group === undefined ? "O": group,
                    concentration: Number(concentration),
                    alignments: Number(record.alignments)
                } satisfies CorrelationData
            }
        }).filter(isCorrelationData);
    }

    const getMaxAlignments = (): number => {
        return correlationData.reduce(function(prev, current) {
            return (prev.alignments > current.alignments) ? prev : current
        }).alignments
    }
    const getMinAlignments = (): number => {
        return correlationData.reduce(function(prev, current) {
            return (prev.alignments < current.alignments) ? prev : current
        }).alignments
    }

    const getMaxConcentration = (): number => {
        return correlationData.reduce(function(prev, current) {
            return (prev.concentration > current.concentration) ? prev : current
        }).concentration
    }
    const getMinConcentration = (): number => {
        return correlationData.reduce(function(prev, current) {
            return (prev.concentration < current.concentration) ? prev : current
        }).concentration
    }

    let tabSet: number = 0;
    let correlationData: CorrelationData[] = getErccCorrelationData(selectedModel.quality.ercc);
    let correlationCoefficient: number = 0;

    let minAlignments: number = 0;
    let maxAlignments: number = 0;
    let minConcentration: number = 0;
    let maxConcentration: number = 0;


    interface GroupCoefficients {
        A: number,
        B: number,
        C: number,
        D: number
    }

    let groupCoefficients: GroupCoefficients;

    const getGroupCoefficients = (correlationData: CorrelationData[]): GroupCoefficients => {

        let groupA = correlationData.filter(x => x.group === "A");
        let pearsonA = sampleCorrelation(groupA.map(x => x.concentration), groupA.map(x => x.alignments));
        let groupB = correlationData.filter(x => x.group === "B");
        let pearsonB = sampleCorrelation(groupB.map(x => x.concentration), groupB.map(x => x.alignments));
        let groupC = correlationData.filter(x => x.group === "C");
        let pearsonC = sampleCorrelation(groupC.map(x => x.concentration), groupC.map(x => x.alignments));
        let groupD = correlationData.filter(x => x.group === "D");
        let pearsonD = sampleCorrelation(groupD.map(x => x.concentration), groupD.map(x => x.alignments));

        return {
            A: pearsonA*pearsonA,
            B: pearsonB*pearsonB,
            C: pearsonC*pearsonC,
            D: pearsonD*pearsonD,
        };
    }

    if (correlationData.length) {
        let pearsonCorrelation = sampleCorrelation(correlationData.map(x => x.concentration), correlationData.map(x => x.alignments))
        correlationCoefficient = pearsonCorrelation*pearsonCorrelation;
        
        maxAlignments = getMaxAlignments();
        minAlignments = getMinAlignments();
        maxConcentration = getMaxConcentration();
        minConcentration = getMinConcentration();

        groupCoefficients = getGroupCoefficients(correlationData);

    }
    
    // Skeleton color preset:
    var style = getComputedStyle(document.body);

    let groupColors: any = {
        "A": `rgb(${style.getPropertyValue('--color-primary-200').split(" ").join(",")})`,
        "B": `rgb(${style.getPropertyValue('--color-primary-400').split(" ").join(",")})`,
        "C": `rgb(${style.getPropertyValue('--color-primary-600').split(" ").join(",")})`,
        "D": `rgb(${style.getPropertyValue('--color-primary-800').split(" ").join(",")})`,
    };

    
    correlationData.sort((a, b) => a.group > b.group ? 1 : -1);
    
    let correlationPlotOptions: any = {
        theme: ChartTheme.G100,
        color: { scale: groupColors},
        title: ``,
        height: '450px',
        grid: {
            y: {
                enabled: false
            },
            x: {
                enabled: false
            }
        },
        axes: {
            left: { title: "Alignments", mapsTo: 'alignments', scaleType: ScaleTypes.LOG, ticks: { min: 0 }, domain: [minAlignments-(0.5*minAlignments), maxAlignments+(0.2*maxAlignments)]},  // extend plot field to fully show highest point
            bottom: { title: "Concentration (attomoles/ul)",  mapsTo: 'concentration', scaleType: ScaleTypes.LOG, domain: [minConcentration-(0.5*minConcentration), maxConcentration+(0.2*maxConcentration)], ticks: {
                formatter: (tick: any, i: number) => {
                    if (typeof(tick) === 'number')
                        return new Intl.NumberFormat('en', { notation: 'scientific' }).format(tick) 
                    return tick.toString()
                },
                number: 5,
                min: 0,
            }}

        }
    }



</script>

<div class="">
    <TabGroup>
        <Tab bind:group={tabSet} name="overview" value={0}>
            <span class="text-sm">Results</span>
        </Tab>
        <Tab bind:group={tabSet} name="correlation" value={1}>
            <span class="text-sm">Correlation</span>
        </Tab>
        <Tab bind:group={tabSet} name="description" value={2}>
            <span class="text-sm">Description</span>
        </Tab>
        <Tab bind:group={tabSet} name="configuration" value={3}>
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
                                                ERCC / EDCC library controls
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
                                <div><span class="opacity-60 mr-4">Reads removed</span></div>
                                <div>{formatAsThousands(selectedQualityControlSummary.ercc_reads)}</div>
                                <div><span class="opacity-60 mr-4">Constructs detected</span></div>
                                <div>{selectedQualityControlSummary.ercc_constructs} / 92</div>
                                <div><span class="opacity-60 mr-4">Input mass</span></div>
                                <div>{selectedQualityControlSummary.ercc_input_mass} pg</div>
                                <div><span class="opacity-60 mr-4">Mass per read</span></div>
                                <div>{selectedQualityControlSummary.ercc_mass_per_read} pg</div>
                            </div>
                        </div>
                    </section>
                    <footer class="card-footer">
                
                        
                    </footer>
                </div>
            
            {:else if tabSet === 1} 
                
                {#if correlationData.length}
                    <p class="text-sm opacity-60">
                        <span>r<sup>2</sup> = {correlationCoefficient.toFixed(4)}</span>
                        <span class="ml-5 text-xs opacity-40">r<sup>2</sup><sub>A</sub> = {groupCoefficients.A.toFixed(4)}</span>
                        <span class="ml-5 text-xs opacity-40">r<sup>2</sup><sub>B</sub> = {groupCoefficients.B.toFixed(4)}</span>
                        <span class="ml-5 text-xs opacity-40">r<sup>2</sup><sub>C</sub> = {groupCoefficients.C.toFixed(4)}</span>
                        <span class="ml-5 text-xs opacity-40">r<sup>2</sup><sub>D</sub> = {groupCoefficients.D.toFixed(4)}</span>
                    </p>
                    <ScatterChart data={correlationData} options={correlationPlotOptions}></ScatterChart>
                {:else}
                    <div class="flex justify-center">
                        <p class="text-sm opacity-60 pt-24">
                            ERCC/EDCC alignment data not found
                        </p>
                    </div>
                {/if}

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
                {:else if tabSet === 3}
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

<!-- 
Overwrite the root variable that is responsible
for setting the background colors of the Carbon  
plot - there is probably a way to adjust these 
to the Skeleton UI / Tailwind color scheme 
variables but it is eluding me...
-->
<style lang="postcss">
   :root {
        --cds-grid-bg: rgb(0, 0, 0, 0);
   }
</style>