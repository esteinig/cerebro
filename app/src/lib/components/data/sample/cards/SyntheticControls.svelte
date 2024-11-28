<script lang="ts">
	import { ERCC_CONCENTRATIONS, ERCC_GROUPS } from "$lib/utils/helpers";
	import type { Cerebro, ErccControl } from "$lib/utils/types";
	import { sampleCorrelation } from "simple-statistics";
	import { ScatterChart, ChartTheme, ScaleTypes } from "@carbon/charts-svelte";

	export let selectedModel: Cerebro;

	interface CorrelationData {
		group: string;
		alignments: number;
		concentration: number;
	}

	const isCorrelationData = (item: CorrelationData | undefined): item is CorrelationData => !!item;

	const getErccCorrelationData = (ercc: ErccControl | null): CorrelationData[] =>
		ercc?.records
			.map((record) => {
				let concentration = ERCC_CONCENTRATIONS.get(record.reference);
				let group = ERCC_GROUPS.get(record.reference);
				if (concentration !== undefined && record.alignments > 0) {
					return {
						group: group || "O",
						concentration: Number(concentration),
						alignments: Number(record.alignments),
					};
				}
			})
			.filter(isCorrelationData) ?? [];

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

	let correlationData: CorrelationData[] = getErccCorrelationData(selectedModel.quality.controls.ercc);

	let correlationCoefficient: number = 0;
    let minAlignments: number = 0;
    let maxAlignments: number = 0;
    let minConcentration: number = 0;
    let maxConcentration: number = 0;

	interface GroupCoefficients {
		A: number;
		B: number;
		C: number;
		D: number;
	}

	let groupCoefficients: GroupCoefficients = { A: 0, B: 0, C: 0, D: 0 };

	const getGroupCoefficients = (correlationData: CorrelationData[]): GroupCoefficients => {
		const calculateGroup = (group: string) => {
			const filteredData = correlationData.filter((x) => x.group === group);
			if (filteredData.length === 0) return 0;
			return sampleCorrelation(
				filteredData.map((x) => x.concentration),
				filteredData.map((x) => x.alignments)
			) ** 2;
		};

		return {
			A: calculateGroup("A"),
			B: calculateGroup("B"),
			C: calculateGroup("C"),
			D: calculateGroup("D"),
		};
	};

	if (correlationData.length > 0) {
		const pearsonCorrelation = sampleCorrelation(
			correlationData.map((x) => x.concentration),
			correlationData.map((x) => x.alignments)
		);
		correlationCoefficient = pearsonCorrelation ** 2;


        maxAlignments = getMaxAlignments();
        minAlignments = getMinAlignments();
        maxConcentration = getMaxConcentration();
        minConcentration = getMinConcentration();

		groupCoefficients = getGroupCoefficients(correlationData);
	}

	let groupColors: Record<string, string> = {};

	if (typeof window !== "undefined") {
		const style = getComputedStyle(document.body);
		groupColors = {
			'A': `rgb(${style.getPropertyValue("--color-primary-200").split(" ").join(",")})`,
			'B': `rgb(${style.getPropertyValue("--color-primary-400").split(" ").join(",")})`,
			'C': `rgb(${style.getPropertyValue("--color-primary-600").split(" ").join(",")})`,
			'D': `rgb(${style.getPropertyValue("--color-primary-800").split(" ").join(",")})`,
		};
	}

	$: correlationPlotOptions = {
        theme: ChartTheme.G100,
        color: { scale: groupColors },
        toolbar: { enabled: false }, 
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
    {#if correlationData.length}
        <p class="text-xs opacity-60 flex justify-start pb-8">
            <span>ERCC</span>
            <span class="ml-5 mr-5">n = {selectedModel.quality.reads.ercc_constructs} / 92</span>
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

<link rel="stylesheet" href="/carbon.css">