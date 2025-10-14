<script lang="ts">
	
    import { ERCC_CONCENTRATIONS, ERCC_GROUPS, getCssVariableAsHex } from "$lib/utils/helpers";
	import type { Cerebro, ErccControl } from "$lib/utils/types";
	import { sampleCorrelation } from "simple-statistics";
	import { ScatterChart, ChartTheme, ScaleTypes, type ScatterChartOptions } from "@carbon/charts-svelte";
	import { storeTheme } from "$lib/stores/stores";

	export let selectedModel: Cerebro;

	interface CorrelationData {
		group: string;
		alignments: number;
		concentration: number;
	}

	const isCorrelationData = (item: CorrelationData | undefined): item is CorrelationData => !!item;

	function getErccCorrelationData(ercc: ErccControl | null): CorrelationData[] {
		if (ercc === null) return [];

		return ercc.records
			.map((record) => {
				const concentration = ERCC_CONCENTRATIONS.get(record.reference);
				const group = ERCC_GROUPS.get(record.reference);
				if (concentration !== undefined && record.alignments > 0) {
					return {
						group: group || "O",
						concentration: Math.max(Number(concentration), 1e-3), // Ensure non-zero
						alignments: Math.max(Number(record.alignments), 1e-3), // Ensure non-zero
					};
				}
			})
			.filter(isCorrelationData) ?? [];
	}

	let correlationData: CorrelationData[] = [];
	$: correlationData = getErccCorrelationData(selectedModel.quality.controls.ercc);

	let correlationCoefficient: number = 0;
	let groupCoefficients = { A: 0, B: 0, C: 0, D: 0 };

	// helpers
	const toNum = (v: unknown) => {
		const n = typeof v === 'string' ? Number(v.trim()) : Number(v);
		return Number.isFinite(n) ? n : null;
	};
	const hasVariance = (arr: number[]) => new Set(arr).size > 1;

	function cleanPairs(data: CorrelationData[]) {
		return (data ?? [])
			.map(d => [toNum(d.concentration), toNum(d.alignments)] as const)
			.filter(([x, y]) => x !== null && y !== null) as Array<[number, number]>;
	}

	const getGroupCoefficients = (data: CorrelationData[]) => {
		const safe = (group: string) => {
			const pairs = cleanPairs(data.filter(d => d.group === group));
			if (pairs.length < 2) return 0;
			const xs = pairs.map(p => p[0]);
			const ys = pairs.map(p => p[1]);
			if (!hasVariance(xs) || !hasVariance(ys)) return 0;
			return sampleCorrelation(xs, ys) ** 2;
		};
		return {
			A: safe('A'),
			B: safe('B'),
			C: safe('C'),
			D: safe('D'),
		};
	};

	$: {
		const pairs = cleanPairs(correlationData);
		if (pairs.length >= 2) {
			const xs = pairs.map(p => p[0]);
			const ys = pairs.map(p => p[1]);
			correlationCoefficient = (hasVariance(xs) && hasVariance(ys))
			? sampleCorrelation(xs, ys) ** 2
			: 0;
		} else {
			correlationCoefficient = 0;
		}
		groupCoefficients = getGroupCoefficients(correlationData);
	}

	let groupColors: Record<string, string> = {};
	
	if (typeof window !== "undefined") {
		groupColors = {
			A: getCssVariableAsHex("--color-primary-200", $storeTheme) ?? "#ffffff",
			B: getCssVariableAsHex("--color-primary-400", $storeTheme) ?? "#ffffff",
			C: getCssVariableAsHex("--color-primary-600", $storeTheme) ?? "#ffffff",
			D: getCssVariableAsHex("--color-primary-800", $storeTheme) ?? "#ffffff",
		};
	}

	let correlationPlotOptions: ScatterChartOptions;

	$: {
		if (correlationData.length > 0) {
			const minAlignments = Math.max(Math.min(...correlationData.map((d) => d.alignments)), 1e-3);
			const maxAlignments = Math.max(...correlationData.map((d) => d.alignments)) * 1.2;

			const minConcentration = Math.max(Math.min(...correlationData.map((d) => d.concentration)), 1e-3);
			const maxConcentration = Math.max(...correlationData.map((d) => d.concentration)) * 1.2;

            // Generate logarithmic ticks
            const generateLogTicks = (min: number, max: number) => {
                const ticks = [];
                let value = Math.pow(10, Math.floor(Math.log10(min)));
                while (value <= max) {
                    ticks.push(value);
                    value *= 10;
                }
                return ticks;
            };

	        const logTicksAlignments = generateLogTicks(minAlignments, maxAlignments);

			correlationPlotOptions = {
				theme: ChartTheme.G100,
				color: { scale: groupColors },
				toolbar: { enabled: false },
				title: '',
				height: '450px',
				grid: { y: { enabled: false }, x: { enabled: false } },
				axes: {
					left: {
						title: "Alignments",
						mapsTo: "alignments",
						scaleType: ScaleTypes.LOG,
						domain: [minAlignments-(0.5*minAlignments), maxAlignments+(0.5*maxAlignments)],
                        ticks: {
                            values: logTicksAlignments, // Explicitly set log ticks
                            formatter: (tick: any) =>
                                typeof tick === "number"
                                    ? new Intl.NumberFormat("en", { notation: "scientific" }).format(tick)
                                    : tick.toString(),
                        },
					},
					bottom: {
						title: "Concentration (attomoles/μL)",
						mapsTo: "concentration",
						scaleType: ScaleTypes.LOG,
						domain: [minConcentration-(0.5*minConcentration), maxConcentration+(0.5*maxConcentration)],
						ticks: {
							formatter: (tick: any) =>
								typeof tick === "number"
									? new Intl.NumberFormat("en", { notation: "scientific" }).format(tick)
									: tick.toString(),
							number: 5,
						},
					},
				},
			};
		}
	}
</script>


<div class="">                
    {#if correlationData.length > 0}
        <p class="text-xs opacity-60 flex justify-start pb-8">
            <span>ERCC</span>
            <span class="ml-5 mr-5">n = {selectedModel.quality.reads.ercc_constructs} / 92</span>
            <span>r<sup>2</sup> = {correlationCoefficient.toFixed(4)}</span>
            <span class="ml-5 text-xs opacity-40">r<sup>2</sup><sub>A</sub> = {groupCoefficients.A.toFixed(4)}</span>
            <span class="ml-5 text-xs opacity-40">r<sup>2</sup><sub>B</sub> = {groupCoefficients.B.toFixed(4)}</span>
            <span class="ml-5 text-xs opacity-40">r<sup>2</sup><sub>C</sub> = {groupCoefficients.C.toFixed(4)}</span>
            <span class="ml-5 text-xs opacity-40">r<sup>2</sup><sub>D</sub> = {groupCoefficients.D.toFixed(4)}</span>
        </p>
        <div class="qc-ercc-chart">
            <ScatterChart data={correlationData} options={correlationPlotOptions}></ScatterChart>
        </div>
        
    {:else}
        <div class="flex justify-center">
            <p class="text-sm opacity-60 pt-24">
                ERCC/EDCC alignment data not found
            </p>
        </div>
    {/if}

</div>


<style lang="postcss">
	/* Scoped globally but applied only within this component */
    :global(.qc-ercc-chart .chart-grid-backdrop) {
        display: none;
    }
</style>