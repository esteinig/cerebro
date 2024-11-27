<script lang="ts">
	import { BarChartSimple, ChartTheme, ScaleTypes } from '@carbon/charts-svelte';
	import type { OrganismControl, Cerebro } from '$lib/utils/types';

	export let selectedModel: Cerebro;

	// Extract data from the OrganismControl
	const getOrganismControlData = (organismControl: OrganismControl | null) => {
		if (!organismControl || !organismControl.records) return [];
		return organismControl.records.map((record) => ({
			group: record.reference,
			coverage: record.coverage,
		}));
	};

	let organismControlData = getOrganismControlData(selectedModel.quality.controls.organism);

	// Sort by alignments for better visualization
	organismControlData.sort((a, b) => a.coverage - b.coverage);

    let lineFillColor: string = "rgb(255, 0, 0)"
    
	let barChartColors: Record<string, string> = {};
    
	if (typeof window !== "undefined") {
		const style = getComputedStyle(document.body);
        const generateColors = (data: { group: string }[]) => {
            const colorGroups = [
                "--color-primary-",
                "--color-secondary-",
                "--color-tertiary-",
            ];
            const steps = [200, 400, 600, 800];

            let colorScale: Record<string, string> = {};
            data.forEach((item, index) => {
                const colorGroup = colorGroups[Math.floor(index / steps.length) % colorGroups.length];
                const step = steps[index % steps.length];
                const colorVar = `${colorGroup}${step}`;
                colorScale[item.group] = `rgb(${style.getPropertyValue(colorVar).split(" ").join(",")})`;
            });

            return colorScale;
        };
        barChartColors = generateColors(organismControlData);
        lineFillColor = `rgba(${style.getPropertyValue('--color-secondary-400').split(" ").join(",")})`
	}
    

	// Define chart options
	$: barChartOptions = {
		theme: ChartTheme.G100,
        toolbar: { enabled: false }, 
		color: {
			scale: barChartColors
		},
		title: "",
		height: "450px",
		grid: {
			y: { enabled: false },
			x: { enabled: false },
		},
		axes: {
			left: {
				title: "",
				mapsTo: "group",
				scaleType: ScaleTypes.LABELS,
			},
			bottom: {
				title: "Genome coverage (%)",
				mapsTo: "coverage",
				scaleType: ScaleTypes.LINEAR,

                thresholds: [
                    {
                        value: 80, // x position for the line
                        label: "Pass",
                        fillColor: lineFillColor
                    },
                ],
			},
		},
		bars: {
			maxWidth: 32, // Increases bar width
		},
	};
</script>

<div class="">   
    {#if organismControlData.length}
        
    <p class="text-xs opacity-60 flex justify-start pb-8">
        <span>Internal control coverage</span>
    </p>
        <BarChartSimple data={organismControlData} options={barChartOptions} />
    {:else}
        <div class="flex justify-center">
            <p class="text-sm opacity-60 pt-24">No organism control data available</p>
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