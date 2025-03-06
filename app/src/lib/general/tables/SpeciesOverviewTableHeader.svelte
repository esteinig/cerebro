<script lang="ts">
	import { getValueDisplay } from "$lib/utils/helpers";
	import { DisplayData, ProfileTool } from "$lib/utils/types";
	import CircleIndicator from "../icons/CircleIndicator.svelte";

	export let sortColumn: string;
	export let displayData: DisplayData;
	export let tool: ProfileTool | string;
	export let circleColor: string | null = null;
	export let sortOrder: "asc" | "desc" | null = null;

	// Mapping from header display names to sort keys.
	const toolToSortKey: Record<string, string> = {
		"Average": "total",
		"Assembly": "contigs",
	};

	// Derive the header's sort key.
	$: headerSortKey = toolToSortKey[tool] ?? tool.toString().toLowerCase();

</script>

<div class="text-right flex flex-col items-end">
	<!-- Top row: tool name and SVG next to each other -->
	<div class="opacity-60 flex items-center">
		<span>{tool === ProfileTool.Blast ? "Assembly" : tool}</span>
		{#if sortColumn === headerSortKey}
			{#if sortOrder === "asc"}
				<svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="ml-1 h-3 w-3">
					<path stroke-linecap="round" stroke-linejoin="round" d="m4.5 15.75 7.5-7.5 7.5 7.5" />
				</svg>
			{:else if sortOrder === "desc"} 
				<svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="ml-1 h-3 w-3">
					<path stroke-linecap="round" stroke-linejoin="round" d="m19.5 8.25-7.5 7.5-7.5-7.5" />
				</svg>
			{/if}
		{:else}     
			<svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="ml-1 h-3 w-3">
				<path stroke-linecap="round" stroke-linejoin="round" d="M8.25 15 12 18.75 15.75 15m-7.5-6L12 5.25 15.75 9" />
			</svg>   
		{/if}
	</div>
	<!-- Bottom row: circle indicator and value display -->
	<div class="flex items-center">
		{#if circleColor}
			<CircleIndicator circleClass="mt-0.5 mr-1" color={circleColor} />
		{/if}
        {#if tool === ProfileTool.Blast}
            <span class="text-xs opacity-40 mr-3">bp</span>
        {:else if tool === "Assembly"}    
		    <span class="text-xs opacity-40 mr-3">contigs</span>
        {:else}
		    <span class="text-xs opacity-40 mr-3">{getValueDisplay(displayData)}</span>
        {/if}
	</div>
</div>