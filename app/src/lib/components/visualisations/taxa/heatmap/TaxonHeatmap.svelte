<script lang="ts">
    import * as d3 from 'd3';
    import { selectedTaxa, selectedServerFilterConfig, storeTheme } from '$lib/stores/stores';
    import { type Taxon, type TaxonOverviewRecord, PathogenDetectionTool, PathogenDetectionMode, DisplayData } from '$lib/utils/types';
    import CerebroApi, { ApiResponse } from "$lib/utils/api";
    import { page } from '$app/stores';
    import { getCssVariableAsHex } from '$lib/utils/helpers';

    const publicApi = new CerebroApi();

    export let selectedIdentifiers: string[] = [];
    export let width: number = 1024;
    export let height: number = 768;
    export let tool: PathogenDetectionTool = PathogenDetectionTool.Ganon2;
    export let mode: PathogenDetectionMode = PathogenDetectionMode.Sequence;
    export let displayData: DisplayData = DisplayData.Rpm;

    let container: HTMLDivElement;
    let svg: SVGSVGElement;
    let g: SVGGElement;

    let taxa: Taxon[] = [];
    let uniqueDetectionIds: string[] = [];
    let dataMatrix: Array<{ row: string; column: string; value: number | null }> = [];
    let hoveredCell = { row: null, column: null };
    let loading = false;

    let rowOrder: string[] = []; // Tracks the rows in the heatmap
    let columnOrder: string[] = [];

    const getAggregatedTaxa = async (selectedIdentifiers: string[], selectedTaxa: TaxonOverviewRecord[]) => {
        if (!selectedTaxa || selectedTaxa.length === 0) return;

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.taxa}?team=${$page.params.team}&db=${$page.params.db}&project=${$page.params.project}&id=${selectedIdentifiers.join(",")}&overview=false&taxid=${selectedTaxa.map(taxon => taxon.taxid).join(",")}`,
            { 
                method: 'POST',  
                mode: 'cors',
                credentials: 'include', 
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify($selectedServerFilterConfig) 
            } as RequestInit,
            $page.data.refreshToken
        );
        
        loading = false;
        
        if (response.ok) {
            taxa = response.json.data.taxa;

            // Update rowOrder to reflect the changes
            rowOrder = [
                ...rowOrder.filter((row) => taxa.some((taxon) => taxon.name === row)), // Keep rows that still exist
                ...taxa.map((taxon) => taxon.name).filter((name) => !rowOrder.includes(name)) // Append new rows at the end
            ];
        }
    };

    $: getAggregatedTaxa(selectedIdentifiers, $selectedTaxa);

    // Derive unique Detection IDs (columns)
    $: uniqueDetectionIds = Array.from(new Set(taxa.flatMap((taxon) => taxon.evidence.records.map((record) => record.id))));

    // Preserve order and append new columns
    $: columnOrder = [...columnOrder, ...uniqueDetectionIds.filter((id) => !columnOrder.includes(id))];

    // Prepare data matrix for heatmap
    $: dataMatrix = rowOrder.flatMap((row) =>
        columnOrder.map((column) => {
            const matchingData = taxa.flatMap((taxon) =>
                taxon.evidence.records
                    .filter((r) => r.id === column)
                    .flatMap((record) =>
                        record.results
                            .filter((result) => result.tool === tool && result.mode === mode)
                            .map((result) => ({
                                row: taxon.name,
                                column: record.id,
                                value: result[displayData] ?? null,
                            }))
                    )
            ).find((d) => d.row === row && d.column === column);
            return matchingData || { row, column, value: null };
        })
    );

    $: rowMaxValues = taxa.reduce((acc, taxon) => {
        acc[taxon.name] = Math.max(...dataMatrix.filter((d) => d.row === taxon.name).map((d) => d.value || 0));
        return acc;
    }, {} as Record<string, number>);

    $: rowColorScales = Object.fromEntries(
        Object.entries(rowMaxValues).map(([row, maxValue]) => {
            const startColor = getCssVariableAsHex('--color-primary-200', $storeTheme);
            const endColor = getCssVariableAsHex('--color-primary-600', $storeTheme);
            return [
                row,
                d3.scaleLinear()
                    .domain([0, maxValue])
                    .range([startColor, endColor]),
            ];
        })
    );

    const margin = { top: 60, right: 20, bottom: 60, left: 150 };

    $: heatmapWidth = width - margin.left - margin.right;
    $: heatmapHeight = height - margin.top - margin.bottom;

    $: xScale = d3.scaleBand().domain(columnOrder).range([0, heatmapWidth]).padding(0.05);

    $: yScale = d3.scaleBand().domain(rowOrder).range([0, heatmapHeight]).padding(0.05);

    const handleMouseOver = (row, column) => {
        hoveredCell = { row, column };
    };

    const handleMouseOut = () => {
        hoveredCell = { row: null, column: null };
    };

    function getNumberPrecision(displayData: DisplayData): number {
        if (displayData == DisplayData.Reads) {
            return 0
        } else if (displayData == DisplayData.Rpm)  {
            return 2
        } else {
            return 4
        }
    }
</script>



<div class="grid grid-cols-4 h-full w-full">
    <!-- First Column: Heatmap -->
    <div id="Heatmap" bind:this={container} class="col-span-3 h-full p-8 relative">
        <svg
            bind:this={svg}
            width={width}
            height={height}
            class="bg-transparent"
        >
            <!-- X-axis Labels -->
            <g transform={`translate(${margin.left}, ${margin.top - 10})`}>
                {#each uniqueDetectionIds as id}
                    <text
                        x={xScale(id) + xScale.bandwidth() / 2}
                        y="-10"
                        transform={`rotate(-45, ${xScale(id) + xScale.bandwidth() / 2}, 0)`}
                        style="fill: rgb(var(--color-primary-400)); text-anchor: start; dominant-baseline: hanging; font-size: 0.7em;"
                    >
                        {id}
                    </text>
                {/each}
            </g>

            <!-- Y-axis Labels -->
            <g transform={`translate(${margin.left - 10}, ${margin.top})`}>
                {#each taxa as taxon}
                    <text
                        x="0"
                        y={yScale(taxon.name) + yScale.bandwidth() / 2}
                        style="fill: rgb(var(--color-primary-400)); text-anchor: end; dominant-baseline: middle; font-size: 0.7em;"
                    >
                        {taxon.name}
                    </text>
                {/each}
            </g>

            <!-- Heatmap Rectangles -->
            <g transform={`translate(${margin.left}, ${margin.top})`}>
                {#each dataMatrix as data}
                    <rect 
                        x={xScale(data.column)} 
                        y={yScale(data.row)} 
                        width={xScale.bandwidth()} 
                        height={yScale.bandwidth()} 
                        fill={data.value !== null ? rowColorScales[data.row](data.value) : "rgba(0, 0, 0, 0)"}
                        on:mouseover={() => handleMouseOver(data.row, data.column)}
                        on:mouseout={handleMouseOut}
                    ></rect>
                    {#if hoveredCell.row === data.row && hoveredCell.column === data.column}
                        <text
                            x={xScale(data.column) + xScale.bandwidth() / 2}
                            y={yScale(data.row) + yScale.bandwidth() / 2}
                            text-anchor="middle"
                            dominant-baseline="middle"
                            font-size={Math.min(xScale.bandwidth(), yScale.bandwidth()) / 3 + "px"}
                            fill="black"
                        >
                            {data.value ? data.value.toFixed(getNumberPrecision(displayData)) : ""}
                        </text>
                    {/if}
                {/each}
            </g>
        </svg>
    </div>

    <!-- Second Column: Controls -->
    <div class="col-span-1 h-full flex flex-col gap-4 p-4 pt-16">
        <label class="label">
            <span class="font-medium mb-2">Tool</span>
            <select bind:value={tool} class="select">
                <option value="{PathogenDetectionTool.Ganon2}">{PathogenDetectionTool.Ganon2}</option>
                <option value="{PathogenDetectionTool.Kraken2}">{PathogenDetectionTool.Kraken2}</option>
                <option value="{PathogenDetectionTool.Bracken}">{PathogenDetectionTool.Bracken}</option>
                <option value="{PathogenDetectionTool.Metabuli}">{PathogenDetectionTool.Metabuli}</option>
                <option value="{PathogenDetectionTool.Kmcp}">{PathogenDetectionTool.Kmcp}</option>
                <option value="{PathogenDetectionTool.Sylph}">{PathogenDetectionTool.Sylph}</option>
            </select>
        </label>

        <label class="label">
            <span class="font-medium mb-2">Mode</span>
            <select bind:value={mode} class="select">
                <option value="{PathogenDetectionMode.Sequence}">{PathogenDetectionMode.Sequence}</option>
                <option value="{PathogenDetectionMode.Profile}">{PathogenDetectionMode.Profile}</option>
            </select>
        </label>

        <label class="label">
            <span class="font-medium mb-2">Data</span>
            <select bind:value={displayData} class="select">
                <option value="{DisplayData.Rpm}">Reads per million (RPM)</option>
                <option value="{DisplayData.Reads}">Reads</option>
                <option value="{DisplayData.Abundance}">Abundance</option>
            </select>
        </label>
    </div>
</div>

<style>
    rect:hover {
        stroke: #000;
        stroke-width: 1px;
    }
</style>
