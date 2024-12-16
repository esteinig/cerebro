<script lang="ts">
    import * as d3 from 'd3';
    import { selectedTaxa, selectedServerFilterConfig, storeTheme, navigationLoading } from '$lib/stores/stores';
    import { type Taxon, type TaxonOverviewRecord, PathogenDetectionTool, AbundanceMode, DisplayData, HeatmapRowOrder, HeatmapColorScheme, DomainName, FileTag, type Cerebro } from '$lib/utils/types';
    import CerebroApi, { ApiResponse } from "$lib/utils/api";
    import { page } from '$app/stores';
    import { getCssVariableAsHex } from '$lib/utils/helpers';

    const publicApi = new CerebroApi();

    export let selectedModels: Cerebro[] = [];
    export let selectedIdentifiers: string[] = [];
    export let width: number = 1024;
    export let height: number = 768;
    export let tool: PathogenDetectionTool = PathogenDetectionTool.Ganon2;
    export let mode: AbundanceMode = AbundanceMode.Sequence;
    export let displayData: DisplayData = DisplayData.Rpm;

    let container: HTMLDivElement;
    let svg: SVGSVGElement;
    let g: SVGGElement;

    let taxa: Taxon[] = [];
    let uniqueDetectionIds: string[] = [];
    let dataMatrix: Array<{ row: string; column: string; value: number | null }> = [];
    let hoveredCell = { row: null, column: null };
    let loading = false;

    let selectedRowOrder: string = HeatmapRowOrder.Domain;
    let colorScheme: string = HeatmapColorScheme.Domain;

    let rowOrder: string[] = []; // Tracks the rows in the heatmap
    let columnOrder: string[] = [];

    let domainColors: Map<string, {start: string, end: string}> = new Map([
        [DomainName.Viruses, { start: '--color-tertiary-400', end: '--color-tertiary-900'}],
        [DomainName.Archaea, { start: '--color-primary-100', end: '--color-primary-900'}],
        [DomainName.Bacteria, { start: '--color-primary-100', end: '--color-primary-900'}],
        [DomainName.Eukaryota, { start: '--color-secondary-100', end: '--color-secondary-900'}],
    ])   

    const getAggregatedTaxa = async (selectedIdentifiers: string[], selectedTaxa: TaxonOverviewRecord[], dateRange: [string, string] | null = null) => {
        if (!selectedTaxa || selectedTaxa.length === 0) return;

        navigationLoading.set(true);

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
        
        navigationLoading.set(false);
        
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

    // Sort taxa by domain name if the toggle is active
    $: {
        if (selectedRowOrder === HeatmapRowOrder.Domain) {
            taxa.sort((a, b) => {
                const domainA = a.level?.domain_name || ""; // Default to an empty string if undefined
                const domainB = b.level?.domain_name || ""; // Default to an empty string if undefined
                return domainA.localeCompare(domainB);
            });
            rowOrder = taxa.map(taxon => taxon.name);
        }
    }



    // Derive unique Detection IDs (columns)
    $: uniqueDetectionIds = Array.from(new Set(taxa.flatMap((taxon) => taxon.evidence.records.map((record) => record.id))));

    // Preserve order and append new columns
    // $: columnOrder = [...columnOrder, ...uniqueDetectionIds.filter((id) => !columnOrder.includes(id))];

    $: {
        
        // Trigger so the order is dynamically computed when taxa are selected

        if ($selectedTaxa && $selectedTaxa.length < 1) {
            columnOrder = columnOrder
        }
        // Map column identifiers to their corresponding tags from the selectedModels
        const columnTags = new Map<string, FileTag>();
        selectedModels.forEach((model: Cerebro) => {
            const tags = model.sample.tags;
            tags.forEach((tag: any) => {
                if (Object.values(FileTag).includes(tag as FileTag)) {
                    columnTags.set(model.name, tag as FileTag);
                }
            });
        });

        // Filter uniqueDetectionIds into groups
        const ntcAndEnvColumns = uniqueDetectionIds.filter(
            (id) => columnTags.get(id) === FileTag.NTC || columnTags.get(id) === FileTag.ENV
        );
        const otherColumns = uniqueDetectionIds.filter(
            (id) => !ntcAndEnvColumns.includes(id)
        );

        // Merge with existing columnOrder, ensuring all columns are included
        const updatedColumnOrder = [
            ...new Set([
                ...ntcAndEnvColumns,
                ...columnOrder.filter((id) => !ntcAndEnvColumns.includes(id) && !otherColumns.includes(id)),
                ...otherColumns,
            ]),
        ];

        // Insert a "gap" placeholder between the groups
        const ntcAndEnvIndex = updatedColumnOrder.findIndex((id) => ntcAndEnvColumns.includes(id));
        const lastNtcEnvIndex = updatedColumnOrder.lastIndexOf(ntcAndEnvColumns[ntcAndEnvColumns.length - 1]);
        if (ntcAndEnvIndex !== -1 && lastNtcEnvIndex !== -1 && lastNtcEnvIndex + 1 < updatedColumnOrder.length) {
            updatedColumnOrder.splice(lastNtcEnvIndex + 1, 0, "gap");
        }

        columnOrder = updatedColumnOrder;
    }

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
        // acc[taxon.name] = Math.max(...dataMatrix.filter((d) => d.row === taxon.name).map((d) => d.value || 0));
        acc[taxon.name] = Math.max(...dataMatrix.map((d) => d.value || 0));
        return acc;
    }, {} as Record<string, number>);

    $: rowColorScales = Object.fromEntries(
        Object.entries(rowMaxValues).map(([row, maxValue]) => {
            
            let startColor = getCssVariableAsHex('--color-primary-200', $storeTheme);
            let endColor = getCssVariableAsHex('--color-primary-600', $storeTheme);

            if (colorScheme === HeatmapColorScheme.Uniform) {
                startColor = getCssVariableAsHex('--color-primary-200', $storeTheme);
                endColor = getCssVariableAsHex('--color-primary-600', $storeTheme);
            } else if (colorScheme === HeatmapColorScheme.Domain) {
                taxa.filter((taxon) => row == taxon.name).map((taxon) => {
                    if (taxon.level.domain_name !== undefined) {
                        startColor = getCssVariableAsHex(domainColors.get(taxon.level.domain_name)?.start, $storeTheme);
                        endColor = getCssVariableAsHex(domainColors.get(taxon.level.domain_name)?.end, $storeTheme);
                    } 
                })

            }

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

    const handleMouseOver = (row: any, column: any) => {
        hoveredCell = { row, column };
    };

    const handleMouseOut = () => {
        hoveredCell = { row: null, column: null };
    };

    function getNumberPrecision(displayData: DisplayData): number {
        if (displayData == DisplayData.Reads || displayData == DisplayData.Bases) {
            return 0
        } else if (displayData == DisplayData.Rpm || displayData == DisplayData.Bpm)  {
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
                {#each rowOrder as row}
                    {#if taxa.find(taxon => taxon.name === row)}
                        <text
                            x="0"
                            y={yScale(row) + yScale.bandwidth() / 2}
                            style="fill: rgb(var(--color-primary-400)); text-anchor: end; dominant-baseline: middle; font-size: 0.7em;"
                        >
                            {row}
                        </text>
                    {/if}
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
                        on:focus={() => handleMouseOver(data.row, data.column)}
                        on:blur={handleMouseOut}
                        role="figure"
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
            <span class="font-medium mb-2">Classifier</span>
            <select bind:value={tool} class="select">
                <option value="{PathogenDetectionTool.Kraken2}">{PathogenDetectionTool.Kraken2}</option>
                <option value="{PathogenDetectionTool.Ganon2}">{PathogenDetectionTool.Ganon2}</option>
                <option value="{PathogenDetectionTool.Bracken}">{PathogenDetectionTool.Bracken}</option>
                <option value="{PathogenDetectionTool.Metabuli}">{PathogenDetectionTool.Metabuli}</option>
                <option value="{PathogenDetectionTool.Kmcp}">{PathogenDetectionTool.Kmcp}</option>
                <option value="{PathogenDetectionTool.Sylph}">{PathogenDetectionTool.Sylph}</option>
                <option value="{PathogenDetectionTool.Vircov}">{PathogenDetectionTool.Vircov}</option>
                <option value="{PathogenDetectionTool.BlastContig}">{PathogenDetectionTool.BlastContig}</option>
            </select>
        </label>

        <label class="label">
            <span class="font-medium mb-2">Abundance</span>
            <select bind:value={mode} class="select">
                <option value="{AbundanceMode.Sequence}">{AbundanceMode.Sequence}</option>
                <option value="{AbundanceMode.Profile}">{AbundanceMode.Profile}</option>
                <option value="{AbundanceMode.Bases}">{AbundanceMode.Bases}</option>
            </select>
        </label>

        <label class="label">
            <span class="font-medium mb-2">Data values</span>
            <select bind:value={displayData} class="select">
                <option value="{DisplayData.Rpm}">Reads per million (RPM)</option>
                <option value="{DisplayData.Bpm}">Bases per million (Bpm)</option>
                <option value="{DisplayData.Reads}">Reads</option>
                <option value="{DisplayData.Bases}">Bases</option>
                <option value="{DisplayData.Abundance}">Abundance</option>
            </select>
        </label>

        <label class="label">
            <span class="font-medium mb-2">Color scheme</span>
            <select bind:value={colorScheme} class="select">
                <option value="{HeatmapColorScheme.Domain}">Rank: Domain (Superkingdom)</option>
                <option value="{HeatmapColorScheme.Uniform}">Uniform</option>
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
