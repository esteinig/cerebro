<script lang="ts">
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import SpeciesOverviewTableHeader from "$lib/general/tables/SpeciesOverviewTableHeader.svelte";
	import { defaultPrefetchClientFilterConfig } from "$lib/stores/stores";
	import { transformTaxonOverview } from "$lib/utils/helpers";
	import { AbundanceMode, DisplayData, DisplayTotal, ProfileTool, type Taxon, type TaxonEvidence, type TaxonOverviewRecord } from "$lib/utils/types";
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";
	import { createEventDispatcher } from "svelte";

    export let taxa: Taxon[] | null = [];
    export let selectedName: string = "";

    export let displayMode: AbundanceMode = AbundanceMode.Mixed;
    export let displayData: DisplayData = DisplayData.Rpm;
    export let displayTotal: DisplayTotal = DisplayTotal.Average;

    export let collapseByGenus: boolean = false;
    export let baseWeight: number = 1.0;

    // Container for filtered table row data
    let tableData: TaxonOverviewRecord[] = [];
    

    let selectedTaxon: string = "";
    $: selectedTaxon = selectedName; 

    // Handler for clicking a header.
    function sortByColumn(column: keyof TaxonOverviewRecord) {
        if (sortColumn === column) {
            // Toggle direction if the same column is clicked.
            sortDirection = sortDirection === "asc" ? "desc" : "asc";
        } else {
            sortColumn = column;
            sortDirection = "desc"; // default direction when changing columns
        }
    }

    function sortTableData(sortColumn: keyof TaxonOverviewRecord, sortDirection: string) {

        return tableData.sort((a, b) => {
            const aVal = a[sortColumn];
            const bVal = b[sortColumn];
            if (typeof aVal === "number" && typeof bVal === "number") {
                return sortDirection === "asc" ? aVal - bVal : bVal - aVal;
            } else if (typeof aVal === "string" && typeof bVal === "string") {
                return sortDirection === "asc"
                    ? aVal.localeCompare(bVal)
                    : bVal.localeCompare(aVal);
            }
            return 0;
        });
    }

    // Number precision to display in table
    function getNumberPrecision(displayData: DisplayData): number {
        if (displayData == DisplayData.Reads || displayData == DisplayData.Bases) {
            return 0
        } else if (displayData == DisplayData.Rpm || displayData == DisplayData.Bpm)  {
            return 1
        } else {
            return 4
        }
    }
    let sortColumn: keyof TaxonOverviewRecord = "total";
    let sortDirection: "asc" | "desc" = "desc";

    $: {
        if (taxa) {
        tableData = transformTaxonOverview(taxa, displayMode, displayData, displayTotal);
        if (collapseByGenus) {
            tableData = collapseGenera(tableData, baseWeight);
        }
        } else {
        tableData = [];
        }
        tableData = sortTableData(sortColumn, sortDirection);
    }


    function selectTaxon(overview: TaxonOverviewRecord) {

        if (selectedTaxon === overview.name) {
            selectedTaxon = "";
            dispatch('select', { index: 0, item: null });
        } else {
            selectedTaxon = overview.name;
            dispatch('select', { index: 0, item: overview });
        }

    }

    const dispatch = createEventDispatcher<{
      select: { index: number; item: TaxonOverviewRecord | null; };
    }>();

    const RPM_FIELDS: (keyof TaxonOverviewRecord)[] = [
        'vircov','kraken2','bracken','metabuli','ganon2','kmcp','sylph'
    ];

    function profileScore(r: TaxonOverviewRecord, weight = 1.0): number {
        const totalRpm = RPM_FIELDS.reduce((s, k) => s + Math.max(0, Number(r[k] ?? 0)), 0);
        const bases     = Math.max(0, Number(r.blast ?? 0));
        const baseScore = bases > 0 ? Math.log10(bases) * weight : 0;
        return totalRpm + baseScore;
    }

    function collapseGenera(rows: TaxonOverviewRecord[], weight = 1.0): TaxonOverviewRecord[] {
        const allowed = new Set(['Eukaryota','Bacteria','Archaea']);
        const byGenus = new Map<string, TaxonOverviewRecord[]>();

        for (const r of rows) {
        if (!r.genus || !r.domain || !allowed.has(r.domain)) continue;
        const list = byGenus.get(r.genus) ?? [];
        list.push(r);
        byGenus.set(r.genus, list);
        }

        // genera with >= 5 species
        const collapseSet = new Set<string>();
        for (const [g, list] of byGenus) {
        if (list.length >= 5) collapseSet.add(g);
        }

        const bestByGenus = new Map<string, TaxonOverviewRecord>();
        for (const g of collapseSet) {
        const list = byGenus.get(g)!;
        let best = list[0];
        let bestScore = profileScore(best, weight);
        for (let i = 1; i < list.length; i++) {
            const s = profileScore(list[i], weight);
            if (s > bestScore) { best = list[i]; bestScore = s; }
        }
        bestByGenus.set(g, best);
        }

        // keep:
        // - all rows in genera not collapsed
        // - for collapsed genera, keep only the best species
        const out: TaxonOverviewRecord[] = [];
        for (const r of rows) {
        if (r.genus && collapseSet.has(r.genus)) {
            if (bestByGenus.get(r.genus) === r) out.push(r);
        } else {
            out.push(r);
        }
        }
        return out;
    }

</script>


<div>

    {#if tableData.length === 0}
        <p class="flex justify-center text-lg pb-4">No taxa available for this category</p>
    {:else}
        <div>
            <ListBox>
                <ListBoxItem group="header" name="header" value="qc" active='variant-soft' hover='hover:cursor-default' rounded='rounded-token'>
                    <div class="grid grid-cols-11 sm:grid-cols-11 md:grid-cols-11 gap-x-1 gap-y-4 w-full text-sm">
                        <div class="col-span-1 flex flex-col items-start">
                            <span class="opacity-60">Domain</span>
                        </div>
                        
                        <div class="col-span-2 flex flex-col items-start">
                            <span class="opacity-60">Species</span>
                        </div>
                                                    
                        <div on:click|preventDefault={() => sortByColumn("total")} role="columnheader">
                            <SpeciesOverviewTableHeader tool={displayTotal} sortColumn={sortColumn} displayData={displayData} circleColor={null} sortOrder={sortDirection}/>
                        </div>
                        
                        {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Vircov)}
                            <div on:click|preventDefault={() => sortByColumn("vircov")} role="columnheader">
                                <SpeciesOverviewTableHeader tool={ProfileTool.Vircov} displayData={displayData} sortColumn={sortColumn} circleColor="bg-primary-500" sortOrder={sortDirection}/>
                            </div>
                        {/if}
                        
                        {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Kraken2)}
                            <div on:click|preventDefault={() => sortByColumn("kraken2")} role="columnheader">
                                <SpeciesOverviewTableHeader tool={ProfileTool.Kraken2} displayData={displayData} sortColumn={sortColumn} circleColor="bg-secondary-400" sortOrder={sortDirection}/>
                            </div>
                        {/if}
                        
                        {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Bracken)}
                            <div on:click|preventDefault={() => sortByColumn("bracken")} role="columnheader">
                                <SpeciesOverviewTableHeader tool={ProfileTool.Bracken} displayData={displayData} sortColumn={sortColumn} circleColor="bg-secondary-500" sortOrder={sortDirection}/>
                            </div>
                        {/if}
                        
                        {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Metabuli)}
                            <div on:click|preventDefault={() => sortByColumn("metabuli")} role="columnheader">
                                <SpeciesOverviewTableHeader tool={ProfileTool.Metabuli} displayData={displayData} sortColumn={sortColumn} circleColor="bg-secondary-600" sortOrder={sortDirection}/>
                            </div>
                        {/if}
                        
                        {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Ganon2)}
                            <div on:click|preventDefault={() => sortByColumn("ganon2")} role="columnheader">
                                <SpeciesOverviewTableHeader tool={ProfileTool.Ganon2} displayData={displayData} sortColumn={sortColumn} circleColor="bg-secondary-700" sortOrder={sortDirection}/>
                            </div>
                        {/if}

                        {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Sylph)}
                            <div on:click|preventDefault={() => sortByColumn("sylph")} role="columnheader">
                                <SpeciesOverviewTableHeader tool={ProfileTool.Sylph} displayData={displayData} sortColumn={sortColumn} circleColor="bg-secondary-800" sortOrder={sortDirection}/>
                            </div>
                        {/if}
                        
                        {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Blast)}
                            <div on:click|preventDefault={() => sortByColumn("blast")} role="columnheader">
                                <SpeciesOverviewTableHeader tool={ProfileTool.Blast} displayData={displayData} sortColumn={sortColumn} circleColor="bg-tertiary-500" sortOrder={sortDirection}/>
                            </div>
                        {/if}


                        {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Blast)}
                            <div on:click|preventDefault={() => sortByColumn("contigs")} role="columnheader">
                                <SpeciesOverviewTableHeader tool={"Assembly"} displayData={displayData} sortColumn={sortColumn} circleColor="bg-tertiary-500" sortOrder={sortDirection}/>
                            </div>
                        {/if}

                        <div class="text-center">
                            <span class="opacity-60">Tools</span>
                        </div>

                        <!-- <div class="col-span-1">Tags</div> -->
                        <!-- <div class="text-right">Kmcp</div> -->
                    </div>
                </ListBoxItem>
                {#each tableData as overview, i}
                    <ListBoxItem bind:group={selectedTaxon} name={overview.name} value={overview.name} active='variant-ghost' rounded='rounded-token'  on:click={() => selectTaxon(overview)}> 
                            
                                <div class="grid grid-cols-11 sm:grid-cols-11 md:grid-cols-11 gap-x-1 gap-y-4 w-full text-sm}">
                                    
                                    <div class="opacity-70">{overview.domain}</div>
                                    <div class="col-span-2 truncate italic">{overview.name}</div>
                                    <div class="text-right">{overview.total > 0 ? overview.total.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                    
                                    {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Vircov)}
                                        <div class="text-right">{overview.vircov > 0 ? overview.vircov.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                    {/if}
                                    {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Kraken2)}
                                        <div class="text-right">{overview.kraken2 > 0 ? overview.kraken2.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                    {/if}
                                    {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Bracken)}
                                        <div class="text-right">{overview.bracken > 0 ? overview.bracken.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                    {/if}
                                    {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Metabuli)}
                                        <div class="text-right">{overview.metabuli > 0 ? overview.metabuli.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                    {/if}
                                    {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Ganon2)}
                                        <div class="text-right">{overview.ganon2 > 0 ? overview.ganon2.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                    {/if}
                                    {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Sylph)}
                                        <div class="text-right">{overview.sylph > 0 ? overview.sylph.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                    {/if}
                                    {#if defaultPrefetchClientFilterConfig.tools.includes(ProfileTool.Blast)}
                                        <div class="text-right">{overview.blast > 0 ? overview.blast.toFixed(getNumberPrecision(DisplayData.Bases)) : "-"}</div>
                                        <div class="text-right">{overview.contigs > 0 ? overview.contigs.toFixed(getNumberPrecision(DisplayData.Bases)) : "-"}</div>
                                    {/if}
                                        <!-- <div class="text-right">{overview.kmcp > 0 ? overview.kmcp.toFixed(getNumberPrecision(displayData)) : "-"}</div> -->
                                    <div class="flex justify-center items-center pt-1">

                                        <div class="grid grid-cols-8 sm:grid-cols-8 md:grid-cols-8 text-sm">

                                            {#if overview.vircov > 0}
                                                <div class="rounded-full bg-primary-400 h-2 w-2"></div>
                                            {:else}
                                                <div></div>
                                            {/if}
                                            {#if overview.kraken2 > 0}
                                                <div class="rounded-full bg-secondary-800 h-2 w-2"></div>
                                            {:else}
                                                <div></div>
                                            {/if}
                                            {#if overview.bracken > 0}
                                                <div class="rounded-full bg-secondary-700 h-2 w-2"></div>
                                            {:else}
                                                <div></div>
                                            {/if}
                                            {#if overview.metabuli > 0}
                                                <div class="rounded-full bg-secondary-600 h-2 w-2"></div>
                                            {:else}
                                                <div></div>
                                            {/if}
                                            {#if overview.ganon2 > 0}
                                                <div class="rounded-full bg-secondary-400 h-2 w-2"></div>
                                            {:else}
                                                <div></div>
                                            {/if}
                                            <!-- {#if overview.kmcp > 0}
                                                <div class="rounded-full bg-secondary-300 h-2 w-2"></div>
                                            {:else}
                                                <div></div>
                                            {/if} -->
                                            {#if overview.sylph > 0}
                                                <div class="rounded-full bg-secondary-200 h-2 w-2"></div>
                                            {:else}
                                                <div></div>
                                            {/if}
                                            {#if overview.blast > 0}
                                                <div class="rounded-full bg-tertiary-500 h-2 w-2"></div>
                                            {:else}
                                                <div></div>
                                            {/if}
                                        </div>
                                    </div>
                                    <div class="info-buttons flex items-center align-center justify-end gap-x-1">                                           
                                    </div>
                                </div>
                            </ListBoxItem>
                {/each}
            </ListBox>
        </div>
    {/if}
</div>

