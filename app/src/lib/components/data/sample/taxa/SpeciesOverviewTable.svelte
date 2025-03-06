<script lang="ts">
	import type { ClientFilterConfig, Taxon, TaxonOverviewRecord, TaxonEvidence } from "$lib/utils/types";
    import { DisplayData, DisplayTotal, ProfileTool } from "$lib/utils/types";
	import { getToastStore, ListBox, ListBoxItem, Paginator, ProgressRadial, type PaginationSettings } from "@skeletonlabs/skeleton";
	import { selectedTaxonHighlightConfig, selectedClientFilterConfig, selectedTaxa, selectedIdentifiers, selectedServerFilterConfig, selectedPrevalenceContamConfig, selectedModels } from "$lib/stores/stores";
    import { AbundanceMode, FileTag } from "$lib/utils/types";
	import TaxonEvidenceOverview from "./evidence/TaxonEvidenceOverview.svelte";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import { page } from "$app/stores";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import CircleIndicator from "$lib/general/icons/CircleIndicator.svelte";
	import { baseTags, extractTaxon, getValueDisplay, TaxRank } from "$lib/utils/helpers";
    import { navigationLoading } from '$lib/stores/stores';
	import SpeciesOverviewTableHeader from "$lib/general/tables/SpeciesOverviewTableHeader.svelte";

    export let displayData: DisplayData = DisplayData.Rpm;
    export let displayMode: AbundanceMode = AbundanceMode.Sequence;

    let pagination: boolean = true;

    // Selected taxonomic identifier
    let selectedTaxid: string;

    // Taxa returned server-side filtered
    let taxa: Taxon[] = [];

    // Taxa filtered client-side
    let filteredData: Taxon[] = [];

    // Contamination taxids
    let contamTaxid: string[] = [];

    // Container for filtered table row data
    let tableData: TaxonOverviewRecord[] = [];

    // Display data modes for table
    let displayTotal: DisplayTotal = DisplayTotal.Average;

    let taxonEvidence: TaxonEvidence | null = null;

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

    

    let publicApi = new CerebroApi();
    let toastStore = getToastStore();

    let loading: boolean = false;

    async function getAggregatedTaxaOverview(selectedIdentifiers: string[]) {

        // Clear the contamination taxonomic identifiers
        contamTaxid = []

        loading = true;
        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.taxa}?team=${$page.params.team}&db=${$page.params.db}&project=${$page.params.project}&id=${selectedIdentifiers.join(",")}&overview=true`,
            { 
                method: 'POST',  
                mode: 'cors',
                credentials: 'include', 
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify($selectedServerFilterConfig) 
            } as RequestInit,
            $page.data.refreshToken, toastStore, "Taxa loaded"
        )

        loading = false;

        if (response.ok){
            taxa = response.json.data.taxa;

            let tags = baseTags($selectedModels.map(model => model.sample.tags), true);

            await getPrevalenceContamination(taxa.map(taxon => taxon.taxid), tags, $selectedPrevalenceContamConfig.threshold, $selectedPrevalenceContamConfig.min_rpm);
        }
    }    

    async function getPrevalenceContamination(taxid: string[], tags: string[], threshold: number, min_rpm: number) {

        $navigationLoading = true; // use the navigation loading for background loading status

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.taxaPrevalenceContamination}?team=${$page.params.team}&db=${$page.params.db}&project=${$page.params.project}`,
            { 
                method: 'POST',  
                mode: 'cors',
                credentials: 'include', 
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify({
                    taxid: taxid,
                    tags: tags,
                    threshold: threshold,
                    min_rpm: min_rpm,
                    sample_type: null
                }) 
            } as RequestInit,
            $page.data.refreshToken, toastStore, "Contaminants identified from project libraries"
        )

        $navigationLoading = false;

        if (response.ok){
            contamTaxid = response.json.data.taxid;
        }
    }    

    $: if ($selectedIdentifiers.length > 0) {
        getAggregatedTaxaOverview($selectedIdentifiers);
    }
    function transformTaxonOverview(
        taxa: Taxon[],
        mode: AbundanceMode,
        field: DisplayData
    ): TaxonOverviewRecord[] {
        // Mapping for Mixed mode: tool -> expected AbundanceMode.
        const mixedToolModeMapping: Record<string, AbundanceMode> = {
            [ProfileTool.Vircov]: AbundanceMode.Sequence,
            [ProfileTool.Kraken2]: AbundanceMode.Sequence,
            [ProfileTool.Metabuli]: AbundanceMode.Sequence,
            [ProfileTool.Ganon2]: AbundanceMode.Sequence,
            [ProfileTool.Kmcp]: AbundanceMode.Sequence,
            [ProfileTool.Bracken]: AbundanceMode.Profile,
            [ProfileTool.Sylph]: AbundanceMode.Sequence,
            [ProfileTool.Blast]: AbundanceMode.Bases,
        };

        return taxa.map((taxon) => {
            const aggregatedResults = {
                vircov: 0,
                kraken2: 0,
                metabuli: 0,
                ganon2: 0,
                kmcp: 0,
                bracken: 0,
                sylph: 0,
                blast: 0,
                contigs: 0, // New aggregated field for BLAST contigs.
            };

            let contributingTools = 0; // Track the number of tools contributing to the average

            taxon.evidence.profile.forEach((record) => {
                // Only include records from selected tools.
                if (!$selectedClientFilterConfig.tools.includes(record.tool)) {
                    return;
                }

                // Determine the data field to use. For Blast, always use Bases.
                let data = field;
                if (record.tool === ProfileTool.Blast) {
                    data = DisplayData.Bases;
                }

                // Mixed mode: include only records with the expected mode for each tool.
                if (mode === AbundanceMode.Mixed) {
                    const expectedMode = mixedToolModeMapping[record.tool];
                    if (!expectedMode) {
                        return;
                    }
                    if (record.mode !== expectedMode) {
                        return;
                    }
                } else {
                    // Non-mixed mode: only include records matching the provided mode.
                    if (record.mode !== mode) {
                        return;
                    }
                }

                // Aggregate data.
                const key = record.tool.toLowerCase() as keyof typeof aggregatedResults;
                if (key in aggregatedResults) {
                    aggregatedResults[key] += record[data];

                    // If the tool is Blast, also aggregate the contigs.
                    if (record.tool === ProfileTool.Blast) {
                        // Ensure that the record has a contigs property.
                        aggregatedResults.contigs += record.contigs;
                    }

                    // Exclude Blast contributions when counting tools for average.
                    if (record[data] > 0 && record.tool !== ProfileTool.Blast) {
                        contributingTools++;
                    }
                }
            });

            // Calculate the total sum, excluding 'blast' and 'contigs'.
            const totalSum = Object.entries(aggregatedResults).reduce((sum, [key, value]) => {
                return (key === "blast" || key === "contigs") ? sum : sum + value;
            }, 0);

            // Calculate the average (if at least one tool contributed).
            const average = contributingTools > 0 ? totalSum / contributingTools : 0;

            return {
                taxid: taxon.taxid,
                name: taxon.name,
                domain: extractTaxon(taxon.lineage, TaxRank.Domain),
                genus: extractTaxon(taxon.lineage, TaxRank.Genus),
                profile: taxon.evidence.profile.length > 0,
                alignment: taxon.evidence.alignment.length > 0,
                assembly: taxon.evidence.assembly.length > 0,
                ...aggregatedResults,
                total: displayTotal === DisplayTotal.Sum ? totalSum : average, // Use the adjusted total
            };
        }).sort((a, b) => b.total - a.total)
    }

    const applyClientSideFilters = (taxa: Taxon[], clientFilterConfig: ClientFilterConfig | null): Taxon[] => {
        
        if (clientFilterConfig === null) {
            return taxa
        }

        return taxa.filter(taxon => {

            // Filter decision
            let taxonomy: boolean = true;
            let modules: boolean = true;
            let contam: boolean = true;
            let syndrome: boolean = true;

            // Explicit taxonomy filters
            if ($selectedClientFilterConfig.domains.length) {
                taxonomy = $selectedClientFilterConfig.domains.includes(extractTaxon(taxon.lineage, TaxRank.Domain));
            }
            if ($selectedClientFilterConfig.genera.length) {
                taxonomy = $selectedClientFilterConfig.genera.includes(extractTaxon(taxon.lineage, TaxRank.Genus));
            }
            if ($selectedClientFilterConfig.species.length) {
                taxonomy = $selectedClientFilterConfig.species.includes(extractTaxon(taxon.lineage, TaxRank.Species));
            }

            // if ($selectedTaxonHighlightConfig.contamination.species.some(species => taxon.name.includes(species))) {
            //     contam = false
            // } else  {
            //     contam = true
            // }

            // if ($selectedTaxonHighlightConfig.syndrome.species.some(species => taxon.name.includes(species))) {
            //     syndrome = true
            // } else {
            //     syndrome = false
            // }
                
            return (taxonomy && modules) // && metrics  -> || syndrome will always show syndromic ones even if they should be client filtered
        })
    }

    let paginationSettings: PaginationSettings = {
        page: 0,
        limit: 500,
        size: taxa.length,
        amounts: [5, 10, 20, 50, 100, 500],
    };

    $: {    
        // Whenever data changes apply client side filters...
        filteredData = applyClientSideFilters(taxa, $selectedClientFilterConfig);

        // ... and transform filtered data into the overview table rows
        tableData = transformTaxonOverview(filteredData, displayMode, displayData);
    }


    // Sorting state variables.
    let sortColumn: keyof TaxonOverviewRecord = "total";
    let sortDirection: "asc" | "desc" = "desc";

    $: { 
        // Reset to first page if paginated


        // Sort the transformed data based on sortColumn and sortDirection.
        tableData = sortTableData(sortColumn, sortDirection);

        // If pagination is enabled, slice the data for the page and limits as configured in table footer
        if (pagination) {
            paginationSettings.size = filteredData.length;
            tableData = tableData.slice(
                paginationSettings.page * paginationSettings.limit,
                paginationSettings.page * paginationSettings.limit + paginationSettings.limit
            );
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

    // Apply pagination to sorted data.
    // $: {
    //     if (pagination) {
    //         paginationSettings.size = tableData.length;
    //         tableData = tableData.slice(
    //             paginationSettings.page * paginationSettings.limit,
    //             paginationSettings.page * paginationSettings.limit + paginationSettings.limit
    //         );
    //     } else {
    //         tableData = tableData;
    //     }
    // }

    const getTaxonBackgroundColor = (overview: TaxonOverviewRecord): string =>  {
        if ($selectedTaxonHighlightConfig.contamination.species.some(species => overview.name.toLowerCase().includes(species.toLowerCase()))) {
            return 'variant-soft-secondary rounded-token py-1.5 px-2'
        } else if ($selectedTaxonHighlightConfig.syndrome.species.some(species => overview.name.toLowerCase().includes(species.toLowerCase()))) {
            return 'variant-soft-tertiary rounded-token py-1.5 px-2'
        } else {
            return 'rounded-token py-1.5 px-2'
        }
    }


    const addSelectedTaxon = (overview: TaxonOverviewRecord) => {
        selectedTaxa.update(currentTaxa => {
            const index = currentTaxa.findIndex(taxon => taxon.taxid === overview.taxid);
            if (index > -1) {
                taxonEvidence = null;
                // If it exists, remove it
                return currentTaxa.filter((_, i) => i !== index);
            } else {
                const taxonRecord = filteredData.find(taxon => taxon.taxid === overview.taxid)
                taxonEvidence = taxonRecord ? taxonRecord.evidence : null;
                // If it doesn't exist, add it
                return [...currentTaxa, overview];
            }
        });
    };


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

</script>

<div>

    {#if loading}
        <div class="flex justify-center py-24">
            <ProgressRadial width="sm:w-12 md:w-24" stroke={20} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
        </div>
    {:else}

        {#if !taxa.length}
            <div class="flex justify-center py-16 "><ErrorAnimation /></div>
            <p class="flex justify-center text-lg pb-4">No taxa available</p>
        {:else}
            <div>
                <ListBox>
                    <ListBoxItem group="header" name="header" value="qc" active='variant-soft' hover='hover:cursor-default' rounded='rounded-token'>
                        <div class="grid grid-cols-12 sm:grid-cols-12 md:grid-cols-12 gap-x-1 gap-y-4 w-full text-sm">
                            <div class="col-span-1 flex flex-col items-start">
                                <span class="opacity-60">Domain</span>
                            </div>
                            
                            <div class="col-span-2 flex flex-col items-start">
                                <span class="opacity-60">Species</span>
                            </div>
                                                        
                            <div on:click|preventDefault={() => sortByColumn("total")} role="columnheader">
                                <SpeciesOverviewTableHeader tool={displayTotal} sortColumn={sortColumn} displayData={displayData} circleColor={null} sortOrder={sortDirection}/>
                            </div>
                            
                            {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Vircov)}
                                <div on:click|preventDefault={() => sortByColumn("vircov")} role="columnheader">
                                    <SpeciesOverviewTableHeader tool={ProfileTool.Vircov} displayData={displayData} sortColumn={sortColumn} circleColor="bg-primary-500" sortOrder={sortDirection}/>
                                </div>
                            {/if}
                            
                            {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Kraken2)}
                                <div on:click|preventDefault={() => sortByColumn("kraken2")} role="columnheader">
                                    <SpeciesOverviewTableHeader tool={ProfileTool.Kraken2} displayData={displayData} sortColumn={sortColumn} circleColor="bg-secondary-400" sortOrder={sortDirection}/>
                                </div>
                            {/if}
                            
                            {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Bracken)}
                                <div on:click|preventDefault={() => sortByColumn("bracken")} role="columnheader">
                                    <SpeciesOverviewTableHeader tool={ProfileTool.Bracken} displayData={displayData} sortColumn={sortColumn} circleColor="bg-secondary-500" sortOrder={sortDirection}/>
                                </div>
                            {/if}
                            
                            {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Metabuli)}
                                <div on:click|preventDefault={() => sortByColumn("metabuli")} role="columnheader">
                                    <SpeciesOverviewTableHeader tool={ProfileTool.Metabuli} displayData={displayData} sortColumn={sortColumn} circleColor="bg-secondary-600" sortOrder={sortDirection}/>
                                </div>
                            {/if}
                            
                            {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Ganon2)}
                                <div on:click|preventDefault={() => sortByColumn("ganon2")} role="columnheader">
                                    <SpeciesOverviewTableHeader tool={ProfileTool.Ganon2} displayData={displayData} sortColumn={sortColumn} circleColor="bg-secondary-700" sortOrder={sortDirection}/>
                                </div>
                            {/if}

                            {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Sylph)}
                                <div on:click|preventDefault={() => sortByColumn("sylph")} role="columnheader">
                                    <SpeciesOverviewTableHeader tool={ProfileTool.Sylph} displayData={displayData} sortColumn={sortColumn} circleColor="bg-secondary-800" sortOrder={sortDirection}/>
                                </div>
                            {/if}
                            
                            {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Blast)}
                                <div on:click|preventDefault={() => sortByColumn("blast")} role="columnheader">
                                    <SpeciesOverviewTableHeader tool={ProfileTool.Blast} displayData={displayData} sortColumn={sortColumn} circleColor="bg-tertiary-500" sortOrder={sortDirection}/>
                                </div>
                            {/if}

                            {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Blast)}
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
                        {#if contamTaxid.includes(overview.taxid) ? $selectedClientFilterConfig.contam.display : true}
                            <ListBoxItem bind:group={selectedTaxid} name={overview.name} value={overview.taxid} active='' hover="hover:variant-soft" regionDefault={getTaxonBackgroundColor(overview)} rounded='rounded-token' on:click={() => addSelectedTaxon(overview)}> 
                                
                                    <div class="grid grid-cols-12 sm:grid-cols-12 md:grid-cols-12 gap-x-1 gap-y-4 w-full text-sm {contamTaxid.includes(overview.taxid) ? 'opacity-20' : ''}">
                                        
                                        <div class="opacity-70">{overview.domain}</div>
                                        <div class="col-span-2 truncate italic">{overview.name}</div>
                                        <div class="text-right">{overview.total > 0 ? overview.total.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                        
                                        {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Vircov)}
                                            <div class="text-right">{overview.vircov > 0 ? overview.vircov.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                        {/if}
                                        {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Kraken2)}
                                            <div class="text-right">{overview.kraken2 > 0 ? overview.kraken2.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                        {/if}
                                        {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Bracken)}
                                            <div class="text-right">{overview.bracken > 0 ? overview.bracken.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                        {/if}
                                        {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Metabuli)}
                                            <div class="text-right">{overview.metabuli > 0 ? overview.metabuli.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                        {/if}
                                        {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Ganon2)}
                                            <div class="text-right">{overview.ganon2 > 0 ? overview.ganon2.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                        {/if}
                                        {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Sylph)}
                                            <div class="text-right">{overview.sylph > 0 ? overview.sylph.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                                        {/if}
                                        {#if $selectedClientFilterConfig.tools.includes(ProfileTool.Blast)}
                                            <div class="text-right">{overview.blast > 0 ? overview.blast.toFixed(getNumberPrecision(DisplayData.Bases)) : "-"}</div>
                                            <div class="text-right">{overview.contigs > 0 ? overview.contigs.toFixed(getNumberPrecision(DisplayData.Bases)) : "-"}</div>
                                        {/if}
                                            <!-- <div class="text-right">{overview.kmcp > 0 ? overview.kmcp.toFixed(getNumberPrecision(displayData)) : "-"}</div> -->
                                        <div class="flex justify-end items-center pt-1">

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
                                    </div>
                                {#if selectedTaxid === overview.taxid && taxonEvidence}
                                    <TaxonEvidenceOverview taxonEvidence={taxonEvidence}></TaxonEvidenceOverview>
                                {/if}
                            </ListBoxItem>
                        {/if}
                    {/each}
                </ListBox>
            </div>
            {#if pagination}
                <div class="mt-8">
                    <Paginator bind:settings={paginationSettings} showFirstLastButtons={false} showPreviousNextButtons={true} class="mt-2" select="paginator-select select text-xs" regionControl="opacity-30 dark:variant-filled-surface dark:opacity-60 text-xs"/>
                </div>
            {/if}
        {/if}
    {/if}
</div>