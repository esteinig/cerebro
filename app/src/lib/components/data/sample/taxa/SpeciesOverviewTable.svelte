<script lang="ts">
	import type { ClientFilterConfig, TaxonHighlightConfig, TaxonOverview, TaxonOverviewRecord } from "$lib/utils/types";
    import { DisplayData, DisplayTotal } from "$lib/utils/types";
	import { ListBox, ListBoxItem, Paginator, type PaginationSettings } from "@skeletonlabs/skeleton";
	import { selectedTaxonHighlightConfig, selectedClientFilterConfig, selectedTaxa } from "$lib/stores/stores";
    import { PathogenDetectionMode } from "$lib/utils/types";

    // export let serverFilterConfig: CerebroFilterConfig | CerebroFilterConfig[];
    
    // export let candidateButton: boolean = true;

    export let pagination: boolean = true;
    export let taxonOverview: TaxonOverview[] = [];
    export let modelNameTags: Map<string, string[]> = new Map();

    let displayMode: PathogenDetectionMode = PathogenDetectionMode.Sequence;
    let displayData: DisplayData = DisplayData.Rpm;
    let displayTotal: DisplayTotal = DisplayTotal.Average;


    let selectedTaxid: string;

    function getNumberPrecision(displayData: DisplayData): number {
        if (displayData == DisplayData.Reads) {
            return 0
        } else if (displayData == DisplayData.Rpm)  {
            return 2
        } else {
            return 4
        }
    }
    
    function transformTaxonOverview(
        overviews: TaxonOverview[],
        mode: PathogenDetectionMode,
        field: DisplayData
    ): TaxonOverviewRecord[] {

        return overviews.map((overview) => {
            const aggregatedResults = {
                vircov: 0,
                kraken2: 0,
                metabuli: 0,
                ganon2: 0,
                kmcp: 0,
                bracken: 0,
                sylph: 0,
            };

            let contributingTools = 0; // Track the number of tools contributing to the average

            overview.evidence.forEach((result) => {
                if (result.mode === mode) {
                    const key = result.tool.toLowerCase() as keyof typeof aggregatedResults;
                    if (key in aggregatedResults) {
                        if (result[field] > 0) {
                            contributingTools++;
                        }
                        aggregatedResults[key] += result[field];
                    }
                }
            });

            // Calculate the average
            const totalSum = Object.values(aggregatedResults).reduce((sum, value) => sum + value, 0);
            const average = contributingTools > 0 ? totalSum / contributingTools : 0;

            return {
                taxid: overview.taxid,
                name: overview.name,
                domain: overview.domain,
                genus: overview.genus,
                sample_names: overview.sample_names,
                kmer: overview.kmer,
                alignment: overview.alignment,
                assembly: overview.assembly,
                ...aggregatedResults,
                total: displayTotal == DisplayTotal.Sum ? totalSum : average,
            };
        }).sort((a, b) => b.total - a.total);
    }

    let filteredData: TaxonOverview[] = [];
    let tableData: TaxonOverviewRecord[] = transformTaxonOverview(
        filteredData, 
        displayMode, 
        displayData
    );

    const isTagData = (item: string[] | undefined): item is string[] => { return !!item }
    
    const getSelectedTags = (names: string[]) => {
        return names.map(name => modelNameTags.get(name)).filter(isTagData);
    }

    const applyClientSideFilters = (taxondOverview: TaxonOverview[], clientFilterConfig: ClientFilterConfig | null): TaxonOverview[] => {
        
        if (clientFilterConfig === null) {
            return taxonOverview
        }

        return taxonOverview.filter(taxonOverview => {

            // Filter decision
            let taxonomy: boolean = true;
            let modules: boolean = true;
            let contam: boolean = true;
            let syndrome: boolean = true;

            // Explicit taxonomy filters
            if ($selectedClientFilterConfig.domains.length) {
                taxonomy = $selectedClientFilterConfig.domains.includes(taxonOverview.domain);
            }
            if ($selectedClientFilterConfig.genera.length) {
                taxonomy = $selectedClientFilterConfig.genera.includes(taxonOverview.genus);
            }
            if ($selectedClientFilterConfig.species.length) {
                taxonomy = $selectedClientFilterConfig.species.includes(taxonOverview.name);
            }

            if ($selectedTaxonHighlightConfig.contamination.species.some(species => taxonOverview.name.includes(species))) {
                contam = false
            } else  {
                contam = true
            }

            if ($selectedTaxonHighlightConfig.syndrome.species.some(species => taxonOverview.name.includes(species))) {
                syndrome = true
            } else {
                syndrome = false
            }
                
            return (taxonomy && modules && contam) || syndrome  // && metrics
        })
    }

    let paginationSettings: PaginationSettings = {
        page: 0,
        limit: 50,
        size: taxonOverview.length,
        amounts: [5, 10, 20, 50, 100, 500],
    };

    $: {    
        filteredData = applyClientSideFilters(filteredData, $selectedClientFilterConfig);
        tableData = transformTaxonOverview(filteredData, displayMode, displayData);

        if (pagination) {
            paginationSettings.size = filteredData.length;
            tableData = tableData.slice(
                paginationSettings.page * paginationSettings.limit,
                paginationSettings.page * paginationSettings.limit + paginationSettings.limit
            );
        }
       
    }

    const getTaxonBackgroundColor = (overview: TaxonOverviewRecord): string =>  {
        if ($selectedTaxonHighlightConfig.contamination.species.some(species => overview.name.includes(species))) {
            return 'variant-soft-secondary rounded-token py-1.5 px-2'
        } else if ($selectedTaxonHighlightConfig.syndrome.species.some(species => overview.name.includes(species))) {
            return 'variant-soft-tertiary rounded-token py-1.5 px-2'
        } else {
            return 'rounded-token py-1.5 px-2'
        }
    }

    const getTaxonHover = (overview: TaxonOverviewRecord): string =>  {
        if ($selectedTaxonHighlightConfig.syndrome.species.some(species => overview.name.includes(species))){
            return 'hover:variant-soft-secondary'
        } else if ($selectedTaxonHighlightConfig.syndrome.species.some(species => overview.name.includes(species))) {
            return 'hover:variant-soft-tertiary'
        } else {
            return 'hover:variant-soft'
        }
    }


    export const addSelectedTaxon = (overview: TaxonOverviewRecord) => {
        selectedTaxa.update(currentTaxa => {
            const index = currentTaxa.findIndex(taxon => taxon.taxid === overview.taxid);
            if (index > -1) {
                // If it exists, remove it
                return currentTaxa.filter((_, i) => i !== index);
            } else {
                // If it doesn't exist, add it
                return [...currentTaxa, overview];
            }
        });
    };

</script>

<div>
    <div>
        <ListBox>
            <ListBoxItem group="header" name="header" value="qc" active='variant-soft' hover='hover:cursor-default' rounded='rounded-token'>
                <div class="grid grid-cols-12 sm:grid-cols-12 md:grid-cols-12 gap-x-1 gap-y-4 w-full text-sm opacity-60">
                    <div class="col-span-1">Domain</div>
                    <div class="col-span-2">Species</div>
                    <div class="col-span-1">Tags</div>
                    <div class="text-right">{displayTotal}</div>
                    <div class="text-right">Kraken2</div>
                    <div class="text-right">Bracken</div>
                    <div class="text-right">Metabuli</div>
                    <div class="text-right">Ganon2</div>
                    <div class="text-right">Vircov</div>
                    
                    <!-- <div class="text-right">Kmcp</div>
                    <div class="text-right">Sylph</div> -->
                    <div class="text-right">Modules</div>
                </div>
            </ListBoxItem>
            {#each tableData as overview, i}
                <ListBoxItem bind:group={selectedTaxid} name={overview.name} value={overview.taxid} active='' hover={getTaxonHover(overview)} regionDefault={getTaxonBackgroundColor(overview)} rounded='rounded-token' on:click={() => addSelectedTaxon(overview)}> 
                    
                    <div class="grid grid-cols-12 sm:grid-cols-12 md:grid-cols-12 gap-x-1 gap-y-4 w-full text-sm">
                        <div class="col-span-1 opacity-70">{overview.domain}</div>
                        <div class="col-span-2 truncate italic">{overview.name}</div>
                        <div class="col-span-1 truncate">
                            {#each getSelectedTags(overview.sample_names) as tags}
                                <span class="code bg-primary-500/30 text-primary-700 dark:bg-primary-500/20 dark:text-primary-400 mr-1">{tags.join("-")}</span>
                            {/each}
                        </div>
                        <div class="text-right">{overview.total > 0 ? overview.total.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                        <div class="text-right">{overview.kraken2 > 0 ? overview.kraken2.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                        <div class="text-right">{overview.bracken > 0 ? overview.bracken.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                        <div class="text-right">{overview.metabuli > 0 ? overview.metabuli.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                        <div class="text-right">{overview.ganon2 > 0 ? overview.ganon2.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                        <div class="text-right">{overview.vircov > 0 ? overview.vircov.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                        <!-- <div class="text-right">{overview.kmcp > 0 ? overview.kmcp.toFixed(getNumberPrecision(displayData)) : "-"}</div>
                        <div class="text-right">{overview.sylph > 0 ? overview.sylph.toFixed(getNumberPrecision(displayData)) : "-"}</div> -->
                        <div class="flex justify-end gap-x-1 items-center align-center pt-1">

                            <div class="grid grid-cols-7 sm:grid-cols-7 md:grid-cols-7 gap-x-1 text-sm">
                                {#if overview.kraken2 > 0}
                                    <div class="rounded-full bg-primary-600 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                {#if overview.bracken > 0}
                                    <div class="rounded-full bg-primary-600 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                {#if overview.metabuli > 0}
                                    <div class="rounded-full bg-secondary-400 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                {#if overview.ganon2 > 0}
                                    <div class="rounded-full bg-secondary-600 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                {#if overview.kmcp > 0}
                                    <div class="rounded-full bg-secondary-800 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                {#if overview.sylph > 0}
                                    <div class="rounded-full bg-tertiary-200 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                {#if overview.vircov > 0}
                                    <div class="rounded-full bg-tertiary-400 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                <!-- {#if overview.assembly}
                                    <div class="rounded-full bg-tertiary-500 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if} -->
                            </div>
                        </div>
                    </div>
                    <!-- {#if taxonEvidence === overview.taxid}
                        <TaxonEvidenceOverview 
                            taxid={overview.taxid} 
                            taxonOverview={overview} 
                            serverFilterConfig={getServerConfig(i)} 
                            selectedIdentifiers={selectedIdentifiers}
                            selectedTags={getSelectedTags(overview.names).map(x => x.join("-"))}
                            candidateButton={candidateButton}
                        ></TaxonEvidenceOverview>
                    {/if} -->
                </ListBoxItem>
            {/each}
        </ListBox>
    </div>
    {#if pagination}
        <div class="mt-8">
            <Paginator bind:settings={paginationSettings} showFirstLastButtons={false} showPreviousNextButtons={true} class="mt-2" select="paginator-select select text-xs" regionControl="opacity-30 dark:variant-filled-surface dark:opacity-60 text-xs"/>
        </div>
    {/if}
</div>