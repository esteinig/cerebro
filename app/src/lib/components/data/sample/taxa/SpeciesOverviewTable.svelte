<script lang="ts">
	import type { CerebroFilterConfig, ClientFilterConfig, TaxonHighlightConfig, TaxonOverview } from "$lib/utils/types";
	import { ListBox, ListBoxItem, Paginator, type PaginationSettings } from "@skeletonlabs/skeleton";
	import TaxonEvidenceOverview from "./evidence/TaxonEvidenceOverview.svelte";

    export let taxonOverview: TaxonOverview[] = [];
    export let modelNameTags: Map<string, string[]> = new Map();
    export let clientFilterConfig: ClientFilterConfig | null = null;
    export let serverFilterConfig: CerebroFilterConfig | CerebroFilterConfig[];
    export let taxonHighlightConfig: TaxonHighlightConfig;
    export let selectedIdentifiers: string[];
    export let candidateButton: boolean = true;
    export let pagination: boolean = true;

    let filteredData: TaxonOverview[] = [];
    let tableData: TaxonOverview[] = [];

    const isTagData = (item: string[] | undefined): item is string[] => { return !!item }
    
    const getSelectedTags = (names: string[]) => {
        return names.map(name => modelNameTags.get(name)).filter(isTagData);
    }

    const applyClientSideFilters = (clientFilterConfig: ClientFilterConfig | null): TaxonOverview[] => {
        
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
            if (clientFilterConfig?.domains.length) {
                taxonomy = clientFilterConfig?.domains.includes(taxonOverview.domain);
            }
            if (clientFilterConfig?.genera.length) {
                taxonomy = clientFilterConfig?.genera.includes(taxonOverview.genus);
            }
            if (clientFilterConfig?.species.length) {
                taxonomy = clientFilterConfig?.species.includes(taxonOverview.name);
            }

            // Explicit minimum aggregate metric filters
            let metrics = taxonOverview.rpm >= clientFilterConfig?.minimum.rpm && 
                taxonOverview.contigs >= clientFilterConfig?.minimum.contigs &&
                taxonOverview.contigs_bases >= clientFilterConfig?.minimum.bases;
            
            // Explicit module filter - if none set retain this taxon overview
           
            if (clientFilterConfig?.modules.alignment && clientFilterConfig?.modules.kmer && clientFilterConfig?.modules.assembly) {
                modules = taxonOverview.alignment && taxonOverview.kmer && taxonOverview.assembly;
            } else if (clientFilterConfig?.modules.alignment && clientFilterConfig?.modules.kmer) {
                modules = taxonOverview.alignment && taxonOverview.kmer;
            } else if (clientFilterConfig?.modules.alignment && clientFilterConfig?.modules.assembly) {
                modules = taxonOverview.alignment && taxonOverview.assembly;
            } else if (clientFilterConfig?.modules.kmer && clientFilterConfig?.modules.assembly) {
                modules = taxonOverview.kmer && taxonOverview.assembly;
            } else if (clientFilterConfig?.modules.alignment) {
                modules = taxonOverview.alignment;
            } else if (clientFilterConfig?.modules.kmer) {
                modules = taxonOverview.kmer;
            } else if (clientFilterConfig?.modules.assembly) {
                modules = taxonOverview.assembly
            } else {
                modules = true
            }

            if (taxonHighlightConfig.contamination.species.some(species => taxonOverview.name.includes(species))) {
                contam = false
            } else  {
                contam = true
            }

            if (taxonHighlightConfig.syndrome.species.some(species => taxonOverview.name.includes(species))) {
                syndrome = true
            } else {
                syndrome = false
            }
                
            return (taxonomy && metrics && modules && contam) || syndrome
        })
    }

    let paginationSettings: PaginationSettings = {
        page: 0,
        limit: 50,
        size: taxonOverview.length,
        amounts: [5, 10, 20, 50, 100, 500],
    };

    $: {    
        filteredData = applyClientSideFilters(clientFilterConfig);
        if (pagination) {
            paginationSettings.size = filteredData.length;
            tableData = filteredData.slice(
                paginationSettings.page * paginationSettings.limit,
                paginationSettings.page * paginationSettings.limit + paginationSettings.limit
            );
        } else {
            tableData = filteredData;
        }
       
    }


    let selectedTaxid: string;
    let taxonEvidence: string | null = null;

    // Can be multiple e.g. saved from candidate taxa
    const getServerConfig = (i: number): CerebroFilterConfig => {
        return Array.isArray(serverFilterConfig) ? serverFilterConfig[i] : serverFilterConfig
    }

    const getTaxonBackgroundColor = (overview: TaxonOverview): string =>  {
        if (taxonHighlightConfig.contamination.species.some(species => overview.name.includes(species))) {
            return 'variant-soft-secondary rounded-token py-1.5 px-2'
        } else if (taxonHighlightConfig.syndrome.species.some(species => overview.name.includes(species))) {
            return 'variant-soft-tertiary rounded-token py-1.5 px-2'
        } else {
            return 'rounded-token py-1.5 px-2'
        }
    }

    const getTaxonHover = (overview: TaxonOverview): string =>  {
        if (taxonHighlightConfig.syndrome.species.some(species => overview.name.includes(species))){
            return 'hover:variant-soft-secondary'
        } else if (taxonHighlightConfig.syndrome.species.some(species => overview.name.includes(species))) {
            return 'hover:variant-soft-tertiary'
        } else {
            return 'hover:variant-soft'
        }
    }

</script>

<div>
    <div>
        <ListBox>
            <ListBoxItem group="header" name="header" value="qc" active='variant-soft' hover='hover:cursor-default' rounded='rounded-token'>
                <div class="grid grid-cols-11 sm:grid-cols-11 md:grid-cols-11 gap-x-1 gap-y-4 w-full text-sm opacity-60">
                    <div class="col-span-1">Domain</div>
                    <div class="col-span-2">Species</div>
                    <div class="col-span-2">Tags</div>
                    <div class="text-right">RPM</div>
                    <div class="text-right opacity-40">K-mer</div>
                    <div class="text-right opacity-40">Alignment</div>
                    <div class="text-right">Contigs</div>
                    <div class="text-right opacity-40">Bases</div>
                    <div class="text-right">Modules</div>
                </div>
            </ListBoxItem>
            {#each tableData as overview, i}
                <ListBoxItem bind:group={selectedTaxid} name={overview.name} value={overview.taxid} active='' hover={getTaxonHover(overview)} regionDefault={getTaxonBackgroundColor(overview)} rounded='rounded-token' on:click={() => taxonEvidence === null ? taxonEvidence = overview.taxid : taxonEvidence = null}>
                    
                    <div class="grid grid-cols-11 sm:grid-cols-11 md:grid-cols-11 gap-x-1 gap-y-4 w-full text-sm">
                        <div class="col-span-1 opacity-70">{overview.domain}</div>
                        <div class="col-span-2 truncate italic">{overview.name}</div>
                        <div class="col-span-2 truncate">
                            {#each getSelectedTags(overview.names) as tags}
                                <span class="code bg-primary-500/30 text-primary-700 dark:bg-primary-500/20 dark:text-primary-400 mr-1">{tags.join("-")}</span>
                            {/each}
                        </div>
                        <div class="text-right">{overview.rpm.toFixed(1)}</div>
                        <div class="text-right opacity-40">{overview.rpm_kmer.toFixed(1)}</div>
                        <div class="text-right opacity-40">{overview.rpm_alignment.toFixed(1)}</div>
                        <div class="text-right">{overview.contigs.toFixed(0)}</div>
                        <div class="text-right opacity-40">{overview.contigs_bases.toFixed(0)}</div>
                        <div class="flex justify-end gap-x-1 items-center align-center pt-1">

                            <div class="grid grid-cols-3 sm:grid-cols-3 md:grid-cols-3 gap-x-1 w-1/3 text-sm">
                                {#if overview.alignment}
                                    <div class="rounded-full bg-primary-500 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                {#if overview.kmer}
                                    <div class="rounded-full bg-secondary-500 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                                {#if overview.assembly}
                                    <div class="rounded-full bg-tertiary-500 h-2 w-2"></div>
                                {:else}
                                    <div></div>
                                {/if}
                            </div>
                        </div>
                    </div>
                    {#if taxonEvidence === overview.taxid}
                        <TaxonEvidenceOverview 
                            taxid={overview.taxid} 
                            taxonOverview={overview} 
                            serverFilterConfig={getServerConfig(i)} 
                            selectedIdentifiers={selectedIdentifiers}
                            selectedTags={getSelectedTags(overview.names).map(x => x.join("-"))}
                            candidateButton={candidateButton}
                        ></TaxonEvidenceOverview>
                    {/if}
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