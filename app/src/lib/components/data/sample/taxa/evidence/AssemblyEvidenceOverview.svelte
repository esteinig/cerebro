<script lang="ts">
	import { formatAsThousands } from "$lib/utils/helpers";
	import type { BlastLcaRecord } from "$lib/utils/types";
	import { Paginator, type PaginationSettings } from "@skeletonlabs/skeleton";

    export let evidence: BlastLcaRecord[];
    
    evidence.sort((a, b) => a.length < b.length ? 1 : -1);

    let tableData: BlastLcaRecord[] = [];
    
    const formatReference = (ref: string): string | undefined => {
        // NCBI NR annotations
        if (ref.includes("gi|")){
            return ref.split("|").at(-2)
        }
        // Virosaurus annotations
        if (ref.includes(";")){
            return ref.split(";").at(0)
        }
        return ref
    }

    let paginationSettings: PaginationSettings = {
        page: 0,
        limit: 50,
        size: evidence.length,
        amounts: [5, 10, 20, 50, 100, 500],
    };

    $: {    
        paginationSettings.size = evidence.length;
		tableData = evidence.slice(
			paginationSettings.page * paginationSettings.limit,
			paginationSettings.page * paginationSettings.limit + paginationSettings.limit
		);
    }
</script>

<div>
    <div class="table-container text-xs">
        {#if evidence.length} 
            <table class="table table-hover table-compact">
                <thead>
                    <tr>
                        <th class="text-center">Taxid</th>
                        <th class="text-left">Reference</th>
                        <th class="text-right">Length</th>
                        <th class="text-right">Coverage</th>
                        <th class="text-right">Identity</th>
                        <th class="text-center">Tool</th>
                        <th class="text-center">Database</th>
                    </tr>
                </thead>
                <tbody>
                    {#each tableData as record}

                        <tr>
                            <td class="opacity-60 text-center">
                                <a href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={record.taxid}" target="_blank" rel="noopener noreferrer" class="text-xs">{record.taxid}</a>           
                            </td>
                            <td class="text-left truncate">
                                <span class="text-xs"><a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM={record.tool === 'diamond' ? 'blastp' : 'blastn'}&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&QUERY={formatReference(record.reference)}" target="_blank" rel="noopener noreferrer">{formatReference(record.reference)}</a></span>
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{formatAsThousands(record.length)}</span>
                                
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{record.coverage.toFixed(1)}</span>
                                
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{record.identity.toFixed(1)}</span>
                                
                            </td>
                            <td class="text-center">
                                <span class="text-xs">{record.tool}</span>
                            </td>
                            <td class="text-center">
                                <span class="text-xs">{record.db}</span>
                            </td>
                        </tr>
                    {/each}
                </tbody>
            </table>
            {:else}
                <p class="flex ml-2 text-sm opacity-30">
                    No assembly data available
                </p>
            {/if}
    </div>
    
    <div class="mt-2">
        <Paginator bind:settings={paginationSettings} showFirstLastButtons={false} showPreviousNextButtons={true} class="mt-2" select="paginator-select select text-xs" regionControl="opacity-30 dark:variant-filled-surface dark:opacity-60 text-xs"/>
    </div>
</div>