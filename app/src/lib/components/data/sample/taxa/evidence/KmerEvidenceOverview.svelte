<script lang="ts">
	import { formatAsThousands } from "$lib/utils/helpers";
	import type { Kraken2UniqRecord } from "$lib/utils/types";

    export let evidence: Kraken2UniqRecord[];
    
    evidence.sort((a, b) => a.reads < b.reads ? 1 : -1);

</script>

<div>
    <div class="table-container text-xs">
        {#if evidence.length} 
            <table class="table table-hover table-compact">
                <thead>
                    <tr>
                        <th class="text-center">Taxid</th>
                        <th class="text-center">Rank</th>
                        <th class="text-right">Reads</th>
                        <th class="text-right">Kmers</th>
                        <th class="text-right">Unique</th>
                        <th class="text-center">Tool</th>
                        <th class="text-center">Database</th>
                    </tr>
                </thead>
                <tbody>
                    {#each evidence as record}

                        <tr>
                            <td class="opacity-60 text-center">
                                <a href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={record.taxid}" target="_blank" rel="noopener noreferrer" class="text-xs">{record.taxid}</a>           
                            </td>
                            <td class="text-center">
                                <span class="text-xs">{record.tax_level}</span>
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{formatAsThousands(record.reads)}</span>
                                
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{formatAsThousands(record.kmers)}</span>
                                
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{formatAsThousands(record.kmers_unique)}</span>
                                
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
                    No k-mer data available
                </p>
            {/if}
    </div>
</div>