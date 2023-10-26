<script lang="ts">
	import { formatAsThousands } from "$lib/utils/helpers";
	import type { VircovScanRemapRecord } from "$lib/utils/types";

    export let evidence: VircovScanRemapRecord[];

    evidence.sort((a, b) => a.scan_reads < b.scan_reads ? 1 : -1);

</script>

<div>
    <div class="table-container text-xs">
        {#if evidence.length} 
            <table class="table table-hover table-compact">
                <thead>
                    <tr>
                        <th class="text-center">Taxid</th>
                        <th>Reference</th>
                        <th class="text-left">Length</th>
                        <th class="text-right">Reads</th>
                        <th class="text-right">Regions</th>
                        <th class="text-right">Coverage</th>
                        <th class="text-right opacity-30">Reads</th>
                        <th class="text-right opacity-30">Regions</th>
                        <th class="text-right opacity-30">Coverage</th>
                        <th class="text-right opacity-30">Depth</th>
                    </tr>
                </thead>
                <tbody>
                    {#each evidence as record}

                        <tr class="text-xs">
                            <td class="text-center opacity-60">
                                <a href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={record.taxid}" target="_blank" rel="noopener noreferrer">{record.taxid}</a>           
                            </td>
                            <td class="truncate">
                                <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&QUERY={record.reference}" target="_blank" rel="noopener noreferrer" class="text-xs">{record.reference}</a>
                            </td>
                            <td class="text-left">
                                <span class="text-xs">{formatAsThousands(record.reference_length)}</span>
                                
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{formatAsThousands(record.scan_reads)}</span>
                                
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{record.scan_regions}</span>
                                
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{record.scan_coverage.toFixed(1)}</span>
                                
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{formatAsThousands(record.remap_reads) ?? "-"}</span>
                                
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{record.remap_regions ?? "-"}</span>
                                
                            </td>
                            <td class="text-right">
                                <span class="text-xs">{record.remap_coverage?.toFixed(1) ?? "-"}</span>
                                
                            </td>

                            <td class="text-right">
                                <span class="text-xs">{record.remap_mean_depth?.toFixed(1) ?? "-"}</span>
                                
                            </td>
                        </tr>
                    {/each}
                </tbody>
            </table>
            {:else}
                <p class="flex ml-2 text-sm opacity-30">
                    No alignment data available
                </p>
            {/if}
    </div>
</div>