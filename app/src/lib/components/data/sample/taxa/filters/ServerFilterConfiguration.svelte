<script lang="ts">
	import type { CerebroFilterConfig } from "$lib/utils/types";

    export let serverFilterConfig: CerebroFilterConfig;

    let numericInputClass: string = "[appearance:textfield] [&::-webkit-outer-spin-button]:appearance-none [&::-webkit-inner-spin-button]:appearance-none";

</script>

<div>
    
    <p class=""><span class="opacity-40">Alignment evidence</span></p>
    <div class="p-4">
        <div class="grid grid-cols-3 sm:grid-cols-2 md:grid-cols-2 gap-y-4 gap-x-4 w-full text-sm">
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Reads</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.alignment_min_reads}/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Regions</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.alignment_min_regions}/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Coverage (%)</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.alignment_min_coverage}/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Aligned Bases</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.alignment_min_bases}/>
                </label>
            </div>

            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Reference Length</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.alignment_min_ref_length}/>
                </label>
            </div>
        </div>
    </div>

    <p class=""><span class="opacity-40">K-mer evidence</span></p>
    <div class="p-4">
        <div class="grid grid-cols-3 sm:grid-cols-2 md:grid-cols-2 gap-y-4 gap-x-4 w-full text-sm">
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Reads</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.kmer_min_reads}/>
                </label>
            </div>
        </div>
    </div>

    <p class=""><span class="opacity-40">Assembly evidence</span></p>
    <div class="p-4">
        <div class="grid grid-cols-3 sm:grid-cols-2 md:grid-cols-2 gap-y-4 gap-x-4 w-full text-sm">
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Contig Length</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.assembly_min_contig_length}/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Contig Coverage (%)</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.assembly_min_contig_coverage}/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Contig Identity (%)</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.assembly_min_contig_identity}/>
                </label>
            </div>
        </div>
    </div>
</div>