<script lang="ts">
    import type { TaxonFilterConfig } from "$lib/utils/types";

    export let serverFilterConfig: TaxonFilterConfig;

    let numericInputClass: string = "[appearance:textfield] [&::-webkit-outer-spin-button]:appearance-none [&::-webkit-inner-spin-button]:appearance-none";
</script>

<div>
    <p class=""><span class="opacity-40">Taxonomic Filters</span></p>
    <div class="p-4">
        <div class="grid grid-cols-3 sm:grid-cols-2 md:grid-cols-2 gap-y-4 gap-x-4 w-full text-sm">
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Rank</span>
                    <select class="input text-xs" bind:value={serverFilterConfig.rank}>
                        <option value={null}>Any</option>
                        <option value="Superkingdom">Superkingdom</option>
                        <option value="Phylum">Phylum</option>
                        <option value="Class">Class</option>
                        <option value="Order">Order</option>
                        <option value="Family">Family</option>
                        <option value="Genus">Genus</option>
                        <option value="Species">Species</option>
                        <option value="Strain">Strain</option>
                        <option value="NoRank">No Rank</option>
                        <option value="Other">Other</option>
                    </select>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Domains</span>
                    <input type="text" class="input text-xs" bind:value={serverFilterConfig.domains} placeholder="Comma-separated domains"/>
                </label>
            </div>
        </div>
    </div>

    <p class=""><span class="opacity-40">Detection Tools and Modes</span></p>
    <div class="p-4">
        <div class="grid grid-cols-3 sm:grid-cols-2 md:grid-cols-2 gap-y-4 gap-x-4 w-full text-sm">
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Tools</span>
                    <input type="text" class="input text-xs" bind:value={serverFilterConfig.tools} placeholder="Comma-separated tools"/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Modes</span>
                    <input type="text" class="input text-xs" bind:value={serverFilterConfig.modes} placeholder="Comma-separated modes"/>
                </label>
            </div>
        </div>
    </div>

    <p class=""><span class="opacity-40">Thresholds</span></p>
    <div class="p-4">
        <div class="grid grid-cols-3 sm:grid-cols-2 md:grid-cols-2 gap-y-4 gap-x-4 w-full text-sm">
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Reads</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.min_reads}/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum RPM</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.min_rpm} step="0.1"/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60">Minimum Abundance (%)</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.min_abundance} step="0.01"/>
                </label>
            </div>

            <div>
                <label class="label">
                    <span class="text-xs opacity-60">NTC Ratio (if selected)</span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.ntc_ratio} step="0.1"/>
                </label>
            </div>
        </div>
    </div>
</div>
