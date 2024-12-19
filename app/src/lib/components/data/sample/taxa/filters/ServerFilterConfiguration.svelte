<script lang="ts">
    import CircleIndicator from "$lib/general/icons/CircleIndicator.svelte";
    import { ProfileTool, type PrevalenceContaminationConfig, type TaxonFilterConfig } from "$lib/utils/types";

    export let serverFilterConfig: TaxonFilterConfig;
    export let prevalenceContamConfig: PrevalenceContaminationConfig;

    const kmerTools: ProfileTool[] = [
        ProfileTool.Bracken,
        ProfileTool.Metabuli,
        ProfileTool.Ganon2,
    ];

    const alignmentTools: string[] = [
        "Alignment"
    ]

    const assemblyTools: string[] = [
        "Assembly"
    ];
    
    let numericInputClass: string = "[appearance:textfield] [&::-webkit-outer-spin-button]:appearance-none [&::-webkit-inner-spin-button]:appearance-none";
</script>

<div>
    <p class=""><span class="opacity-40">Evidence thresholds</span></p>
    <div class="p-4">
        <div class="grid grid-cols-3 sm:grid-cols-2 md:grid-cols-2 gap-y-4 gap-x-4 w-full text-sm">
            <div>
                <label class="label">
                    <span class="text-xs opacity-60 flex items-center">Minimum reads<CircleIndicator circleClass="ml-auto" /><CircleIndicator circleClass="ml-1 mr-2" color="bg-secondary-500"/></span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.min_reads}/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60 flex items-center">Minimum reads per million (rpm)<CircleIndicator circleClass="ml-auto" /><CircleIndicator circleClass="ml-1 mr-2" color="bg-secondary-500"/></span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.min_rpm} step="0.1"/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60 flex items-center">Minimum contig length (bp)<CircleIndicator circleClass="ml-auto mr-2" color="bg-tertiary-500"/></span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.min_bases}/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60 flex items-center">Minimum contig length per million (bpm)<CircleIndicator circleClass="ml-auto mr-2" color="bg-tertiary-500"/></span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.min_bpm} step="0.1"/>
                </label>
            </div>
            <div>
                <label class="label">
                    <span class="text-xs opacity-60 flex items-center">Minimum abundance (%)<CircleIndicator circleClass="ml-auto" /><CircleIndicator circleClass="ml-1" color="bg-secondary-500"/><CircleIndicator circleClass="ml-1 mr-2" color="bg-tertiary-500"/></span>
                    <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.min_abundance} step="0.01"/>
                </label>
            </div>
        </div>
    </div>
    <p class="pt-4"><span class="opacity-40">Negative template control</span></p>
    <div class="p-4">
        <div class="grid grid-cols-3 sm:grid-cols-2 md:grid-cols-2 gap-y-4 gap-x-4 w-full text-sm">
            <label class="label">
                <span class="text-xs opacity-60 flex items-center">NTC : Library ratio ({">"})<CircleIndicator circleClass="ml-auto" /><CircleIndicator circleClass="ml-1 mr-2" color="bg-secondary-500"/></span>
                <input type="number" class="input text-xs {numericInputClass}" bind:value={serverFilterConfig.ntc_ratio} step="0.1"/>
            </label>
            <div class="pt-6">
                <span class="text-xs opacity-40">More than {serverFilterConfig.ntc_ratio} x RPM must be present in NTC libraries compared to the sample library to remove tool-specific evidence for a taxon. NTC RPM are summed if multiple are selected. Ratios are computed separately for DNA and RNA.</span>
            </div>
        </div>
    </div>

    <p class="pt-4"><span class="opacity-40">Prevalence contamination</span></p>
    <div class="p-4">
        <div class="grid grid-cols-3 sm:grid-cols-2 md:grid-cols-2 gap-y-4 gap-x-4 w-full text-sm">
            <label class="label">
                <span class="text-xs opacity-60 flex items-center">Prevalence contamination threshold (fraction)<CircleIndicator circleClass="ml-auto" /><CircleIndicator circleClass="ml-1" color="bg-secondary-500"/><CircleIndicator circleClass="ml-1 mr-2" color="bg-tertiary-500"/></span>
                <input type="number" class="input text-xs {numericInputClass}" bind:value={prevalenceContamConfig.threshold} step="0.01" min="0" max="1"/>
            </label>
            <label class="label">
                <span class="text-xs opacity-60 flex items-center">Minimum reads per million (rpm)<CircleIndicator circleClass="ml-auto" /><CircleIndicator circleClass="ml-1" color="bg-secondary-500"/><CircleIndicator circleClass="ml-1 mr-2" color="bg-tertiary-500"/></span>
                <input type="number" class="input text-xs {numericInputClass}" bind:value={prevalenceContamConfig.min_rpm} step="0.01" min="0"/>
            </label>
            <div class="">
                <span class="text-xs opacity-40">Prevalence contamination identifies any taxa present in a percentage of other libraries ({"> "}{prevalenceContamConfig.threshold*100}% of libraries in the current project) with the same tag as the selected library (DNA/RNA)</span>
            </div>
        </div>
    </div>

    <!-- <p class=""><span class="opacity-40">Taxonomy</span></p>
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
        </div>
    </div> -->
    
</div>
