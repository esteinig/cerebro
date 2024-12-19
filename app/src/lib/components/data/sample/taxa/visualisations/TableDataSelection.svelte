<script lang="ts">
	import { DisplayData } from "$lib/utils/types";

    export let selectedTableData: DisplayData;
</script>


<div class="flex gap-x-4">
    {#each [DisplayData.Rpm, DisplayData.Reads, DisplayData.Abundance] as viz}
        <button
            class="chip {selectedTableData === viz ? 'variant-filled' : 'variant-soft'}"
            on:click={() => { selectedTableData = viz }}
            on:keypress
        >
            <span>{viz}</span>
        </button>
    {/each}
</div>