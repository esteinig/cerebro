<script lang="ts">
	import { AbundanceMode } from "$lib/utils/types";

    export let selectedTableMode: AbundanceMode;
</script>


<div class="flex gap-x-4">
    {#each [ AbundanceMode.Mixed, AbundanceMode.Sequence, AbundanceMode.Profile, AbundanceMode.Bases] as viz}
        <button
            class="chip {selectedTableMode === viz ? 'variant-filled' : 'variant-soft'}"
            on:click={() => {selectedTableMode = viz }}
            on:keypress
        >
            <span>{viz.toLowerCase()}</span>
        </button>
    {/each}
</div>