<script lang="ts">
	import { DisplayVisualisation } from "$lib/utils/types";

    export let selectedVisualisation: DisplayVisualisation;
    export let numberModels: number = 0;

</script>


<div class="flex gap-x-4">
    {#each [ DisplayVisualisation.Threshold, DisplayVisualisation.Table, DisplayVisualisation.Run, DisplayVisualisation.History,] as viz}
        <button
            class="chip {selectedVisualisation === viz ? 'variant-filled' : 'variant-soft'}"
            on:click={() => { selectedVisualisation = viz }}
            on:keypress
        >
            <span>{viz}</span>
        </button>
    {/each}
</div>