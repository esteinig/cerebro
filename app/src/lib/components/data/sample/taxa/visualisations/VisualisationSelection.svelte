<script lang="ts">
    export let selectedVisualisation: string;
    export let numberModels: number = 0;
</script>


<div class="flex gap-x-4">
    {#each ['Table', 'History', 'Heatmap'] as viz}
        <button
            class="chip {selectedVisualisation === viz ? 'variant-filled' : 'variant-soft'}"
            on:click={() => { selectedVisualisation = viz }}
            on:keypress
            disabled={numberModels <= 1 && viz === "Heatmap"}
        >
            <span>{viz}</span>
        </button>
    {/each}
</div>