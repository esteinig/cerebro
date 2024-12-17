<script lang="ts">
    export let selectedVisualisation: string;
    
</script>


<div class="flex gap-x-4">
    {#each ['Pathogens', 'Heatmap'] as viz}
        <button
            class="chip {selectedVisualisation === viz ? 'variant-filled' : 'variant-soft'}"
            on:click={() => { selectedVisualisation = viz }}
            on:keypress
        >
            <span>{viz}</span>
        </button>
    {/each}
</div>