<script lang="ts">
	import type { Cerebro } from "$lib/utils/types";
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";
	import FileTagChip from "$lib/general/icons/FileTagChip.svelte";
    
    export let models: Cerebro[];
    export let identifiers: string[];
    export let variant: string = "sample";
    export let variantColor: string = "secondary";
</script>

<div>
    {#if models.length}
        <ListBox multiple>
            {#each models as model}
                <ListBoxItem bind:group={identifiers} name="sample" value={model.id} active='variant-ghost' hover='hover:variant-soft' rounded='rounded-token'>
                    <div class="">
                        <div>
                            {#if variant === "sample"}
                                {model.sample.id} <span class="text-xs opacity-60 ml-1">{model.run.id}</span>
                            {:else if variant == "control"}
                                {model.sample.id} <span class="text-xs opacity-60 ml-1">{model.run.id}</span>
                            {:else}
                                {model.name}
                            {/if}
                        </div>
                    </div>
                    <div class=""> 
                        <FileTagChip tags={model.sample.tags} join={true} separator={"-"}/>
                        <span class="code bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90">Specimen</span> 
                        <span class="code bg-{variantColor}-500/30 text-{variantColor}-700 dark:bg-{variantColor}-500/20 dark:text-{variantColor}-400">{model.sample.sample_type}</span>
                        <span class="code bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90">Group</span> 
                        <span class="code bg-{variantColor}-500/30 text-{variantColor}-700 dark:bg-{variantColor}-500/20 dark:text-{variantColor}-400">{model.sample.sample_group}</span>
                    </div>
                </ListBoxItem>
            {/each}
        </ListBox>
    {:else}
        <div class="pl-5 opacity-60">Data not available</div>
    {/if}
</div>