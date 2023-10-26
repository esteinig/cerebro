<script lang="ts">
	import type { Cerebro } from "$lib/utils/types";
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    export let selectedIdentifiers: string[];
    export let models: Cerebro[] = [];
    export let variant: string = "sample";
    export let variantColor: string = "secondary";

</script>

<div>
    <ListBox multiple>
        {#each models as model}
            <ListBoxItem bind:group={selectedIdentifiers} name="sample" value={model.id} active='variant-ghost' hover='hover:variant-soft' rounded='rounded-token'>
                <div class="ml-5">
                    <div>
                        {#if variant === "sample"}
                            {model.sample.id} <span class="text-xs opacity-60 ml-1">{model.run.id}</span>
                        {:else if variant == "control"}
                            {model.sample.id} <span class="text-xs opacity-60 ml-1">{model.run.id}</span>
                        {:else}
                            {model.name}
                        {/if}
                </div>
                    <div> 
                        <span class="code bg-{variantColor}-500/30 text-{variantColor}-700 dark:bg-{variantColor}-500/20 dark:text-{variantColor}-400">{model.sample.tags?.join("-")}</span>
                        <span class="code bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90">Specimen</span> 
                        <span class="code bg-{variantColor}-500/30 text-{variantColor}-700 dark:bg-{variantColor}-500/20 dark:text-{variantColor}-400">{model.sample.sample_type}</span>
                        <span class="code bg-gray-600 text-gray-300 dark:bg-gray-500/50 dark:text-gray-300/90">Group</span> 
                        <span class="code bg-{variantColor}-500/30 text-{variantColor}-700 dark:bg-{variantColor}-500/20 dark:text-{variantColor}-400">{model.sample.sample_group}</span>
                    </div>
            </ListBoxItem>
        {/each}
    </ListBox>
</div>