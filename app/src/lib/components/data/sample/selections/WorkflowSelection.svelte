<script lang="ts">
	import { getDateTimeString } from "$lib/utils/helpers";
	import type { WorkflowConfig } from "$lib/utils/types";
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    export let selectedWorkflowIdentifier: string;
    export let workflows: WorkflowConfig[] = [];
    

</script>

<ListBox>
    {#each workflows as workflow}
	    <ListBoxItem bind:group={selectedWorkflowIdentifier} name="workflow" value={workflow.id} active='variant-ghost' hover='hover:variant-soft'>
            <div class="ml-5">
                <p class="text-sm mb-1">{getDateTimeString(workflow.completed)} - {workflow.name}</p>
                <p class="text-xs opacity-60">
                    {#if workflow.description}
                        {workflow.description}
                    {:else}
                        No description available for this workflow
                    {/if}
                </p>   
            </div>
        </ListBoxItem>
    {/each}
</ListBox>