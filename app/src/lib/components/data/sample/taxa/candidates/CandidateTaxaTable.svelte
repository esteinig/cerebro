<script lang="ts">
	import CandidateIcon from "$lib/general/icons/CandidateIcon.svelte";
	import { getDateTimeStringUtc } from "$lib/utils/helpers";
    import { DecisionType, type PriorityTaxon, type PriorityTaxonDecision } from "$lib/utils/types";

    // Unique priority taxa for the selected models
    export let priorityTaxa: PriorityTaxon[] = [];
    export let selectedCandidate: PriorityTaxon | null = null;

    const getAcceptCount = (decisions: PriorityTaxonDecision[]): number => {
        return decisions.filter(decision => decision.decision === DecisionType.Accept).length
    }
    const getRejectCount = (decisions: PriorityTaxonDecision[]): number => {
        return decisions.filter(decision => decision.decision === DecisionType.Reject).length
    }

</script>


<div>
    <div class="table-container text-xs">
        {#if priorityTaxa.length} 
            <table class="table table-hover table-compact">
                <thead>
                    <tr>
                        <th class="text-center">Category</th>
                        <th class="text-left">Taxid</th>
                        <th class="text-left">Domain</th>
                        <th class="text-left">Scientific Name</th>
                        <th class="text-left">User</th>
                        <th class="text-left">Date</th>
                        <th class="text-center">Accept<span class="mx-2">:</span>Reject</th>
                    </tr>
                </thead>
                <tbody>
                    {#each priorityTaxa as priorityTaxon}

                        <tr class="text-xs hover:cursor-pointer" on:click={() => { selectedCandidate = priorityTaxon; }}>
                            <td class="text-center opacity-60">
                                <div class="flex justify-center">
                                    <CandidateIcon priorityTaxonType={priorityTaxon.taxon_type}></CandidateIcon>
                                </div>
                            </td>
                            <td class="text-left">
                                <span class="ml-1">{priorityTaxon.taxon_overview.taxid}</span>
                            </td>
                            <td class="text-left">
                                <span class="ml-1">{priorityTaxon.taxon_overview.domain}</span>
                            </td>
                            <td class="text-left">
                                <span class="ml-1 italic">{priorityTaxon.taxon_overview.name}</span>
                            </td>
                            <td class="text-left">
                                <span class="ml-1">{priorityTaxon.user_name}</span>
                            </td>
                            <td class="text-left">
                                <span class="ml-1">{getDateTimeStringUtc(priorityTaxon.date)}</span>
                            </td>
                            <td class="text-center">
                                <div class="flex justify-center align-center items-center">
                                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-3 h-3 mt-0.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                                    </svg>
                                    <span class="ml-1">{getAcceptCount(priorityTaxon.decisions)}</span>
                                    <span class="mx-2"> : </span>
                                    <span class="mr-1">{getRejectCount(priorityTaxon.decisions)}</span>
                                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5"  class="w-3 h-3 mt-0.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M6 18L18 6M6 6l12 12" stroke-linecap="round" stroke-linejoin="round"></path>
                                    </svg>
                                </div>
                            </td>
                        </tr>
                    {/each}
                </tbody>
            </table>
            {:else}
                <p class="flex ml-2 text-sm opacity-30">
                    No candidate taxa available for this sample
                </p>
            {/if}
    </div>
</div>