<script lang="ts">
    import { Role } from "$lib/utils/types";
    import type { Cerebro, PriorityTaxon, WorkflowConfig } from "$lib/utils/types";
	import ReportAuditTable from "./report/ReportAuditTable.svelte";
	import NewReportForm from "./report/NewReportForm.svelte";
	import { page } from "$app/stores";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import CandidateTaxaTable from "./taxa/candidates/CandidateTaxaTable.svelte";
	import { storeSettings } from "$lib/stores/stores";

    export let selectedModels: Cerebro[];
    export let selectedIdentifiers: string[];
    export let selectedWorkflowConfiguration: WorkflowConfig;

    let priorityTaxa: PriorityTaxon[];

    $: {
        priorityTaxa = selectedModels.flatMap(cerebro => cerebro.sample.priority).filter((value, index, self) => 
            self.findIndex(v => v.id === value.id) === index
        )
    }


    let selectedCandidate: PriorityTaxon | null = null;

</script>

<div class="p-4">
    <div>

        <div class="">
            <p class="pb-2 opacity-60 text-2xl">Reports</p>
        </div>
        <ReportAuditTable selectedModels={selectedModels}></ReportAuditTable>
        
        <div class="pt-4">
            <p class="pb-2 opacity-60 text-2xl">Candidates</p>
        </div>
        {#if priorityTaxa}
            <CandidateTaxaTable priorityTaxa={priorityTaxa} bind:selectedCandidate={selectedCandidate}></CandidateTaxaTable>
        {:else}

             <p class="flex justify-center text-sm pb-8">No candidate taxa selected for this sample</p>
            
        {/if}

        {#if $page.data.userData.roles.includes(Role.Report)}
            <div class="pt-8 pb-3">
                <p class="opacity-60 text-2xl">New Report</p>
            </div>
            {#if !$storeSettings.acceptedCommentRisk}
                <div class="flex justify-center p-16">
                    <aside class="alert variant-ghost-warning text-sm w-3/4">
                        <svg class="w-36" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M12 9v3.75m-9.303 3.376c-.866 1.5.217 3.374 1.948 3.374h14.71c1.73 0 2.813-1.874 1.948-3.374L13.949 3.378c-.866-1.5-3.032-1.5-3.898 0L2.697 16.126zM12 15.75h.007v.008H12v-.008z" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg>
                        <div class="alert-message">
                            <h3 class="h3 font-bold">Warning</h3>
                            <p>
                                Submitting reports carries a risk of exposing protected patient information 
                                in the comment fields. Reports are generated on the server. Completed reports are stored in the database.
                            </p>
                            <p> Please consider compliance with local regulatory and ethics frameworks before submitting a report. </p>				
                            <button type="button" class="btn-icon-lg mr-5"  on:click={() => $storeSettings.acceptedCommentRisk = true}>
                                Accept
                            </button>
                            <button type="button" class="btn-icon-lg"  on:click={() => $storeSettings.acceptedCommentRisk = false}>
                                Reject
                            </button>
                        </div>
                    </aside>
                </div>
            {:else}
                <div class="p-4 border border-surface-500 rounded-xl bg-surface-500/5">
                    <NewReportForm selectedIdentifiers={selectedIdentifiers} selectedModels={selectedModels} selectedWorkflowConfiguration={selectedWorkflowConfiguration} selectedPriorityTaxon={selectedCandidate}></NewReportForm>
                </div>
            {/if}
        {:else}
            
            <div class="flex justify-center py-16 ">
                <ErrorAnimation></ErrorAnimation>
            </div>
            <p class="flex justify-center text-sm pb-8">You do not have permissions to generate a report.</p>
            
        {/if}
    </div> 

</div>