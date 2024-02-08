<script lang="ts">
	import { invalidate } from "$app/navigation";
	import { page } from "$app/stores";
	import CandidateIcon from "$lib/general/icons/CandidateIcon.svelte";
	import { storeSettings } from "$lib/stores/stores";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import { getDateTimeStringUtc } from "$lib/utils/helpers";
    import { DecisionType, type PriorityTaxon, type PriorityTaxonDecision, type PriorityTaxonDecisionSchema } from "$lib/utils/types";
	import { ProgressRadial, getToastStore } from "@skeletonlabs/skeleton";

    export let priorityTaxon: PriorityTaxon;

    const publicApi = new CerebroApi();
    const toastStore = getToastStore();

    const getUserDecision = (priorityTaxon: PriorityTaxon): PriorityTaxonDecision | null => {
        let userDecisions = priorityTaxon.decisions.filter(decision => decision.user_id === $page.data.userData.id);
        if (!userDecisions.length) {
            return null
        }
        return userDecisions[0]
    }

    
    let loading: boolean = false;
    let userDecision: PriorityTaxonDecision | null; 
    let decisionComment: string = "";

    $: {
        userDecision = getUserDecision(priorityTaxon);
    }

    const submitDecision = async(decision: DecisionType) => {
        loading = true;

		let decisionSchema = {
			id: priorityTaxon.id,
            decision: decision,
            decisionComment: decisionComment,
            taxonName: priorityTaxon.taxon_overview.name,
            taxonTaxid: priorityTaxon.taxon_overview.taxid,
            taxonType: priorityTaxon.taxon_type
		} satisfies PriorityTaxonDecisionSchema;


		let response: ApiResponse = await publicApi.fetchWithRefresh(
			`${publicApi.routes.cerebro.updatePriorityTaxaDecision}?db=${$page.params.db}&project=${$page.params.project}&id=${priorityTaxon.cerebro_identifiers.join(",")}`, {
				method: 'PATCH',
				mode: 'cors',
				credentials: 'include',
				body: JSON.stringify(decisionSchema),
				headers:  { 'Content-Type': 'application/json' }
			} satisfies RequestInit,
			$page.data.refreshToken, toastStore, `Decision submitted`
		);

		if (response.ok) {
			// Updating page data with reload
			await invalidate("sample:data")
		}

		loading = false;
    }

</script>

<div class="card max-w-xl2">
    <header class="card-header">
        <div class="p-4 space-y-4">
           
            <div>
                <ol class="breadcrumb justify-start text-lg">
                    <li class="crumb opacity-60">{priorityTaxon.taxon_overview.domain}</li>
                    <li class="crumb-separator" aria-hidden>&rsaquo;</li>
                    <li class="crumb italic">{priorityTaxon.taxon_overview.name}</li>
                </ol>
            </div>

            <div class="flex items-center align-center">

                <div class="flex align-center items-center text-sm gap-x-2 mr-5">
                    {#each priorityTaxon.evidence_tags as tag}
                        <span class="code text-xs mt-1 mr-1">{tag}</span>
                    {/each}
                    <CandidateIcon priorityTaxonType={priorityTaxon.taxon_type} size={4}></CandidateIcon>
                    {priorityTaxon.taxon_type}

                    
                </div>
            </div>

            
        </div>
    </header>
	<section class="p-8 space-y-2">      
        {#if loading}
            <div class="flex justify-center py-24">
                <ProgressRadial width="sm:w-12 md:w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
            </div>
        {:else}
            <p class="opacity-80">Submission</p>
            <div class="p-4">
                <p class="opacity-60 mb-1"><span class="font-semibold">{priorityTaxon.user_name}</span> ({getDateTimeStringUtc(priorityTaxon.date)})</p>
                
                {#if priorityTaxon.comment}
                    <p class="opacity-60 text-xs w-3/4">{priorityTaxon.comment}</p>
                {:else}
                    <p class="opacity-60 text-xs">No submission comment</p>
                {/if}
            </div>

            <p class="opacity-80">Decision</p>
            <div class="p-4">
                {#if userDecision?.decision === DecisionType.Accept}
                    <p class="opacity-60 text-xs">You <span class="font-semibold">accepted</span> this candidate</p>
                {:else if userDecision?.decision === DecisionType.Reject}
                        <p class="opacity-60 text-xs">You <span class='font-semibold'>rejected</span> this candidate</p>
                {:else}
                    <p class="opacity-60 text-xs">You have not made a decision</p>
                {/if}
                <p class="opacity-60 text-xs mt-4">
                    You can submit or update a decision on this candidate selection as <span class="font-semibold">{priorityTaxon.taxon_type}</span>.
                    Your submission may include a comment on the reason for your decision.
                </p>
                

                <div class="pt-5 bg-surface-500/0 {$storeSettings.acceptedCommentRisk ? "blur-none" : "blur-xl"}">
                    <div class="grid-cols-[auto_1fr_auto] rounded-container-token">
                        <textarea
                            bind:value={decisionComment}
                            class="textarea"
                            name="prompt"
                            id="prompt"
                            placeholder="Comment on your decision..."
                            rows="4"
                            disabled={!$storeSettings.acceptedCommentRisk}
                        />
                    </div>
                </div>
                <div class="text-end pt-5"> 

                    <button class="btn variant-outline-tertiary" on:click={() => submitDecision(DecisionType.Reject)}>
                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5"  class="w-5 h-5 mr-1" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M6 18L18 6M6 6l12 12" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg>
                        Reject
                    </button>
                    <button class="btn variant-outline-primary ml-2" on:click={() => submitDecision(DecisionType.Accept)}>
                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-5 h-5 mr-1" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg>
                        Accept
                    </button>
                </div>
            </div>
        {/if}
        
    </section>
	<footer class="card-footer">

    </footer>
</div>