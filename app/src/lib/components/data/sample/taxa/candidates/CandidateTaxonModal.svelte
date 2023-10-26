<script lang="ts">
	import { invalidate } from '$app/navigation';
    import { page } from '$app/stores';
	import { CerebroApi, type ApiResponse } from '$lib/utils/api';
	import { PriorityTaxonType, type CandidateTaxonFormData, type Taxon, type PriorityTaxonSchema, type CerebroFilterConfig, type TaxonOverview } from '$lib/utils/types';


	/** Exposes parent props to this component. */
	export let parent: any;

    export let taxonData: Taxon;
    export let taxonOverview: TaxonOverview;
    export let selectedIdentifiers: string[];
    export let selectedTags: string[];
    export let serverFilterConfig: CerebroFilterConfig;


	// Stores
	import { ProgressRadial, getModalStore, getToastStore } from '@skeletonlabs/skeleton';
	const modalStore = getModalStore();
    const toastStore = getToastStore();
    
	// Form Data
	const formData: CandidateTaxonFormData = {
        userIdentifier: $page.data.userData.id,
        userName: $page.data.userData.name,
		comment: "",
        taxonType: PriorityTaxonType.Unknown
	};

    const publicApi = new CerebroApi();
    let loading: boolean = false;

    const submitPriorityTaxon = async() => {
        loading = true;

        let priorityTaxonSchema: PriorityTaxonSchema = {
            user_id: formData.userIdentifier,
            user_name: formData.userName,
            comment: formData.comment,
            taxon_type: formData.taxonType,
            evidence_tags: selectedTags,
            cerebro_identifiers: selectedIdentifiers,
            filter_config: serverFilterConfig,
            taxon_overview: taxonOverview,
            decisions: []
        };

        let _: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.addPriorityTaxon}?db=${$page.params.db}&project=${$page.params.project}&id=${selectedIdentifiers.join(",")}`,
            { 
                method: 'POST',  
                mode: 'cors',
                credentials: 'include', 
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify(priorityTaxonSchema) 
            } as RequestInit,
            $page.data.refreshToken, toastStore, "Candidate taxon submitted"
        )

        // Updating page data with reload
        await invalidate("sample:data")

        loading = false;

        // Return the form data as response and close modal
		if ($modalStore[0].response) $modalStore[0].response(formData);
		modalStore.close();

    } 


</script>

<!-- @component This example creates a simple form modal. -->

{#if $modalStore[0]}
	<div class="modal-example-form card p-4 w-modal shadow-xl space-y-4">
		<header class="text-2xl font-bold">Candidate taxon</header>
		
        <article class="text-sm opacity-60">
            <slot></slot>
        </article>
        
        {#if loading}
            <div class="flex justify-center items-center py-24">
                <ProgressRadial width="sm:w-12 md:w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
            </div>
        {:else}
            <form class="modal-form border border-surface-500 p-4 space-y-4 rounded-container-token">
                <label class="label">
                    <span>Scientific name</span>
                    <input class="input" type="text" disabled value={taxonData.name} />
                </label>
                <label class="label">
                    <span>Taxonomy</span>
                    <input class="input" type="text" disabled value={`${taxonData.rank} (${taxonData.taxid})`}/>
                </label>

                <label class="label">
                    <span>Submitted by</span>
                    <input class="input" type="text" disabled bind:value={formData.userName}/>
                </label>
                <label class="label">
                    <span>Submission date</span>
                    <input class="input" type="text" disabled value={new Date()}/>
                </label>
                <label class="label">
                    <span>Comment</span>
                    <textarea class="input" rows="2" bind:value={formData.comment} required placeholder="Reason for this selection..." />
                </label>

                <label class="label">
                    <span>Taxon category</span>
                    <p class="opacity-60 text-xs">Select a broad classification for this candidate</p>
                    <div class="pt-1">
                    {#each [PriorityTaxonType.Unknown, PriorityTaxonType.Pathogen, PriorityTaxonType.Contaminant] as c}
                        <span
                            class="chip {formData.taxonType === c ? 'variant-filled' : 'variant-soft'} mr-2"
                            on:click={() => { formData.taxonType = c }}
                            on:keypress
                        >
                            {#if formData.taxonType === c}
                                <span>
                                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-3 h-3" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                                    </svg>
                                </span>
                            {/if}
                            <span>{c}</span>
                        </span>
                    {/each}
                    </div>
                </label>
            </form>
        {/if}
		<!-- prettier-ignore -->
		<footer class="modal-footer {parent.regionFooter}">
        <button class="btn {parent.buttonNeutral}" on:click={parent.onClose}>{parent.buttonTextCancel}</button>
        <button class="btn {parent.buttonPositive}" on:click={submitPriorityTaxon}>Submit</button>
    </footer>
	</div>
{/if}