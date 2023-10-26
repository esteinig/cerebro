<script lang="ts">
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
    import { page } from "$app/stores";
	import type { CandidateTaxonFormData, CerebroFilterConfig, PriorityTaxonSchema, Taxon, TaxonOverview } from "$lib/utils/types";
	import { ProgressRadial, getToastStore } from "@skeletonlabs/skeleton";
	import { onMount } from "svelte";
    import type { ModalComponent, ModalSettings, ToastSettings } from "@skeletonlabs/skeleton";
	import AlignmentEvidenceOverview from "./AlignmentEvidenceOverview.svelte";
	import KmerEvidenceOverview from "./KmerEvidenceOverview.svelte";
	import AssemblyEvidenceOverview from "./AssemblyEvidenceOverview.svelte";
    import { getModalStore } from '@skeletonlabs/skeleton';
	import CandidateTaxon from "../candidates/CandidateTaxonModal.svelte";

    export let taxid: string;
    export let taxonOverview: TaxonOverview;
    export let serverFilterConfig: CerebroFilterConfig;
    export let selectedIdentifiers: string[];
    export let selectedTags: string[];
    export let candidateButton: boolean = true;

    const publicApi: CerebroApi = new CerebroApi();
    const toastStore = getToastStore();
    const modalStore = getModalStore();

    let loading: boolean = false;

    let taxonData: Taxon;

    let modalComponent: ModalComponent;
    let modal: ModalSettings;

    const getTaxonEvidence = async() => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.taxa}?db=${$page.params.db}&project=${$page.params.project}&id=${selectedIdentifiers.join(",")}&overview=false&taxid=${taxid}`,
            { 
                method: 'POST',  
                mode: 'cors',
                credentials: 'include', 
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify(serverFilterConfig) 
            } as RequestInit,
            $page.data.refreshToken, toastStore, "Taxon evidence loaded"
        )

        loading = false;

        if (response.ok){
            // Should only have one response entry 
            // from taxonomic identifier query
            if (!response.json.data || !response.json.data.taxa.length){
                toastStore.trigger({
                    message: `Error: no taxa found for taxonomic identifier: ${taxid}`,
                    background: "variant-filled-tertiary"
                } satisfies ToastSettings)
            } else if (response.json.data.taxa.length > 1) {
                toastStore.trigger({
                    message: `Error: multiple taxa found for taxonomic identifier ${taxid}`,
                    background: "variant-filled-tertiary"
                } satisfies ToastSettings)
            } else {
                taxonData = response.json.data.taxa[0];
            }
        }

        modalComponent = {
            // Pass a reference to your custom component
            ref: CandidateTaxon,
            // Add the component properties as key/value pairs
            props: { 
                taxonData: taxonData,
                taxonOverview: taxonOverview,
                selectedTags: selectedTags,
                selectedIdentifiers: selectedIdentifiers,
                serverFilterConfig: serverFilterConfig
            },
            // Provide a template literal for the default component slot
            slot: '<p>Select this taxon as a candidate for reporting</p>'
        } satisfies ModalComponent;

        modal = {
            type: 'component',
            // Pass the component directly:
            component: modalComponent,
        } satisfies ModalSettings;
    }
    
    onMount(async() => {
        await getTaxonEvidence();
    })

</script>

<div class="p-8 mt-4 border border-primary-500 border-opacity-30 rounded-lg">
    {#if loading}
        <div class="flex justify-center">
            <ProgressRadial width="sm:w-12 md:w-12" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
        </div>
    {:else}
        
        {#if taxonData}
            <div class="flex align-center">
                <p class="text-2xl opacity-60 italic pb-5">{taxonData.name}<span class="text-sm ml-5 not-italic"> {taxonData.rank} <a href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={taxonData.taxid}" target="_blank" rel="noopener noreferrer">{taxonData.taxid}</a></span></p>
                
                <div class="ml-auto">
                    <!-- <button class="btn variant-outline-primary mr-2" type="button">
                        <div class="h-5 w-5 mr-2">
                            <svg xmlns="http://www.w3.org/2000/svg" height="1em" class="mask-icon mt-0.5" viewBox="0 0 576 512">
                                <path d="M288 64C64 64 0 160 0 272S80 448 176 448h8.4c24.2 0 46.4-13.7 57.2-35.4l23.2-46.3c4.4-8.8 13.3-14.3 23.2-14.3s18.8 5.5 23.2 14.3l23.2 46.3c10.8 21.7 33 35.4 57.2 35.4H400c96 0 176-64 176-176s-64-208-288-208zM96 256a64 64 0 1 1 128 0A64 64 0 1 1 96 256zm320-64a64 64 0 1 1 0 128 64 64 0 1 1 0-128z"/>
                            </svg>
                        </div>
                        Mask
                    </button> -->
                    {#if candidateButton}
                        <button class="btn variant-outline-primary mr-2" type="button" on:click={() => modalStore.trigger(modal)}>
                            <div class="h-5 w-5 mr-2">
                                <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                    <path d="M7.864 4.243A7.5 7.5 0 0119.5 10.5c0 2.92-.556 5.709-1.568 8.268M5.742 6.364A7.465 7.465 0 004.5 10.5a7.464 7.464 0 01-1.15 3.993m1.989 3.559A11.209 11.209 0 008.25 10.5a3.75 3.75 0 117.5 0c0 .527-.021 1.049-.064 1.565M12 10.5a14.94 14.94 0 01-3.6 9.75m6.633-4.596a18.666 18.666 0 01-2.485 5.33" stroke-linecap="round" stroke-linejoin="round"></path>
                                </svg>
                            </div>
                            Candidate
                        </button>
                    {/if}
                    <a class="btn variant-outline-primary" href="{$page.url}/{taxonData.taxid}">
                        <div class="h-5 w-5 mr-2">
                            <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                <path d="M13.5 4.5L21 12m0 0l-7.5 7.5M21 12H3" stroke-linecap="round" stroke-linejoin="round"></path>
                            </svg>
                        </div>
                        Details
                    </a>
                </div>
            </div>
            
            {#if taxonData.evidence.alignment.length} 
                <div class="mb-8">
                    <p class="text-primary-500 opacity-90 mb-2">Alignment</p>
                    <AlignmentEvidenceOverview evidence={taxonData.evidence.alignment}></AlignmentEvidenceOverview>
                </div>
            {/if}
            {#if taxonData.evidence.kmer.length} 
                <div class="mb-8">
                    <p class="text-secondary-500 opacity-90 mb-2">K-mer</p>
                    <KmerEvidenceOverview evidence={taxonData.evidence.kmer}></KmerEvidenceOverview>
                </div>
            {/if}
            {#if taxonData.evidence.assembly.length} 
            <div class="mb-8">
                    <p class="text-tertiary-500 opacity-90 mb-2">Assembly</p>
                    <AssemblyEvidenceOverview evidence={taxonData.evidence.assembly}></AssemblyEvidenceOverview>
                </div>
            {/if}
        {/if}
    {/if}
</div>




<style lang="postcss">
    .mask-icon {
        -webkit-filter: invert(100%);
        filter: invert(100%);
    }
</style>