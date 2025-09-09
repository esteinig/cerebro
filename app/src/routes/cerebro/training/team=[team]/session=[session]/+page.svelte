<script lang="ts">
    import { onMount } from "svelte";
    import { getModalStore, getToastStore, ProgressRadial, type ModalSettings } from "@skeletonlabs/skeleton";
    import type { TrainingPrefetchData } from "$lib/utils/types";
    import { page } from "$app/stores";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import PrefetchSlides from "$lib/components/data/sample/taxa/PrefetchSlides.svelte";

    export let data;
    export let collection: string = "Test";

    const modalStore = getModalStore();

    let prefetch: TrainingPrefetchData[] = [];

    let publicApi = new CerebroApi();
    let toastStore = getToastStore();


    let loading = false;

    async function getPrefetchDataCollection(collection: string) {
        loading = true;
        let response: ApiResponse = await publicApi.fetchWithRefresh(
        `${publicApi.routes.training.getCollection}?team=${data.selectedTeam.name}&collection=${collection}`,
        {
            method: "GET",
            mode: "cors",
            credentials: "include",
        } as RequestInit,
        $page.data.refreshToken,
        toastStore,
        "Training data loaded"
        );

        loading = false;

        if (response.ok) {
            prefetch = response.json.data as TrainingPrefetchData[];
        }
    }


    async function startTrainingSession(collection: string) {
        loading = true;
        let response: ApiResponse = await publicApi.fetchWithRefresh(
        `${publicApi.routes.training.startSession}?team=${data.selectedTeam.name}&collection=${collection}`,
        {
            method: "POST",
            mode: "cors",
            credentials: "include",
        } as RequestInit,
        $page.data.refreshToken,
        toastStore,
        "Training data loaded"
        );

        loading = false;

        if (response.ok) {
            prefetch = response.json.data as TrainingPrefetchData[];
        }
    }

    function onClassify(e: any) {
        const { index, item, label } = e.detail;
        // Persist label or move item, etc.

        console.log(index, item, label)
    }

    function getRandomOneToThree(): number {
        return Math.floor(Math.random() * 3) + 1;
    }

    async function openInfo() {

        let body = `

        <p class="opacity-60 mt-2">Your task is to select a (single) pathogen candidate species from the provided taxonomic profiling data.</p>

        <p class="opacity-60 mt-2"> Click on the threshold categories to see a filtered selection of organisms in this training dataset. Click on any organism to select a positive identification;
        if you believe the sample is negative, simply click 'Next'. Once you move on to the next slide, you are unable to change your decisions.</p> 

        <p class="opacity-60 mt-4">Good luck, bug hunter! 𖢥 </p>
        `

        const modal: ModalSettings = {
            type: 'confirm',
            image: `/bughuntress_${getRandomOneToThree()}.png`,
            title: 'Welcome',
            body: body,
            response: async(confirmed: boolean) => {
                if (confirmed) {
                    
                }
            }
        };
        modalStore.trigger(modal);
    }

    onMount(() => {
        getPrefetchDataCollection(collection);
        openInfo();
        
    });

</script>

{#if loading}
<div class="flex justify-center pt-36">
    <ProgressRadial width="sm:w-12 md:w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
</div>
{:else}
    <PrefetchSlides items={prefetch}
    on:classify={onClassify}
    />
{/if}