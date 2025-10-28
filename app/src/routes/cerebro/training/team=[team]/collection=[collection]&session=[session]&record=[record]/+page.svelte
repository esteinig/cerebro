<script lang="ts">
    import { onMount } from "svelte";
    import { getModalStore, getToastStore, ProgressRadial, type ToastSettings, type ModalSettings } from "@skeletonlabs/skeleton";
    import type { PatchTrainingRecord, TestResult, TrainingPrefetchData, TrainingResult } from "$lib/utils/types";
    import { page } from "$app/stores";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import PrefetchSlide from "$lib/components/data/sample/taxa/PrefetchSlide.svelte";
	import { navigationLoading } from "$lib/stores/stores.js";
	import { goto, invalidate } from "$app/navigation";
	import { formatAsPercentage } from "$lib/utils/helpers.js";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";

    export let data;

    const publicApi = new CerebroApi();

    const toastStore = getToastStore();
    const modalStore = getModalStore();

    async function updateTrainingRecord(sample_index: number, record_id: string, test_result: TestResult, candidate: string | null) {

        $navigationLoading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.training.updateRecord}?team=${$page.params.team}`,
            {
                method: "PATCH",
                mode: "cors",
                credentials: "include",
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify({
                    "session_id": $page.params.session,
                    "record_id": record_id,
                    "result": test_result,
                    "candidates": candidate ? [candidate] : null
                } as PatchTrainingRecord) 
            } as RequestInit,
            $page.data.refreshToken,
            toastStore,
            null
        );
        
        $navigationLoading = false;

    }


    async function completeTrainingSession() {

        $navigationLoading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.training.completeSession}/${$page.params.session}?team=${$page.params.team}`,
            {
                method: "PATCH",
                mode: "cors",
                credentials: "include"
            } as RequestInit,
            $page.data.refreshToken,
            toastStore,
            null
        );

        $navigationLoading = false;

        if (response.ok) {
            let trainingResult: TrainingResult | null = response.json.data;
            openCompletion(trainingResult)
        }
    }

    interface ClassifyDetail {
        index: number
        item: TrainingPrefetchData
        label: TestResult
        candidate: string | null
    }

    function onNext(e: CustomEvent<ClassifyDetail>)  {
        let { index, item, label, candidate } = e.detail;
        updateTrainingRecord(index, item.id, label, candidate)
        goto(`/cerebro/training/team=${$page.params.team}/collection=${$page.params.collection}&session=${$page.params.session}&record=${index+1}`)
    }

    function onPrevious(e: CustomEvent<ClassifyDetail>)  {
        let { index, item, label, candidate } = e.detail;
        updateTrainingRecord(index, item.id, label, candidate)
        goto(`/cerebro/training/team=${$page.params.team}/collection=${$page.params.collection}&session=${$page.params.session}&record=${index-1}`)
    }

    function getRandomOneToThree(): number {
        return Math.floor(Math.random() * 3) + 1;
    }


    function getRandomeOneToTwo(): number {
        return Math.floor(Math.random() * 2) + 1;
    }

    async function openInfo() {

        let body = `
        <div class="text-large">
        <p class="opacity-60 mt-2">Your task is to select a (single) pathogen candidate species from the provided taxonomic profiling data.</p>
        <p class="opacity-60 mt-2">Click on the information button on the left hand side of the header to see more information about the data</p>
        </div>
        `

        const modal: ModalSettings = {
            type: 'confirm',
            image: `/bughuntress_${getRandomOneToThree()}.png`,
            title: 'Good luck, bug hunter!',
            body: body,
            response: async(confirmed: boolean) => {
                if (confirmed) {
                    // Confirmation response if necessary
                }
            }
        };
        modalStore.trigger(modal);
    }

    async function openCompletion(trainingResult: TrainingResult | null) {

        let results: string = ""; 

        if (trainingResult) {
            results = `
            <p class="mt-2">
                <span class="opacity-80 mr-4">
                    Sensitivity: ${formatAsPercentage(trainingResult.sensitivity)}
                </span> 
                <span class="opacity-80 mr-4">
                    Specificity: ${formatAsPercentage(trainingResult.specificity)}
                </span> 
                <span class="opacity-80 mr-4">
                    Positives: ${trainingResult.true_positive} / ${trainingResult.true_positive+trainingResult.false_negative} 
                </span>  
                <span class="opacity-80 mr-4">
                    Negatives: ${trainingResult.true_negative} / ${trainingResult.true_negative+trainingResult.false_positive} 
                </span> 
            </p>   
            `
        }

        let body = `
        <div class="text-large">  
        ${results}        
        <p class="opacity-60 mt-4">You can download a certificate of completion for your last training session in the overview table.</p> 
        <p class="opacity-60 mt-4"></p>    
        </div>
        `

        const modal: ModalSettings = {
            type: 'confirm',
            image: `/bug_huntress_slam_${getRandomeOneToTwo()}.png`,
            title: 'Well done, bug slayer!',
            body: body,
            response: async(confirmed: boolean) => {
                goto(`/cerebro/training/team=${$page.params.team}`, { invalidateAll: true })
            }
        };
        modalStore.trigger(modal);
    }

    function getSampleNumber() {
        return data.trainingOverview.find(o => o.collection === $page.params.collection)?.samples ?? 0;
    }

    onMount(() => {
        openInfo();
    });

</script>

{#if data.trainingData && data.trainingRecord}
    <PrefetchSlide trainingData={data.trainingData} trainingRecord={data.trainingRecord} trainingSlides={getSampleNumber()} on:next={onNext} on:previous={onPrevious} on:complete={completeTrainingSession}/>

{:else}
    <div class="mt-[10%] flex justify-center items-center">
        <div class="text-center space-y-4">
            <div class="flex justify-center">
                <ErrorAnimation />
            </div>
            <h3 class="h3">
                Failed to fetch training data from server
            </h3>
        </div>
    </div>
{/if}