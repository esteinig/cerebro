<script lang="ts">
	import { goto } from "$app/navigation";
	import { page } from "$app/stores";
	import { type TrainingPrefetchOverview } from "$lib/utils/types";
	import { navigationLoading } from "$lib/stores/stores";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import { getToastStore } from "@skeletonlabs/skeleton";

    export let data: TrainingPrefetchOverview[] = [];
    export let userName: string = "CerebroUser";

    const publicApi = new CerebroApi();
    const toastStore = getToastStore();

    async function downloadCertificate(session_id: string | null) {

        if (!session_id) {
            return
        }

        $navigationLoading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.training.downloadCertificate}/${session_id}/certificate?team=${$page.params.team}`,
            {
                method: "GET",
                mode: "cors",
                credentials: "include"
            } as RequestInit,
            $page.data.refreshToken,
            toastStore,
            null
        );

        $navigationLoading = false;

        if (response.ok) {
            let pdfData: number[] = response.json.data;
            downloadPDF(pdfData)
        }
    }

    export function makeCertificateFilename(userName: string, now: Date = new Date()): string {
        const pad = (n: number) => String(n).padStart(2, "0");
        const y = now.getFullYear();
        const m = pad(now.getMonth() + 1);
        const d = pad(now.getDate());
        const clean = userName.replace(/\s+/g, "");
        return `${y}${m}${d}_${clean}_CerebroTraining_Certificate.pdf`;
    }

    const downloadPDF = (pdfData: number[]) => {

        const blob = new Blob([new Uint8Array(pdfData)], { type: 'application/pdf' });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');

        link.href = url;
        link.download = makeCertificateFilename(userName);
        link.click();
        URL.revokeObjectURL(url);
    };

</script>

<div>
    <div class="table-container pt-8 text-xl">
        <table class="table table-hover table-compact">
            <thead>
                <tr>
                    <th>Training Dataset</th>
                    <th>Training Samples</th>
                    <th>Dataset Description</th>
                    <th>Operations</th>
                </tr>
            </thead>
            <tbody>
                {#if data.length === 0}
					<tr>
						<td colspan="4" class="text-center opacity-60">No datasets available</td>
					</tr>
				{:else}
					{#each data as overview}
						<tr>
							<td class="!align-middle !text-lg">{overview.collection}</td>
							<td class="!align-middle !text-lg">{overview.samples}</td>
							<td class="!align-middle !text-lg opacity-80">{overview.description}</td>
							<td class="grid grid-cols-2 gap-x-2 gap-y-2">
                                <button type="button" class="btn btn-md variant-outline-primary" on:click={() => goto(`${$page.url}/collection=${overview.collection}&session=0&record=0`)}>

                                    <div class="w-5 h-5 mr-2">
                                        <svg aria-hidden="true" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg">
                                            <path clip-rule="evenodd" d="M3 10a.75.75 0 01.75-.75h10.638L10.23 5.29a.75.75 0 111.04-1.08l5.5 5.25a.75.75 0 010 1.08l-5.5 5.25a.75.75 0 11-1.04-1.08l4.158-3.96H3.75A.75.75 0 013 10z" fill-rule="evenodd"></path>
                                        </svg>
                                    </div>           
                                    New Session     
                                </button>
                                <button type="button" class="btn btn-md variant-outline-secondary ml-4" on:click={() => downloadCertificate(overview.session_id)} disabled={overview.session_id ? false : true}>

                                    <div class="w-5 h-5 mr-2">
                                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                            <path d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" stroke-linecap="round" stroke-linejoin="round"></path>
                                        </svg>
                                    </div>            
                                    Certificate    
                                </button>
							</td>
						</tr>
					{/each}
				{/if}
            </tbody>
        </table>
    </div>
</div>