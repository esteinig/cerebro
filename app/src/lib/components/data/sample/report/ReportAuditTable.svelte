<script lang="ts">
	import { page } from "$app/stores";
	import { CerebroApi, type ApiResponse } from "$lib/utils/api";
	import { getDateTimeString, getDateTimeStringUtc } from "$lib/utils/helpers";
	import type { Cerebro, ReportEntry } from "$lib/utils/types";
	import { getToastStore, type ToastSettings } from "@skeletonlabs/skeleton";

    export let selectedModels: Cerebro[];
    let selectedModelReports: ReportEntry[];

    $: selectedModelReports = selectedModels.flatMap(cerebro => cerebro.sample.reports).filter((value, index, self) => {
		return self.findIndex(v => v.id === value.id) === index;
	})

    const publicApi: CerebroApi = new CerebroApi();
    const toastStore = getToastStore();

    let loading: boolean = false;

    const downloadReportCopy = async(id: string) => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.getReport}/${id}?db=${$page.params.db}&project=${$page.params.project}`,
            { 
                method: 'GET',  
                mode: 'cors',
                credentials: 'include',
            } as RequestInit,
            $page.data.refreshToken, toastStore, null
        )

        loading = false;

        if (response.ok) {
            let data: Blob;
            let fileName: string;
            let now = getDateTimeString(new Date().toISOString(), false);
            let report_entry: ReportEntry | undefined = response.json.data.report_entry;

            if (report_entry && report_entry.report_pdf) {
                // PDF base64-encoded string to blob
                const byteArray = Uint8Array.from(
                    atob(response.json.data.report_entry.report_text)
                    .split('')
                    .map(char => char.charCodeAt(0))
                );
                data = new Blob([byteArray], { type: 'application/pdf' });
                fileName = `${now}_${$page.params.sample}_ClinicalReport.pdf`;
            } else if (report_entry && !report_entry.report_pdf) {
                // LaTeX report as string
                data = new Blob([response.json.data.tex], {type: 'text/plain'});
                fileName = `${$page.params.sample}_Report_${now}_ClinicalReport.tex`;
            } else {
                toastStore.trigger({
                    message: response.json.message ?? "Failed to extract report from response",
                    background: "variant-filled-tertiary"
                } satisfies ToastSettings);
                return undefined
            };
            
            toastStore.trigger({
                message: "Retrieved report",
                background: "variant-filled-primary"
            } satisfies ToastSettings);

            let link = document.createElement("a");
            if (link.download !== undefined) { // feature detection
                // Browsers that support HTML5 download attribute
                let url = URL.createObjectURL(data);
                link.setAttribute("href", url);
                link.setAttribute("download", fileName);
                link.style.visibility = 'hidden';
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
            }
        }
    }


</script>

<div class="table-container text- pb-4">
    {#if selectedModelReports.length} 
        <table class="table table-hover table-compact text-sm">
            <thead>
                <tr>
                    <th class="text-center">Date</th>
                    <th class="text-center">Report ID</th>
                    <th class="text-center">Submission</th>
                    <th class="text-center">Negative</th>
                    <th class="text-center">Organism</th>
                    <th class="text-center">Review date</th>
                    <th class="text-center">Download</th>
                </tr>
            </thead>
            <tbody>
                {#each selectedModelReports as reportEntry}

                    <tr class="text-xs hover:cursor-pointer align-center">
                        <td class="text-center opacity-60">
                            {getDateTimeStringUtc(reportEntry.date)}
                        </td>
                        <td class="text-center opacity-60">
                            {reportEntry.id.substring(0, 8)}
                        </td>
                        <td class="text-center opacity-60">
                            {reportEntry.user_name}
                        </td>
                        <td class="text-center opacity-60">
                            {reportEntry.negative}
                        </td>
                        <td class="text-center opacity-60 {reportEntry.organism === "NEGATIVE" ? "" : "italic"}">
                            {reportEntry.organism}
                        </td>
                        <td class="text-center opacity-60">
                            {reportEntry.review_date}
                        </td>
                        <td class="text-center">
                            <button class="btn variant-outline-primary text-xs p-1" type="button" on:click={() => downloadReportCopy(reportEntry.id)}>
                                <div class="h-3 w-3">
                                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" stroke-linecap="round" stroke-linejoin="round"></path>
                                      </svg>
                                </div>
                            </button>
                        </td>
                    </tr>

                {/each}
            </tbody>
        </table>
        {:else}
            <p class="flex p-8 text-sm opacity-30">
                No reports have been generated for selected libraries of sample <span class="font-semibold ml-1">{$page.params.sample}</span>
            </p>
        {/if}
</div>