<script lang="ts">
	import { invalidate } from "$app/navigation";
	import { page } from "$app/stores";
	import { CerebroApi, type ApiResponse } from "$lib/utils/api";
	import { getDateTimeString, getDateTimeStringUtc } from "$lib/utils/helpers";
	import type { Cerebro, ReportEntry } from "$lib/utils/types";
    import { Role } from "$lib/utils/types";
	import { getToastStore, type PopupSettings, type ToastSettings } from "@skeletonlabs/skeleton";
    import { popup } from '@skeletonlabs/skeleton';

    export let selectedModels: Cerebro[];

    let selectedModelReports: ReportEntry[];

    $: selectedModelReports = selectedModels.flatMap(cerebro => cerebro.sample.reports).filter((value, index, self) => {
		return self.findIndex(v => v.id === value.id) === index;
	})

    const publicApi: CerebroApi = new CerebroApi();
    const toastStore = getToastStore();

    let loading: boolean = false;

    const deleteReportCopy = async(report_id: string) => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.deleteReport}/${report_id}?db=${$page.params.db}&project=${$page.params.project}`,
            { 
                method: 'DELETE',  
                mode: 'cors',
                credentials: 'include',
            } as RequestInit,
            $page.data.refreshToken, toastStore, null
        )
        loading = false;

        if (response.ok) {
            toastStore.trigger({
                message: response.json.message ?? "Failed to delete report",
                background: "variant-filled-tertiary"
            } satisfies ToastSettings);

            invalidate("sample:data");
        }

    }

    const downloadReportCopy = async(report_id: string) => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.getReport}/${report_id}?db=${$page.params.db}&project=${$page.params.project}`,
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

    const deletePopup: PopupSettings = {
        event: 'hover',
        target: 'deletePopup',
        placement: 'right'
    };
    const downloadPopup: PopupSettings = {
        event: 'hover',
        target: 'downloadPopup',
        placement: 'left'
    };

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
                            <button class="btn variant-outline-primary text-xs p-1 [&>*]:pointer-events-none" type="button" id="downloadPopup" on:click={() => downloadReportCopy(reportEntry.id)} use:popup={downloadPopup}>
                                <div class="h-3 w-3">
                                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" stroke-linecap="round" stroke-linejoin="round"></path>
                                    </svg>
                                </div>
                            </button>
                            {#if $page.data.userData.roles.includes(Role.Data)}
                                <button class="btn variant-outline-tertiary text-xs p-1 ml-2  [&>*]:pointer-events-none" type="button" id="deletePopup" on:click={() => deleteReportCopy(reportEntry.id)} use:popup={deletePopup}>
                                    <div class="h-3 w-3">
                                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                            <path d="M14.74 9l-.346 9m-4.788 0L9.26 9m9.968-3.21c.342.052.682.107 1.022.166m-1.022-.165L18.16 19.673a2.25 2.25 0 01-2.244 2.077H8.084a2.25 2.25 0 01-2.244-2.077L4.772 5.79m14.456 0a48.108 48.108 0 00-3.478-.397m-12 .562c.34-.059.68-.114 1.022-.165m0 0a48.11 48.11 0 013.478-.397m7.5 0v-.916c0-1.18-.91-2.164-2.09-2.201a51.964 51.964 0 00-3.32 0c-1.18.037-2.09 1.022-2.09 2.201v.916m7.5 0a48.667 48.667 0 00-7.5 0" stroke-linecap="round" stroke-linejoin="round"></path>
                                        </svg>                               
                                    </div>
                                </button>
                            {/if}
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