<script lang="ts">
    import { page } from "$app/stores";
	import LogCard from "./cards/LogCard.svelte";
	import LogHeaderCard from "./cards/LogHeaderCard.svelte";
    import LogTable from "./tables/LogTable.svelte";
    import { getToastStore } from '@skeletonlabs/skeleton';
	import type { RequestLog } from "$lib/utils/types";
    import type { ToastSettings } from '@skeletonlabs/skeleton';

    const toastStore = getToastStore();

    let selectedRow: string[] = [];
    let selectedRequestLog: RequestLog | null = null;
    let showHeaders: string[] = ["referer", "x-forwarded-for", "x-real-ip", "user-agent"]
    
    const getSelectedRequestLog = (selectedRow: string[]): RequestLog | null => {

        if (selectedRow.length > 0){
            // LogTable selects only the `id` field of the RequestLog
            let requestLogId: string = selectedRow[0];

            let match = $page.data.logs.filter(
                (log: RequestLog) => log.id === requestLogId
            )[0]
            if (match === undefined) {
                const errorToast: ToastSettings = {
                    message: "Failed to select a RequestLog entry!",
                    background: 'variant-filled-warning',
                };
                toastStore.trigger(errorToast);
                return null
            }
            return match
        }
        return null
    }   

    $: {
        selectedRequestLog = getSelectedRequestLog(selectedRow);
    }


</script>

<div class="grid grid-rows-2">
    <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-16 pt-10">
        <div>
            <LogTable bind:selectedRow={selectedRow} />
        </div>
        <div class="pt-16">
            {#if selectedRequestLog !== null}
                <LogCard bind:requestLog={selectedRequestLog}/>
            {/if}
        </div>
    </div>
    <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 gap-16 pt-10">
        <div>
            {#if selectedRequestLog !== null}
                <LogHeaderCard bind:requestLog={selectedRequestLog} showHeaders={showHeaders}></LogHeaderCard>
            {/if}
        </div>
    </div>
</div>