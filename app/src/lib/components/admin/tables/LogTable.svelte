<script lang="ts">
	import { page } from "$app/stores";
    import { getDateTime } from "$lib/utils/helpers";
    import { InputChip, SlideToggle } from '@skeletonlabs/skeleton';
    import ClientPaginatedTable from "$lib/general/tables/ClientPaginatedTable.svelte";

	import type { RequestLog } from "$lib/utils/types";

    export let selectedRow: string[] = [];

    let showCriticalLogs: boolean = false;

    interface LogData {
        id: string,
        date: string,
        time: string,
        module: string,
        action: string
    }

    let data: LogData[] = [];
    let filteredData: LogData[] = data;
    let searchTerms: string[] = [];

    $: {
        let baseData = showCriticalLogs ? $page.data.criticalLogs : $page.data.logs;

        data = baseData.map((log: RequestLog) => {
            let [date, time] = getDateTime(log.date);
            return {
                id: log.id,
                date: date,
                time: time,
                module: log.module,
                action: log.action,
                email: log.access_details.user_email
            } as LogData
        });
    }

    $: {
        if (searchTerms.length == 0) {
            filteredData = data;
        } else {
            filteredData = data.filter(logData => {
                let tableFields: string[] = [
                    logData.time, 
                    logData.date,  
                    logData.module,
                    logData.action,
                ]
                let decision = searchTerms.every(term => tableFields.includes(term));
                return decision
            });
        }
        selectedRow = filteredData.length > 0 ? [filteredData[0].id] : [];
    }

    $: {
        if (showCriticalLogs) {
            selectedRow = filteredData.length > 0 ? [filteredData[0].id] : [];
        }
    }

    let columns: string[] = [
        "date", 
        "time",
        "module", 
        "action",
        "email",
    ];
    let header: string[] = [
        "Date", 
        "Time", 
        "Module", 
        "Action",
        "email"
    ];
    let showEntries: number[] = [
        10, 
        20, 
        30, 
        50, 
        100
    ];
    
    let selectColumns: string[] = [
        "id"
    ];
    
</script>
<div class="flex items-center"> 
    
    <InputChip name="logSearch" placeholder="Search request logs..." class="mb-5 w-1/2" bind:value={searchTerms} allowUpperCase/> 
    <button class="btn btn-md variant-outline-primary mb-5 ml-5">
        <div class="flex items-center gap-2">
            <svg aria-hidden="true" fill="none" class="h-5 w-5" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
            <path d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" stroke-linecap="round" stroke-linejoin="round"></path>
        </svg>
        Download
        </div>
    </button>
    <div class="ml-auto mb-2">
        <SlideToggle name="critical-logs" bind:checked={showCriticalLogs} active="variant-filled-error dark:variant-filled-error">Critical Logs</SlideToggle>
    </div>
</div>
<ClientPaginatedTable data={filteredData} columns={columns} header={header} showEntries={showEntries} show={10} interactive={true} selectColumns={selectColumns} bind:selectedRow={selectedRow}></ClientPaginatedTable>
