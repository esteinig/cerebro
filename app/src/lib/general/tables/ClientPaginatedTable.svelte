<script lang="ts">
import { Table, tableMapperValues } from '@skeletonlabs/skeleton';
import type { TableSource } from '@skeletonlabs/skeleton';

import { Paginator } from '@skeletonlabs/skeleton';
import type { PaginationSettings } from '@skeletonlabs/skeleton';

export let data: any[] = [];
export let columns: string[] = [];
export let header: string[]  = [];

export let interactive: boolean = false;
export let selectedRow: any[] = [];
export let selectColumns: string[] = [];

export let showEntries: number[] = [1, 2, 5, 10];
export let show: number = 10;

export let pagination: boolean = true;
export let regionHead: string | undefined = undefined;
export let regionBody: string| undefined = undefined;

let tableSimple: TableSource = { head: header, body: [], meta: [] }

let paginationSettings: PaginationSettings = {
	page: 0,
	limit: show,
	size: data.length,
	amounts: showEntries,
} satisfies PaginationSettings;

$: {
	let tableData = data;
	if (pagination) {
		paginationSettings.size = data.length;
		tableData = data.slice(
			paginationSettings.page * paginationSettings.limit,
			paginationSettings.page * paginationSettings.limit + paginationSettings.limit
		);
	}
	
	tableSimple.body = tableMapperValues(tableData, columns);
	tableSimple.meta = tableMapperValues(tableData, selectColumns);
}

const selectionHandler = (element: any) => {
	selectedRow = element.detail
}

</script>


<Table source={tableSimple} interactive={interactive} on:selected={selectionHandler} regionHead={regionHead} regionBody={regionBody}/>
{#if pagination}
<Paginator bind:settings={paginationSettings} showFirstLastButtons={false} showPreviousNextButtons={true} class="mt-2" select="paginator-select select text-xs" regionControl="opacity-30 dark:variant-filled-surface dark:opacity-60 text-xs"/>
{/if}