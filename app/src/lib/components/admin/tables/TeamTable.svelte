<script lang="ts">
	import { page } from "$app/stores";
    import ClientPaginatedTable from "$lib/general/tables/ClientPaginatedTable.svelte";
    
    // Use this for binding data when updating 
    // page data on successful calls

    export let selectedRow: string[] = [];
    
    let columns: string[] = ["name", "description"];
    let header: string[] = ["Team", "Description"];

    // Columns to select from the underlying team data
    let selectColumns: string[] = ["id"];    
    let showEntries: number[] = [5, 10, 20, 50]

</script>
    
<ClientPaginatedTable 
    data={$page.data.teams} 
    columns={columns} 
    header={header} 
    showEntries={showEntries} 
    show={10} 
    interactive={true} 
    selectColumns={selectColumns} 
    bind:selectedRow={selectedRow}
></ClientPaginatedTable>
