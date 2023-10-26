<script lang="ts">

    import { page } from "$app/stores";
    import ClientPaginatedTable from "$lib/general/tables/ClientPaginatedTable.svelte";

    export let selectedRow: string[] = [];

    let columns: string[] = ["name", "positions", "roles", "email", "verified"];
    let header: string[] = ["User", "Positions", "Permissions", "Email", "Verified"];

    let selectColumns: string[] = ["id"];
    
</script>
    
<ClientPaginatedTable 
    data={$page.data.users} 
    columns={columns} 
    header={header} 
    interactive={true} 
    selectColumns={selectColumns} 
    bind:selectedRow={selectedRow}
></ClientPaginatedTable>