
<script lang="ts">

    import Search from '$lib/general/datatable/components/Search.svelte';
    import ThFilter from '$lib/general/datatable/components/ThFilter.svelte';
    import ThSort from '$lib/general/datatable/components/ThSort.svelte';
    import RowCount from '$lib/general/datatable/components/RowCount.svelte';
    import RowsPerPage from '$lib/general/datatable/components/RowsPerPage.svelte';
    import Pagination from '$lib/general/datatable/components/Pagination.svelte';

    import { DataHandler } from '@vincjo/datatables';

    // Init data handler - CLIENT
    let data = [
        { id: 1, first_name: 'Tobie', last_name: 'Vint', email: 'tvint0@fotki.com' },
        { id: 2, first_name: 'Zacharias', last_name: 'Cerman', email: 'zcerman1@sciencedirect.com' },
        { id: 3, first_name: 'Gérianna', last_name: 'Bunn', email: 'gbunn2@foxnews.com' },
        { id: 4, first_name: 'Bee', last_name: 'Saurin', email: 'bsaurin3@live.com' },
        { id: 5, first_name: 'Méyère', last_name: 'Granulette', email: 'mgranul4@yellowbook.com' }
    ];
    
    const handler = new DataHandler(data, { rowsPerPage: 5 });
    const rows = handler.getRows();
</script>

<div class=" overflow-x-auto space-y-4">
<!-- Header -->
<header class="flex justify-between gap-4">
    <Search {handler} />
    <RowsPerPage {handler} />
</header>
<!-- Table -->
<table class="table table-hover table-compact w-full table-auto">
    <thead>
        <tr>
            <ThSort {handler} orderBy="first_name">First name</ThSort>
            <ThSort {handler} orderBy="last_name">Last name</ThSort>
            <ThSort {handler} orderBy="email">Email</ThSort>
        </tr>
        <tr>
            <ThFilter {handler} filterBy="first_name" />
            <ThFilter {handler} filterBy="last_name" />
            <ThFilter {handler} filterBy="email" />
        </tr>
    </thead>
    <tbody>
        {#each $rows as row}
            <tr>
                <td>{row.first_name}</td>
                <td>{row.last_name}</td>
                <td>{row.email}</td>
            </tr>
        {/each}
    </tbody>
</table>
<!-- Footer -->
<footer class="flex justify-between">
    <RowCount {handler} />
    <Pagination {handler} />
</footer>
</div>