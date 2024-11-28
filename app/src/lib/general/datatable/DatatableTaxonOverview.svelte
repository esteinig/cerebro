
<script lang="ts">

    import Search from '$lib/general/datatable/components/Search.svelte';
    import ThFilter from '$lib/general/datatable/components/ThFilter.svelte';
    import ThSort from '$lib/general/datatable/components/ThSort.svelte';
    import RowCount from '$lib/general/datatable/components/RowCount.svelte';
    import RowsPerPage from '$lib/general/datatable/components/RowsPerPage.svelte';
    import Pagination from '$lib/general/datatable/components/Pagination.svelte';

    import { DataHandler } from '@vincjo/datatables';
	import { type TaxonOverviewRecord, type TaxonOverview } from '$lib/utils/types';
    
    export let data: TaxonOverview[] = [];

    let taxonOverviewRecords = transformTaxonOverview(data, "rpm");
    let handler = new DataHandler(taxonOverviewRecords, { rowsPerPage: 5 });
    let rows = handler.getRows();


    function transformTaxonOverview(
        overviews: TaxonOverview[],
        field: 'reads' | 'rpm' | 'abundance'
    ): TaxonOverviewRecord[] {
        return overviews.map((overview) => {
            const aggregatedResults = {
                kraken2_sequence: 0,
                kraken2_profile: 0,
                metabuli_sequence: 0,
                metabuli_profile: 0,
                ganon2_sequence: 0,
                ganon2_profile: 0,
                kmcp_sequence: 0,
                kmcp_profile: 0,
                sylph_sequence: 0,
                sylph_profile: 0,
            };

            
            overview.evidence.forEach((result) => {
                const key = `${result.tool.toLowerCase()}_${result.mode.toLowerCase()}` as keyof typeof aggregatedResults;
                if (key in aggregatedResults) {
                    aggregatedResults[key] += result[field]; // Safe to add now since `field` is numeric
                }
            });

            return {
                taxid: overview.taxid,
                name: overview.name,
                domain: overview.domain,
                genus: overview.genus,
                sample_names: overview.sample_names,
                kmer: overview.kmer,
                alignment: overview.alignment,
                assembly: overview.assembly,
                ...aggregatedResults,
            };
        });
    }

    $: {
        taxonOverviewRecords = transformTaxonOverview(data, "rpm");
    }

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
            <ThSort {handler} orderBy="domain">Domain</ThSort>
            <ThSort {handler} orderBy="genus">Genus</ThSort>
            <ThSort {handler} orderBy="name">Taxon</ThSort>
            <ThSort {handler} orderBy="kraken2">Kraken2</ThSort>
            <ThSort {handler} orderBy="bracken">Bracken</ThSort>
            <ThSort {handler} orderBy="metabuli">Metabuli</ThSort>
            <ThSort {handler} orderBy="ganon2">Ganon2</ThSort>
            <ThSort {handler} orderBy="kmcp">Kmcp</ThSort>
            <ThSort {handler} orderBy="sylph">Sylph</ThSort>
        </tr>
        <tr>
            <ThFilter {handler} filterBy="domain" />
            <ThFilter {handler} filterBy="genus" />
            <ThFilter {handler} filterBy="name" />
            <ThFilter {handler} filterBy="kraken2_sequence" />
            <ThFilter {handler} filterBy="bracken_sequence" />
            <ThFilter {handler} filterBy="metabuli_sequence" />
            <ThFilter {handler} filterBy="ganon2_sequence" />
            <ThFilter {handler} filterBy="kmcp_sequence" />
            <ThFilter {handler} filterBy="sylph_sequence" />
        </tr>
    </thead>
    <tbody>
        {#each $rows as row}
            <tr>
                <td>{row.domain}</td>
                <td>{row.genus}</td>
                <td>{row.name}</td>
                <td>{row.kraken2_sequence}</td>
                <td>{row.bracken_sequence}</td>
                <td>{row.metabuli_sequence}</td>
                <td>{row.ganon2_sequence}</td>
                <td>{row.kmcp_sequence}</td>
                <td>{row.sylph_sequence}</td>
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