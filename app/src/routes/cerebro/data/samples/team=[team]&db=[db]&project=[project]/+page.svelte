<script lang="ts">
    import SampleOverview from "$lib/components/data/samples/SampleOverview.svelte";

    export let data;
    
</script>

<SampleOverview sampleOverviewData={data.sampleOverview}></SampleOverview>