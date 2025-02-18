<script lang="ts">
	import Report from "$lib/components/production/reports/Report.svelte";
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    let selectedView: string = "report";

</script>

<div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-4 gap-16 pt-10">
    <div class="gap-y-4 md:mt-5 col-span-1 text-lg">
        <div class="mb-4">
            <p class="opacity-80 p-2">Clinical Report</p>
                <ListBox class="border-token p-4 border-primary-500">
                    <ListBoxItem bind:group={selectedView} name="medium" value="report" active='variant-soft' rounded='rounded-token'>New Report</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="configuration", active='variant-soft' rounded='rounded-token' disabled>Configuration</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="appendixA", active='variant-soft' rounded='rounded-token' disabled>Appendix A: Laboratory</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="appendixB", active='variant-soft' rounded='rounded-token' disabled>Appendix B: Bioinformatics</ListBoxItem>
                </ListBox>
        </div>
    </div>
    <div class="col-span-3 pt-12">
        <Report selectedView={selectedView}></Report>
    </div>
</div>