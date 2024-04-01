<script lang="ts">
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    export let selectedView: string;
 
</script>
<p class="opacity-80 p-2">Database collections</p>
<ListBox class="border p-4 border-primary-500">
	<ListBoxItem bind:group={selectedView} name="medium" value="samples" active='variant-soft' rounded='rounded-token'>Samples & Libraries</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="runs", active='variant-soft' rounded='rounded-token'>Sequencing Runs</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="reports", active='variant-soft' rounded='rounded-token'>Clinical Reports</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="reports", active='variant-soft' rounded='rounded-token'>Pathogens</ListBoxItem>
</ListBox>