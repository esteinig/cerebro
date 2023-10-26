<script lang="ts">
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    export let selectedView: string;
 
</script>

<ListBox>
	<ListBoxItem bind:group={selectedView} name="medium" value="qc" active='variant-soft' rounded='rounded-token'>Quality Control</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="classification", active='variant-soft' rounded='rounded-token'>Classification</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="candidates", active='variant-soft' rounded='rounded-token'>Candidates</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="comments", active='variant-soft' rounded='rounded-token'>Comments</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="reports", active='variant-soft' rounded='rounded-token'>Reports</ListBoxItem>
</ListBox>