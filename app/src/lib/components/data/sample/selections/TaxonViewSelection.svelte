<script lang="ts">
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    export let selectedView: string;
 
</script>

<ListBox>
	<ListBoxItem bind:group={selectedView} name="medium" value="qc" active='variant-soft' rounded='rounded-token'>Taxonomy</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="classification", active='variant-soft' rounded='rounded-token'>Comparison</ListBoxItem>
</ListBox>