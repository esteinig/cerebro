<script lang="ts">
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    export let selectedView: string;
 
</script>

<ListBox>
	<ListBoxItem bind:group={selectedView} name="medium" value="samples" active='variant-soft' rounded='rounded-token'>Samples & Libraries</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="runs", active='variant-soft' rounded='rounded-token' disabled>Sequencing Runs</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="reports", active='variant-soft' rounded='rounded-token' disabled>Clinical Reports</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="reports", active='variant-soft' rounded='rounded-token' disabled>Pathogens</ListBoxItem>
</ListBox>