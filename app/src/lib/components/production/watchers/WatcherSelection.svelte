<script lang="ts">
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    export let selectedView: string;
 
</script>

<ListBox class="border p-4 border-primary-500">
	<ListBoxItem bind:group={selectedView} name="medium" value="team" active='variant-soft' rounded='rounded-token'>Team Settings</ListBoxItem>
	<ListBoxItem bind:group={selectedView} name="medium" value="user", active='variant-soft' rounded='rounded-token'>User Settings</ListBoxItem>
</ListBox>