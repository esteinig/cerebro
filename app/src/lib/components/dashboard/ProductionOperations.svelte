<script lang="ts">

    import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    let selectedView: string = "samples";

</script>

<p class="opacity-60">Production interface</p>
<div class="card border p-4 border-primary-500">
    <ListBox class="">
        <ListBoxItem bind:group={selectedView} name="medium" value="samples" active='variant-soft' rounded='rounded-token'>Samples & Reports</ListBoxItem>
        <ListBoxItem bind:group={selectedView} name="medium" value="runs", active='variant-soft' rounded='rounded-token'>Sequencing Runs </ListBoxItem>
        <ListBoxItem bind:group={selectedView} name="medium" value="taxdb", active='variant-soft' rounded='rounded-token'>Contamination</ListBoxItem>
        <ListBoxItem bind:group={selectedView} name="medium" value="taxdb", active='variant-soft' rounded='rounded-token'>Taxonomicon</ListBoxItem>
    </ListBox>
</div>