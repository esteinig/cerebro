<script lang="ts">
	import { getInitials } from "$lib/utils/helpers";
    import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    let selectedView: string = "samples";
</script>

<div>
    <p class="opacity-60">Production tools</p>
    <div class="card w-full border border-primary-500">
        <section>
            <div class="p-4">
                <ListBox class="opacity-20">
                    <ListBoxItem bind:group={selectedView} name="medium" value="sheets", active='variant-soft' rounded='rounded-token'>Sample Sheet Editor</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="watchers", active='variant-soft' rounded='rounded-token'>Input Watchers</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="watchers", active='variant-soft' rounded='rounded-token'>File System and Backup</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="taxdb", active='variant-soft' rounded='rounded-token'>Databases and Taxonomy</ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="qa", active='variant-soft' rounded='rounded-token'>Quality Assurance Reports</ListBoxItem>
                </ListBox>
            </div>
        </section>
        <footer class="card-footer">
        </footer>
    </div>
</div>