<script lang="ts">
	import Announcements from "$lib/components/dashboard/Announcements.svelte";
	import CerebroRelease from "$lib/components/dashboard/CerebroRelease.svelte";
	import TeamActivity from "$lib/components/dashboard/TeamActivity.svelte";
    import UserProfileSummary from "$lib/components/dashboard/UserProfileSummary.svelte";
	import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    let selectedView: string = "samples";
</script>


 
<div class="container mt-[2.5%] mx-auto max-w-[90%] max-h-[90%] h-full">
    
    <div class="grid sm:grid-cols-1 md:grid-cols-8 gap-8 h-full">
        <div class="md:col-span-2">
            <div class="grid grid-rows-2 md:grid-rows-2 gap-16">
                <div class="row-span-1">
                    <UserProfileSummary></UserProfileSummary>
                </div>
                <div class="row-span-1">
                    <p class="opacity-60">Navigation</p>
                    <div class="card border p-4 border-primary-500">
                        <ListBox class="">
                            <ListBoxItem bind:group={selectedView} name="medium" value="samples" active='variant-soft' rounded='rounded-token'>Samples & Clinical Reports</ListBoxItem>
                            <ListBoxItem bind:group={selectedView} name="medium" value="runs", active='variant-soft' rounded='rounded-token'>Sequence Runs & Data </ListBoxItem>
                            <ListBoxItem bind:group={selectedView} name="medium" value="taxdb", active='variant-soft' rounded='rounded-token'>Taxonomy & Databases</ListBoxItem>
                            <ListBoxItem bind:group={selectedView} name="medium" value="qa", active='variant-soft' rounded='rounded-token'>Quality Assurance</ListBoxItem>
                            <ListBoxItem bind:group={selectedView} name="medium" value="tools", active='variant-soft' rounded='rounded-token'>Other Tools</ListBoxItem>
                        
                        </ListBox>
                    </div>                    
                </div>
            </div>
        </div>
        <div class="md:col-span-4">
            <div class="grid grid-rows-1 md:grid-rows-1 gap-y-8 h-full">
                <div class="row-span-1">
                    <TeamActivity></TeamActivity>
                </div>
            </div>
        </div>

        <div class="md:col-span-2">
            <div class="grid grid-rows-1 md:grid-rows-2 gap-y-8 h-full">
                <div class="row-span-1">
                    <Announcements></Announcements>
                </div>
                <div class="row-span-1">
                    <CerebroRelease></CerebroRelease>
                </div>
            </div>
        </div>
    </div>
</div>

