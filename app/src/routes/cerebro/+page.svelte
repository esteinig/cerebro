<script lang="ts">
    import { page } from "$app/stores";
	import Announcements from "$lib/components/dashboard/Announcements.svelte";
	import CerebroRelease from "$lib/components/dashboard/CerebroRelease.svelte";
	import TeamActivity from "$lib/components/dashboard/TeamActivity.svelte";
    import UserProfileSummary from "$lib/components/dashboard/UserProfileSummary.svelte";

    console.log($page.data.userTeams)
</script>


 
<div class="container mt-[2.5%] mx-auto max-w-[90%] max-h-[90%] h-full">
    
    <div class="grid sm:grid-cols-1 md:grid-cols-8 gap-8 h-full">
        <div class="md:col-span-2">
            <div class="grid grid-rows-1 md:grid-rows-1 gap-y-8">
                <div class="row-span-1">
                    <UserProfileSummary></UserProfileSummary>
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

