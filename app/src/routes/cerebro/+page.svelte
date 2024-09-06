<script lang="ts">
	import Announcements from "$lib/components/dashboard/Announcements.svelte";
	import ChangeLog from "$lib/components/dashboard/ChangeLog.svelte";
	import PipelineDeployment from "$lib/components/dashboard/PipelineDeployment.svelte";
	import ProductionOperations from "$lib/components/dashboard/ProductionOperations.svelte";
	import TeamActivity from "$lib/components/dashboard/TeamActivity.svelte";
    import UserProfileSummary from "$lib/components/dashboard/UserProfileSummary.svelte";
</script>


 
<div class="container mt-[2.5%] mx-auto max-w-[90%] max-h-[90%] h-full">
    
    <div class="grid sm:grid-cols-1 md:grid-cols-8 gap-8 h-full">
        <div class="md:col-span-2">
            <div class="grid grid-rows-2 md:grid-rows-2 gap-y-12">
                <div class="row-span-1">
                    <UserProfileSummary></UserProfileSummary>
                </div>
                <div class="row-span-1">
                    <div class="mb-8">
                        <ProductionOperations></ProductionOperations>   
                    </div>
                    <div>
                        <PipelineDeployment></PipelineDeployment>      
                    </div>   
                </div>
            </div>
        </div>
        <div class="md:col-span-4">
            <div class="grid grid-rows-1 md:grid-rows-1 gap-y-8 h-full">
                <div class="row-span-1 blur-sm">
                    <TeamActivity></TeamActivity>
                </div>
            </div>
        </div>

        <div class="md:col-span-2">
            <div class="grid grid-rows-2 md:grid-rows-2 gap-y-8 ">
                <div class="row-span-1 blur-sm">
                    <Announcements></Announcements>
                </div>
                <div class="row-span-1 blur-sm">
                    <ChangeLog></ChangeLog>
                </div>
            </div>
        </div>
    </div>
</div>

