<script lang="ts">
    import { page } from "$app/stores";
	import UserSettingsSelection from "$lib/components/user/UserSettingsSelection.svelte";
    
    let selectedView: string = "team";

    let selectedTeamId: string = $page.data.userTeams[0];

    function changeTeam() {
    }
</script>


<div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-4 gap-16 pt-10 ">
    <div class="col-span-1">

        <div class="grid grid-rows-2 gap-8">
            <div class="row-span-1">
                <UserSettingsSelection selectedView={selectedView} />
            </div>
            <div class="row-span-1">
                <div class="border p-4 border-primary-500 flex gap-4 rounded-2xl">
                    <div class="w-1/2">
                        <p class="mb-1"><span class="opacity-60">Team</span></p>
                        <select id="teamSelect" class="select" bind:value={selectedTeamId} on:change={changeTeam}>
                            {#each $page.data.userTeams as team}
                                <option value={team.id}>{team.name}</option>
                            {/each}
                        </select>
                    </div>
                    <div class="w-1/2 mt-6">
                        <p class="opacity-40 text-sm">Change settings and configurations shared with your team</p>
                    </div>
                </div>
                
            </div>
        </div>
    </div>
    <div class="col-span-1 md:col-span-3">

        <div class="border p-4 border-primary-500 flex gap-4 rounded-2xl">
        
        </div>
    </div>
</div>
