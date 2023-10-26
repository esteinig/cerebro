<script lang="ts">
    import TeamTable from "./tables/TeamTable.svelte";
    import TeamCard from "./cards/TeamCard.svelte";
    import NewTeamCard from "./cards/NewTeamCard.svelte";
	import { page } from "$app/stores";
	import type { Team } from "$lib/utils/types";
	import { getToastStore, type ToastSettings } from "@skeletonlabs/skeleton";

    const toastStore = getToastStore();

    let newTeamCard: boolean = true;
    let selectedRow: string[] = [];

    let selectedTeam: Team | null = null;

    const getSelectedTeam = (selectedRow: string[]): Team | null => {

        if (selectedRow.length > 0){
            let teamId: string = selectedRow[0];

            let match = $page.data.teams.filter(
                (team: Team) => team.id === teamId
            )[0]
            if (match === undefined) {
                const errorToast: ToastSettings = {
                    message: "Failed to select a RequestLog entry!",
                    background: 'variant-filled-warning',
                };
                toastStore.trigger(errorToast);
                return null
            }
            return match
        }
        return null
    }   

    $: {
        selectedTeam = getSelectedTeam(selectedRow);
        if (selectedTeam !== null){
            newTeamCard = false;
        }
    }



</script>

<div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-16 pt-10">
    <div class="">
        <TeamTable bind:selectedRow={selectedRow}/>
        <div class="flex justify-center gap-4 mt-2">
            <button type="button" class="btn variant-outline-success" on:click={() => newTeamCard = true}>New Team</button>
        </div>
    </div>
    <div class="">
        {#if newTeamCard}
            <NewTeamCard />
        {:else}
            {#if selectedTeam !== null}
                <TeamCard bind:team={selectedTeam} bind:selectedRow/>
            {/if}
        {/if}
    </div>
</div>