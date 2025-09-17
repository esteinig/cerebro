<script lang="ts">
	import { goto } from "$app/navigation";
    import { page } from "$app/stores";

	import TrainingOverview from "$lib/components/training/TrainingOverview.svelte";

    $: selectedTeamName = $page.data.selectedTeam.name;

    async function changeTeam() {
        await goto(`/cerebro/training/team=${selectedTeamName}`, { invalidateAll: true })
    }

</script>

<div>
	<div class="w-1/4">
		<p class="mb-1"><span class="opacity-60">Team</span></p>
		<select id="teamSelect" class="select" bind:value={selectedTeamName} on:change={changeTeam}>
			{#each $page.data.userTeams as team}
				<option value={team.name}>{team.name}</option>
			{/each}
		</select>
	</div>
	<TrainingOverview data={$page.data.trainingOverview}></TrainingOverview>
</div>