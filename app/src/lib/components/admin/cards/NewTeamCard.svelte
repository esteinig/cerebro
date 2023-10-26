
<script lang="ts">
	import { invalidate } from "$app/navigation";
    import { page } from "$app/stores";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
    import { sanitizeMongoDatabaseName } from "$lib/utils/helpers";
	import type { RegisterTeamSchema } from "$lib/utils/types";
    import { ProgressRadial, getToastStore } from '@skeletonlabs/skeleton';

    let publicApi: CerebroApi = new CerebroApi();
    const toastStore = getToastStore();

    let name: string;
    let lead: string;
    let description: string;


    // Use loading so we can invalidate and refresh
    // page data on successful team creation
    let loading: boolean = false;

    const registerTeam = async() => {

        loading = true;

        let teamRegistration: RegisterTeamSchema = {
            teamLead: lead,
            teamName: name,
            teamDescription: description,
            databaseName: name,
            databaseMongoName: sanitizeMongoDatabaseName(name),
            databaseDescription: 'Default team database',
            projectName: `Default`,
            projectMongoName: `default`,
            projectDescription: `Default sample collection`
        };

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            publicApi.routes.teams.register, {
                method: 'POST',
                mode: 'cors',
                credentials: 'include',
                body: JSON.stringify(teamRegistration),
                headers:  { 'Content-Type': 'application/json' }
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, "New team created"
        );

        if (response.ok){
            await invalidate("admin:data"); // Reload page data
            await invalidate("cerebro:data"); // Reload user base data
        }
        loading = false;

    }  

</script>

<div class="card max-w-xl2">
	<header class="card-header">
        <div class="p-4 space-y-4">
            <h4 class="h4">Create a new team</h4>
            <p class="text-sm opacity-60">
                Teams are assigned their own database. Team members have shared access to data collections (projects)
                and can collaborate via comments and taxon selections in the sample view. Additional team members can be 
                added later.
            </p>
        </div>
    </header>
	<section class="p-4">
        {#if loading}
        <div class="flex justify-center p-4 space-y-4">
            <ProgressRadial width="w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
        </div>
        {:else}
            <form class="p-4 space-y-4" on:submit|preventDefault={registerTeam}>
                <label class="label">
                    <span>Name</span>
                    <input type="text" placeholder="Enter team name" class="input" bind:value={name} required/>
                </label>
                <label class="label">
                    <span>Lead</span>
                    <select class="select" bind:value={lead}>
                        {#each $page.data.users as user}
                            <option value={user.id}>{user.name}</option>
                        {/each}
                    </select>
                </label>
                <label class="label">
                    <span>Description</span>
                    <textarea class="textarea" rows="4" placeholder="Enter team description" bind:value={description} required/>
                </label>
                <div class="flex justify-center gap-4 mt-2">
                    <button type="submit" class="btn btn-lg variant-outline-tertiary">Create</button>
                </div>
            </form>
        {/if}    

    </section>
	<footer class="card-footer">
        
    </footer>
</div>