
<script lang="ts">
	import { invalidate } from "$app/navigation";
	import { page } from "$app/stores";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import type { Team, UpdateTeamSchema, User } from "$lib/utils/types";
	import { getToastStore, ProgressRadial, type ModalSettings, getModalStore } from "@skeletonlabs/skeleton";
	import MembersTable from "../tables/MembersTable.svelte";

    export let team: Team;
    export let selectedRow: string[] = [];

    const toastStore = getToastStore();
    const modalStore = getModalStore();

    let publicApi: CerebroApi = new CerebroApi();
    let loading: boolean = false;

    let update: boolean = false;
    let updateDescription: string = "";
    let updateMembers: string[] = []


    const deleteTeam = async() => {

        loading = true;
        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.teams.delete}/${team.id}`, {
                method: 'DELETE',
                mode: 'cors',
                credentials: 'include'
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, "Team deleted"
        );

        if (response.ok){
            await invalidate("admin:data"); // Reload admin route data
            await invalidate("cerebro:data"); // Reload user base data
            selectedRow = [];  // Deselect upstream after deletion
        }
        loading = false;

    }  

    const updateTeam = async() => {

        loading = true;

        let updateTeamSchema: UpdateTeamSchema = {
            teamName: null, // do not allow team name updates
            teamDescription: updateDescription,
        }

        await publicApi.fetchWithRefresh(
            `${publicApi.routes.teams.update}/${team.id}`, {
                method: 'PATCH',
                mode: 'cors',
                credentials: 'include',
                body: JSON.stringify(updateTeamSchema),
                headers:  { 'Content-Type': 'application/json' }
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, "Team updated"
        );
        
        await Promise.all<Promise<ApiResponse>>(
            updateMembers.map(userId => addTeamMember(userId))
        );

        await invalidate("admin:data"); // Reload page data
        await invalidate("cerebro:data"); // Reload user base data

        selectedRow = [team.id]; // Reload of card data
        loading = false;
        
        
    }


    const addTeamMember = async(userId: string): Promise<ApiResponse> => {

        return publicApi.fetchWithRefresh(
            `${publicApi.routes.teams.removeUser}?user_id=${userId}&team_id=${team.id}&update=add`, {
                method: 'PATCH',
                mode: 'cors',
                credentials: 'include',
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, "Team member added"
        );

    }



    // Associate the user identifiers with names and actual users
    let members: User[] = []

    $: {
        selectedRow;
        update = false;
        updateDescription = team.description;

    }

    $: {
        members = team.users.map(id => {
            return $page.data.users.filter((user: User) => user.id === id)
        }).flatMap((x: User) => x ? x : [])
    }

    const openTeamDeleteVerification = async() => {

        const modal: ModalSettings = {
            type: 'confirm',
            title: 'Confirm Deletion',
            body: `<p class="opacity-60 my-2">Delete the selected team?</p><p><span class="font-semibold text-red-500">All databases and projects associated with this team will be deleted.</span></p>`,
            response: async(confirmed: boolean) => {
                if (confirmed) {
                    await deleteTeam()
                }
            }
        };
        modalStore.trigger(modal);
    }

</script>

<div class="card max-w-xl2">
	<header class="card-header">
        <div class="p-4 space-y-4">
            <div class="flex">
                <h2 class="h2">{team.name}</h2>
                <div class="ml-auto">
                    {#if update}
                        <button type="submit" class="btn btn-sm variant-outline-success m-2 " on:click={() => { updateTeam() }}>Submit</button>
                        <button type="submit" class="btn btn-sm variant-outline-warning m-2 " on:click={() => update = false}>Cancel</button>
                    {:else}
                        <button type="submit" class="btn btn-sm variant-outline-warning m-2 " on:click={() => update = true}>Update</button>
                        <button type="submit" class="btn btn-sm variant-outline-error m-2 " on:click={openTeamDeleteVerification}>Delete</button>
                    {/if}
                </div>
            </div>
        </div>
    </header>
	<section class="p-4">
        {#if loading}
            <div class="flex justify-center p-4 space-y-4">
                <ProgressRadial width="w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
            </div>
        {:else}
            <div class="px-4 space-y-2">
                <p><span class="opacity-60">Description</span></p>
                {#if update}
                    <textarea class="textarea" rows="4" placeholder="{team.description}" bind:value={updateDescription} required/>
                {:else}
                    <p> {team.description}</p>
                {/if}
            </div>
            
            <div class="px-4 pt-4 space-y-2">
                <p><span class="opacity-60">Databases</span></p>
                {#each team.databases as db}
                    <div class="flex mt-2 text-xs">
                        <svg class="w-5 h-5" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M20.25 6.375c0 2.278-3.694 4.125-8.25 4.125S3.75 8.653 3.75 6.375m16.5 0c0-2.278-3.694-4.125-8.25-4.125S3.75 4.097 3.75 6.375m16.5 0v11.25c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125V6.375m16.5 0v3.75m-16.5-3.75v3.75m16.5 0v3.75C20.25 16.153 16.556 18 12 18s-8.25-1.847-8.25-4.125v-3.75m16.5 0c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg>
                        <span class="code text-xs ml-2">{db.name}</span><span class="ml-2 text-xs">{db.description}</span>
                    </div>
                {/each}
            </div>

            <div class="px-4 pt-10 space-y-2">
                <p><span class="opacity-60">Members</span></p>
                <MembersTable bind:users={members} bind:update={update} team={team} bind:selectedRow={selectedRow}></MembersTable>
            </div>

            {#if update}
                <div class="px-4 pt-10 space-y-2">
                    <p><span class="opacity-60">New member selection </span><span class="ml-5 opacity-30 text-sm"><kbd class="kbd text-xs">CTRL</kbd> / <kbd class="kbd text-xs">⌘</kbd></span></p>
                    <select class="select" multiple bind:value={updateMembers}>
                        {#each $page.data.users as user}
                            <option value={user.id}>{user.name} ({user.email})</option>
                        {/each}
                    </select>
                </div>

            {/if}
        {/if}
    </section>
	<footer class="card-footer">

        <div class="flex justify-center gap-4 mt-2">
        </div>
    </footer>
</div>