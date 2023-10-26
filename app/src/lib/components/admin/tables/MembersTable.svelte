<script lang="ts">
	import { invalidate } from "$app/navigation";
	import { page } from "$app/stores";
	import { CerebroApi, ApiResponse } from "$lib/utils/api";
	import type { Team, User } from "$lib/utils/types";
	import { getToastStore, ProgressRadial } from "@skeletonlabs/skeleton";

    export let users: User[];
    export let team: Team;
    export let update: boolean = false;
    export let selectedRow: string[] = [];

    let loading: boolean = false

    let publicApi = new CerebroApi();
        const toastStore = getToastStore();

    const removeTeamMember = async(userId: string) => {

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.teams.removeUser}?user_id=${userId}&team_id=${team.id}&update=remove`, {
                method: 'PATCH',
                mode: 'cors',
                credentials: 'include',
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, "Removed team member"
        );

        if (response.ok){
            await invalidate("admin:data"); // reload page data
            selectedRow = [team.id]; // refresh current card
        }
        loading = false;

    }

</script>

<div class="table-container">
    {#if loading}
        <div class="flex justify-center">
            <ProgressRadial width="sm:w-12 md:w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
        </div>
    {:else}
        <table class="table table-hover table-compact" >
            <tbody>
                {#each users as user, i}
                    <tr class="">
                        <td class="">{user.title ? user.title: ""} {user.name}</td>
                        <td>{user.positions.join(", ")}</td>
                        <td>{user.roles.join(", ")}</td>
                        <td>{user.email}</td>
                        <td>
                            {#if update}
                                <div class="flex justify-end">
                                    <button type="button" class="btn-xs variant-outline-none p-1 text-xs" on:click={() => removeTeamMember(user.id)}>Remove</button>
                                </div>
                            {:else}
                                {#if user.verified}
                                    <svg aria-hidden="true" class="h-4 w-4" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                                    </svg>
                                {:else}
                                    <svg aria-hidden="true" class="h-4 w-4" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M6 18L18 6M6 6l12 12" stroke-linecap="round" stroke-linejoin="round"></path>
                                    </svg>
                                {/if}

                            {/if}
                            
                        </td>
                    </tr>
                {/each}
            </tbody>
        </table>
    {/if}
	
</div>