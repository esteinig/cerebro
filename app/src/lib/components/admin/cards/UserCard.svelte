<script lang="ts">

	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import type { User, UpdateUserScheme } from "$lib/utils/types";
	import { Avatar, getToastStore, ProgressRadial, type ModalSettings } from "@skeletonlabs/skeleton";
    import { getModalStore } from '@skeletonlabs/skeleton';
	import { invalidate } from "$app/navigation";
	import { getDateTime, getUserRoles, getInitials } from "$lib/utils/helpers";
	import { page } from "$app/stores";

    export let user: User;
    export let selectedRow: string[] = [];

    const toastStore = getToastStore();
    const modalStore = getModalStore();

    let publicApi: CerebroApi = new CerebroApi();
    let loading: boolean = false;
    let update: boolean = false;

    let [dateCreated, timeCreated]: string[] = [];
    let [dateUpdated, timeUpdated]: string[] = [];


    $: {
        [dateCreated, timeCreated] = getDateTime(user.created);
        [dateUpdated, timeUpdated] = getDateTime(user.updated);
    }

    let adminRole: boolean = false,
        reportRole: boolean = false,
        dataRole: boolean = false,
        botRole: boolean = false;

    let verified: boolean = false;
    let email: string = "";
    let name: string = "";
    let title: string | null = null;
    let position: string = "";


    let initials = getInitials(user.name);

    const deleteUser = async() => {

        loading = true;
        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.users.delete}/${user.id}`, {
                method: 'DELETE',
                mode: 'cors',
                credentials: 'include'
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, "User deleted"
        );

        if (response.ok){
            await invalidate("admin:data"); // Reload page data
            await invalidate("cerebro:data"); // Reload user base data
            selectedRow = [];  // Deselect upstream after deletion
        }
        loading = false; 

    }  
    const updateUser = async() => {

        let updateUserData: UpdateUserScheme = {
            name: name,
            email: email,
            title: title,
            password: null,
            positions: [position],
            verified: verified,
            roles: getUserRoles(adminRole, reportRole, dataRole, botRole),
        }

        loading = true;
        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.users.update}/${user.id}`, {
                method: 'PATCH',
                mode: 'cors',
                credentials: 'include',
                body: JSON.stringify(updateUserData),
                headers:  { 'Content-Type': 'application/json' }
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, "User updated"
        );

        if (response.ok){
            await invalidate("admin:data"); // Reload page data
            await invalidate("cerebro:data"); // Reload user base data
            selectedRow = [];  // Deselect upstream after deletion
        }
        loading = false; 
    }  

    const verifyUser = async() => {

        loading = true;
        let _: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.auth.verificationEmail}/${user.id}`, {
                method: 'GET',
                mode: 'cors',
                credentials: 'include'
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, "Verification email sent"
        );
        loading = false; 

    }  
    
    const openVerificationModal = async() => {

        const modal: ModalSettings = {
            type: 'confirm',
            title: 'Please Confirm',
            body: `Send verification and password email to: <p class="mt-2 font-bold">${user.email}</p>`,
            response: async(confirmed: boolean) => {
                if (confirmed) {
                    await verifyUser()
                }
            }
        };
        modalStore.trigger(modal);
    }

    $: {
        selectedRow;
        update = false;

        // Reassign changed user update view values

        verified = user.verified;
        email = user.email;
        name = user.name;
        title = user.title;
        position = user.positions.join(", ");

        user.roles.sort();
    }

</script>


<div class="card max-w-xl2">
	<header class="card-header">
        <div class="p-4">
            <div class="flex">
                <div class="flex items-center">
                    <Avatar src="" width="w-24" rounded="rounded-full" initials={initials} />
                    <div class="flex-none ml-5">
                        <div class="md:h4 lg:h3">{user.title ? user.title : ""} {user.name}</div>
                        <div class="flex items-center">
                            <h4 class="h4 opacity-60 pr-1">{user.positions.join(", ")}</h4>
                            {#if user.verified}
                                <svg aria-hidden="true" class="h-4 w-4 text-green-500" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                    <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                                </svg>
                            {:else}
                                <svg aria-hidden="true" class="h-4 w-4 text-red-500 pt-0.5" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                    <path d="M6 18L18 6M6 6l12 12" stroke-linecap="round" stroke-linejoin="round"></path>
                                </svg>
                            {/if}
                        </div>
                    </div>
                    

                </div>
            </div>
        </div>
    </header>
	<section class="p-4">
        {#if loading}
            <div class="flex justify-center p-4 space-y-4">
                <ProgressRadial width="sm:w-12 md:w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
            </div>
        {:else}
            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 lg:grid-cols-2 gap-x-2">
                <div class="px-10 pt-4 space-y-2">
                    <p><span class="opacity-60">Email</span></p>
                    {#if update}
                        <input type="text" placeholder="Email" class="input" value={email} required/>
                    {:else}
                        <p>{user.email}</p>
                    {/if}
                </div>
                <div class="px-10 pt-4 space-y-2">
                    {#if update}
                        <p class="opacity-60 text-sm">Changing email address forces re-verification by email or manual verification check by an administrator. Please note the security implications of manually verifying users.</p>
                        <p class="opacity-60 text-sm">Passwords can be changed by users on their profile pages or manually by the stack administrator in the database itself.</p>
                    {:else}
                        <p><span class="opacity-60">Registered</span></p>
                        <p class="text-sm">{dateCreated} @ {timeCreated} (UTC)</p>
                    {/if}
                </div>
                
            </div>
            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 lg:grid-cols-2 gap-x-2">
                <div class="px-10 pt-4 space-y-2">
                    <p><span class="opacity-60">Permissions</span></p>
                    {#if update}
                        <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 lg:grid-cols-2 gap-x-2">
                            <div class="p-2 space-y-2 bg-transparent dark:bg-transparent">
                                <label class="flex items-center space-x-2">
                                    <input class="checkbox" type="checkbox" checked={adminRole}/>
                                    <p>Admin</p>
                                </label>
                                <label class="flex items-center space-x-2">
                                    <input class="checkbox" type="checkbox" checked={botRole} />
                                    <p>Bot</p>
                                </label>
                                <label class="flex items-center space-x-2">
                                    <input class="checkbox" type="checkbox" checked={dataRole} />
                                    <p>Data</p>
                                </label>
                                <label class="flex items-center space-x-2">
                                    <input class="checkbox" type="checkbox" checked={reportRole} />
                                    <p>Report</p>
                                </label>
                            </div>
                            <div class="p-2">
                                <label class="flex items-center space-x-2 pt-2">
                                    <input class="checkbox" type="checkbox" checked={verified} />
                                    <p>Verified</p>
                                </label>
                                <div class="opacity-30 text-xs pt-2">
                                    Please note the security implications of updating verification
                                    status without verification email.
                                </div>
                            </div>
                        </div>
                    {:else}
                        <div class="flex gap-4">
                            {#each user.roles as role}
                                <div class="code text-cs">{role}</div>
                            {/each}
                        </div>
                    {/if}
                </div>
                <div class="px-10 pt-4 space-y-2">
                    {#if !update}
                        <p><span class="opacity-60">Last updated</span></p>
                        <p class="text-sm">{dateUpdated} @ {timeUpdated} (UTC)</p>
                    {/if}
                </div>
            </div>
        {/if}
        <div class="pl-8 pt-8 flex pl">
            {#if update}
                <button type="submit" class="btn btn-sm variant-outline-success m-2" on:click={updateUser}>Submit</button>
                <button type="submit" class="btn btn-sm variant-outline-warning m-2" on:click={() =>  { update = false; loading = false }}>Cancel</button>
            {:else}
                <button type="submit" class="btn btn-sm variant-outline-success m-2" on:click={openVerificationModal}>Verify</button>
                <button type="submit" class="btn btn-sm variant-outline-warning m-2" on:click={() => update = true }>Update</button>
                <button type="submit" class="btn btn-sm variant-outline-error m-2" on:click={deleteUser}>Delete</button>
            {/if}
        </div>
    </section>
	<footer class="card-footer">

        
    </footer>
</div>