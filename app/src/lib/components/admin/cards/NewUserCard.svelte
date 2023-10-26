
<script lang="ts">
	import { invalidate } from "$app/navigation";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import type { RegisterUserSchema } from "$lib/utils/types";
    import { getUserRoles } from "$lib/utils/helpers";
    import { ProgressRadial, getToastStore } from '@skeletonlabs/skeleton';
	import { page } from "$app/stores";

    const toastStore = getToastStore();

    let loading: boolean = false;

    let publicApi = new CerebroApi();

    let setUserPassword: boolean = false;

    let name: string = "";
    let email: string = "";
    let title: string = "";
    let position: string = "";
    let verified: boolean = false;
    let password: string = "";
    let passwordConfirmation: string = "";

    let adminRole: boolean = false;
    let dataRole: boolean = false;
    let reportRole: boolean = false;
    let botRole: boolean = false;
  

    const registerUser = async() => {

        let newUser: RegisterUserSchema = {
            name: name,
            email: email,
            title: title,
            password: password,
            verified: verified,
            roles: getUserRoles(adminRole, reportRole, dataRole, botRole),
            positions: [position],
        };

        loading = true;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.users.register}`, {
                method: 'POST',
                mode: 'cors',
                credentials: 'include',
                body: JSON.stringify(newUser),
                headers:  { 'Content-Type': 'application/json' }
            } satisfies RequestInit,
            $page.data.refreshToken, toastStore, "User added"
        );

        if (response.ok){
            await invalidate("admin:data"); // Reload page data
            await invalidate("cerebro:data"); // Reload user base data
        }
        loading = false;

    }

    let passwordConfirmationClass: string = "";

    $: {
        if (password) {
            passwordConfirmationClass = password === passwordConfirmation ? "input-success" : "input-error";
        }
    }



</script>

<div class="card max-w-xl2">
	<header class="card-header">
        <div class="p-4 space-y-4">
            <h4 class="h4">Register a new user</h4>
            <p class="text-sm opacity-60">
                Create a new user with credentials and permissions. Users must be <span class="font-semibold">verified</span> to authenticate with the application
                programming interface and login to the application. 
                Verification emails can be sent with the SMTP server configured for this stack. Emails will let the user confirm their 
                email address and set their password. You can read more about this process in the <a href="" class="text-blue-500 dark:text-blue-300">documentation</a>.
             </p>
            <div class="flex">
                <div class="text-xs opacity-30 w-[70%] mr-3">
                    Manual password and user verification breaks the secure registration flow as the password
                    has to be communicated to the user and is not private until it is changed. Please consider the deployment context of your stack before using
                    manual settings.
                </div>
                <div class="ml-auto mb-1 space-y-2">
                    <label class="flex items-center space-x-2">
                        <input class="checkbox" type="checkbox" bind:checked={setUserPassword} />
                        <p class="text-xs">Manual password</p>
                    </label>
                    <label class="flex items-center space-x-2">
                        <input class="checkbox" type="checkbox" bind:checked={verified} />
                        <p class="text-xs">User is verified</p>
                    </label>
                  </div>
            </div>
            
        </div>
    </header>
	<section class="p-4">
        {#if loading}
        <div class="flex justify-center">
            <ProgressRadial width="w-24 p-4 space-y-4" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
        </div>
        {:else}
            <form class="p-4 space-y-4" on:submit|preventDefault={registerUser}>
                           
                <label class="label">
                    <span>Email</span>
                    <input type="text" placeholder="Email" class="input" bind:value={email} required/>
                </label>

                <label class="label">
                    <span>Name</span>
                    <input type="text" placeholder="Full name" class="input" bind:value={name} required/>
                </label>
                <label class="label">
                    <span>Title</span>
                    <input type="text" placeholder="Title (optional)" class="input" bind:value={title}/>
                </label>
                <label class="label">
                    <span>Position</span>
                    <input type="text" placeholder="Position" class="input" bind:value={position} required/>
                </label>
                                   
                {#if setUserPassword}
                    <label class="label">
                        <span>Password</span>
                        <input type="password" placeholder="Enter a  password..." class="input {passwordConfirmationClass}" bind:value={password}  required/>
                    </label>
                    <label class="label">
                        <span>Password confirmation</span>
                        <input type="password" placeholder="Confirm password ..." class="input {passwordConfirmationClass}" bind:value={passwordConfirmation} required/>
                    </label>
                {/if}
                
                <label class="label">
                    <span>Permissions</span>
                    <div class="input-group input-group-divider grid-cols-4 p-3 bg-transparent dark:bg-transparent">
                        <label class="flex items-center space-x-2">
                            <input class="checkbox" type="checkbox" bind:checked={adminRole} />
                            <p>Admin</p>
                        </label>
                        <label class="flex items-center space-x-2">
                            <input class="checkbox" type="checkbox"  bind:checked={dataRole} />
                            <p>Data</p>
                        </label>
                        <label class="flex items-center space-x-2">
                            <input class="checkbox" type="checkbox"  bind:checked={reportRole} />
                            <p>Report</p>
                        </label>
                        <label class="flex items-center space-x-2">
                            <input class="checkbox" type="checkbox"  bind:checked={botRole} />
                            <p>Bot</p>
                        </label>
                    </div>
                </label>
                <div class="flex justify-center">
                    <button type="submit" class="btn btn-lg variant-outline-tertiary">Submit</button>
                </div>
            </form>
        {/if}    

    </section>
	<footer class="card-footer">
        
    </footer>
</div>