<script lang="ts">
   import { goto } from '$app/navigation';
   import { env } from '$env/dynamic/public';
   import { ProgressRadial } from '@skeletonlabs/skeleton';
   import type { AuthLoginSchema } from '$lib/utils/types';
    
    let userLogin: AuthLoginSchema = { 
        email: "", 
        password: ""
    };

    let loading: boolean = false;
    let failure: boolean = false;
    
    const login = async() => {
        
        loading = true;
        
        const requestInitOptions: RequestInit = {
            headers:  { 'Content-Type': 'application/json' },
            method: 'POST',
            mode: 'cors',
            credentials: 'include',
            body: JSON.stringify(userLogin)
        };

        try {
            let requestResponse: Response = await fetch(env.PUBLIC_CEREBRO_API_URL+"/auth/login", requestInitOptions);
            if (requestResponse.ok) {
                setTimeout(() => goto('/cerebro'), 800);
            } else {
                failedLogin();
            }
        } catch (requestError) {
            failedLogin();
        }       
    };

    const failedLogin = () => {
        // Make user experience more expected and 
        // show loading radial on failed login,
        // then timeout with generic error
        failure = true;
        setTimeout(() => loading = false, 800);
        setTimeout(() => failure = false, 3000);
    };

</script>
    <div class="card p-6 space-y-6 shadow-xl text-left"  data-sveltekit-preload-data="off">
        <form class="space-y-4">
            
            <label class="label">
                <span>Email</span>
                <input type="email" placeholder="" class="input" bind:value={userLogin.email} />
            </label>
            <label class="label">
                <span>Password</span>
                <input type="password" placeholder="" class="input" bind:value={userLogin.password} />
            </label>

            {#if loading}
                <div class="flex justify-center">
                    <ProgressRadial width="w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
                </div>
            {:else if failure}
                <div class="space-y-2">
                    <div class="flex justify-center">
                        <svg aria-hidden="true" fill="none" stroke="currentColor" class="h-16 w-16 text-secondary-500" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M12 9v3.75m9-.75a9 9 0 11-18 0 9 9 0 0118 0zm-9 3.75h.008v.008H12v-.008z" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg> 
                    </div>
                    <p class="text-center">Unauthorized</p>  
                </div> 
            {:else}
                <button class="btn variant-filled-primary w-full" on:click={login}>            
                    Login
                </button>
            {/if}
        </form>
    </div>