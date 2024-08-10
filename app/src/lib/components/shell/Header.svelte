<script lang="ts">
	import { enhance } from '$app/forms';
	import { goto, invalidateAll } from '$app/navigation';
	import { navigating, page } from '$app/stores';
    import { env } from '$env/dynamic/public';
	import { storeTheme } from '$lib/stores/stores';

    import { AppBar, popup, ProgressBar } from '@skeletonlabs/skeleton';
    import { LightSwitch } from '@skeletonlabs/skeleton';
    import { getToastStore } from '@skeletonlabs/skeleton';
    import type { ToastSettings } from '@skeletonlabs/skeleton';
	import type { SubmitFunction } from '@sveltejs/kit';

    const toastStore = getToastStore();

    const logout = async() => {
                
        const requestInitOptions: RequestInit = {
            method: 'GET',
            mode: 'cors',
            credentials: 'include'
        };

        try {
            let requestResponse: Response = await fetch(env.PUBLIC_CEREBRO_API_URL+"/auth/logout", requestInitOptions);
            if (requestResponse.ok) {
                goto('/login');
            } else {
                failedLogout(`Logout failed (${requestResponse.status} - ${requestResponse.statusText})`);
            }
        } catch (requestError) {
            failedLogout("Logout failed, server is unavailable");
        }       
    };

    const failedLogout = (message: string) => {
        const errorToast: ToastSettings = {
            message: message,
            background: 'variant-filled-warning',
        };
        toastStore.trigger(errorToast);
    }

	const themes = [
		{ type: 'hamlindigo', name: 'Hamlindigo', icon: '👔' },
		{ type: 'wintry', name: 'Wintry', icon: '🌨️', badge: '' },
		{ type: 'rocket', name: 'Rocket', icon: '🚀' },
		{ type: 'vintage', name: 'Vintage', icon: '📺' },
		{ type: 'crimson', name: 'Crimson', icon: '⭕' },
		{ type: 'gold-nouveau', name: 'Gold Nouveau', icon: '💫' },
		{ type: 'modern', name: 'Modern', icon: '🤖' },
        { type: 'skeleton', name: 'Skeleton', icon: '💀' },
		{ type: 'seafoam', name: 'Seafoam', icon: '🧜‍♀️' },
		{ type: 'sahara', name: 'Sahara', icon: '🏜️' },
		{ type: 'dali', name: 'Dali', icon: '🎨' },
	];

	const setTheme: SubmitFunction = ({ formData }) => {
		const theme = formData.get('theme')?.toString();

		if (theme) {
			$storeTheme = theme;
		}
        invalidateAll()
	};

</script>

<AppBar>
	<svelte:fragment slot="lead">
        <h2 class="h2"><a href="/cerebro">Cerebro</a></h2>
    </svelte:fragment>

	<svelte:fragment slot="trail">

    <!-- Theme -->
    <div>
        <!-- trigger -->
        <button class="btn hover:variant-soft-primary" use:popup={{ event: 'click', target: 'theme', closeQuery: 'a[href]' }}>
            <i class="fa-solid fa-palette text-lg md:!hidden" />
            <span class="hidden md:inline-block">Theme</span>
            <i class="fa-solid fa-caret-down opacity-50" />
        </button>
        <!-- popup -->
        <div class="card p-4 w-60 shadow-xl" data-popup="theme">
            <div class="space-y-4">
                <section class="flex justify-between items-center">
                    <h6 class="h6">Mode</h6>
                    <LightSwitch />
                </section>
                <hr />
                <nav class="list-nav p-4 -m-4 max-h-64 lg:max-h-[500px] overflow-y-auto">
                    <form action="/?/setTheme" method="POST" use:enhance={setTheme}>
                        <ul>
                            {#each themes as { icon, name, type, badge }}
                                <li>
                                    <button
                                        class="option w-full h-full"
                                        type="submit"
                                        name="theme"
                                        value={type}
                                        class:bg-primary-active-token={$storeTheme === type}
                                    >
                                        <span>{icon}</span>
                                        <span class="flex-auto text-left">{name}</span>
                                        {#if badge}<span class="badge variant-filled-secondary">{badge}</span>{/if}
                                    </button>
                                </li>
                            {/each}
                        </ul>
                    </form>
                </nav>
            </div>
            <!-- <div class="arrow bg-surface-100-800-token" /> -->
        </div>
    </div>

    <!-- Data Interface -->
    <a href="/cerebro/data/samples" class="btn-icon btn-icon-sm variant-filled">
        <div class="h-5 w-5">
            <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                <path d="M3 13.125C3 12.504 3.504 12 4.125 12h2.25c.621 0 1.125.504 1.125 1.125v6.75C7.5 20.496 6.996 21 6.375 21h-2.25A1.125 1.125 0 013 19.875v-6.75zM9.75 8.625c0-.621.504-1.125 1.125-1.125h2.25c.621 0 1.125.504 1.125 1.125v11.25c0 .621-.504 1.125-1.125 1.125h-2.25a1.125 1.125 0 01-1.125-1.125V8.625zM16.5 4.125c0-.621.504-1.125 1.125-1.125h2.25C20.496 3 21 3.504 21 4.125v15.75c0 .621-.504 1.125-1.125 1.125h-2.25a1.125 1.125 0 01-1.125-1.125V4.125z" stroke-linecap="round" stroke-linejoin="round"></path>
            </svg>
        </div>
    </a>
    
    <!-- Admin Section -->
    {#if $page.data.admin }
        <a href="/cerebro/admin/logs" class="btn-icon btn-icon-sm variant-filled">
            <div class="h-5 w-5">
                <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                    <path d="M4.5 12a7.5 7.5 0 0015 0m-15 0a7.5 7.5 0 1115 0m-15 0H3m16.5 0H21m-1.5 0H12m-8.457 3.077l1.41-.513m14.095-5.13l1.41-.513M5.106 17.785l1.15-.964m11.49-9.642l1.149-.964M7.501 19.795l.75-1.3m7.5-12.99l.75-1.3m-6.063 16.658l.26-1.477m2.605-14.772l.26-1.477m0 17.726l-.26-1.477M10.698 4.614l-.26-1.477M16.5 19.794l-.75-1.299M7.5 4.205L12 12m6.894 5.785l-1.149-.964M6.256 7.178l-1.15-.964m15.352 8.864l-1.41-.513M4.954 9.435l-1.41-.514M12.002 12l-3.75 6.495" stroke-linecap="round" stroke-linejoin="round"></path>
                </svg>
            </div>
        </a>
    {/if}


    <!-- User Profile -->
    <a href="/cerebro/user" class="btn-icon btn-icon-sm variant-filled">
        <div class="h-5 w-5">
            <svg aria-hidden="true" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg">
                <path d="M10 8a3 3 0 100-6 3 3 0 000 6zM3.465 14.493a1.23 1.23 0 00.41 1.412A9.957 9.957 0 0010 18c2.31 0 4.438-.784 6.131-2.1.43-.333.604-.903.408-1.41a7.002 7.002 0 00-13.074.003z"></path>
            </svg>
        </div>
    </a>

    <!-- Logout -->
    <button type="button" class="btn-icon btn-icon-sm variant-filled" on:click={logout}>
        <div class="h-5 w-5">
            <svg aria-hidden="true" fill="currentColor" viewBox="0 0 20 20" xmlns="http://www.w3.org/2000/svg">
                <path clip-rule="evenodd" d="M3 4.25A2.25 2.25 0 015.25 2h5.5A2.25 2.25 0 0113 4.25v2a.75.75 0 01-1.5 0v-2a.75.75 0 00-.75-.75h-5.5a.75.75 0 00-.75.75v11.5c0 .414.336.75.75.75h5.5a.75.75 0 00.75-.75v-2a.75.75 0 011.5 0v2A2.25 2.25 0 0110.75 18h-5.5A2.25 2.25 0 013 15.75V4.25z" fill-rule="evenodd"></path>
                <path clip-rule="evenodd" d="M6 10a.75.75 0 01.75-.75h9.546l-1.048-.943a.75.75 0 111.004-1.114l2.5 2.25a.75.75 0 010 1.114l-2.5 2.25a.75.75 0 11-1.004-1.114l1.048-.943H6.75A.75.75 0 016 10z" fill-rule="evenodd"></path>
            </svg>
        </div>
    </button>
        
    </svelte:fragment>
</AppBar>