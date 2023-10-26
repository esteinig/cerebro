<script lang="ts">
	import { getInitials } from "$lib/utils/helpers";
import { Avatar } from "@skeletonlabs/skeleton";

    export let timestamp: string;
    export let user: string;
	export let activity: string;
	export let target: string = '';
	export let description: string = '';
    export let link: string = '';
	export let badgeText: string = '';
    export let targetPath: string[] = [];

    let userInitials = getInitials(user);

</script>

<li class="mb-10 ml-10 rounded-tr-lg rounded-br-lg rounded-bl-lg bg-surface-500/10">
	
    <!-- Slot: Icon -->
	<span class="absolute flex items-center justify-center w-8 h-8 variant-filled-primary rounded-full -left-3 ring-8 ring-surface-200-700-token">
        <slot name="icon">
            <Avatar initials={userInitials}></Avatar>
        </slot>
    </span>
	<div class="p-4">
        <h3 class="flex items-center mb-1 space-x-2 -translate-y-2 justify-between cursor-pointer">
            <div class="group">
                <div class="group-hover:hidden">
                    <div class="text-sm">
                        <span class="font-semibold">{user} </span> <span class="opacity-60">{activity} </span> <span class="font-semibold">{target}</span>    
                        {#if badgeText}
                            <span class="badge variant-soft-primary">{badgeText}</span>
                        {/if}
                    </div>
                </div>
                <div class="hidden group-hover:flex">
                    {#if targetPath.length}
                        <div class="text-xs opacity-75">
                            <ol class="breadcrumb justify-center">
                                {#each targetPath as path, i}
                                    <li class="crumb">{path}</li>
                                    {#if !(i === targetPath.length-1)}
                                        <li class="crumb-separator" aria-hidden>&rsaquo;</li>
                                    {/if}
                                {/each}
                            </ol>
                        </div>
                    {:else}
                        <div class="text-sm">
                            <span class="font-semibold">{user} </span> <span class="opacity-60">{activity} </span> <span class="font-semibold">{target}</span>    
                            {#if badgeText}
                                <span class="badge variant-soft-primary">{badgeText}</span>
                            {/if}
                        </div>
                    {/if}
                </div>
            </div>
            

            <time class="block mb-2 opacity-30 text-sm">{timestamp}</time>
        </h3>
        <div class="flex items-start justify-between cursor-pointer">

            <p class="opacity-75 text-sm max-w-[90%]">{@html description}</p>
            {#if link}
                <a class="btn-sm variant-ghost rounded-lg hover:" href="{link}" target="_blank" rel="noopener noreferrer">
                    <svg class="w-3" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                        <path d="M13.5 6H5.25A2.25 2.25 0 003 8.25v10.5A2.25 2.25 0 005.25 21h10.5A2.25 2.25 0 0018 18.75V10.5m-10.5 6L21 3m0 0h-5.25M21 3v5.25" stroke-linecap="round" stroke-linejoin="round"></path>
                    </svg>
                </a>
            {/if}
        </div>

        <!-- Slot: Default -->
        {#if $$slots.default}
            <div><slot /></div>
        {/if}

        
    </div>
</li>
