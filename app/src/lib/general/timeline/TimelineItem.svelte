<script lang="ts">
	export let title: string;
	export let timestamp: string;
	export let author: string;
	export let paragraphs: string[];
	export let badgeText: string = '';
    export let link: string = '';
    export let linkText: string = 'Link';
</script>

<li class="mb-10 ml-10">
	
    <!-- Slot: Icon -->
	<span class="absolute flex items-center justify-center w-6 h-6 variant-filled-primary rounded-full -left-3 ring-8 ring-surface-200-700-token">
        <slot name="icon">&hearts;</slot>
    </span>
		
	<h3 class="flex items-center mb-1 space-x-2 -translate-y-2 text-sm">
        <span>{title}</span>
        {#if badgeText}
            <span class="badge variant-soft-primary">{badgeText}</span>
        {/if}
    </h3>
    <time class="block mb-2 opacity-30 text-sm">{timestamp}<span class="ml-2 opacity-60">by {author}</span></time>
    {#each paragraphs as description}
	    <p class="mb-4 opacity-75 text-sm">{description}</p>
    {/each}

    {#if link}
        <div>
            <a class="mt-4 btn text-xs variant-ghost" href="{link}" target="_blank" rel="noopener noreferrer">
                <svg class="w-4 mr-2" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                    <path d="M13.5 6H5.25A2.25 2.25 0 003 8.25v10.5A2.25 2.25 0 005.25 21h10.5A2.25 2.25 0 0018 18.75V10.5m-10.5 6L21 3m0 0h-5.25M21 3v5.25" stroke-linecap="round" stroke-linejoin="round"></path>
                </svg>
                {linkText}
            </a>
        </div>
    {/if}
	<!-- Slot: Default -->
	{#if $$slots.default}
        <div><slot /></div>
    {/if}
</li>
