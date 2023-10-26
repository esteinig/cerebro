<script lang="ts">
	import { page } from '$app/stores';
	import ErrorAnimation from '$lib/general/error/ErrorAnimation.svelte';
</script>

{#if $page}
<div class="mt-[10%] flex justify-center items-center">
    <div class="text-center space-y-4">
        <div class="flex justify-center">
            <ErrorAnimation />
        </div>
        <h3 class="h3">
            Something went wrong
        </h3>
        
        <p>{#if $page.error} {$page.error.message} {/if}</p>
    </div>
</div>
{/if}