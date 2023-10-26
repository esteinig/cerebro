<script lang="ts">
	import { page } from '$app/stores';
	import ErrorAnimation from '$lib/general/error/ErrorAnimation.svelte';
</script>

{#if $page}
<div class="mt-[10%] flex justify-center items-center">
    <div class="text-center space-y-8">
        <div class="flex justify-center">
            <ErrorAnimation />
        </div>
        <h4 class="h4">
            This page is under construction
        </h4>
    </div>
</div>
{/if}