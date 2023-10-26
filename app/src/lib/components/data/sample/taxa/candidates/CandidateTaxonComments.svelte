<script lang="ts">
	import type { CommentBubble, PriorityTaxon } from "$lib/utils/types";
	import { onMount } from "svelte";
	import { Avatar, ProgressRadial } from "@skeletonlabs/skeleton";
	import { page } from "$app/stores";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import { getDateTimeStringUtc } from "$lib/utils/helpers";
	import { storeSettings } from "$lib/stores/stores";

    export let priorityTaxon: PriorityTaxon;
    
    let commentFeed: CommentBubble[] = [];

	$: {
        commentFeed = priorityTaxon.decisions.map(decision => {
            return {
                id: decision.id,
                host: decision.user_id === $page.data.userData.id,
                name: decision.user_name,
                timestamp: getDateTimeStringUtc(decision.date),
                message: decision.comment.length ? decision.comment : "No comment submitted",
                color: 'variant-soft-primary',
                title: '',
                position: `${decision.decision}ed`
            } satisfies CommentBubble
        });
    }
    
	let loading: boolean = false;

	let elemChat: HTMLElement;

	const scrollChatBottom = (behavior?: ScrollBehavior): void => {
		elemChat.scrollTo({ top: elemChat.scrollHeight, behavior });
	}

	onMount(() => {
		if ($storeSettings.acceptedCommentRisk && commentFeed.length){
			scrollChatBottom();
		}
	})

	$: {		
		commentFeed.sort((a, b) => a.timestamp.localeCompare(b.timestamp));
		commentFeed = commentFeed.filter((value, index, self) => {
			return self.findIndex(v => v.id === value.id) === index;
		})
	}

</script>

<div>
    <div class="w-full grid grid-cols-1 gap-1">
        <div class="h-full grid grid-rows-[1fr_auto] gap-1">
			{#if loading}
				<div class="flex justify-center py-24">
					<ProgressRadial width="sm:w-12 md:w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
				</div>
			{:else}
                {#if !$storeSettings.acceptedCommentRisk}
                    <div class="flex justify-center">
                        <aside class="alert variant-ghost-warning w-3/4">
                            <svg class="w-48" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                <path d="M12 9v3.75m-9.303 3.376c-.866 1.5.217 3.374 1.948 3.374h14.71c1.73 0 2.813-1.874 1.948-3.374L13.949 3.378c-.866-1.5-3.032-1.5-3.898 0L2.697 16.126zM12 15.75h.007v.008H12v-.008z" stroke-linecap="round" stroke-linejoin="round"></path>
                            </svg>
                            <div class="alert-message">
                                <h3 class="h3 font-bold">Warning</h3>
                                <p>
                                    Submitting comments carries a risk of exposing protected patient information in your 
                                    discussion with other team members. All comments are logged on the server.
                                </p>
                                <p> Please consider compliance with local regulatory and ethics frameworks before submitting comments.</p>				
                                <button type="button" class="btn-icon-lg mr-5"  on:click={() => $storeSettings.acceptedCommentRisk = true}>
                                    Accept
                                </button>
                                <button type="button" class="btn-icon-lg"  on:click={() => $storeSettings.acceptedCommentRisk = false}>
                                    Reject
                                </button>
                            </div>
                        </aside>
                    </div>
                {:else}
                    {#if commentFeed.length}
                        <div class="bg-surface-500/0 p-4 overflow-y-auto {$storeSettings.acceptedCommentRisk ? "blur-none" : "blur-xl"}">
                            <section bind:this={elemChat} class="w-full max-h-[600px] p-4 overflow-y-auto space-y-4">
                                {#each commentFeed as bubble, i}
                                    {#if bubble.host === true}
                                        <!-- User comment -->
                                        <div class="grid grid-cols-[auto_1fr] gap-2">
                                            <Avatar initials="TM" width="w-12" />
                                            <div class="card p-4 variant-soft-primary rounded-tl-none space-y-2">
                                                <header class="flex justify-between items-center">
                                                    <p class="font-bold">{bubble.title ?? ""}{bubble.name}<span class="opacity-60 ml-2">{bubble.position}</span></p>
                                                    <small class="opacity-50 text-xs">{getDateTimeStringUtc(bubble.timestamp, false)}</small>
                                                </header>
                                                <div class="flex justify-between items-center">
                                                    <div>{bubble.message}</div>
                                                </div>
                                            </div>
                                        </div>
                                    {:else}
                                        <!-- Team member comment -->
                                        <div class="grid grid-cols-[1fr_auto] gap-2">
                                            <div class="card p-4 rounded-tr-none space-y-2 variant-soft">
                                                <header class="flex justify-between items-center">
                                                    <p class="font-bold">{bubble.title ?? ""}{bubble.name}<span class="opacity-60 ml-2">{bubble.position}</span></p>
                                                    <small class="opacity-50">{getDateTimeStringUtc(bubble.timestamp, false)}</small>
                                                </header>
                                                <p>{bubble.message}</p>
                                            </div>
                                            <Avatar initials="TM" width="w-12" />
                                        </div>
                                    {/if}
                                {/each}
                            </section>
                        </div>
                    {:else}
                        <div class="{$storeSettings.acceptedCommentRisk ? "blur-none" : "blur-xl"}">

                            <div class="flex justify-center py-16 ">
                                <ErrorAnimation></ErrorAnimation>
                            </div>
                            <p class="flex justify-center text-sm pb-8">No decision comments for this candidate</p>
                        </div>
                    {/if}

                {/if}
				
			{/if}
        </div>
    </div>
</div>