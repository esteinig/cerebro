<script lang="ts">
	import type { Cerebro, CommentBubble, SampleCommentSchema } from "$lib/utils/types";
	import { onMount } from "svelte";
	import { CerebroApi, type ApiResponse } from "$lib/utils/api";
	import { Avatar, ProgressRadial, getToastStore } from "@skeletonlabs/skeleton";
	import { page } from "$app/stores";
	import { invalidate } from "$app/navigation";
	import ErrorAnimation from "$lib/general/error/ErrorAnimation.svelte";
	import { getDateTimeStringUtc } from "$lib/utils/helpers";
	import { storeSettings } from "$lib/stores/stores";

	export let selectedDatabaseId: string;
	export let selectedProjectId: string;
	
    const toastStore = getToastStore();
    const publicApi = new CerebroApi();

	let commentFeed: CommentBubble[] = [];
	let loading: boolean = false;

	const refreshComments = async() => {

		loading = true;

		await invalidate("sample:data");

		loading = false;

	}

	const submitComment = async() => {

		loading = true;

		let commentData = {
			user_id: $page.data.userData.id,
			user_name: $page.data.userData.name,
			user_title: $page.data.userData.title,
			user_positions: $page.data.userData.positions,
			comment_text: comment
		} satisfies SampleCommentSchema;

		let updateSamples = $page.data.sampleCerebro.map((model: Cerebro) => model.id).join(",");

		let response: ApiResponse = await publicApi.fetchWithRefresh(
			`${publicApi.routes.cerebro.addSampleComment}?db=${selectedDatabaseId}&project=${selectedProjectId}&id=${updateSamples}`, {
				method: 'POST',
				mode: 'cors',
				credentials: 'include',
				body: JSON.stringify(commentData),
				headers:  { 'Content-Type': 'application/json' }
			} satisfies RequestInit,
			$page.data.refreshToken, toastStore, `Comment submitted`
		);

		if (response.ok) {
			// Updating page data with reload
			await invalidate("sample:data")
		}

		loading = false;
	}

	const deleteComment = async(commentIdentifier: string) => {

		loading = true;

		let updateSamples = $page.data.sampleCerebro.map((model: Cerebro) => model.id).join(",");

		let response: ApiResponse = await publicApi.fetchWithRefresh(
			`${publicApi.routes.cerebro.deleteSampleComment}?db=${selectedDatabaseId}&project=${selectedProjectId}&id=${updateSamples}&comment_id=${commentIdentifier}`, {
				method: 'DELETE',
				mode: 'cors',
				credentials: 'include',
			} satisfies RequestInit,
			$page.data.refreshToken, toastStore, `Comment deleted`
		);

		if (response.ok) {
			// Updating page data with reload
			await invalidate("sample:data")
		}

		loading = false;
	}

	let comment = '';
	let elemChat: HTMLElement;

	const scrollChatBottom = (behavior?: ScrollBehavior): void => {
		elemChat.scrollTo({ top: elemChat.scrollHeight, behavior });
	}

	onMount(() => {
		if (commentFeed.length && elemChat){
			scrollChatBottom();
		}
	})

	$: {
		let comments: CommentBubble[] = $page.data.sampleCerebro.map((model: Cerebro) => {
			return model.sample.comments.map(comment => {
				return {
					id: comment.comment_id,
					host: comment.user_id === $page.data.userData.id,
					name: comment.user_name,
					timestamp: getDateTimeStringUtc(comment.comment_date),
					message: comment.comment_text,
					color: 'variant-soft-primary',
					title: comment.user_title,
					position: comment.user_positions.join(", ")
				} satisfies CommentBubble
			})
		}).flat()
		
		comments.sort((a, b) => a.timestamp.localeCompare(b.timestamp));
		commentFeed = comments.filter((value, index, self) => {
			return self.findIndex(v => v.id === value.id) === index;
		})

	}

</script>

<div>
    <div class="flex items-center justify-center {$storeSettings.acceptedCommentRisk ? "blur-none" : "blur-xl"}">
        <button type="button" class="btn-icon-sm"  on:click={refreshComments}>
            <div class="w-5 h-5">
                <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                    <path d="M16.023 9.348h4.992v-.001M2.985 19.644v-4.992m0 0h4.992m-4.993 0l3.181 3.183a8.25 8.25 0 0013.803-3.7M4.031 9.865a8.25 8.25 0 0113.803-3.7l3.181 3.182m0-4.991v4.99" stroke-linecap="round" stroke-linejoin="round"></path>
                </svg>
            </div>
        </button>
    </div>
   
    <div class="w-full grid grid-cols-1 gap-1">
        <div class="h-full grid grid-rows-[1fr_auto] gap-1">
			{#if loading}
				<div class="flex justify-center py-24">
					<ProgressRadial width="sm:w-12 md:w-24" stroke={100} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
				</div>
			{:else}
				{#if !$storeSettings.acceptedCommentRisk}
					<div class="flex justify-center">
						<aside class="alert variant-ghost-warning w-1/2">
							<svg class="w-48" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
								<path d="M12 9v3.75m-9.303 3.376c-.866 1.5.217 3.374 1.948 3.374h14.71c1.73 0 2.813-1.874 1.948-3.374L13.949 3.378c-.866-1.5-3.032-1.5-3.898 0L2.697 16.126zM12 15.75h.007v.008H12v-.008z" stroke-linecap="round" stroke-linejoin="round"></path>
							</svg>
							<div class="alert-message">
								<h3 class="h3 font-bold">Warning</h3>
								<p>Submitting comments carries a risk of exposing protected patient information in your 
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
													<small class="opacity-50">{bubble.timestamp}</small>
												</header>
												<div class="flex justify-between items-center">
													<div>{bubble.message}</div>
													<button class="btn-icon-sm w-4 h-4" on:click={() => deleteComment(bubble.id)}>
														<svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
															<path d="M14.74 9l-.346 9m-4.788 0L9.26 9m9.968-3.21c.342.052.682.107 1.022.166m-1.022-.165L18.16 19.673a2.25 2.25 0 01-2.244 2.077H8.084a2.25 2.25 0 01-2.244-2.077L4.772 5.79m14.456 0a48.108 48.108 0 00-3.478-.397m-12 .562c.34-.059.68-.114 1.022-.165m0 0a48.11 48.11 0 013.478-.397m7.5 0v-.916c0-1.18-.91-2.164-2.09-2.201a51.964 51.964 0 00-3.32 0c-1.18.037-2.09 1.022-2.09 2.201v.916m7.5 0a48.667 48.667 0 00-7.5 0" stroke-linecap="round" stroke-linejoin="round"></path>
														</svg>
													</button>
												</div>
											</div>
										</div>
									{:else}
										<!-- Team member comment -->
										<div class="grid grid-cols-[1fr_auto] gap-2">
											<div class="card p-4 rounded-tr-none space-y-2 variant-soft">
												<header class="flex justify-between items-center">
													<p class="font-bold">{bubble.title ?? ""}{bubble.name}<span class="opacity-60 ml-2">{bubble.position}</span></p>
													<small class="opacity-50">{bubble.timestamp}</small>
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
							<p class="flex justify-center text-sm pb-8">Submit the first comment for this sample</p>
							
						</div>
					{/if}
					<div class="bg-surface-500/0 p-4 {$storeSettings.acceptedCommentRisk ? "blur-none" : "blur-xl"}">
						<div class="input-group input-group-divider grid-cols-[auto_1fr_auto] rounded-container-token">
							<button class="input-group-shim">+</button>
							<textarea
								bind:value={comment}
								class="bg-transparent border-0 ring-0"
								name="prompt"
								id="prompt"
								placeholder="Write a comment..."
								rows="1"
								required
							/>
							<button class="bg-primary-500/30" on:click={submitComment}>Send</button>
						</div>
					</div>
				{/if}
			{/if}
				
        </div>
    </div>
</div>