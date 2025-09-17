<script lang="ts">
	import { goto } from "$app/navigation";
	import { page } from "$app/stores";
	import { type TrainingPrefetchOverview } from "$lib/utils/types";

    export let data: TrainingPrefetchOverview[] = [];

</script>

<div>
    <div class="table-container pt-8 text-xl align-center">
        <table class="table table-hover table-compact">
            <thead>
                <tr>
                    <th>Training Dataset</th>
                    <th>Training Samples</th>
                    <th>Dataset Description</th>
                    <th>Operations</th>
                </tr>
            </thead>
            <tbody>
                {#if data.length === 0}
					<tr>
						<td colspan="4" class="text-center opacity-60">No datasets available</td>
					</tr>
				{:else}
					{#each data as overview}
						<tr>
							<td class="">{overview.collection}</td>
							<td>{overview.samples}</td>
							<td class="opacity-80">{overview.description}</td>
							<td class="">
                                <button type="button" class="btn btn-md variant-outline-primary" on:click={() => goto(`${$page.url}/collection=${overview.collection}&session=0&record=0`)}>

                                    <div class="w-5 h-5 mr-2">
                                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                            <path d="M3 16.5v2.25A2.25 2.25 0 005.25 21h13.5A2.25 2.25 0 0021 18.75V16.5M16.5 12L12 16.5m0 0L7.5 12m4.5 4.5V3" stroke-linecap="round" stroke-linejoin="round"></path>
                                        </svg>
                                    </div>            
                                    New Session     
                                </button>
							</td>
						</tr>
					{/each}
				{/if}
            </tbody>
        </table>
    </div>
</div>