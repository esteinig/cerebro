<script lang="ts">
  import { onMount } from "svelte";
  import { getToastStore, ProgressRadial } from "@skeletonlabs/skeleton";
  import type { TrainingPrefetchData } from "$lib/utils/types";
  import { page } from "$app/stores";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";

  export let collection: string = "Test";

  let prefetch: TrainingPrefetchData[] = [];
  let publicApi = new CerebroApi();
  let toastStore = getToastStore();
  let loading = false;

  async function getPrefetchDataCollection(collection: string) {
    loading = true;
    let response: ApiResponse = await publicApi.fetchWithRefresh(
      `${publicApi.routes.training.getCollection}?team=${$page.params.team}&collection=${collection}`,
      {
        method: "GET",
        mode: "cors",
        credentials: "include",
      } as RequestInit,
      $page.data.refreshToken,
      toastStore,
      "Training data loaded"
    );

    loading = false;

    if (response.ok) {
      prefetch = response.json.data as TrainingPrefetchData[];
    }
  }

  onMount(() => {
    getPrefetchDataCollection(collection);
  });
</script>


<div>


  {#if loading}
    <div class="flex justify-center py-24">
        <ProgressRadial width="sm:w-12 md:w-24" stroke={20} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
    </div>
  {:else}

  {/if}

</div>