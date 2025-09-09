<script lang="ts">
    import { createEventDispatcher, onMount } from 'svelte';
    import PrefetchTable from './PrefetchTable.svelte';
	import { type TrainingPrefetchData } from '$lib/utils/types';
	import { getUuidShort } from '$lib/utils/helpers';
	import { navigationLoading } from '$lib/stores/stores';
  
  
    export let items: TrainingPrefetchData[] = [];
    export let startIndex = 0;

    let showPrimary = false;
    let showSecondary = false;
    let showTarget = false;
  
    const dispatch = createEventDispatcher<{
      change: { index: number; item: TrainingPrefetchData };
      classify: { index: number; item: TrainingPrefetchData; label: 'positive' | 'negative' };
    }>();
  
    let index = Math.min(Math.max(startIndex, 0), Math.max(items.length - 1, 0));
  
    $: hasData = items && items.length > 0;
    $: current = hasData ? items[index] : null;
  
    function next() {
      if (!hasData) return;
      selectedCandidateName = null;
      showPrimary = true;
      showSecondary = true;
      showTarget = true;
      hasSelection = false;
      index = Math.min(index + 1, items.length - 1);
      dispatch('change', { index, item: items[index] });
    }
  
    function prev() {
      if (!hasData) return;
      selectedCandidateName = null;
      showPrimary = true;
      showSecondary = true;
      showTarget = true;
      hasSelection = false;
      index = Math.max(index - 1, 0);
      dispatch('change', { index, item: items[index] });
    }
  
    function classify(label: 'positive' | 'negative') {
      if (!hasData) return;
      dispatch('classify', { index, item: items[index], label });
    }
  
    // Optional: arrow-key navigation
    function onKey(e: KeyboardEvent) {
      if (e.key === 'ArrowRight') next();
      if (e.key === 'ArrowLeft') prev();
    }
    
    // run only in browser
    onMount(() => {
        window.addEventListener('keydown', onKey);
        return () => window.removeEventListener('keydown', onKey);
    });

    let selectedCandidateName: string | null = null;

    function onSelect(e: CustomEvent<{ index: number; item: { name: string } | null }>) {
        selectedCandidateName = e.detail.item ? e.detail.item.name : null;
        classify(selectedCandidateName ? 'positive' : 'negative');
    }

    $: hasSelection = !!selectedCandidateName;

  </script>
  
  {#if hasData}
    <div class="w-full flex items-center justify-center gap-3 py-3 sticky top-0 bg-base-100/80 backdrop-blur z-10">
      
      <!-- <button class="px-3 py-1 border rounded disabled:opacity-40"
              on:click={prev}
              disabled={index === 0}>
        ◀ Prev
      </button> -->
  
      <div class="text-sm opacity-70">
        Sample {index + 1} / {items.length}
        {#if current?.id}<span class="ml-2">• ID: {getUuidShort(current.id)}</span>{/if}
      </div>
    
      <span class="mx-2 h-5 w-px bg-gray-300" />


      <div class="text-sm opacity-70">
        Sample type: {current?.prefetch.config.sample_type}
      </div>
  
      <span class="mx-2 h-5 w-px bg-gray-300" />
      <button
        class={`px-3 py-1 border rounded
                ${hasSelection
                    ? 'bg-primary-500 text-white hover:bg-primary-600'
                    : 'hover:bg-primary-50'}`}
        on:click={() => classify('positive')}
        aria-label="Mark sample as positive">
        Positive
      </button>

      <button
        class={`px-3 py-1 border rounded
                ${!hasSelection
                    ? 'bg-tertiary-500 text-white hover:bg-tertiary-600'
                    : 'hover:bg-tertiary-50'}`}
        on:click={() => classify('negative')}>
        Negative
      </button>
    
      
      {#if selectedCandidateName}
        <p class="mx-5">{selectedCandidateName}</p>
      {:else}
        <p class="mx-5">No pathogen candidate selected</p>
      {/if}
      
      <button class="px-3 py-1 border rounded disabled:opacity-40"
              on:click={next}
              disabled={index === items.length - 1}>
        Next<span class="ml-3">▶</span> 
      </button>

    </div>
  
    <div class="space-y-6">
        <section>
          <p class="flex items-center justify-between text-lg opacity-60 py-4" >
            <span on:click={() => (showPrimary = !showPrimary)} class="cursor-pointer">Above-threshold taxonomic profile</span>
            <button class="text-xs px-2 py-1 border rounded"
              on:click={() => (showPrimary = !showPrimary)}>
              {showPrimary ? 'Hide' : 'Show'}
            </button>
          </p>
          {#if showPrimary}
            <PrefetchTable taxa={current?.prefetch.primary} on:select={onSelect} />
          {/if}
        </section>
      
        <section>
          <p class="flex items-center justify-between text-lg opacity-60 py-4">
            <span on:click={() => {  (showSecondary = !showSecondary)}} class="cursor-pointer">Below-threshold taxonomic profile</span>
           
            <button class="text-xs px-2 py-1 border rounded"
              on:click={() => (showSecondary = !showSecondary)}>
              {showSecondary ? 'Hide' : 'Show'}
            </button>
          </p>
          {#if showSecondary}
            <PrefetchTable taxa={current?.prefetch.secondary} on:select={onSelect} />
          {/if}
        </section>
      
        <section>
          <p class="flex items-center justify-between text-lg opacity-60 py-4">
            <span on:click={() => (showTarget = !showTarget)} class="cursor-pointer">Target-threshold taxonomic profile</span>
            <button class="text-xs px-2 py-1 border rounded"
              on:click={() => (showTarget = !showTarget)}>
              {showTarget ? 'Hide' : 'Show'}
            </button>
          </p>
          {#if showTarget}
            <PrefetchTable taxa={current?.prefetch.target} on:select={onSelect} />
          {/if}
        </section>
      </div>
  {:else}
    <div class="py-12 opacity-60 text-center">No samples.</div>
  {/if}
  
  <style>
    /* Tailwind assumed; keep minimal if not using Tailwind */
  </style>
  