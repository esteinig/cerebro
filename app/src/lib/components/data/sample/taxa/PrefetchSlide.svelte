<script lang="ts">
  import { createEventDispatcher, onMount } from 'svelte';
  import PrefetchTable from './PrefetchTable.svelte';
  import { CiqaSampleType, type TrainingRecord, type TestResult, type TrainingPrefetchData } from '$lib/utils/types';
  import { getUuidShort } from '$lib/utils/helpers';
	import { page } from '$app/stores';

  export let trainingRecord: TrainingRecord;
  export let trainingData: TrainingPrefetchData;
  export let trainingSlides: number = 0;

  $: trainingIndex = Number.parseInt($page.params.record, 10);
  if (Number.isNaN(trainingIndex)) trainingIndex = 0;

  let showPrimary = false;
  let showSecondary = false;
  let showTarget = false;

  const dispatch = createEventDispatcher<{
    next: { index: number; item: TrainingPrefetchData; label: TestResult; candidate: string | null };
    previous: { index: number; item: TrainingPrefetchData; label: TestResult; candidate: string | null };
    complete: {  };
  }>();

  $: selectedCandidateName = trainingRecord.candidates ? trainingRecord.candidates[0] ?? null : null; // first candidate from multiple (not implemented yet)
  $: selectedLabel = trainingRecord.result;
  
  function next() {
    dispatch('next', { index: trainingIndex, item: trainingData, label: selectedLabel, candidate: selectedCandidateName });
  }

  function prev() {
    dispatch('previous', { index: trainingIndex, item: trainingData, label: selectedLabel, candidate: selectedCandidateName });
  }

  function submit() {
    dispatch('complete', {})
  }

  function getSampleTypeName(sample_type: CiqaSampleType): string {
    if (sample_type == CiqaSampleType.Eye) return 'Vitreous fluid';
    if (sample_type == CiqaSampleType.Csf) return 'Cerebrospinal fluid';
    return sample_type as string;
  }

  function onKey(e: KeyboardEvent) {
    if (e.key === 'ArrowRight') next();
    if (e.key === 'ArrowLeft') prev();
  }

  onMount(() => {     
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  });

  function onSelect(e: CustomEvent<{ index: number; item: { name: string } | null }>) {
    const name = e.detail.item ? e.detail.item.name : null;
    selectedCandidateName = name;
    selectedLabel = name ? 'positive' : 'negative';
  }
</script>
  
    <div class="w-full flex items-center justify-center gap-3 py-3 sticky top-0 bg-base-100/80 backdrop-blur z-10">
      
      <button class="px-3 py-1 border rounded disabled:opacity-40"
              on:click={prev}
              disabled={trainingIndex === 0}>
        ◀ Prev
      </button>
  
      <div class="text-xl opacity-70">
        Sample {trainingIndex + 1} / {trainingSlides} • 
        <span class="ml-1 chip variant-filled-primary">{getSampleTypeName(trainingData.prefetch.config.sample_type)}</span>
        <span class="ml-1">• {getUuidShort(trainingData.id)}</span>
      </div>

  
      <span class="mx-2 h-5 w-px bg-gray-300" />
      <button
        class={`px-3 py-1 border rounded
                ${selectedCandidateName
                    ? 'bg-primary-500 text-white'
                    : ''}`}
        disabled>
        Positive
      </button>

      <button
        class={`px-3 py-1 border rounded
                ${!selectedCandidateName
                    ? 'bg-tertiary-500 text-white'
                    : ''}`}
        disabled>
        Negative
      </button>
    
      
      {#if selectedCandidateName}
        <p class="mx-5">{selectedCandidateName}</p>
      {:else}
        <p class="mx-5">No pathogen candidate selected</p>
      {/if}
      
      {#if trainingIndex === trainingSlides - 1}
        <button class="px-3 py-1 border rounded disabled:opacity-40"
              on:click={submit}>
        Submit<span class="ml-3">▶</span> 
        </button>
      {:else}
        <button class="px-3 py-1 border rounded disabled:opacity-40"
                on:click={next}
                disabled={trainingIndex === trainingSlides - 1}>
          Next<span class="ml-3">▶</span> 
        </button>
      {/if}

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
            <PrefetchTable taxa={trainingData.prefetch.primary}  selectedName={selectedCandidateName ?? ""}  on:select={onSelect} />
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
            <PrefetchTable taxa={trainingData.prefetch.secondary}   selectedName={selectedCandidateName ?? ""}  on:select={onSelect} />
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
            <PrefetchTable taxa={trainingData.prefetch.target}   selectedName={selectedCandidateName ?? ""}  on:select={onSelect} />
          {/if}
        </section>
      </div>  
  <style>
    /* Tailwind assumed; keep minimal if not using Tailwind */
  </style>
  