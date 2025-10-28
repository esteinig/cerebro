<script lang="ts">
  import { createEventDispatcher, onMount } from 'svelte';
  import PrefetchTable from './PrefetchTable.svelte';
  import { CiqaSampleType, type TrainingRecord, type TestResult, type TrainingPrefetchData } from '$lib/utils/types';
  import { getUuidShort } from '$lib/utils/helpers';
	import { page } from '$app/stores';
	import { getModalStore, type ModalSettings, SlideToggle } from '@skeletonlabs/skeleton';
	import InfoModal from '$lib/components/training/InfoModal.svelte';

  const modalStore = getModalStore();

  export let trainingRecord: TrainingRecord;
  export let trainingData: TrainingPrefetchData;
  export let trainingSlides: number = 0;

  // optional weight control if we want to tweak later
  let baseWeight: number = 1.0;
  let collapseByGenus = false;
  
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
    showPrimary = false;
    showSecondary = false;
    showTarget = false;
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

  async function openInfoModal() {

    const modal: ModalSettings = {
      type: 'component',
      component: {
        ref: InfoModal,
        // optional props
        props: { title: 'Wanted: Bugs!', confirmLabel: 'Ok' },
        // default slot HTML
        slot: `
          <div class="text-base">
            <p class="mt-2">You are a diagnostic assistant trained to support the interpretation of metagenomic sequencing results for infectious disease diagnosis. Use your expertise in microbiology, pathology, metagenomics, clinical diagnostics, and infectious diseases to help determine whether a case is infectious (positive) or non-infectious (negative).</p>
            <p class="mt-2">We conducted metagenomic sequencing for pathogen detection and diagnosis (Illumina PE, RNA and DNA libraries). Filtering the taxonomic profiling data from the bioinformatics pipeline produced three subsets of the same dataset: 
              primary threshold (present at high abundance), secondary threshold (present at moderate to low abundance), and a target threshold for high priority pathogens (present at low abundance). Our pipeline uses multiple methods for 
              pathogen detection - alignment (reads per million, RPM), k-mer classifiers (read per million, RPM) and metagenome assembly (contigs, bases). RPM and contigs/bases provided for each detected species are the primary evidence.</p>
            <p class="mt-2">Species names are taxonomic species names (genus name and species name). If you do not know a species, assume that the provided species name is correct - do not interpret unknown species names as another species you know. 
              You must make your considerations and determinations based on the species, not the genus. For viral species, traditional abbreviations and alternative names may be provided, which you can assume are valid alternatives to 
              the output species name.</p>
              <p class="mt-2">You can use the 'Best Species' slider to collapse genera with more than four species represented in the data and retain only the best species for each genus (Archaea/Bacteria/Eukaryota)</p>
          </div>` 
      },
      modalClasses: 'w-[80vw] max-w-[80vw] h-[80vh] !max-h-[80vh] flex', 
      backdropClasses: 'backdrop-blur-md bg-black/40 items-center',
    };
    modalStore.trigger(modal);
  }

</script>
  
    <div class="w-full flex items-center justify-center gap-3 py-3 sticky top-0 bg-base-100/80 backdrop-blur z-10">
      
      <button class="px-3 py-1 rounded disabled:opacity-40 flex items-center" on:click={openInfoModal}>
              <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-8 h-8 stroke-success-500 mr-2">
                <path stroke-linecap="round" stroke-linejoin="round" d="m11.25 11.25.041-.02a.75.75 0 0 1 1.063.852l-.708 2.836a.75.75 0 0 0 1.063.853l.041-.021M21 12a9 9 0 1 1-18 0 9 9 0 0 1 18 0Zm-9-3.75h.008v.008H12V8.25Z" />
              </svg>
              <p class="opacity-80">Information</p>
      </button>
      
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


      <div class="flex items-center gap-2 ml-5">
        <SlideToggle name="collapseByGenus" bind:checked={collapseByGenus} />
        <span class="text-sm opacity-70 select-none">Best Species</span>
      </div>

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
            <PrefetchTable taxa={trainingData.prefetch.primary}  selectedName={selectedCandidateName ?? ""}  
            collapseByGenus={collapseByGenus} on:select={onSelect} />
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
            <PrefetchTable taxa={trainingData.prefetch.secondary}   selectedName={selectedCandidateName ?? ""} 
            collapseByGenus={collapseByGenus} on:select={onSelect} />
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
            <PrefetchTable taxa={trainingData.prefetch.target}   selectedName={selectedCandidateName ?? ""} 
            collapseByGenus={collapseByGenus} on:select={onSelect} />
          {/if}
        </section>
      </div>  
  <style>
    /* Tailwind assumed; keep minimal if not using Tailwind */
  </style>
  