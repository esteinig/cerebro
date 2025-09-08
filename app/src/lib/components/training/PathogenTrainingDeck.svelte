<!-- src/lib/components/training/PathogenTrainingDeck.svelte -->
<script lang="ts">
    import SpeciesOverviewTable from "../data/sample/taxa/SpeciesOverviewTable.svelte";
    import TaxonHistory from "$lib/components/visualisations/taxa/history/TaxonHistory.svelte";
    import { DisplayData, DisplayVisualisation, AbundanceMode } from "$lib/utils/types";
    import { selectedTaxa as selectedTaxaStore } from "$lib/stores/stores";
    import type { DeckConfig, SlideConfig, SlideResponse, UUID } from "./types";
    import { SlideDecision } from "./types";
    import { onMount } from "svelte";
  
    export let deck: DeckConfig;             // provided by page load
    export let userId: UUID;                 // provided by page (validated on server)
  
    let current = 0;
    let sessionId: UUID | null = null;
  
    // selection snapshot per slide
    const perSlideSelections = new Map<string, { taxid: string; name: string }[]>();
    const perSlideComments = new Map<string, string>(); // optional
  
    // timing
    let deckStartMs = 0;
    let slideStartMs = 0;
    const perSlideMs = new Map<string, number>();
  
    // queue of emitted events (for robust retry you could persist this)
    const responses: SlideResponse[] = [];
  
    $: slide = deck.slides[current];
  
    function now() { return performance.now(); }
  
    function resetSelectedStore() {
      // wipe UI selection in the shared store when we switch slides
      selectedTaxaStore.set([]);
    }
    
    let commentText = "";

    function loadSelectionsForSlide(sc: SlideConfig) {
      resetSelectedStore();
      const prev = perSlideSelections.get(sc.id) ?? [];
      // restore per-slide comment
      commentText = perSlideComments.get(sc.id) ?? "";
    }

    function onCommentInput(e: Event) {
      const val = (e.target as HTMLTextAreaElement).value;
      commentText = val;
      perSlideComments.set(slide.id, val);
    }

    function captureCurrentSelection(sc: SlideConfig) {
      const snapshot = $selectedTaxaStore.map(t => ({ taxid: t.taxid, name: t.name }));
      perSlideSelections.set(sc.id, snapshot);
      return snapshot;
    }
  
    async function emitDecision(decision: SlideDecision) {
      if (!sessionId) return;
  
      // capture selection
      const selected = captureCurrentSelection(slide);
      const minSel = deck.minSelection ?? 1;
      if (selected.length < minSel) {
        // TODO: replace with toast store if preferred
        alert(`Please select at least ${minSel} taxon row(s) before submitting.`);
        return;
      }
  
      // accumulate slide timing
      const spent = now() - slideStartMs;
      perSlideMs.set(slide.id, (perSlideMs.get(slide.id) ?? 0) + spent);
      slideStartMs = now(); // reset for next slide entry
  
      const event: SlideResponse = {
        slideId: slide.id,
        decision,
        selectedTaxa: selected,
        msSpent: perSlideMs.get(slide.id)!,
        comment: perSlideComments.get(slide.id),
        clientEventId: crypto.randomUUID()
      };
  
      responses.push(event);
  
      // INSERT POINT: send to API (best-effort)
      try {
        // await api.postEvent(sessionId, event);
        console.error("postEvent failed; will rely on completeSession to resend aggregate");
      } catch (e) {
        console.error("postEvent failed; will rely on completeSession to resend aggregate", e);
      }
  
      // move next if not at end
      if (current < deck.slides.length - 1) {
        current += 1;
      } else {
        await finishDeck();
      }
    }
  
    async function finishDeck() {
      if (!sessionId) return;
      const totalMs = now() - deckStartMs;
  
      try {
        // await api.completeSession(sessionId, totalMs);
        alert("Training deck complete. Thanks!");
      } finally {
        // UX: navigate to a “done” state; keep it simple here
        alert("Training deck complete. Thanks!");
      }
    }
  
    function goNext() {
      captureCurrentSelection(slide);
      if (current < deck.slides.length - 1) {
        // accumulate time for the slide before switching
        const spent = now() - slideStartMs;
        perSlideMs.set(slide.id, (perSlideMs.get(slide.id) ?? 0) + spent);
        current += 1;
        slideStartMs = now();
        loadSelectionsForSlide(deck.slides[current]);
      }
    }
  
    function goBack() {
      captureCurrentSelection(slide);
      if (current > 0) {
        const spent = now() - slideStartMs;
        perSlideMs.set(slide.id, (perSlideMs.get(slide.id) ?? 0) + spent);
        current -= 1;
        slideStartMs = now();
        loadSelectionsForSlide(deck.slides[current]);
      }
    }
  
    function clearSelection() {
      selectedTaxaStore.set([]);
      perSlideSelections.delete(slide.id);
    }
  
    onMount(async () => {
      deckStartMs = now();
      slideStartMs = now();

      // // INSERT POINT: start session
      // const s = await api.startSession(deck.id, userId);
      // sessionId = s.sessionId;
  
      // load initial selections (none on first run)
      loadSelectionsForSlide(deck.slides[0]);
    });
</script>
  
<!-- Controls header -->
<div class="flex items-center justify-between gap-3 py-2">
<div class="flex items-center gap-2">
    <span class="opacity-70 text-sm">{deck.name}</span>
    <span class="opacity-50 text-xs">slide {current+1} / {deck.slides.length}</span>
</div>
<div class="flex items-center gap-2">
    <button class="btn btn-sm variant-ghost" on:click={goBack} disabled={current===0}>Back</button>
    <button class="btn btn-sm variant-ghost" on:click={goNext} disabled={current===deck.slides.length-1}>Next</button>
    <button class="btn btn-sm variant-soft" on:click={clearSelection}>Clear selection</button>
    <button class="btn btn-sm variant-filled-primary" on:click={() => emitDecision(SlideDecision.Positive)}>Mark Positive</button>
    <button class="btn btn-sm variant-filled-secondary" on:click={() => emitDecision(SlideDecision.Negative)}>Mark Negative</button>
</div>
</div>

<!-- Slide intro -->
<div class="my-3 p-3 rounded-token variant-soft">
<div class="text-base font-medium">{slide.title}</div>
{#if slide.description}
    <div class="opacity-70 text-sm mt-1">{slide.description}</div>
{/if}
</div>

<!-- Optionally show history for selected taxa -->
{#if $selectedTaxaStore.length > 0}
<div class="mb-4">
    {#each $selectedTaxaStore as t}
    <TaxonHistory selectedTaxon={t}/>
    {/each}
</div>
{/if}

<!-- The table itself with per-slide configuration -->
<SpeciesOverviewTable
  displayData={slide.forceRpmForValidationPlate ? DisplayData.Rpm : slide.displayData}
  displayMode={slide.displayMode}
  selectedVisualisation={slide.selectedVisualisation}
  serverFilterConfig={slide.serverFilterConfig ?? null}
  disablePrevalenceContamination={slide.disablePrevalenceContamination ?? false}
/>

<!-- Optional comment box -->
<div class="mt-4">
<label class="text-sm opacity-70">Comment (optional)
    <textarea
    class="textarea w-full mt-1"
    rows="3"
    bind:value={commentText}
    on:input={onCommentInput}
    placeholder="Any notes about this decision…"
    />
</label>
</div>
  