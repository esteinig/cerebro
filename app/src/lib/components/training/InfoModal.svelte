<script lang="ts">
  import { getModalStore } from '@skeletonlabs/skeleton';
  export let title = 'Good luck, bug hunter!';
  export let confirmLabel = 'Ok';
  export let cancelLabel = 'Close';

  const modalStore = getModalStore();

</script>

<div class="flex flex-col h-full">
  <header class="px-4 py-3 text-lg font-medium border-b">{title}</header>

  <!-- scrollable body -->
  <section class="flex-1 min-h-0 overflow-auto px-4 py-3 max-w-none">
    <slot />
  </section>

  <footer class="px-4 py-3 border-t flex justify-end gap-2">
    <button class="btn variant-outline" on:click={() => modalStore.close()}>
      {cancelLabel}
    </button>
    <button class="btn variant-filled" on:click={() => modalStore.close()}>
      {confirmLabel}
    </button>
  </footer>
</div>