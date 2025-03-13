<script lang="ts">
	import { FileTag } from "$lib/utils/types";

    export let tags: string[] | FileTag[];
    export let join: boolean = false;
    export let separator: string = "-";
    export let variantColor: "primary" | "secondary" | "tertiary" = "primary";

    $: displayTags = join ? [tags.join(separator)] : tags;

</script>

{#each displayTags as tag}
    <span class="code bg-{variantColor}-500/30 text-{variantColor}-700 dark:bg-{variantColor}-500/20 dark:text-{variantColor}-400 ml-1">{tag}</span>
{/each}