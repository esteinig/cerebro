<script lang="ts">
	import type { HighlightConfig } from "$lib/utils/types";
	import { Autocomplete, InputChip } from "@skeletonlabs/skeleton";
    import type { AutocompleteOption } from "@skeletonlabs/skeleton";

    export let title: string;
    export let highlightConfig: HighlightConfig;

    let speciesOptions: AutocompleteOption<string, null>[] = highlightConfig.species.map(species => {
        return { label: species, value: species, keywords: species }
    });

    const speciesWhitelist: string[] = speciesOptions.map(opt => opt.value)

    function onSpeciesSelection(event: CustomEvent<AutocompleteOption>): void {
        if (highlightConfig.species.includes(event.detail.value as string) === false) {
			highlightConfig.species = [...highlightConfig.species, event.detail.value as string];
			speciesChip = '';
		}
    }
    let speciesChip: string = '';

</script>

<div>
    <p class=""><span class="opacity-40">{title}</span></p>
    <div class="p-4 w-full">
        <p class="text-xs opacity-40"></p>
        <div class="pb-2">
            <InputChip bind:input={speciesChip} bind:value={highlightConfig.species} name="species" placeholder="Species" whitelist={speciesWhitelist} allowUpperCase/>
        </div>
        {#if speciesChip.length}
            <div class="my-3">
                <div class="card w-full max-w-3/4 p-4 overflow-y-auto text-sm" tabindex="-1">
                    <Autocomplete
                        bind:input={speciesChip}
                        options={speciesOptions}
                        denylist={highlightConfig.species}
                        on:selection={onSpeciesSelection}
                    />
                </div>
            </div>
        {/if}
    </div>
</div>