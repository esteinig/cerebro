<script lang="ts">
	import type { Cerebro, ClientFilterConfig } from "$lib/utils/types";
	import { Autocomplete, InputChip } from "@skeletonlabs/skeleton";
    import type { AutocompleteOption } from "@skeletonlabs/skeleton";

    export let clientFilterConfig: ClientFilterConfig;
    export let selectedModels: Cerebro[] = [];

    let numericInputClass: string = "[appearance:textfield] [&::-webkit-outer-spin-button]:appearance-none [&::-webkit-inner-spin-button]:appearance-none";

    const domainOptions: AutocompleteOption<string, null>[] = [
        { label: 'Eukaryota', value: 'Eukaryota', keywords: 'eukaryotes, protozoan, fungus, amoeba, worms' },
        { label: 'Bacteria', value: 'Bacteria', keywords: 'bacteria' },
        { label: 'Viruses', value: 'Viruses', keywords: 'viruses' },
        { label: 'Archaea', value: 'Archaea', keywords: 'archaea' }
    ];

    const domainWhitelist: string[] = domainOptions.map(opt => opt.value)


    function onDomainSelection(event: CustomEvent<AutocompleteOption>): void {
        if (clientFilterConfig.domains.includes(event.detail.value as string) === false) {
			clientFilterConfig.domains = [...clientFilterConfig.domains, event.detail.value as string];
			domainChip = '';
		}
    }

    function onGenusSelection(event: CustomEvent<AutocompleteOption>): void {
        if (clientFilterConfig.genera.includes(event.detail.value as string) === false) {
			clientFilterConfig.genera = [...clientFilterConfig.genera, event.detail.value as string];
			genusChip = '';
		}
    }

    function onSpeciesSelection(event: CustomEvent<AutocompleteOption>): void {
        if (clientFilterConfig.species.includes(event.detail.value as string) === false) {
			clientFilterConfig.species = [...clientFilterConfig.species, event.detail.value as string];
			speciesChip = '';
		}
    }

    let domainChip: string = '';
    let genusChip: string = '';
    let speciesChip: string = '';

</script>

<div>
    

    <p class=""><span class="opacity-40">Taxonomy filters</span></p>
    <div class="p-4 w-full">
        <p class="text-xs opacity-40"></p>
        <div class="pb-2">
            <InputChip bind:input={domainChip} bind:value={clientFilterConfig.domains} name="domains" placeholder="Domains" whitelist={domainWhitelist} allowUpperCase/></div>
            {#if domainChip.length}
                <div class="my-3">
                    <div class="card w-full max-w-3/4 p-4 overflow-y-auto text-sm" tabindex="-1">
                        <Autocomplete
                            bind:input={domainChip}
                            options={domainOptions}
                            denylist={clientFilterConfig.domains}
                            on:selection={onDomainSelection}
                        />
                    </div>
                </div>
            {/if}
        <div class="pb-2">
        <InputChip bind:input={genusChip} bind:value={clientFilterConfig.genera} name="genus" placeholder="Genera" allowUpperCase /></div>
        {#if genusChip.length}
            <div class="my-3">
                <div class="card w-full max-w-3/4 p-4 overflow-y-auto text-sm" tabindex="-1">
                    <Autocomplete
                        bind:input={genusChip}
                        options={domainOptions}
                        denylist={clientFilterConfig.genera}
                        on:selection={onGenusSelection}
                    />
                </div>
            </div>
        {/if}
        <div class="pb-2">
        
            <InputChip bind:input={speciesChip} bind:value={clientFilterConfig.species} name="species" placeholder="Species" allowUpperCase /></div>
        {#if speciesChip.length}
            <div class="my-3">
                <div class="card w-full max-w-3/4 p-4 overflow-y-auto text-sm" tabindex="-1">
                    <Autocomplete
                        bind:input={speciesChip}
                        options={domainOptions}
                        denylist={clientFilterConfig.species}
                        on:selection={onSpeciesSelection}
                    />
                </div>
            </div>
        {/if}
    </div>

    <p class=""><span class="opacity-40">Overview metrics</span></p>
    <div class="p-4">

        <div class="flex gap-x-8 pb-4 w-3/4">
            <label class="label">
                <span class="text-xs opacity-60">Minimum RPM</span>
                <input type="number" class="input text-xs {numericInputClass}" bind:value={clientFilterConfig.minimum.rpm}/>
            </label>
            <label class="label">
                <span class="text-xs opacity-60">Minimum Contigs</span>
                <input type="number" class="input text-xs {numericInputClass}" bind:value={clientFilterConfig.minimum.contigs}/>
            </label>
        </div>
    </div>

    <p class=""><span class="opacity-40">Module classification</span></p>
    <div class="p-4">
        <span class="chip {clientFilterConfig.modules.alignment ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.modules.alignment ? clientFilterConfig.modules.alignment = false : clientFilterConfig.modules.alignment = true }} on:keypress>
            {#if clientFilterConfig.modules.alignment}
                <span>
                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-3 h-3" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                        <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                    </svg>
                </span>
            {/if}
            <span>Alignment</span>
            <div class="rounded-full bg-primary-500 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.modules.kmer ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.modules.kmer ? clientFilterConfig.modules.kmer = false : clientFilterConfig.modules.kmer = true }} on:keypress>
            {#if clientFilterConfig.modules.kmer}
                <span>
                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-3 h-3" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                        <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                    </svg>
                </span>
            {/if}
            <span>K-mer</span>
            <div class="rounded-full bg-secondary-500 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.modules.assembly ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.modules.assembly ? clientFilterConfig.modules.assembly = false : clientFilterConfig.modules.assembly = true }} on:keypress>
            {#if clientFilterConfig.modules.assembly}
                <span>
                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-3 h-3" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                        <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                    </svg>
                </span>
            {/if}
            <span>Assembly</span>
            <div class="rounded-full bg-tertiary-500 h-2 w-2 ml-2"></div>
        </span>
    </div>
</div>