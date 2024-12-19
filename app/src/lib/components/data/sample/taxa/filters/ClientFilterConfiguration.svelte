<script lang="ts">
	import type { ClientFilterConfig } from "$lib/utils/types";
    import { DomainName } from "$lib/utils/types";
	import { Autocomplete, InputChip } from "@skeletonlabs/skeleton";
    import type { AutocompleteOption } from "@skeletonlabs/skeleton";

    export let clientFilterConfig: ClientFilterConfig;

    let numericInputClass: string = "[appearance:textfield] [&::-webkit-outer-spin-button]:appearance-none [&::-webkit-inner-spin-button]:appearance-none";


    const domainOptions: AutocompleteOption<string, null>[] = [
        { label: 'Eukaryota', value: 'Eukaryota', keywords: 'eukaryotes, protozoan, fungus, amoeba, worms' },
        { label: 'Bacteria', value: 'Bacteria', keywords: 'bacteria' },
        { label: 'Viruses', value: 'Viruses', keywords: 'viruses' },
        { label: 'Archaea', value: 'Archaea', keywords: 'archaea' }
    ];


    function modifyDomainFilter(domain: DomainName): void {
        if (clientFilterConfig.domains.includes(domain) === false) {
			clientFilterConfig.domains = [...clientFilterConfig.domains, domain];
		} else {
            clientFilterConfig.domains = clientFilterConfig.domains.filter(d => d !== domain)
        }
        console.log(clientFilterConfig.domains)
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

    let genusChip: string = '';
    let speciesChip: string = '';

</script>

<div>
    <p class=""><span class="opacity-40">Filter by taxonomy</span></p>
    <div class="p-4 w-full">
        <p class="text-xs opacity-40"></p>
        <div class="py-4">
            <span class="chip {clientFilterConfig.domains.includes(DomainName.Archaea) ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { modifyDomainFilter(DomainName.Archaea) }} on:keypress aria-hidden>
                {#if clientFilterConfig.domains.includes(DomainName.Archaea)}
                    <span>
                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-3 h-3" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg>
                    </span>
                {/if}
                <span>{DomainName.Archaea}</span>
            </span>
            <span class="chip {clientFilterConfig.domains.includes(DomainName.Bacteria) ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { modifyDomainFilter(DomainName.Bacteria) }} on:keypress  aria-hidden>
                {#if clientFilterConfig.domains.includes(DomainName.Bacteria)}
                    <span>
                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-3 h-3" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg>
                    </span>
                {/if}
                <span>{DomainName.Bacteria}</span>
            </span>
            <span class="chip {clientFilterConfig.domains.includes(DomainName.Eukaryota) ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { modifyDomainFilter(DomainName.Eukaryota) }} on:keypress  aria-hidden>
                {#if clientFilterConfig.domains.includes(DomainName.Eukaryota)}
                    <span>
                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-3 h-3" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg>
                    </span>
                {/if}
                <span>{DomainName.Eukaryota}</span>
            </span>
            <span class="chip {clientFilterConfig.domains.includes(DomainName.Viruses) ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { modifyDomainFilter(DomainName.Viruses) }} on:keypress  aria-hidden>
                {#if clientFilterConfig.domains.includes(DomainName.Viruses)}
                    <span>
                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-3 h-3" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg>
                    </span>
                {/if}
                <span>{DomainName.Viruses}</span>
            </span>
        </div>
        <div class="pb-2">
            <InputChip bind:input={genusChip} bind:value={clientFilterConfig.genera} name="genus" placeholder="Genus" allowUpperCase class="w-3/4" />
        </div>
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
            <InputChip bind:input={speciesChip} bind:value={clientFilterConfig.species} name="species" placeholder="Species" allowUpperCase class="w-3/4"/>
        </div>
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

    <p class=""><span class="opacity-40">Modules</span></p>
    <div class="p-4">
        <span class="chip {clientFilterConfig.modules.alignment ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.modules.alignment ? clientFilterConfig.modules.alignment = false : clientFilterConfig.modules.alignment = true }} on:keypress aria-hidden>
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
        <span class="chip {clientFilterConfig.modules.profile ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.modules.profile ? clientFilterConfig.modules.profile = false : clientFilterConfig.modules.profile = true }} on:keypress  aria-hidden>
            {#if clientFilterConfig.modules.profile}
                <span>
                    <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" class="w-3 h-3" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                        <path d="M4.5 12.75l6 6 9-13.5" stroke-linecap="round" stroke-linejoin="round"></path>
                    </svg>
                </span>
            {/if}
            <span>K-mer</span>
            <div class="rounded-full bg-secondary-500 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.modules.assembly ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.modules.assembly ? clientFilterConfig.modules.assembly = false : clientFilterConfig.modules.assembly = true }} on:keypress  aria-hidden>
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