<script lang="ts">
	import CheckmarkIcon from "$lib/general/icons/CheckmarkIcon.svelte";
	import type { ClientFilterConfig } from "$lib/utils/types";
    import { DomainName, ProfileTool } from "$lib/utils/types";
	import { Autocomplete, InputChip } from "@skeletonlabs/skeleton";
    import type { AutocompleteOption } from "@skeletonlabs/skeleton";

    export let clientFilterConfig: ClientFilterConfig;

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

    function updateToolSelection(tool: ProfileTool) {
        if (clientFilterConfig.tools.includes(tool)) {
            clientFilterConfig.tools = clientFilterConfig.tools.filter(t => t !== tool)
        } else {
            clientFilterConfig.tools = [...clientFilterConfig.tools, tool]
        }
    }

    let genusChip: string = '';
    let speciesChip: string = '';

</script>

<div>
    <p class=""><span class="opacity-40">Taxonomy</span></p>
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

    <p class="pt-4"><span class="opacity-40">Profiling modules</span></p>
    <div class="p-4">
        <span class="chip {clientFilterConfig.modules.alignment ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.modules.alignment ? clientFilterConfig.modules.alignment = false : clientFilterConfig.modules.alignment = true }} on:keypress aria-hidden>
            {#if clientFilterConfig.modules.alignment}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>Alignment</span>
            <div class="rounded-full bg-primary-500 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.modules.profile ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.modules.profile ? clientFilterConfig.modules.profile = false : clientFilterConfig.modules.profile = true }} on:keypress  aria-hidden>
            {#if clientFilterConfig.modules.profile}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>K-mer</span>
            <div class="rounded-full bg-secondary-500 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.modules.assembly ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.modules.assembly ? clientFilterConfig.modules.assembly = false : clientFilterConfig.modules.assembly = true }} on:keypress  aria-hidden>
            {#if clientFilterConfig.modules.assembly}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>Assembly</span>
            <div class="rounded-full bg-tertiary-500 h-2 w-2 ml-2"></div>
        </span>
    </div>
    <p class="pt-4"><span class="opacity-40">Profiling tools</span></p>
    <div class="p-4 grid grid-cols-4 gap-2">
        <span class="chip {clientFilterConfig.tools.includes(ProfileTool.Vircov) ? 'variant-ghost-primary' : 'variant-soft'} mr-2" on:click={() => { updateToolSelection(ProfileTool.Vircov) }} on:keypress aria-hidden>
            {#if clientFilterConfig.tools.includes(ProfileTool.Vircov)}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>{ProfileTool.Vircov}</span>
            <div class="rounded-full bg-primary-400 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.tools.includes(ProfileTool.Kraken2) ? 'variant-ghost-secondary' : 'variant-soft'} mr-2" on:click={() => { updateToolSelection(ProfileTool.Kraken2) }} on:keypress aria-hidden>
            {#if clientFilterConfig.tools.includes(ProfileTool.Kraken2)}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>{ProfileTool.Kraken2}</span>
            <div class="rounded-full bg-secondary-800 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.tools.includes(ProfileTool.Bracken) ? 'variant-ghost-secondary' : 'variant-soft'} mr-2" on:click={() => { updateToolSelection(ProfileTool.Bracken) }} on:keypress aria-hidden>
            {#if clientFilterConfig.tools.includes(ProfileTool.Bracken)}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>{ProfileTool.Bracken}</span>
            <div class="rounded-full bg-secondary-700 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.tools.includes(ProfileTool.Metabuli) ? 'variant-ghost-secondary' : 'variant-soft'} mr-2" on:click={() => { updateToolSelection(ProfileTool.Metabuli) }} on:keypress aria-hidden>
            {#if clientFilterConfig.tools.includes(ProfileTool.Metabuli)}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>{ProfileTool.Metabuli}</span>
            <div class="rounded-full bg-secondary-600 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.tools.includes(ProfileTool.Ganon2) ? 'variant-ghost-secondary' : 'variant-soft'} mr-2" on:click={() => { updateToolSelection(ProfileTool.Ganon2) }} on:keypress aria-hidden>
            {#if clientFilterConfig.tools.includes(ProfileTool.Ganon2)}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>{ProfileTool.Ganon2}</span>
            <div class="rounded-full bg-secondary-500 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.tools.includes(ProfileTool.Sylph) ? 'variant-ghost-secondary' : 'variant-soft'} mr-2" on:click={() => { updateToolSelection(ProfileTool.Sylph) }} on:keypress aria-hidden>
            {#if clientFilterConfig.tools.includes(ProfileTool.Sylph)}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>{ProfileTool.Sylph}</span>
            <div class="rounded-full bg-secondary-400 h-2 w-2 ml-2"></div>
        </span>
        <span class="chip {clientFilterConfig.tools.includes(ProfileTool.Blast) ? 'variant-ghost-tertiary' : 'variant-soft'} mr-2" on:click={() => { updateToolSelection(ProfileTool.Blast) }} on:keypress aria-hidden>
            {#if clientFilterConfig.tools.includes(ProfileTool.Blast)}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>{ProfileTool.Blast}</span>
            <div class="rounded-full bg-tertiary-500 h-2 w-2 ml-2"></div>
        </span>
    </div>
    <p class="pt-4"><span class="opacity-40">Display filtered</span></p>
    <div class="p-4">
        <span class="chip {clientFilterConfig.contam.display ? 'variant-ghost-surface' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.contam.display ? clientFilterConfig.contam.display = false : clientFilterConfig.contam.display = true }} on:keypress aria-hidden>
            {#if clientFilterConfig.contam.display}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>Prevalence Contamination</span>
        </span>
        <span class="chip {clientFilterConfig.evidence.display ? 'variant-ghost-surface' : 'variant-soft'} mr-2" on:click={() => { clientFilterConfig.evidence.display ? clientFilterConfig.evidence.display = false : clientFilterConfig.evidence.display = true }} on:keypress aria-hidden>
            {#if clientFilterConfig.evidence.display}
                <span><CheckmarkIcon /></span>
            {/if}
            <span>NTC Comparison (no evidence)</span>
        </span>
    </div>
</div>