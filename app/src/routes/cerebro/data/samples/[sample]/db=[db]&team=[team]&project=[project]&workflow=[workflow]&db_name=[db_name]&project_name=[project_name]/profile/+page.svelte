<script lang="ts">
	import TaxonomicProfile from '$lib/components/data/sample/TaxonomicProfile.svelte';
    import { selectedIdentifiers, selectedModels, selectedClientFilterConfig, selectedServerFilterConfig, selectedTaxonHighlightConfig } from '$lib/stores/stores';
</script>

<div class="w-full mb-5">
    <p class="mb-1"><span class="opacity-60">Taxonomic Profile</span></p>
    <div class="mb-4 border border-primary-500 p-8 rounded-md">
            <TaxonomicProfile
                selectedIdentifiers={$selectedIdentifiers}
                selectedModels={$selectedModels}
                clientFilterConfig={$selectedClientFilterConfig}
                serverFilterConfig={$selectedServerFilterConfig}
                taxonHighlightConfig={$selectedTaxonHighlightConfig}
            />
    </div>
</div>

