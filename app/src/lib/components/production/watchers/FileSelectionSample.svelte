<script lang="ts">
	import { LibraryControl, LibraryNucleicAcid, LibrarySpecimen, type SeaweedFile } from "$lib/utils/types";

    export let id: string;

    let nucleicAcidTag: LibraryNucleicAcid | null;
    let specimenTag: LibrarySpecimen | null;
    let controlTag: LibraryControl | null;

    let nucleicAcidTags: LibraryNucleicAcid[] = [
        LibraryNucleicAcid.DNA, 
        LibraryNucleicAcid.RNA
    ];
    let specimenTags: LibrarySpecimen[] = [
        LibrarySpecimen.CSF,
        LibrarySpecimen.BLOOD,
        LibrarySpecimen.RESPIRATORY
    ];
    let controlTags: LibraryControl[] = [
        LibraryControl.POSITIVE,
        LibraryControl.NEGATIVE,
        LibraryControl.TEMPLATE,
        LibraryControl.ENVIRONMENT
    ];


</script>

<div class="grid grid-cols-6 gap-16 items-center">
    <div class="opacity-60 col-span-2 truncate">{id}</div>
    <div class="col-span-1 flex justify-start">
        {#each nucleicAcidTags as tag}
            <button
                class="chip {nucleicAcidTag === tag ? 'variant-filled-primary' : 'variant-soft'} ml-2"
                on:click={() => { nucleicAcidTag = nucleicAcidTag === tag ? null : tag }}
                on:keypress
            >
                <span>{tag}</span>
            </button>
        {/each}
    </div>
    <div class="col-span-1 flex justify-start">
        {#each specimenTags as tag}
            <button
                class="chip {specimenTag === tag ? 'variant-filled-primary' : 'variant-soft'} ml-2"
                on:click={() => { specimenTag = specimenTag === tag ? null : tag  }}
                on:keypress
            >
                <span>{tag}</span>
            </button>
        {/each}
    </div>
    <div class="col-span-2 flex justify-start">
        {#each controlTags as tag}
            <button
                class="chip {controlTag === tag ? 'variant-filled-primary' : 'variant-soft'} ml-2"
                on:click={() => { controlTag = controlTag === tag ? null : tag }}
                on:keypress
            >
                <span>{tag}</span>
            </button>
        {/each}
    </div>
</div>