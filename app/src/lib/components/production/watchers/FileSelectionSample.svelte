<script lang="ts">
	import { page } from "$app/stores";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import { FileTag, type FileTagUpdateSchema, type SeaweedFile } from "$lib/utils/types";

    const publicApi = new CerebroApi();

    export let id: string;
    export let files: SeaweedFile[];

    let nucleicAcidTag: FileTag | null = null;
    let controlTag: FileTag | null = null;

    let nucleicAcidTags: FileTag[] = [
        FileTag.DNA, 
        FileTag.RNA
    ];
    let controlTags: FileTag[] = [
        FileTag.POS,
        FileTag.NEG,
        FileTag.TMP,
        FileTag.ENV
    ];

    const updateFileTags = async() => {

        let fileTagUpdateSchema: FileTagUpdateSchema = {
            ids: files.map(file => file.id), tags: []
        }

        if (nucleicAcidTag !== null) {
            fileTagUpdateSchema.tags.push(nucleicAcidTag)
        }
        if (controlTag !== null) {
            fileTagUpdateSchema.tags.push(controlTag)
        }

        let _: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.files.updateTags}?team=${$page.params.team}`,
            { 
                method: 'PATCH',  
                mode: 'cors',
                credentials: 'include', 
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify(fileTagUpdateSchema) 
            } as RequestInit,
            $page.data.refreshToken
        )
    }

    // Iterate over the files
    for (const file of files) {
        // Check for nucleic acid tags
        for (const tag of file.tags) {
            if (nucleicAcidTags.includes(tag)) {
                nucleicAcidTag = tag;
                break; // Stop checking if a match is found
            }
        }
        // Check for control tags
        for (const tag of file.tags) {
            if (controlTags.includes(tag)) {
                controlTag = tag;
                break; // Stop checking if a match is found
            }
        }

        // If both tags are found, no need to continue
        if (nucleicAcidTag && controlTag) {
            break;
        }
    }

</script>

<div class="grid grid-cols-6 gap-8 items-center">
    <div class="opacity-60 col-span-2 truncate">{id}</div>
    <div class="col-span-1 flex justify-start">
        {#each nucleicAcidTags as tag}
            <button
                class="chip {nucleicAcidTag === tag ? 'variant-filled-primary' : 'variant-soft'} ml-2"
                on:click={async() => { nucleicAcidTag = nucleicAcidTag === tag ? null : tag; updateFileTags() }}
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
                on:click={async() =>  { controlTag = controlTag === tag ? null : tag; updateFileTags() }}
                on:keypress
            >
                <span>{tag}</span>
            </button>
        {/each}
    </div>
</div>