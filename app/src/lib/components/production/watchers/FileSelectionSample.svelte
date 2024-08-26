<script lang="ts">
	import { page } from "$app/stores";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import { FileTag, SeaweedFileType, type FileTagUpdateSchema, type SeaweedFile } from "$lib/utils/types";

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

    // Check if the files are of ftype ReadsPaired and if they are not a multiple of 2
    let disableDiv: boolean = files.filter(file => file.ftype === SeaweedFileType.ReadsPaired).length % 2 !== 0;



</script>

<div class="grid grid-cols-6 gap-8 items-center">
    <div class="opacity-60 col-span-2 truncate">{id}</div>
   
    {#if disableDiv}
        <div class="col-span-1 flex justify-start opacity-60 text-red-500">
            <div class="flex items-center">
                <div class="w-5 h-5 mr-2">
                    <svg data-slot="icon" aria-hidden="true" fill="none" stroke-width="1.5" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                        <path d="M12 9v3.75m-9.303 3.376c-.866 1.5.217 3.374 1.948 3.374h14.71c1.73 0 2.813-1.874 1.948-3.374L13.949 3.378c-.866-1.5-3.032-1.5-3.898 0L2.697 16.126ZM12 15.75h.007v.008H12v-.008Z" stroke-linecap="round" stroke-linejoin="round"></path>
                    </svg>
                </div>
                Unpaired Reads
            </div>
        </div>
    {:else}
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
    {/if}
</div>