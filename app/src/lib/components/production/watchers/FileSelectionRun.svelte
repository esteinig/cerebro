<script lang="ts">
	import { getDateTimeStringUtc } from "$lib/utils/helpers";
	import type { SeaweedFile } from "$lib/utils/types";

    export let id: string;
    export let files: SeaweedFile[];

    let numberLibraries = new Set(files.map(file => file.sample_id)).size;
    let numberFiles =  files.length;

    const dates = files.map(file => file.date);
    const minDate = dates.reduce((a, b) => a < b ? a : b);
</script>


<div class="grid grid-cols-7 gap-8 items-center">
    <div class="opacity-60 col-span-2  truncate">{id}</div> 
    <div class="opacity-60 col-span-1">{numberLibraries} {numberLibraries == 1 ? "library" : "libraries"}</div>
    <div class="opacity-60 col-span-1">{numberFiles} {numberFiles == 1 ? "file" : "files"}</div>
    <div class="opacity-60 col-span-3"><span class="mr-2">{getDateTimeStringUtc(minDate, true, true)}</span></div>
</div>
    