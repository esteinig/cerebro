
<script lang="ts">
	import type { RequestLog } from "$lib/utils/types";
    
    export let requestLog: RequestLog;
    export let showHeaders: string[] = [];

    let requestHeaders: string[][]

    $: {
        requestHeaders = requestLog.access_details.request_details.headers.map(header => {
            let [key, value] = header.split(": ");
            if (key === "cookie") {
                return null
            }
            if (showHeaders.length > 0 && showHeaders.includes(key)) {
                return [key, value]
            } else {
                return [key, value]
            }
        }).flatMap(x => x ? [x] : []);

        requestHeaders.sort();
    }

</script>

<div class="card max-w-xl2">
	<section class="p-8 space-y-2">
       <div class="space-y-8">
           <p><span class="opacity-60">Request headers</span></p>

            <div class="grid grid-cols-1 sm:grid-cols-3 md:grid-cols-6 gap-4" >
            {#each requestHeaders as [header, value]}
                <div>
                    <p class="text-xs code py-2">{header}</p>
                    <p class="text-sm pb-2">{value}</p>
                </div>
            {/each}
           </div>
       </div>
       


    </section>
	<footer class="card-footer">

    </footer>
</div>