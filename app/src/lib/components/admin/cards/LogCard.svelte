
<script lang="ts">
    import { getDateTime } from "$lib/utils/helpers";
	import type { RequestLog } from "$lib/utils/types";
    
    export let requestLog: RequestLog;

    let date: string;
    let time: string;

    $: {
        [date, time] = getDateTime(requestLog.date);
    }
    

</script>

<div class="card max-w-xl2">
	<header class="card-header">
        <div class="p-4 space-y-4">
            <h6 class="h6">{date} @ {time}</h6>
            <h4 class="h4">{requestLog.module} - {requestLog.action}</h4>
            <hr>
        </div>
    </header>
	<section class="p-8 space-y-4">
        <div class="space-y-2">
            <p><span class="opacity-60 ">Description</span></p>
            <p> {requestLog.description}</p>
        </div>
        <div  class="space-y-2">
            <p><span class="opacity-60">User email</span></p>
            <p> {requestLog.access_details.user_email}</p>
        </div>
        <div>
             <hr>
        </div>
        <div class="space-y-2">
            <p><span class="opacity-60">Request route</span></p>
            <p> {requestLog.access_details.request_details.route}</p>
        </div>
        <div class="space-y-2">
            <p><span class="opacity-60">Request method</span></p>
            <p> {requestLog.access_details.request_details.method}</p>
        </div>
        <div>
             <hr>
        </div>
    </section>
</div>