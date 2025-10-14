<script lang="ts">
	import { goto } from "$app/navigation";


    import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    let selectedView: string = "runs";

</script>

<p class="opacity-60">Production interface</p>
<div class="card border p-4 border-primary-500">
    <ListBox class="">
        <a href="/cerebro/production/watchers" class="block pointer-events-none opacity-90" tabindex="-1">
            <ListBoxItem bind:group={selectedView} name="medium" value="runs", active='variant-soft' rounded='rounded-token' disabled>
            
                <svelte:fragment slot="lead">
                    <div class="w-4 h-4">
                        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-test-tubes"><path d="M9 2v17.5A2.5 2.5 0 0 1 6.5 22v0A2.5 2.5 0 0 1 4 19.5V2"/><path d="M20 2v17.5a2.5 2.5 0 0 1-2.5 2.5v0a2.5 2.5 0 0 1-2.5-2.5V2"/><path d="M3 2h7"/><path d="M14 2h7"/><path d="M9 16H4"/><path d="M20 16h-5"/></svg>
                    </div>
                </svelte:fragment> 
                Sequencing Runs 
            </ListBoxItem>
        </a>
        <a href="/cerebro/production/reports">
        <ListBoxItem bind:group={selectedView} name="medium" value="samples" active='variant-soft' rounded='rounded-token'>
            <svelte:fragment slot="lead">
                <div class="w-4 h-4">
                    <svg xmlns="http://www.w3.org/2000/svg"  viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-microscope"><path d="M6 18h8"/><path d="M3 22h18"/><path d="M14 22a7 7 0 1 0 0-14h-1"/><path d="M9 14h2"/><path d="M9 12a2 2 0 0 1-2-2V6h6v4a2 2 0 0 1-2 2Z"/><path d="M12 6V3a1 1 0 0 0-1-1H9a1 1 0 0 0-1 1v3"/></svg>
                </div>
            </svelte:fragment>
            Clinical Report
        </ListBoxItem>
        </a>
        <ListBoxItem bind:group={selectedView} name="medium" value="contam", active='variant-soft' rounded='rounded-token' disabled>
            <svelte:fragment slot="lead">
                <div class="w-4 h-4">
                    <svg data-slot="icon" aria-hidden="true" fill="none" stroke-width="1.5" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                        <path d="M12 12.75c1.148 0 2.278.08 3.383.237 1.037.146 1.866.966 1.866 2.013 0 3.728-2.35 6.75-5.25 6.75S6.75 18.728 6.75 15c0-1.046.83-1.867 1.866-2.013A24.204 24.204 0 0 1 12 12.75Zm0 0c2.883 0 5.647.508 8.207 1.44a23.91 23.91 0 0 1-1.152 6.06M12 12.75c-2.883 0-5.647.508-8.208 1.44.125 2.104.52 4.136 1.153 6.06M12 12.75a2.25 2.25 0 0 0 2.248-2.354M12 12.75a2.25 2.25 0 0 1-2.248-2.354M12 8.25c.995 0 1.971-.08 2.922-.236.403-.066.74-.358.795-.762a3.778 3.778 0 0 0-.399-2.25M12 8.25c-.995 0-1.97-.08-2.922-.236-.402-.066-.74-.358-.795-.762a3.734 3.734 0 0 1 .4-2.253M12 8.25a2.25 2.25 0 0 0-2.248 2.146M12 8.25a2.25 2.25 0 0 1 2.248 2.146M8.683 5a6.032 6.032 0 0 1-1.155-1.002c.07-.63.27-1.222.574-1.747m.581 2.749A3.75 3.75 0 0 1 15.318 5m0 0c.427-.283.815-.62 1.155-.999a4.471 4.471 0 0 0-.575-1.752M4.921 6a24.048 24.048 0 0 0-.392 3.314c1.668.546 3.416.914 5.223 1.082M19.08 6c.205 1.08.337 2.187.392 3.314a23.882 23.882 0 0 1-5.223 1.082" stroke-linecap="round" stroke-linejoin="round"></path>
                    </svg>
                </div>
            </svelte:fragment>
            Contamination
        </ListBoxItem>
        <ListBoxItem bind:group={selectedView} name="medium" value="assay", active='variant-soft' rounded='rounded-token' disabled>
	        <svelte:fragment slot="lead">
                <div class="w-4 h-4">
                    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-activity"><path d="M22 12h-4l-3 9L9 3l-3 9H2"/></svg>
                </div>
            </svelte:fragment>
            Assay Monitor
        </ListBoxItem>
    </ListBox>
</div>