<script lang="ts">
	import { page } from "$app/stores";
    import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    let selectedTeam: string = $page.data.userTeams[0]?.name ?? "None";
        
    let selectedView: string = "samples";
</script>

<div>
    <p class="opacity-60">Production operations</p>
    <div class="card w-full border border-primary-500">
        <section>
            <div class="p-4">
                <ListBox>
                    <a href='/cerebro/production/watchers/team={selectedTeam}&watcher=0' class="block pointer-events-none opacity-90" tabindex="-1">
                        <ListBoxItem bind:group={selectedView} name="medium" value="watchers" active='variant-soft' rounded='rounded-token' disabled>
                            <svelte:fragment slot="lead">
                                <div class="w-4 h-4">
                                    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-satellite-dish"><path d="M4 10a7.31 7.31 0 0 0 10 10Z"/><path d="m9 15 3-3"/><path d="M17 13a6 6 0 0 0-6-6"/><path d="M21 13A10 10 0 0 0 11 3"/></svg>
                                </div>
                            </svelte:fragment> 
                                Pipeline Deployment
                        </ListBoxItem>
                    </a>
                    <ListBoxItem bind:group={selectedView} name="medium" value="qa" active='variant-soft' rounded='rounded-token' disabled>
                        <svelte:fragment slot="lead">
                            <div class="w-4 h-4">
                                 <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-book-open-check"><path d="M8 3H2v15h7c1.7 0 3 1.3 3 3V7c0-2.2-1.8-4-4-4Z"/><path d="m16 12 2 2 4-4"/><path d="M22 6V3h-6c-2.2 0-4 1.8-4 4v14c0-1.7 1.3-3 3-3h7v-2.3"/></svg>
                            </div>
                        </svelte:fragment> 
                        Quality Assurance
                    </ListBoxItem>
                    <a href='/cerebro/review/team={selectedTeam}' >
                        <ListBoxItem bind:group={selectedView} name="medium" value="review" active='variant-soft' rounded='rounded-token'>
                            <svelte:fragment slot="lead">
                                <div class="w-4 h-4">
                                    <svg data-slot="icon" aria-hidden="true" fill="none" stroke-width="1.5" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                        <path d="M12 12.75c1.148 0 2.278.08 3.383.237 1.037.146 1.866.966 1.866 2.013 0 3.728-2.35 6.75-5.25 6.75S6.75 18.728 6.75 15c0-1.046.83-1.867 1.866-2.013A24.204 24.204 0 0 1 12 12.75Zm0 0c2.883 0 5.647.508 8.207 1.44a23.91 23.91 0 0 1-1.152 6.06M12 12.75c-2.883 0-5.647.508-8.208 1.44.125 2.104.52 4.136 1.153 6.06M12 12.75a2.25 2.25 0 0 0 2.248-2.354M12 12.75a2.25 2.25 0 0 1-2.248-2.354M12 8.25c.995 0 1.971-.08 2.922-.236.403-.066.74-.358.795-.762a3.778 3.778 0 0 0-.399-2.25M12 8.25c-.995 0-1.97-.08-2.922-.236-.402-.066-.74-.358-.795-.762a3.734 3.734 0 0 1 .4-2.253M12 8.25a2.25 2.25 0 0 0-2.248 2.146M12 8.25a2.25 2.25 0 0 1 2.248 2.146M8.683 5a6.032 6.032 0 0 1-1.155-1.002c.07-.63.27-1.222.574-1.747m.581 2.749A3.75 3.75 0 0 1 15.318 5m0 0c.427-.283.815-.62 1.155-.999a4.471 4.471 0 0 0-.575-1.752M4.921 6a24.048 24.048 0 0 0-.392 3.314c1.668.546 3.416.914 5.223 1.082M19.08 6c.205 1.08.337 2.187.392 3.314a23.882 23.882 0 0 1-5.223 1.082" stroke-linecap="round" stroke-linejoin="round"></path>
                                    </svg>
                                </div>
                            </svelte:fragment> 
                            Bug Board
                        </ListBoxItem>
                    </a>
                    <a href='/cerebro/training/team={selectedTeam}' >
                        <ListBoxItem bind:group={selectedView} name="medium" value="training" active='variant-soft' rounded='rounded-token'>
                            <svelte:fragment slot="lead">
                                <div class="w-4 h-4">
                                    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-book-open-check"><path d="M8 3H2v15h7c1.7 0 3 1.3 3 3V7c0-2.2-1.8-4-4-4Z"/><path d="m16 12 2 2 4-4"/><path d="M22 6V3h-6c-2.2 0-4 1.8-4 4v14c0-1.7 1.3-3 3-3h7v-2.3"/></svg>
                                </div>
                            </svelte:fragment> 
                            Training
                        </ListBoxItem>
                    </a>
                </ListBox>
            </div>
        </section>
        <footer class="card-footer">
        </footer>
    </div>
</div>