<script lang="ts">
	import { goto } from "$app/navigation";
	import { page } from "$app/stores";
    import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";

    let selectedTeam: string = $page.data.userTeams[0]?.name ?? "Default";
        
    let selectedView: string = "samples";
</script>

<div>
    <p class="opacity-60">Production operations</p>
    <div class="card w-full border border-primary-500">
        <section>
            <div class="p-4">
                <ListBox>
                    <ListBoxItem bind:group={selectedView} name="medium" value="watchers" active='variant-soft' rounded='rounded-token' on:click={() => goto(`./cerebro/production/watchers/team=${selectedTeam}&watcher=0`)}>
                        <svelte:fragment slot="lead">
                            <div class="w-4 h-4">
                                <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-satellite-dish"><path d="M4 10a7.31 7.31 0 0 0 10 10Z"/><path d="m9 15 3-3"/><path d="M17 13a6 6 0 0 0-6-6"/><path d="M21 13A10 10 0 0 0 11 3"/></svg>
                             </div>
                        </svelte:fragment> 
                            Pipeline Deployment
                    </ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="taxdb" active='variant-soft' rounded='rounded-token' disabled>
                        <svelte:fragment slot="lead">
                            <div class="w-4 h-4">
                                <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-database"><ellipse cx="12" cy="5" rx="9" ry="3"/><path d="M3 5V19A9 3 0 0 0 21 19V5"/><path d="M3 12A9 3 0 0 0 21 12"/></svg>
                            </div>
                        </svelte:fragment> 
                          Databases and Taxonomy
                    </ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="qa" active='variant-soft' rounded='rounded-token' disabled>
                        <svelte:fragment slot="lead">
                            <div class="w-4 h-4">
                                 <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-book-open-check"><path d="M8 3H2v15h7c1.7 0 3 1.3 3 3V7c0-2.2-1.8-4-4-4Z"/><path d="m16 12 2 2 4-4"/><path d="M22 6V3h-6c-2.2 0-4 1.8-4 4v14c0-1.7 1.3-3 3-3h7v-2.3"/></svg>
                            </div>
                        </svelte:fragment> 
                        Quality Assurance
                    </ListBoxItem>
                    <ListBoxItem bind:group={selectedView} name="medium" value="train" active='variant-soft' rounded='rounded-token'  on:click={() => goto(`./cerebro/training/team=${selectedTeam}`)}>
                        <svelte:fragment slot="lead">
                            <div class="w-4 h-4">
                                 <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" class="lucide lucide-book-open-check"><path d="M8 3H2v15h7c1.7 0 3 1.3 3 3V7c0-2.2-1.8-4-4-4Z"/><path d="m16 12 2 2 4-4"/><path d="M22 6V3h-6c-2.2 0-4 1.8-4 4v14c0-1.7 1.3-3 3-3h7v-2.3"/></svg>
                            </div>
                        </svelte:fragment> 
                        Training Interface
                    </ListBoxItem>
                </ListBox>
            </div>
        </section>
        <footer class="card-footer">
        </footer>
    </div>
</div>