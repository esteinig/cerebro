<script lang="ts">
    import { page } from "$app/stores";
    import { Role } from "$lib/utils/types";
	import { getInitials } from "$lib/utils/helpers";
    import { ListBox, ListBoxItem } from "@skeletonlabs/skeleton";
	import { Accordion, AccordionItem, Avatar } from "@skeletonlabs/skeleton";

    let selectedView: string;
    let initials = getInitials($page.data.userData.name);

    $: console.log(selectedView)


</script>

<p class="opacity-60">User profile</p>
<div class="card border border-primary-500 h-full w-full">
    <header class="card-header">
        <div class="p-4">
            <div class="flex">
                <div class="flex items-center">
                    <Avatar src="" width="w-16" rounded="rounded-full" initials={initials} />
                    <div class="flex-none px-5">
                        <div class="md:h4 lg:h3">{$page.data.userData.title ?? ""} {$page.data.userData.name}</div>
                        <div class="flex items-center">
                            <h4 class="h4 opacity-60 pr-1">{$page.data.userData.positions.join(", ")}</h4>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </header>
    <section>
        <div class="p-4">
            <div class="px-8">
                <p class="text-sm opacity-60 pb-4">Last login </p>
                <p class="text-sm px-4 pb-4">{new Date().toDateString()} {new Date().toLocaleTimeString()}</p>
                
                <p class="text-sm opacity-60 py-4">Permissions</p>
                
                <Accordion>
                    {#each $page.data.userData.roles.sort() as role}
                        <AccordionItem>
                            <svelte:fragment slot="lead">                       
                                {#if role === Role.Admin}
                                    <div class="w-5">
                                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                            <path d="M12 21a9.004 9.004 0 008.716-6.747M12 21a9.004 9.004 0 01-8.716-6.747M12 21c2.485 0 4.5-4.03 4.5-9S14.485 3 12 3m0 18c-2.485 0-4.5-4.03-4.5-9S9.515 3 12 3m0 0a8.997 8.997 0 017.843 4.582M12 3a8.997 8.997 0 00-7.843 4.582m15.686 0A11.953 11.953 0 0112 10.5c-2.998 0-5.74-1.1-7.843-2.918m15.686 0A8.959 8.959 0 0121 12c0 .778-.099 1.533-.284 2.253m0 0A17.919 17.919 0 0112 16.5c-3.162 0-6.133-.815-8.716-2.247m0 0A9.015 9.015 0 013 12c0-1.605.42-3.113 1.157-4.418" stroke-linecap="round" stroke-linejoin="round"></path>
                                        </svg>
                                    </div>
                                {:else if role === Role.Bot}
                                    <div class="w-5">
                                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                            <path d="M6.75 7.5l3 2.25-3 2.25m4.5 0h3m-9 8.25h13.5A2.25 2.25 0 0021 18V6a2.25 2.25 0 00-2.25-2.25H5.25A2.25 2.25 0 003 6v12a2.25 2.25 0 002.25 2.25z" stroke-linecap="round" stroke-linejoin="round"></path>
                                        </svg>
                                    </div>
                                {:else if role === Role.Data}
                                    <div class="w-5">
                                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                            <path d="M20.25 6.375c0 2.278-3.694 4.125-8.25 4.125S3.75 8.653 3.75 6.375m16.5 0c0-2.278-3.694-4.125-8.25-4.125S3.75 4.097 3.75 6.375m16.5 0v11.25c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125V6.375m16.5 0v3.75m-16.5-3.75v3.75m16.5 0v3.75C20.25 16.153 16.556 18 12 18s-8.25-1.847-8.25-4.125v-3.75m16.5 0c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125" stroke-linecap="round" stroke-linejoin="round"></path>
                                        </svg>
                                    </div>
                                {:else if role === Role.Report}
                                    <div class="w-5">
                                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                            <path d="M19.5 14.25v-2.625a3.375 3.375 0 00-3.375-3.375h-1.5A1.125 1.125 0 0113.5 7.125v-1.5a3.375 3.375 0 00-3.375-3.375H8.25m0 12.75h7.5m-7.5 3H12M10.5 2.25H5.625c-.621 0-1.125.504-1.125 1.125v17.25c0 .621.504 1.125 1.125 1.125h12.75c.621 0 1.125-.504 1.125-1.125V11.25a9 9 0 00-9-9z" stroke-linecap="round" stroke-linejoin="round"></path>
                                        </svg>
                                    </div>
                                {:else if role === Role.User}
                                    <div class="w-5">
                                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                            <path d="M15.75 6a3.75 3.75 0 11-7.5 0 3.75 3.75 0 017.5 0zM4.501 20.118a7.5 7.5 0 0114.998 0A17.933 17.933 0 0112 21.75c-2.676 0-5.216-.584-7.499-1.632z" stroke-linecap="round" stroke-linejoin="round"></path>
                                        </svg>
                                    </div>
                                {/if}
                            </svelte:fragment>
                            <svelte:fragment slot="summary">{role}</svelte:fragment>
                            <svelte:fragment slot="content">(content)</svelte:fragment>
                        </AccordionItem>
                    {/each}
                </Accordion>


                <p class="text-sm opacity-60 py-4">Teams</p>

                <Accordion>
                    {#each $page.data.userTeams.sort() as team}
                        <AccordionItem>
                            <svelte:fragment slot="lead"></svelte:fragment>
                            <svelte:fragment slot="summary">{team.name}</svelte:fragment>
                            <svelte:fragment slot="content">
                            <ListBox>
                                {#each team.databases as db}
                                    <ListBoxItem bind:group={selectedView} name="medium" value="{db.name}" active='variant-soft' rounded='rounded-token'>
                                        <svelte:fragment slot="lead">
                                            <div class="w-5">
                                            <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                                <path d="M20.25 6.375c0 2.278-3.694 4.125-8.25 4.125S3.75 8.653 3.75 6.375m16.5 0c0-2.278-3.694-4.125-8.25-4.125S3.75 4.097 3.75 6.375m16.5 0v11.25c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125V6.375m16.5 0v3.75m-16.5-3.75v3.75m16.5 0v3.75C20.25 16.153 16.556 18 12 18s-8.25-1.847-8.25-4.125v-3.75m16.5 0c0 2.278-3.694 4.125-8.25 4.125s-8.25-1.847-8.25-4.125" stroke-linecap="round" stroke-linejoin="round"></path>
                                            </svg>
                                        </div>
                                        </svelte:fragment>
                                        <div class="opacity-60 pl-2 cursor-pointer">{db.name}</div>
                                    </ListBoxItem>
                                {/each} 
                            </ListBox>
                            </svelte:fragment>
                        </AccordionItem>
                    {/each}
                </Accordion>
            </div>
        </div>
        

    </section>
    <footer class="card-footer">

        
    </footer>
</div>