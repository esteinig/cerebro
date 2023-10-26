<script lang="ts">
    import { page } from "$app/stores";
    import { Role } from "$lib/utils/types";
	import { getInitials } from "$lib/utils/helpers";
	import { Accordion, AccordionItem, Avatar } from "@skeletonlabs/skeleton";
	import AnnouncementTimeline from "./AnnouncementTimeline.svelte";

    let initials = getInitials($page.data.userData.name);

</script>


<p class="opacity-60">News</p>
<div class="card w-full border border-primary-500 overflow-y-scroll md:h-[32rem]">
    <header class="card-header">
        
    </header>
    <section>
        <div class="p-4">
            <div class="px-8 pb-4">
                <AnnouncementTimeline></AnnouncementTimeline>
            </div>
            
        </div>
        

    </section>
    <footer class="card-footer">

        
    </footer>
</div>