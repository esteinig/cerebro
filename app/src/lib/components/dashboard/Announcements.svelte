<script lang="ts">
    import { page } from "$app/stores";
	import { getInitials } from "$lib/utils/helpers";
	import AnnouncementTimeline from "./AnnouncementTimeline.svelte";

</script>


<p class="opacity-60">Release notes</p>
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