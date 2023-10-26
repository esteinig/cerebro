<script lang="ts">
	import Timeline from "$lib/general/timeline/Timeline.svelte";
	import TimelineItem from "$lib/general/timeline/TimelineItem.svelte";
</script>

<Timeline>
    <TimelineItem
        title="Development Release"
        timestamp="October 26th, 2023"
        author="Eike Steinig"
        paragraphs={[
            "Welcome to the development version of Cerebro!", 
            "Please note that the current deployment may be slightly less performant than expected due to security consideration on our server. You can deploy the stack locally for a snappier experience.",
            "Please also note that you will likely encounter bugs and unfinished sections. If you see something egregious, you can open an issue in the repository using the bug issue template."
        ]}
        badgeText="Latest"
        link="https://github.com/esteinig/cerebro/issues"
        linkText="GitHub Issues"
    />
</Timeline>