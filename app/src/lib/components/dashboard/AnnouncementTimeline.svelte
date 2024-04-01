<script lang="ts">
	import Timeline from "$lib/general/timeline/Timeline.svelte";
	import TimelineItem from "$lib/general/timeline/TimelineItem.svelte";
</script>

<Timeline>
    <TimelineItem
        title="Cerebro v0.8.0"
        timestamp="30th March 20243"
        author="Eike Steinig"
        paragraphs={[
            "Latest Cerebro release and implementations for NATA accreditation.",
            "This release adds a nanopore implementation of the primary Cerebro short-read, meta-metagenomic taxonomic profiling pipeline including for the viral re-alignment and consensus assembly sub-workflow we use with the pan-viral and other specific viral panels from Twist Biosciences for high coverage viral genome recovery from human samples with higher range titres (e.g. clinical with 30 < Ct < 40 or rapid antigen tests).",
            "It also adds a bacterial assembly and species-typing workflow for ohenotypically challenging isolates with nanopore or hybrid sequencing using automated assembly and GTDB-informed typing."
           ]}
        badgeText="Latest"
        link="https://github.com/esteinig/cerebro/issues"
        linkText="GitHub Releases"
    />
    <TimelineItem
        title="Development Release"
        timestamp="October 26th, 2023"
        author="Eike Steinig"
        paragraphs={[
            "Welcome to the internal development version of Cerebro!", 
            "Please note that the current deployment may be slightly less performant than expected due to security consideration on our server. You can deploy the stack locally for a snappier experience.",
            "Please also note that you will likely encounter bugs and unfinished sections. If you see something egregious, you can open an issue in the repository using the bug issue template."
        ]}
        badgeText=""
        link="https://github.com/esteinig/cerebro/issues"
        linkText="GitHub Issues"
    />
</Timeline>