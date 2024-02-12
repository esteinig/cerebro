<script lang="ts">
	import ActivityItem from "$lib/general/timeline/ActivityItem.svelte";
	import Timeline from "$lib/general/timeline/Timeline.svelte";
</script>

<Timeline>
    <ActivityItem
        user="Chhay Lim"
        activity="created report for sample"
        target="DW-63-76"
        timestamp="October 26th, 2023"      
        description="Clinical report was created and stored in database."  
        link="Test"
    />
    <ActivityItem
        user="Prashanth Ramachandran"
        activity="accepted candidate for sample"
        target="DW-63-76"
        timestamp="October 26th, 2023"
        description="<blockquote class='blockquote text-sm'>Evidence looks reasonable, please report this back to treating clinician.<blockquote/>"
        link="Test"
        targetPath={["VIDRL", "CNS Metagenomics", "Phage Validation", "DW-63-76"]}
    />
    <ActivityItem
        user="Eike Steinig"
        activity="commented on candidate for sample"
        target="DW-63-76"
        timestamp="October 26th, 2023"
        description="<blockquote class='blockquote text-sm'>Evidence clearly distinguishes abundance from usual background contamination. 
        DNA and RNA indicate elevated reads in k-mer and alignment profiling, with some assembled fragments typed at species-level with NCBI NT and NCBI NR.<blockquote/>"
        link="Test"
    />
    <ActivityItem
        user="Eike Steinig"
        activity="selected candidate for sample"
        target="DW-63-76"
        timestamp="October 26th, 2023"
        description="<i>Staphylococcus aureus</i> was selected as candidate pathogen."
        link="Test"
        targetPath={["VIDRL", "CNS Metagenomics", "Phage Validation", "DW-63-76"]}
    />
</Timeline>