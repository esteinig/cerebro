<script lang="ts">
    import { page } from "$app/stores";
    import PathogenTrainingDeck from "$lib/components/training/PathogenTrainingDeck.svelte";
    
    import type { DeckConfig, SlideConfig } from "$lib/components/training/types";
	import { AbundanceMode, DisplayData, DisplayVisualisation } from "$lib/utils/types";


    let testDeck = {
        id: "test",
        name: "Deck Testing",
        slides: [
            {
                id: "test1",
                title: "Test Slide",
                description: "This is a test slide deck for your convenience",
                displayData: DisplayData.Rpm,
                displayMode: AbundanceMode.Mixed,
                selectedVisualisation: DisplayVisualisation.Threshold,
            } as SlideConfig
        ]
    } as DeckConfig;

</script>

<PathogenTrainingDeck deck={testDeck} userId={$page.data.userData.id}></PathogenTrainingDeck>