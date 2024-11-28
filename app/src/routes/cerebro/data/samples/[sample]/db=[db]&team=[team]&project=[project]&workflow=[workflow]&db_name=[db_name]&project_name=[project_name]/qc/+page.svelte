<script lang="ts">
    import QualityControl from '$lib/components/data/sample/QualityControl.svelte';
    import { selectedWorkflowConfiguration, selectedModels } from '$lib/stores/stores';

</script>

<div class="w-full mb-5">
    <p class="mb-1"><span class="opacity-60">Quality Control</span></p>
    <div class="mb-4 border border-primary-500 p-8 rounded-md">
            <QualityControl
                selectedWorkflowConfiguration={$selectedWorkflowConfiguration}
                selectedModels={$selectedModels}
            />
    </div>
</div>

