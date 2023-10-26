<script lang="ts">
	import { page } from "$app/stores";
	import TaxonViewSelection from "$lib/components/data/sample/selections/TaxonViewSelection.svelte";

    let selectedView: string = "";

    let prompt: string = `Provide a brief summary of the clinical significance of E.coli (NCBI taxid: 1111)
     in human patients. Finish in one sentence if this organism is a common lab- or sample-contaminant.`;

    const queryChatGpt = async() => {
        
        const res = await fetch(`https://api.openai.com/v1/chat/completions`, {
            method: 'POST',
            body: JSON.stringify({
                "model": "gpt-3.5-turbo",
                "messages": [{
                    "role": "user",
                    "content":  prompt
                    }]
            }),
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer `
            }
        });

        if (res.ok) {
            return await res.json();
        } else {
            throw new Error(`${res.status} - ${res.statusText}`);
        }
    } 

</script>

<div class="flex justify-center">
    <div>
        <ol class="breadcrumb justify-start">
            <li class="crumb opacity-60">{$page.params.db_name}</li>
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            <li class="crumb opacity-60">{$page.params.project_name}</li>
            <li class="crumb-separator" aria-hidden>&rsaquo;</li>
            <li class="crumb">{$page.params.sample}</li>
        </ol>
    </div>
</div>
<div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-8 gap-16 pt-10">
    <div class="col-span-3">
        <div class="w-full">
            <div class="mb-4 p-4 border border-primary-500 rounded-md">

                <ol class="breadcrumb justify-start p-4">
                    <li class="crumb opacity-60">Domain</li>
                    <li class="crumb-separator" aria-hidden>&rsaquo;</li>
                    <li class="crumb opacity-60">Genus</li>
                    <li class="crumb-separator" aria-hidden>&rsaquo;</li>
                    <li class="crumb">Species</li>
                </ol>
                <TaxonViewSelection bind:selectedView={selectedView}/>
            </div>
        </div>
    </div>
</div>