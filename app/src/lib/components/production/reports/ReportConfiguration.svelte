<script lang="ts">
	import { type ReportTemplateSchema } from "$lib/utils/types";

    let reportTemplateSchema: ReportTemplateSchema = {
        template_id: "",
        template_name: "",
        contact_name: "",
        contact_email: "",
        legal_disclaimer: "",
        legal_disclosure: "",
        legal_liability: "",
        signatures: []
    }

    let reportSettingsTemplate: string | null = null;
  

</script>

<form>
    <div class="p-4">

    <div id="reportConfigurationAuthorisation" class="py-8"> 
        <p class="opacity-40 text-regular">Contact details</p>
        <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-6 w-3/4 pb-12 pt-4">
        <label class="label">
            <span class="text-xs opacity-60">Name</span>
            <input class="input text-xs" bind:value={reportTemplateSchema.contact_name}  placeholder="Contact name" required={true} />
        </label>
        <label class="label">
            <span class="text-xs opacity-60">Email</span>
            <input class="input text-xs" bind:value={reportTemplateSchema.contact_email}  placeholder="Contact email" required={true} />
        </label>
        </div>
    </div>


    <div id="reportConfigurationLegal" class="py-8">

        <p class="opacity-40 text-regular pb-4">Legal statements</p>
        <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 gap-6 w-3/4 pb-12">
        <label class="label">
            <span class="text-xs opacity-60">Disclaimer</span>
            <textarea class="textarea text-xs" rows="8" bind:value={reportTemplateSchema.legal_disclaimer} placeholder="Disclaimer" required={false} />
        </label>
        <label class="label">
            <span class="text-xs opacity-60">Disclosure</span>
            <textarea class="textarea text-xs" rows="8" bind:value={reportTemplateSchema.legal_disclosure} placeholder="Disclosure statement" required={false} />
        </label>
        <label class="label">
            <span class="text-xs opacity-60">Liability</span>
            <textarea class="textarea text-xs" rows="8" bind:value={reportTemplateSchema.legal_liability} placeholder="Liability statement" required={false} />
        </label>
        </div>


    <div id="reportSettingsTemplateSelection" class="py-8"> 
        <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 gap-6 w-3/4">
        <label class="label">
            <span class="text-xs opacity-60">Template</span>
            <select class="input text-xs" bind:value={reportSettingsTemplate}>
                <option value={null}>New report template</option>
            </select>
        </label>
        </div>
    </div>
    </div>
    </div>
</form>