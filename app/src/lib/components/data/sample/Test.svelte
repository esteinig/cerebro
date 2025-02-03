<script lang="ts">

  import { onMount } from 'svelte';
	import { selectedReportSchema } from '$lib/stores/stores';
	import { type PathogenDetectionReport, type ReportTemplateSchema } from '$lib/utils/types';

  import init, { ReportCompiler } from '$lib/wasm/cerebro_report_wasm';

  let compiler: any;
  

  let report: string = "";
  let svgData: string[] = [];
  let pdfData: Uint8Array | null = null;

  let virtualPDF: string = "pdf.typ";
  let virtualSVG: string = "svg.typ";

  let displayConfigurationPage: boolean = false;
  let reportSettingsTemplate: string | null = null;
 
  enum ReportTemplate {
    PathogenDetection = "PathogenDetection"
  }



  const addSignature = () => {
    reportTemplateSchema.signatures = [
      ...reportTemplateSchema.signatures,
      { name: "", position: "", institution: "" }
    ];
  };

  const removeSignature = (index: number) => {
    reportTemplateSchema.signatures = reportTemplateSchema.signatures.filter((_, i) => i !== index);
  };

  const loadCompiler = async () => {
    try {
      await init();
      compiler = new ReportCompiler('/', ReportTemplate.PathogenDetection); // Initialize with root and empty callback
    } catch (error) {
      console.error('Failed to load WASM module:', error);
    }

    if (!report) {
      compileSVG(reportSchema); // Compile on page initilization
    }
  };

  const compilePDF = async () => {

    if (!compiler) {
      console.error('Compiler not loaded yet!');
      return;
    }

    try {
      const result = await compiler.pdf(report, virtualPDF);
      pdfData = new Uint8Array(result);
      downloadPDF();
    } catch (error) {
      console.error('Failed to compile PDF:', error);
    }

  };

  const compileSVG = async (reportSchema: PathogenDetectionReport) => {

    if (!compiler) {
      console.error('Compiler not loaded yet!');
      return;
    }

    try {
      report = compiler.report(reportSchema);
      svgData = await compiler.svg(report, virtualSVG); 
    } catch (error) {
      console.error('Failed to compile SVG:', error);
    }
  };


  const downloadPDF = () => {
    if (pdfData) {
      
      const blob = new Blob([pdfData], { type: 'application/pdf' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');

      link.href = url;
      link.download = 'output.pdf';
      link.click();
      URL.revokeObjectURL(url);
    }
  };

  const handleSaveAction = () => {
    compileSVG(reportSchema);
  };

  const handleKeydown = (event: KeyboardEvent) => {
    if ((event.ctrlKey || event.metaKey) && event.key === 's') {
      event.preventDefault(); // Prevent the default browser "Save Page" action
      handleSaveAction();
    }
  };

  onMount(() => {
    loadCompiler();
    document.addEventListener('keydown', handleKeydown);
    return () => {
      document.removeEventListener('keydown', handleKeydown);
    };
  });
  
  let reportSchema: PathogenDetectionReport = $selectedReportSchema;
  
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

  $: {
    if (compiler) {
      compileSVG(reportSchema);
    }
  }
  
</script>

<div> 
    <div class="grid grid-cols-2">
        <!-- Data input and compiler -->
        <div>
          {#if displayConfigurationPage}

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

              <div id="reportConfigurationAuthorisation" class="py-8"> 
                  
                <div class="grid grid-cols-4 gap-6 w-full pb-8">
                  <p class="opacity-40 text-regular col-span-3">Authorisation signatures</p>
                  <div class="justify-end">

                    <button class="btn-icon btn-icon-sm variant-ghost-primary" on:click={() => addSignature()}>
                      <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-4 h-4">
                        <path stroke-linecap="round" stroke-linejoin="round" d="M12 4.5v15m7.5-7.5h-15" />
                      </svg>
                    </button>
                  </div>
                </div>
                {#if reportTemplateSchema.signatures.length > 0}
                  {#each reportTemplateSchema.signatures as signature, index}
                    <div class="grid grid-cols-4 gap-6 w-full py-2">
                      <label class="label">
                        <span class="text-xs opacity-60">Name</span>
                          <input class="input text-xs" bind:value={signature.name}  placeholder="Name" required={false} />
                      </label>
                      <label class="label">
                        <span class="text-xs opacity-60">Position</span>
                          <input class="input text-xs" bind:value={signature.position}  placeholder="Position" required={false} />
                      </label>
                      <label class="label">
                        <span class="text-xs opacity-60">Institution</span>
                          <input class="input text-xs" bind:value={signature.institution}  placeholder="Institution" required={false} />
                      </label>
                      <button class="btn-icon btn-icon-sm variant-ghost-primary mt-8" on:click={() => removeSignature(index)}>
                        <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-4 h-4">
                          <path stroke-linecap="round" stroke-linejoin="round" d="M5 12h14" />
                        </svg>              
                      </button>
                    </div>
                  {/each}
                {:else}
                    
                    <p class="opacity-60 text-sm">No signatures have been added</p>
                {/if}
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

                <div class="grid grid-cols-3 gap-6 w-3/4 pb-12 pt-12"> 
                  <button class="btn variant-ghost-primary" on:click={() => displayConfigurationPage ? displayConfigurationPage = false : displayConfigurationPage = true }>
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="h-5 w-5 mr-2">
                      <path stroke-linecap="round" stroke-linejoin="round" d="M6.75 15.75 3 12m0 0 3.75-3.75M3 12h18" />
                    </svg>                
                    Back
                  </button> 
                  <button class="btn variant-ghost-primary" on:click={() => displayConfigurationPage ? displayConfigurationPage = false : displayConfigurationPage = true }>
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="h-5 w-5 mr-2">
                      <path stroke-linecap="round" stroke-linejoin="round" d="M7.5 7.5h-.75A2.25 2.25 0 0 0 4.5 9.75v7.5a2.25 2.25 0 0 0 2.25 2.25h7.5a2.25 2.25 0 0 0 2.25-2.25v-7.5a2.25 2.25 0 0 0-2.25-2.25h-.75m-6 3.75 3 3m0 0 3-3m-3 3V1.5m6 9h.75a2.25 2.25 0 0 1 2.25 2.25v7.5a2.25 2.25 0 0 1-2.25 2.25h-7.5a2.25 2.25 0 0 1-2.25-2.25v-.75" />
                    </svg>
                                 
                    Load
                  </button> 
                  <button class="btn variant-ghost-primary" on:click={() => displayConfigurationPage ? displayConfigurationPage = false : displayConfigurationPage = true }>
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="h-5 w-5 mr-2">
                      <path stroke-linecap="round" stroke-linejoin="round" d="M7.5 7.5h-.75A2.25 2.25 0 0 0 4.5 9.75v7.5a2.25 2.25 0 0 0 2.25 2.25h7.5a2.25 2.25 0 0 0 2.25-2.25v-7.5a2.25 2.25 0 0 0-2.25-2.25h-.75m0-3-3-3m0 0-3 3m3-3v11.25m6-2.25h.75a2.25 2.25 0 0 1 2.25 2.25v7.5a2.25 2.25 0 0 1-2.25 2.25h-7.5a2.25 2.25 0 0 1-2.25-2.25v-.75" />
                    </svg>
                                 
                    Save
                  </button> 
                </div>
              </div>

            </div>

          </form>
          {:else}
            <form>
              <div class="p-4">
                <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-3 gap-12 w-5/6 pb-16">
                    <label class="label">
                        <span class="text-xs opacity-60">Laboratory</span>
                        <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.reporting_laboratory} required={true} />
                    </label>
                    <label class="label">
                        <span class="text-xs opacity-60">Date</span>
                        <input type="text" class="input text-xs"  bind:value={reportSchema.patient_header.reporting_date} required={true} />
                    </label>
                    <label class="label">
                        <span class="text-xs opacity-60">Location</span>
                        <input type="text" class="input text-xs" bind:value={reportSchema.footer.reporting_location} required={true} />
                    </label>
                </div>

                <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-3 gap-12 w-5/6 pb-8">
                  <label class="label">
                      <span class="text-xs opacity-60">Patient name</span>
                      <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.patient_name} required={true} />
                  </label>
                  <label class="label">
                      <span class="text-xs opacity-60">Patient DOB</span>
                      <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.patient_dob} required={true} />
                  </label>
                  <label class="label">
                      <span class="text-xs opacity-60">Patient URN</span>
                      <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.patient_urn} required={true} />
                  </label>
                </div>
                <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-3 gap-12 w-5/6 pb-16">
                  <label class="label">
                      <span class="text-xs opacity-60">Hospital site</span>
                      <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.hospital_site} required={true} />
                  </label>

                  <label class="label">
                      <span class="text-xs opacity-60">Laboratory number</span>
                      <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.laboratory_number} required={true} />
                  </label>
                  <label class="label">
                      <span class="text-xs opacity-60">Requested Doctor</span>
                      <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.requested_doctor} required={true} />
                  </label>
                </div>

                <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-12 w-3/4 pb-8">
                  <label class="label">
                      <span class="text-xs opacity-60">Specimen ID</span>
                      <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.specimen_id} required={true} />
                  </label>
                  <label class="label">
                      <span class="text-xs opacity-60">Specimen type</span>
                      <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.specimen_type} required={true} />
                  </label>
                </div>

                <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-12 w-3/4 pb-16">
                  <label class="label">
                      <span class="text-xs opacity-60">Date collected</span>
                  
                      <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.date_collected} required={true} />
                  </label>
                  <label class="label">
                      <span class="text-xs opacity-60">Date received</span>
                      <input type="text" class="input text-xs" bind:value={reportSchema.patient_header.date_received} required={true} />
                  </label>
                </div>

                <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 gap-6 w-3/4 pb-12">
                  <label class="label">
                      <span class="text-xs opacity-60">Comments</span>
                      <textarea class="textarea text-xs" bind:value={reportSchema.patient_result.comments} placeholder="Comments on results" required={true} />
                  </label>
                  <label class="label">
                    <span class="text-xs opacity-60">Orthogonal Tests</span>
                    <textarea class="textarea text-xs" bind:value={reportSchema.patient_result.orthogonal_tests} placeholder="Orthogonal diagnostic test results" required={true} />
                  </label>
                  <label class="label">
                    <span class="text-xs opacity-60">Clinical Notes</span>
                    <textarea class="textarea text-xs" bind:value={reportSchema.patient_result.clinical_notes} placeholder="Clinical notes supporting results" required={true} />
                  </label>
                  <label class="label">
                    <span class="text-xs opacity-60">Actions</span>
                    <textarea class="textarea text-xs" bind:value={reportSchema.patient_result.actions} placeholder="Actions taken based on results" required={true} />
                  </label>
                </div>
                

                <div class="grid grid-cols-3 gap-6 w-3/4 pb-12 pt-12"> 
                  <button class="btn variant-ghost-primary" on:click={compilePDF}>
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-5 h-5 mr-2">
                      <path stroke-linecap="round" stroke-linejoin="round" d="M3 16.5v2.25A2.25 2.25 0 0 0 5.25 21h13.5A2.25 2.25 0 0 0 21 18.75V16.5M16.5 12 12 16.5m0 0L7.5 12m4.5 4.5V3" />
                    </svg>
                    Download
                  </button>
                  <button class="btn variant-ghost-primary" on:click={() => displayConfigurationPage ? displayConfigurationPage = false : displayConfigurationPage = true }>
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-5 h-5 mr-2">
                      <path stroke-linecap="round" stroke-linejoin="round" d="M10.343 3.94c.09-.542.56-.94 1.11-.94h1.093c.55 0 1.02.398 1.11.94l.149.894c.07.424.384.764.78.93.398.164.855.142 1.205-.108l.737-.527a1.125 1.125 0 0 1 1.45.12l.773.774c.39.389.44 1.002.12 1.45l-.527.737c-.25.35-.272.806-.107 1.204.165.397.505.71.93.78l.893.15c.543.09.94.559.94 1.109v1.094c0 .55-.397 1.02-.94 1.11l-.894.149c-.424.07-.764.383-.929.78-.165.398-.143.854.107 1.204l.527.738c.32.447.269 1.06-.12 1.45l-.774.773a1.125 1.125 0 0 1-1.449.12l-.738-.527c-.35-.25-.806-.272-1.203-.107-.398.165-.71.505-.781.929l-.149.894c-.09.542-.56.94-1.11.94h-1.094c-.55 0-1.019-.398-1.11-.94l-.148-.894c-.071-.424-.384-.764-.781-.93-.398-.164-.854-.142-1.204.108l-.738.527c-.447.32-1.06.269-1.45-.12l-.773-.774a1.125 1.125 0 0 1-.12-1.45l.527-.737c.25-.35.272-.806.108-1.204-.165-.397-.506-.71-.93-.78l-.894-.15c-.542-.09-.94-.56-.94-1.109v-1.094c0-.55.398-1.02.94-1.11l.894-.149c.424-.07.765-.383.93-.78.165-.398.143-.854-.108-1.204l-.526-.738a1.125 1.125 0 0 1 .12-1.45l.773-.773a1.125 1.125 0 0 1 1.45-.12l.737.527c.35.25.807.272 1.204.107.397-.165.71-.505.78-.929l.15-.894Z" />
                      <path stroke-linecap="round" stroke-linejoin="round" d="M15 12a3 3 0 1 1-6 0 3 3 0 0 1 6 0Z" />
                    </svg>
                    Configure
                  </button> 
                </div>
                <p class="opacity-40 text-regular mt-4"><kbd class="kbd text-sm">Ctrl + S</kbd> or <kbd class="kbd">Command + S</kbd> (MacOS) to update preview</p>
                
              </div>
            </form>
          {/if}
          
        </div>
        
        <!-- Live preview -->

        <div class="grid grid-rows-12 p-1 gap-y-2 bg-surface-500/5">
          {#await Promise.resolve(svgData)}
            <p>Loading preview...</p>
          {:then data}
            {#each data as svg}
              <div class="svg-preview">
                {@html svg}
              </div>
            {/each}
          {:catch error}
            <p>Error loading preview: {error.message}</p>
          {/await}
        </div>
    </div>
</div>

<style>

    /* Container for the SVG preview */
    .svg-preview {
      display: flex;
      justify-content: center;
      align-items: center;
      width: 100%;
      height: calc(100vh - 18vh); /* Adjust height based on layout */
      overflow: auto;
    }
</style>