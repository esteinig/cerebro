<script lang="ts">

    import { onMount } from 'svelte';
    import { selectedReportSchema } from '$lib/stores/stores';
    import { type PathogenDetectionReport } from '$lib/utils/types';

    import init, { ReportCompiler } from '$lib/wasm/cerebro_report_wasm';

	import ReportConfiguration from './ReportConfiguration.svelte';
	import ReportAppendixLaboratory from './ReportAppendixLaboratory.svelte';
	import ReportAppendixBioinformatics from './ReportAppendixBioinformatics.svelte';
	import { SlideToggle } from '@skeletonlabs/skeleton';
  
    export let selectedView: string = "report";

    let compiler: any;

    let report: string = "";
    let svgData: string[] = [];
    let pdfData: Uint8Array | null = null;
  
    let virtualPDF: string = "pdf.typ";
    let virtualSVG: string = "svg.typ";
    

    enum ReportTemplate {
      PathogenDetection = "PathogenDetection"
    }

    let reportSchema: PathogenDetectionReport = $selectedReportSchema;
  
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

  
    $: if (compiler) compileSVG(reportSchema);


    const addSignature = () => {
        reportSchema.authorisation.signatures = [
        ...reportSchema.authorisation.signatures,
        { name: "", position: "", institution: "" }
      ];
    };
  
    const removeSignature = (index: number) => {
        reportSchema.authorisation.signatures = reportSchema.authorisation.signatures.filter((_, i) => i !== index);
    };
    
  </script>
  
  <div> 
    <div class="grid grid-cols-2">
          <!-- Data input and compiler -->
          <div>
            {#if selectedView === "report"}
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

                    <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 gap-6 w-5/6 pb-6">
                        <div class="flex items-center pt-6">
                        <SlideToggle name="sliderPathogenDetected" bind:checked={reportSchema.patient_result.pathogen_detected} active="bg-primary-500" size="sm"><span class="opacity-60 text-sm">Pathogen detected</span></SlideToggle>
                        </div>
                    </div>
                    <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-6 w-3/4 pb-12">

                        <label class="label">
                        <span class="text-xs opacity-60">Pathogen</span>
                        <input type="text" class="input text-xs" bind:value={reportSchema.patient_result.pathogen_reported} required={true} disabled={!reportSchema.patient_result.pathogen_detected}/>
                        </label>
                        <label class="label">
                        <span class="text-xs opacity-60">Review date</span>
                        <input type="text" class="input text-xs" bind:value={reportSchema.patient_result.review_date} required={true}/>
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
                        {#if reportSchema.authorisation.signatures.length > 0}
                        {#each reportSchema.authorisation.signatures as signature, index}
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

                    </div>
                </form>
            {:else if selectedView === "configuration"}
                <ReportConfiguration />
            {:else if selectedView === "appendixA"}
                <ReportAppendixLaboratory />
            {:else if selectedView === "appendixB"}
                <ReportAppendixBioinformatics />
            {/if}
            <div class="grid grid-cols-3 gap-6 w-3/4 pb-12 pt-12"> 
                <button class="btn variant-ghost-primary" on:click={compilePDF}>
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-5 h-5 mr-2">
                    <path stroke-linecap="round" stroke-linejoin="round" d="M3 16.5v2.25A2.25 2.25 0 0 0 5.25 21h13.5A2.25 2.25 0 0 0 21 18.75V16.5M16.5 12 12 16.5m0 0L7.5 12m4.5 4.5V3" />
                    </svg>
                    Download
                </button>
            </div>
            <p class="opacity-40 text-regular mt-4"><kbd class="kbd text-sm">Ctrl + S</kbd> or <kbd class="kbd">Command + S</kbd> (MacOS) to update preview</p>
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
        height: 100%; /* Adjust height based on layout */
        overflow: auto;
      }
  </style>