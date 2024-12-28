<script lang="ts">

  import { onMount } from 'svelte';
	import { selectedReportSchema } from '$lib/stores/stores';
  import init, { ReportCompiler } from '$lib/wasm/cerebro_report_wasm';

  let compiler: any;
  
  let pdfData: Uint8Array | null = null; // Stores the compiled PDF
  let svgData: string[] = []; // Stores the compiled SVG

  let virtualPdf: string = "pdf.typ";
  let virtualSvg: string = "svg.typ";

  let report: string = "";

  enum ReportTemplate {
    PathogenDetection = "PathogenDetection"
  }

  const loadCompiler = async () => {
    try {
      await init();
      compiler = new ReportCompiler('/', () => null, ReportTemplate.PathogenDetection); // Initialize with root and empty callback
    } catch (error) {
      console.error('Failed to load WASM module:', error);
    }
  };


  const compilePDF = async () => {

    if (!compiler) {
      console.error('Compiler not loaded yet!');
      return;
    }

    try {
      const result = await compiler.pdf(report, virtualPdf);
      pdfData = new Uint8Array(result);
    } catch (error) {
      console.error('Failed to compile PDF:', error);
    }

  };


  const compileSVG = async () => {
        
    if (!compiler) {
      console.error('Compiler not loaded yet!');
      return;
    }

    try {
      report = compiler.report(reportSchema);
    } catch (error) {
      console.error('Failed to template report:', error);
    }

    try {
      const svgResult = await compiler.svg(report, virtualSvg);
      svgData = svgResult;
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

  onMount(async () => {
      loadCompiler();
  });
  
  let reportSchema = $selectedReportSchema;

  $: console.log(reportSchema);

</script>

<div> 
    <div class="grid grid-cols-2">
        <!-- Data input and compiler -->
        <div>
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

              <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-6 w-3/4 pb-12">
                <label class="label">
                  <span class="text-xs opacity-60">Contact</span>
                  <input class="input text-xs" bind:value={reportSchema.patient_result.contact_name}  placeholder="Contact name" required={true} />
                </label>
                <label class="label">
                  <span class="text-xs opacity-60">Email</span>
                  <input class="input text-xs" bind:value={reportSchema.patient_result.contact_email}  placeholder="Contact email" required={true} />
                </label>
              </div>

              <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 gap-6 w-3/4 pb-12">
                <label class="label">
                    <span class="text-xs opacity-60">Comments</span>
                    <textarea class="textarea text-xs" bind:value={reportSchema.patient_result.comments} placeholder="Comments on results" required={true} />
                </label>
                <label class="label">
                  <span class="text-xs opacity-60">Actions</span>
                  <textarea class="textarea text-xs" bind:value={reportSchema.patient_result.comments} placeholder="Actions taken based on results" required={true} />
              </label>
              </div>

              <button class="btn variant-ghost-primary mt-4" on:click={compilePDF}>Compile report</button>
              {#if pdfData}
                  <button  class="btn variant-ghost-primary mt-4"on:click={downloadPDF}>Download PDF</button>
              {/if}
  
              <button class="btn variant-ghost-primary mt-4" on:click={compileSVG}>Update preview</button></div>

          </form>

        </div>
        
        <!-- Live preview -->

        <div class="grid grid-rows-12 p-1 gap-y-2 bg-surface-500/5">
            {#each svgData as svg}
                <div class="svg-preview">
                    {@html svg}
                </div>
            {/each}
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
  
    /* Style applied to dynamically inserted SVGs */
    .svg-preview svg {
        max-width: 100%;
        max-height: 100%;
        width: auto;
        height: auto;
    }
</style>