<script lang="ts">

    import { onMount } from 'svelte';

    import init, { ReportCompiler } from '$lib/wasm/cerebro_report_wasm';


    let compiler: any;
    let pdfData: Uint8Array | null = null;
    let typstText = `

#set text(font: "IBM Plex Sans")

#let first_page_header = {
  show table.cell.where(x: 0): set text(weight: 500)
  show table.cell.where(x: 2): set text(weight: 500)
  table(
    columns: 4,
    [Patient Name], [Patient X], [Specimen ID], [DW-63-103],
    [URN Number], [3421341241], [Date Collected], [01/01/2000],
    [Date of Birth], [01/01/2000], [Date Received], [01/01/2000],
    [Requested Doctor], [Dr. Y], [Specimen Type], [CSF],
    [Hospital Site], [St. Mangione HHS], [Reporting Laboratory], [DoJ],
    [Laboratory Number], [3421341241], [Reporting Date], [01/01/2000],
  )
}

#first_page_header

#outlinebox(title: "Test Result")[#first_page_results]

    
    `;

    const loadCompiler = async () => {
        try {
            await init();
            compiler = new ReportCompiler('/', () => null); // Initialize with root and empty callback
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
            // Compile the Typst text into a PDF
            const result = await compiler.pdf(typstText, 'example.typ');
            pdfData = new Uint8Array(result);
            console.log('PDF compiled successfully!', pdfData);
        } catch (error) {
            console.error('Failed to compile PDF:', error);
        }
    };

    const addFont = async (fontPath: string) => {
        try {
        // Fetch the font file as an ArrayBuffer
        const response = await fetch(fontPath);
        const fontData = new Uint8Array(await response.arrayBuffer());
        compiler.add_font(fontData); // Add the font to the compiler
        console.log(`Font added: ${fontPath}`);
        } catch (error) {
        console.error(`Failed to add font: ${fontPath}`, error);
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

</script>

<textarea class="textarea" bind:value={typstText} rows="10" cols="50"></textarea>
<button class="btn variant-ghost-primary mt-4" on:click={compilePDF}>Compile PDF</button>
{#if pdfData}
  <button  class="btn variant-ghost-primary mt-4"on:click={downloadPDF}>Download PDF</button>
{/if}

