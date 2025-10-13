<script lang="ts">
  
    import { get } from 'svelte/store';
    import { FileDropzone, ProgressRadial } from "@skeletonlabs/skeleton";
    import { selectedReportSchema } from '$lib/stores/stores';
    import type { PathogenDetectionReport, ReportSignature } from '$lib/utils/types';
	  import { page } from '$app/stores';
	  import JSZip from 'jszip';

    export let compiler: any;
    export let reportSchema: PathogenDetectionReport = $selectedReportSchema;

    export let batch: PathogenDetectionReport[] = [];
    let processing: boolean = false;

    let busy = false;
    let progress = 0;

    async function compileOne(r: PathogenDetectionReport): Promise<Uint8Array> {
      const typ = compiler.report(r);
      const pdf = await compiler.pdf(typ, 'pdf.typ');
      return new Uint8Array(pdf);
    }


    export function makeReportFilename(userName: string, i: number, now: Date = new Date()): string {
        const pad = (n: number) => String(n).padStart(2, "0");
        const y = now.getFullYear();
        const m = pad(now.getMonth() + 1);
        const d = pad(now.getDate());
        const clean = userName.replace(/\s+/g, "");
        return `${y}${m}${d}_${clean}_ClinicalReport_${i}.pdf`;
    }

    async function downloadBatchZip() {
      
      if (!compiler || !batch?.length) return;

      busy = true;
      progress = 0;

      const zip = new JSZip();
      const total = batch.length;

      for (let i = 0; i < total; i++) {
        try {
          const pdf = await compileOne(batch[i]);
          const name = makeReportFilename($page.data.userData.name, i+1);
          zip.file(name, pdf);
        } catch (e) {
          console.error('Error compiling', e);
        }
        progress = Math.round(((i + 1) / total) * 100);
      }

      const blob = await zip.generateAsync({ type: 'blob' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'batch_reports.zip';
      a.click();
      URL.revokeObjectURL(url);

      busy = false;
    }
  
    const FIELD_MAP: Record<string, string> = {
      "Laboratory": "patient_header.reporting_laboratory",
      "Date": "patient_header.reporting_date",
      "Location": "footer.reporting_location",
      "Patient name": "patient_header.patient_name",
      "Patient DOB": "patient_header.patient_dob",
      "Patient URN": "patient_header.patient_urn",
      "Hospital site": "patient_header.hospital_site",
      "Laboratory number": "patient_header.laboratory_number",
      "Requested Doctor": "patient_header.requested_doctor",
      "Specimen ID": "patient_header.specimen_id",
      "Specimen type": "patient_header.specimen_type",
      "Date collected": "patient_header.date_collected",
      "Date received": "patient_header.date_received",
      "Pathogen detected": "patient_result.pathogen_detected",
      "Pathogen": "patient_result.pathogen_reported",
      "Review date": "patient_result.review_date",
      "Comments": "patient_result.comments",
      "Orthogonal Tests": "patient_result.orthogonal_tests",
      "Clinical Notes": "patient_result.clinical_notes",
      "Actions": "patient_result.actions"
    };
  
    function onChangeHandler(e: CustomEvent | Event) {
      const file = (e as any)?.detail?.files?.[0] ?? (e.target as HTMLInputElement)?.files?.[0];
      if (!file) return;

      processing = true;

      file.text().then((text: any) => {

        const rows = parseCSV(text);
        if (!rows.length) return;

        const [header, ...data] = rows;
        const idx = indexMap(header);

        batch = data.filter(r => r && r.length).map(r => {
          let base = structuredClone(get(selectedReportSchema)) as PathogenDetectionReport;
          console.log(base);
          // Set configured signatures from new report interface view
          base.authorisation.signatures = reportSchema.authorisation.signatures;
          applyToSchema(base, rowToObj(r, idx));
          return base;
        });

      }).catch((e: any) => {
        
        console.error(e);

      }).finally(() => processing = false);
    }
  
    // ---- CSV template download ----
    function downloadTemplate(includeExampleRow = false) {
      const header = Object.keys(FIELD_MAP);
      const rows: string[][] = [header];
  
      if (includeExampleRow) {
        // blank example row; edit to inject defaults if you like
        rows.push(header.map(() => ""));
      }
  
      const csv = toCSV(rows);
      // BOM for Excel
      const bom = "\uFEFF";
      const blob = new Blob([bom + csv], { type: "text/csv;charset=utf-8" });
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "pathogen_report_template.csv";
      a.click();
      URL.revokeObjectURL(url);
    }
  
    function toCSV(rows: string[][]): string {
      return rows.map(r => r.map(escapeCsv).join(",")).join("\r\n");
    }
    function escapeCsv(v: string): string {
      if (v == null) return "";
      const s = String(v);
      // quote if contains comma, quote, or newline; escape quotes
      return /[",\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
    }
  
    // ---- helpers you already had ----
    function indexMap(header: string[]): Record<string, number> {
      const map: Record<string, number> = {};
      header.forEach((h, i) => map[normalize(h)] = i);
      return map;
    }
    function rowToObj(row: string[], idx: Record<string, number>): Record<string, string> {
      const out: Record<string, string> = {};
      for (const [csvName, path] of Object.entries(FIELD_MAP)) {
        const i = idx[normalize(csvName)];
        if (i != null) out[path] = row[i] ?? "";
      }
      return out;
    }
    function normalize(s: string) { return s.trim().toLowerCase(); }
    function applyToSchema(schema: any, kv: Record<string, string>) {
      for (const [path, raw] of Object.entries(kv)) setDeep(schema, path, coerce(path, raw));
    }
    function setDeep(obj: any, path: string, value: any) {
      const parts = path.split('.');
      let cur = obj;
      for (let i = 0; i < parts.length - 1; i++) cur = cur[parts[i]] ??= {};
      cur[parts.at(-1)!] = value;
    }
    function coerce(path: string, v: string): any {
      if (path === "patient_result.pathogen_detected") {
        const s = v.trim().toLowerCase();
        return s === "true" || s === "1" || s === "yes" || s === "y";
      }
      return v;
    }

    const cleanQuotes = (s: string) => {
        // remove surrounding quotes and unescape double quotes
        if (s.startsWith('"') && s.endsWith('"')) {
          return s.slice(1, -1).replace(/""/g, '"');
        }
        return s;
      };

    // Minimal RFC4180 CSV parser
    function parseCSV(input: string): string[][] {
      const rows: string[][] = [];
      let row: string[] = [];
      let field = '';
      let i = 0, q = false;


      while (i < input.length) {
        const c = input[i];
        if (q) {
          if (c === '"') {
            const n = input[i + 1];
            if (n === '"') { field += '"'; i += 2; continue; }
            q = false; i++; continue;
          }
          field += c; i++; continue;
        }
        if (c === '"') { q = true; i++; continue; }
        if (c === ',') { row.push(cleanQuotes(field)); field = ''; i++; continue; }
        if (c === '\n') { row.push(cleanQuotes(field)); rows.push(row); row = []; field = ''; i++; continue; }
        if (c === '\r') { row.push(cleanQuotes(field)); rows.push(row); row = []; field = ''; i += input[i + 1] === '\n' ? 2 : 1; continue; }
        field += c; i++;
      }
      if (q) throw new Error('Unclosed quote');
      row.push(cleanQuotes(field));
      if (!(row.length === 1 && row[0] === '' && rows.length > 0)) rows.push(row);
      return rows;
    }

  </script>
  
  <div class="grid grid-cols-4 pt-4">
    <FileDropzone name="files" on:change={onChangeHandler} class="col-span-3">
      <svelte:fragment slot="lead">
        <div class="flex justify-center items-center">
          <div class="md:w-12 text-secondary-500">
            <svg data-slot="icon" aria-hidden="true" fill="none" stroke-width="1.5" stroke="currentColor" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                <path d="M9 8.25H7.5a2.25 2.25 0 0 0-2.25 2.25v9a2.25 2.25 0 0 0 2.25 2.25h9a2.25 2.25 0 0 0 2.25-2.25v-9a2.25 2.25 0 0 0-2.25-2.25H15m0-3-3-3m0 0-3 3m3-3V15" stroke-linecap="round" stroke-linejoin="round"></path>
            </svg>                    
          </div>    
        </div>
      </svelte:fragment>
      <svelte:fragment slot="message">CSV File</svelte:fragment>
      <svelte:fragment slot="meta">.csv</svelte:fragment>
    </FileDropzone>
  </div>

  <div class="flex gap-6 w-3/4 pb-12 pt-12 justify-end">


    {#if processing}
      <div class="flex justify-center align-center items-center">
          <ProgressRadial width="w-12" stroke={30} meter="stroke-primary-500" track="stroke-primary-500/30" />
      </div>
    {/if}
    <button class="btn variant-ghost-primary flex items-center" on:click={downloadBatchZip} disabled={batch.length == 0 || busy}>
      <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-5 h-5 mr-2">
        <path stroke-linecap="round" stroke-linejoin="round" d="M3 16.5v2.25A2.25 2.25 0 0 0 5.25 21h13.5A2.25 2.25 0 0 0 21 18.75V16.5M16.5 12 12 16.5m0 0L7.5 12m4.5 4.5V3" />
      </svg>
      
      {busy ? `Building ${progress}%` : 'Download (ZIP)'}
    </button>
  
    <button class="btn variant-ghost-secondary flex items-center" on:click={() => downloadTemplate(false)}>
      <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-5 h-5 mr-2">
        <path stroke-linecap="round" stroke-linejoin="round" d="M3 16.5v2.25A2.25 2.25 0 0 0 5.25 21h13.5A2.25 2.25 0 0 0 21 18.75V16.5M16.5 12 12 16.5m0 0L7.5 12m4.5 4.5V3" />
      </svg>
      CSV Template
    </button>
  </div>