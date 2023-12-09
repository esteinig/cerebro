<script lang="ts">

    import { page } from "$app/stores";
	import CerebroApi, { ApiResponse } from "$lib/utils/api";
	import { baseTags, getDateTimeString } from "$lib/utils/helpers";
    import type { Cerebro, PatientHeaderSchema, PatientResultSchema, PriorityTaxon, ReportSchema, WorkflowConfig } from "$lib/utils/types";
	import { getToastStore, ProgressBar, RadioGroup, RadioItem, SlideToggle, type ToastSettings } from "@skeletonlabs/skeleton";
    
    export let selectedModels: Cerebro[];
    export let selectedIdentifiers: string[];
    export let selectedWorkflowConfiguration: WorkflowConfig;

    const publicApi: CerebroApi = new CerebroApi();
    const toastStore = getToastStore();
    
    let requireFields: boolean = false;
    let loading: boolean = false;
    let pdf: number = 0;

    export let selectedPriorityTaxon: PriorityTaxon | null = null;

    let secureDataSubmission: boolean = true;
    let noPermissionMessage: string = "Not permitted";

    let sampleIsNegative: boolean = false;
    let taxonNegativeControl: boolean = false;

    // Report parameters
    let schemaSampleId: string = $page.params.sample;
    let schemaRunId: string = "";

    let schemaPatientHeaderName: string = "";
    let schemaPatientHeaderUrn: string = "";
    let schemaPatientHeaderDob: string = "";


    let schemaPatientHeaderSpecimenId: string = "";
    let schemaPatientHeaderSpecimenType: string = "";
    let schemaPatientHeaderDateCollected: string = "";
    let schemaPatientHeaderDateReceived: string = "";

    let schemaPatientHeaderRequestedDoctor: string = "";
    let schemaPatientHeaderHospitalSite: string = "";
    let schemaPatientHeaderLaboratoryNumber: string = "";


    let schemaPatientResultNegative: boolean = true;
    let schemaPatientResultOrganism: string = "NEGATIVE";
    let schemaPriorityTaxon: PriorityTaxon | null = null;

    let schemaPatientResultComments: string = "No further comments";
    let schemaPatientResultAction: string = "No further actions taken.";
    let schemaPatientResultContact: string = "";
    let schemaPatientResultReviewDate: string = "";
    let schemaLaboratoryComments: string = "No further comments.";
    let schemaBioinformaticsComments: string = "No further comments.";
    let schemaOtherEvidence: string | null = null;

    let schemaReportingLaboratory: string = ""
    let schemaReportingDate: string = ""
    let schemaReportingLocation: string = ""
    
    $: {
        schemaPatientResultNegative = selectedPriorityTaxon === null ? true : false;
        selectedPriorityTaxon === null ? sampleIsNegative = true : sampleIsNegative = false;
        schemaPatientResultOrganism = selectedPriorityTaxon === null ? "NEGATIVE" : selectedPriorityTaxon.taxon_overview.name;
        schemaPriorityTaxon = selectedPriorityTaxon;
    }
    

    $: {
        schemaRunId = selectedModels.map(cerebro => cerebro.run.id).filter((value, index, self) => self.findIndex(v => v === value) === index).join(", ");
    }


    // Undersores break Latex engine
    const sanitizeLatexString = (value: string) => {
        return value.replace(/_/g, "\\_")
    }

    const createReport = async() => {

        loading = true;


        let reportSchema: ReportSchema = {
            ids: selectedIdentifiers,
            user_id: sanitizeLatexString($page.data.userData.id),
            user_name: sanitizeLatexString($page.data.userData.name),
            sample_id: sanitizeLatexString($page.params.sample),
            run_id: sanitizeLatexString(schemaRunId),
            negative: sampleIsNegative,
            workflow: structuredClone(selectedWorkflowConfiguration),
            patient_header: {
                patient_name: sanitizeLatexString(schemaPatientHeaderName),
                patient_urn: sanitizeLatexString(schemaPatientHeaderUrn),
                patient_dob: sanitizeLatexString(schemaPatientHeaderDob),
                requested_doctor: sanitizeLatexString(schemaPatientHeaderRequestedDoctor),
                hospital_site: sanitizeLatexString(schemaPatientHeaderHospitalSite),
                laboratory_number: sanitizeLatexString(schemaPatientHeaderLaboratoryNumber),
                specimen_id: sanitizeLatexString(schemaPatientHeaderSpecimenId),
                date_collected: sanitizeLatexString(schemaPatientHeaderDateCollected),
                date_received: sanitizeLatexString(schemaPatientHeaderDateReceived),
                specimen_type: sanitizeLatexString(schemaPatientHeaderSpecimenType),
                reporting_laboratory: sanitizeLatexString(schemaReportingLaboratory),
                reporting_date: sanitizeLatexString(schemaReportingDate),
                reporting_location: sanitizeLatexString(schemaReportingLocation),
            } satisfies PatientHeaderSchema,
            patient_result: {
                review_date: sanitizeLatexString(schemaPatientResultReviewDate),
                negative: sampleIsNegative,
                organism: sanitizeLatexString(schemaPatientResultOrganism),
                contact: sanitizeLatexString(schemaPatientResultContact),
                comments: sanitizeLatexString(schemaPatientResultComments),
                actions: sanitizeLatexString(schemaPatientResultAction),
            } satisfies PatientResultSchema,
            laboratory_comments: sanitizeLatexString(schemaLaboratoryComments),
            bioinformatics_comments: sanitizeLatexString(schemaBioinformaticsComments),
            priority_taxon: schemaPriorityTaxon,
            priority_taxon_negative_control: taxonNegativeControl,
            priority_taxon_other_evidence: schemaOtherEvidence === null ? null : sanitizeLatexString(schemaOtherEvidence),
        }
        
        // Format some of the dates and other workflow header parameters
        reportSchema.workflow.id = reportSchema.workflow.id.substring(0, 8)
        reportSchema.workflow.started = getDateTimeString(reportSchema.workflow.started, true)
        reportSchema.workflow.completed = getDateTimeString(reportSchema.workflow.completed, true)
        
        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.createReport}/${pdf === 0 ? "pdf" : "tex"}?db=${$page.params.db}&project=${$page.params.project}`,
            { 
                method: 'POST',  
                mode: 'cors',
                credentials: 'include', 
                headers: { 'Content-Type': 'application/json' }, 
                body: JSON.stringify(reportSchema) 
            } as RequestInit,
            $page.data.refreshToken, toastStore, "Report generated"
        )

        loading = false;

        if (response.ok) {
            let data: Blob;
            let fileName: string;
            let now = getDateTimeString(new Date().toISOString(), false);
            if (response.json.data?.pdf) {
                // PDF base64-encoded string to blob
                const byteArray = Uint8Array.from(
                    atob(response.json.data.report)
                    .split('')
                    .map(char => char.charCodeAt(0))
                );
                data = new Blob([byteArray], { type: 'application/pdf' });
                fileName = `${now}_${$page.params.sample}_ClinicalReport.pdf`;
            } else if (!response.json.data?.pdf) {
                // LaTeX report as string
                data = new Blob([response.json.data.report], {type: 'text/plain'});
                fileName = `${$page.params.sample}_Report_${now}_ClinicalReport.tex`;
            } else {
                toastStore.trigger({
                    message: "Failed to extract report data in response",
                    background: "variant-filled-tertiary"
                } satisfies ToastSettings);
                return undefined
            };

            let link = document.createElement("a");
            if (link.download !== undefined) { // feature detection
                // Browsers that support HTML5 download attribute
                let url = URL.createObjectURL(data);
                link.setAttribute("href", url);
                link.setAttribute("download", fileName);
                link.style.visibility = 'hidden';
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
            }
        }
        
    }

</script>

<div>
    <form>
        <p class="pb-3">Submission</p>
        <p class="text-xs opacity-60 w-3/4 pb-5">
            Submission details for this report
        </p>
            <div class="p-4">
            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-3 gap-12 w-3/4 pb-12">
                <label class="label">
                    <span class="text-sm opacity-60">Reporting laboratory</span>
                    <input type="text" class="input text-sm" bind:value={schemaReportingLaboratory} required={requireFields} />
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Reporting date</span>
                    <input type="text" class="input text-sm"  bind:value={schemaReportingDate} required={requireFields} />
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Reporting location</span>
                    <input type="text" class="input text-sm" bind:value={schemaReportingLocation} required={requireFields} />
                </label>
            </div>
        </div>
        <p class="pb-3">Sample</p>
        <p class="text-xs opacity-60 w-3/4 pb-5">
            Report is generated for a biological sample analysed with a specific workflow run. 
            If multiple samples have been analysed for this patient, or multiple workflows have been run for this sample, they must be reported separately.
        </p>
        <div class="p-4">
            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-4 gap-12 w-full pb-10">
                <label class="label">
                    <span class="text-sm opacity-60">Sample</span>
                    <input type="text" class="input text-sm" value={schemaSampleId} disabled/>
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Runs</span>
                    <input type="text" class="input text-sm" bind:value="{schemaRunId}" disabled/>
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Workflow</span>
                    <input type="text" class="input text-sm" value="{selectedWorkflowConfiguration.id.substring(0,8)} ({selectedWorkflowConfiguration.name})" disabled/>
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Libraries</span>
                    <input type="text" class="input text-sm" value="{baseTags(selectedModels.map(cerebro => cerebro.sample.tags), true).join(", ")}" disabled/>
                </label>
            </div>
        </div>

        <p class="pb-3">Results</p>
        <p class="text-xs opacity-60 w-3/4 pb-5">
            Diagnostic result to be reported. Pathogens must be selected and recorded as candidate taxa in the pathogen category before they can be reported using this form. Multiple detected pathogens must be reported separately.
        </p>
        <div class="p-4">
            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-3 gap-12 w-full pb-10 align-center items-center">
                <label class="label">
                    <span class="text-sm opacity-60">Organism</span>
                    <input type="text" class="input text-sm {schemaPatientResultOrganism === "NEGATIVE" ? "" : "italic"}" bind:value="{schemaPatientResultOrganism}" disabled />
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Results reviewed on</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientResultReviewDate}  placeholder="DD/MM/YYYY" required={requireFields}  />
                </label>
                <label class="label flex align-center items-center">
                    <div>
                        <input class="checkbox" type="checkbox" checked={sampleIsNegative} />
                    </div>
                    <div class="opacity-60 text-sm ml-3 w-1/2">Sample is negative</div>
                </label>
            </div>

            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-3 gap-12 w-full pb-12">
                <label class="label">
                    <span class="text-sm opacity-60">Contact details</span>
                    <textarea class="textarea text-sm" bind:value={schemaPatientResultContact}  placeholder="Contact details for results" required={requireFields} />
                </label>

                <label class="label">
                    <span class="text-sm opacity-60">General comments</span>
                    <textarea class="textarea text-sm" bind:value={schemaPatientResultComments} placeholder="Additional comments on results" required={requireFields} />
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Actions taken</span>
                    <textarea class="textarea text-sm" bind:value={schemaPatientResultAction} placeholder="Actions taken based on results" required={requireFields} />
                </label>
            </div>
        </div>

        <p class="pb-3">Details</p> 
        <p class="text-xs opacity-60 w-3/4 pb-5">
            Request details included in the header of the report. Please note that some of these fields may be protected information. Submission may be disabled by the system administrator. Fields can be filled in manually in the generated report.
        </p>
        <div class="p-4">
            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-3 gap-12 w-full pb-12">
                <label class="label">
                    <span class="text-sm opacity-60">Patient name</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientHeaderName}  placeholder={secureDataSubmission ? noPermissionMessage: ""}  disabled={secureDataSubmission} />
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Patient URN</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientHeaderUrn}  placeholder={secureDataSubmission ? noPermissionMessage: ""}  disabled={secureDataSubmission} />
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Patient DOB</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientHeaderDob}  placeholder={secureDataSubmission ? noPermissionMessage: ""}   disabled={secureDataSubmission} />
                </label>
            </div>
            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-4 gap-12 w-full pb-16">
                <label class="label">
                    <span class="text-sm opacity-60">Specimen ID</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientHeaderSpecimenId} placeholder={secureDataSubmission ? noPermissionMessage: ""} disabled={secureDataSubmission} />
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Specimen type</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientHeaderSpecimenType} placeholder={secureDataSubmission ? noPermissionMessage: ""} disabled={secureDataSubmission} />
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Date collected</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientHeaderDateCollected}  placeholder={secureDataSubmission ? noPermissionMessage: "DD/MM/YYYY"} disabled={secureDataSubmission} />
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Date received</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientHeaderDateReceived} placeholder={secureDataSubmission ? noPermissionMessage: "DD/MM/YYYY"} disabled={secureDataSubmission} />
                </label>
            </div>
            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-3 gap-12 w-3/4 pb-16">
                <label class="label">
                    <span class="text-sm opacity-60">Requested Doctor</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientHeaderRequestedDoctor} placeholder={secureDataSubmission ? noPermissionMessage: ""} disabled={secureDataSubmission} />
                </label>

                <label class="label">
                    <span class="text-sm opacity-60">Hospital site</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientHeaderHospitalSite} placeholder={secureDataSubmission ? noPermissionMessage: ""} disabled={secureDataSubmission} />
                </label>
                <label class="label">
                    <span class="text-sm opacity-60">Laboratory number</span>
                    <input type="text" class="input text-sm" bind:value={schemaPatientHeaderLaboratoryNumber} placeholder={secureDataSubmission ? noPermissionMessage: ""} disabled={secureDataSubmission}  />
                </label>
            </div>
        </div>


        <p class="pb-2">Processing</p>
        <p class="text-xs opacity-60 w-3/4 pb-5">
            Additional comments on sequencing assay and bioinformatics analysis
        </p>
        <div class="p-4">
            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 w-3/4 gap-12 pb-12">
                <label class="label">
                    <span class="text-sm opacity-60">Assay comments</span>
                    <textarea class="textarea text-sm" bind:value={schemaLaboratoryComments} placeholder="No additional comments"/>
                </label>

                <label class="label">
                    <span class="text-sm opacity-60">Bioinformatics comments</span>
                    <textarea class="textarea text-sm" bind:value={schemaBioinformaticsComments} placeholder="No additional comments"/>
                </label>
            </div>
        </div>
        {#if !sampleIsNegative}
            <p class="pb-2">Evidence</p>
            <p class="text-xs opacity-60 w-3/4 pb-5">
                Evidence for a detection is populated automatically from the selected candidate. 
                Additional evidence comments are not considered under the accreditation scheme and 
                are highlighted as such in the report. This field allows you to record additional evidence
                features or results from manual analyses.
            </p>
            <div class="p-4">
                <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-2 gap-12 pb-12">
                    <label class="label">
                        <span class="text-sm opacity-60">Additional evidence</span>
                        <textarea class="textarea text-sm" bind:value={schemaOtherEvidence} placeholder="No additional comments"/>
                    </label>
                    <label class="label flex align-center items-center">
                        <div>
                            <input class="checkbox" type="checkbox" checked={taxonNegativeControl} />
                        </div>
                        <div class="opacity-60 text-sm ml-3 w-1/2">Organism occurs at significant levels in negative controls</div>
                    </label>
                </div>
            </div>
        {/if}
        <p class="pb-2">Format</p>
        <p class="text-xs opacity-60 w-3/4 pb-5">
            Available output formats are templated LaTeX or compiled PDF
        </p>
        <div class="p-4">
            <div class="grid grid-cols-1 sm:grid-cols-1 md:grid-cols-1 gap-12 pb-12">
                <div class="flex items-center align-center space-x-4 text-sm">
                    <RadioGroup>
                        <RadioItem bind:group={pdf} name="justify" value={0}>PDF</RadioItem>
                        <RadioItem bind:group={pdf} name="justify" value={1}>LaTeX</RadioItem>
                    </RadioGroup>
                </div>
            </div>
        </div>

        {#if loading}
        <div>
            <label class="label">
                <span class="mb-1">Requesting report ...</span><span class="ml-2 text-xs opacity-60">this may take a few seconds ...</span>
                <div><ProgressBar /></div>
            </label>
        </div>
        {:else}
            <div class="flex items-center justify-center">
                <button class="btn variant-outline-primary text-lg p-4" type="submit" on:click={createReport}>
                    <div class="h-7 w-7 mr-2">
                        <svg aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                            <path d="M19.5 14.25v-2.625a3.375 3.375 0 00-3.375-3.375h-1.5A1.125 1.125 0 0113.5 7.125v-1.5a3.375 3.375 0 00-3.375-3.375H8.25M9 16.5v.75m3-3v3M15 12v5.25m-4.5-15H5.625c-.621 0-1.125.504-1.125 1.125v17.25c0 .621.504 1.125 1.125 1.125h12.75c.621 0 1.125-.504 1.125-1.125V11.25a9 9 0 00-9-9z" stroke-linecap="round" stroke-linejoin="round"></path>
                        </svg>
                    </div>
                        Generate Report
                </button>
            </div>
        {/if}
    </form>
</div>