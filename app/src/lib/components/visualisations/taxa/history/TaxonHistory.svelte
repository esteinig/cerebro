<script lang="ts">
	
    import { getCssVariableAsHex, TaxRank, transformTaxonOverview } from "$lib/utils/helpers";
	import { type TaxonOverviewRecord, type TaxonHistory, AbundanceMode, DisplayData, DisplayTotal, FileTag } from "$lib/utils/types";
	import { ScatterChart, ChartTheme, ScaleTypes, type ScatterChartOptions } from "@carbon/charts-svelte";
	import { selectedClientFilterConfig, selectedIdentifiers, selectedTaxa, storeTheme } from "$lib/stores/stores";
	import {CerebroApi, ApiResponse } from "$lib/utils/api";
	import { getToastStore, ProgressRadial } from "@skeletonlabs/skeleton";
	import { page } from "$app/stores";


    export let selectedTaxon: TaxonOverviewRecord;
    export let taxRank: TaxRank = TaxRank.Species;
    export let hostTaxLabel: string = "s__Homo sapiens";

	interface TaxonHistoryData {
        id: string;
        sample_id: string,
        sample_tags: string[],
        sample_type: string,
        run_id: string,
        run_date: string,
		group: string;
		taxon_rpm: number;
		host_rpm: number;
	}

    let loading: boolean = false;

    let publicApi = new CerebroApi();
    let toastStore = getToastStore();

    let plotLogScatter: boolean = true;
    let taxonHistoryData: TaxonHistoryData[] = [];

    getTaxonHistory(selectedTaxon);

    function removeSelectedTaxon() {
        $selectedTaxa = $selectedTaxa.filter(taxon => taxon.taxid !== selectedTaxon.taxid)
    }

    async function getTaxonHistory(selectedTaxon: TaxonOverviewRecord) {

        loading = true;

        let taxon_label: string = `${taxRank}${selectedTaxon.name}`;

        let response: ApiResponse = await publicApi.fetchWithRefresh(
            `${publicApi.routes.cerebro.taxonHistory}?team=${$page.params.team}&db=${$page.params.db}&project=${$page.params.project}&taxon_label=${taxon_label}&host_label=${hostTaxLabel}`,
            { 
                method: 'GET',  
                mode: 'cors',
                credentials: 'include', 
            } as RequestInit,
            $page.data.refreshToken, toastStore, `Taxon history retrieved for: ${taxon_label}`
        )

        loading = false;

        if (response.ok){

            taxonHistoryData = response.json.data.map((taxonHistory: TaxonHistory) => {
                
                let taxonOverviewRecords: TaxonOverviewRecord[] = transformTaxonOverview(taxonHistory.taxa, AbundanceMode.Mixed, DisplayData.Rpm, DisplayTotal.Average, $selectedClientFilterConfig);

                 // Avoid division by zero and handle null host_reads.
                let host_rpm = (taxonHistory.input_reads > 0 && taxonHistory.host_reads !== null)
                    ? (taxonHistory.host_reads / taxonHistory.input_reads) * 1000000
                    : 0;


                // Lookup additional host record by checking if the lineage contains hostTaxLabel.
                const hostRecord = taxonOverviewRecords.find(record => record.lineage.includes(hostTaxLabel));
                if (hostRecord) {
                    host_rpm += hostRecord.total;
                }

                // Lookup the record for the selected taxon and use its total as taxon_rpm
                const selectedRecord = taxonOverviewRecords.find(
                    record => record.taxid === selectedTaxon.taxid
                );


                const taxon_rpm = selectedRecord ? selectedRecord.total : 0;
                
                let group = "Other";
                if ($selectedIdentifiers.includes(taxonHistory.id)) {
                    group = "Selected"
                } else if (taxonHistory.sample_tags.includes(FileTag.DNA)) {
                    if (taxonHistory.sample_tags.includes(FileTag.NTC) || taxonHistory.sample_tags.includes(FileTag.ENV) || taxonHistory.sample_tags.includes(FileTag.POS)) {
                        group = "Control"
                    } else {
                        group = "Library (DNA)"
                    }
                } else if (taxonHistory.sample_tags.includes(FileTag.RNA)) {
                    if (taxonHistory.sample_tags.includes(FileTag.NTC) || taxonHistory.sample_tags.includes(FileTag.ENV) || taxonHistory.sample_tags.includes(FileTag.POS)) {
                        group = "Control"
                    } else {
                        group = "Library (RNA)"
                    }
                }
                
                return {
                    id: taxonHistory.id,
                    sample_id: taxonHistory.sample_id,
                    sample_tags: taxonHistory.sample_tags,
                    sample_type: taxonHistory.sample_type,
                    run_id: taxonHistory.run_id,
                    run_date: taxonHistory.run_date,
                    group: group,
                    taxon_rpm: taxon_rpm,
                    host_rpm: host_rpm
                } as TaxonHistoryData
            })

        }
    }    

    
	let groupColors: Record<string, string> = {};
	
	if (typeof window !== "undefined") {
		groupColors = {
			"Control": getCssVariableAsHex("--color-primary-100", $storeTheme) ?? "#ffffff",
			"Library (DNA)": getCssVariableAsHex("--color-primary-600", $storeTheme) ?? "#ffffff",
			"Library (RNA)": getCssVariableAsHex("--color-secondary-600", $storeTheme) ?? "#ffffff",
			"Selected": getCssVariableAsHex("--color-tertiary-600", $storeTheme) ?? "#ffffff",
			"Other": getCssVariableAsHex("--color-primary-100", $storeTheme) ?? "#ffffff",
		};
	}

	let historyPlotOptions: ScatterChartOptions;

    function getLogHistoryPlotOptions(): ScatterChartOptions {

        const minHost = Math.max(Math.min(...taxonHistoryData.map((d) => d.host_rpm)), 1e-1);
        const maxHost = Math.max(...taxonHistoryData.map((d) => d.host_rpm)) * 1.2;

        const minTaxon = Math.max(Math.min(...taxonHistoryData.map((d) => d.taxon_rpm)), 1e-1);
        const maxTaxon = Math.max(...taxonHistoryData.map((d) => d.taxon_rpm)) * 1.2;

        const generateLogTicks = (min: number, max: number) => {
            const ticks = [];
            let value = Math.pow(10, Math.floor(Math.log10(min)));
            while (value <= max) {
                ticks.push(value);
                value *= 10;
            }
            return ticks;
        };

        const logTicksTaxon = generateLogTicks(minTaxon, maxTaxon);
        const logTicksHost = generateLogTicks(minHost, maxHost);

        return {
            theme: ChartTheme.G100,
            color: { scale: groupColors },
            toolbar: { enabled: false },
            title: `${selectedTaxon.name} (log10)`,
            height: '450px',
            grid: { y: { enabled: false }, x: { enabled: false } },
            axes: {
                left: {
                    title: "Taxon RPM",
                    mapsTo: "taxon_rpm",
                    scaleType: ScaleTypes.LOG,
                    domain: [minTaxon-(0.5*minTaxon), maxTaxon+(0.5*maxTaxon)],
                    ticks: {
                        values: logTicksTaxon, // Explicitly set log ticks
                    }
                },
                bottom: {
                    title: "Host RPM",
                    mapsTo: "host_rpm",
                    scaleType: ScaleTypes.LOG,
                    domain: [minHost-(0.5*minHost), maxHost+(0.5*maxHost)],
                    ticks: {
                        values: logTicksHost, // Explicitly set log ticks
                    }
                },
            },
            tooltip: {
                enabled: true,
                customHTML: (data: any, defaultHTML: string): string => {
                    if (data.length > 0) {
                        return `
                            <div class="custom-tooltip">
                                <p><strong>Run ID:</strong> ${data[0].run_id}</p>
                                <p><strong>Run Date:</strong> ${data[0].run_date}</p>
                                <p><strong>Sample Name:</strong> ${data[0].sample_id}</p>
                                <p><strong>Sample Tags:</strong> ${data[0].sample_tags.join(" ,")}</p>
                                <p><strong>Sample Type:</strong> ${data[0].sample_type}</p>
                                <p><strong>Group:</strong> ${data[0].group}</p>
                                <p><strong>Taxon RPM:</strong> ${data[0].taxon_rpm.toFixed(2)}</p>
                                <p><strong>Host RPM:</strong> ${data[0].host_rpm.toFixed(2)}</p>
                            </div>
                        `;
                    } else {
                        return `
                            <div class="custom-tooltip">
                                <h3>Failed to retrieve data for this library</h3>
                            </div>
                        `
                    }
                }
            }
        };
    }

    function getHistoryPlotOptions(): ScatterChartOptions {

        return {
            theme: ChartTheme.G100,
            color: { scale: groupColors },
            toolbar: { enabled: false },
            title: `${selectedTaxon.name}`,
            height: '450px',
            grid: { y: { enabled: false }, x: { enabled: false } },
            axes: {
                left: {
                    title: "Taxon RPM",
                    mapsTo: "taxon_rpm",
                },
                bottom: {
                    title: "Host RPM",
                    mapsTo: "host_rpm",
                },
            },
            tooltip: {
                enabled: true,
                customHTML: (data: any, defaultHTML: string): string => {
                    if (data.length > 0) {
                        return `
                            <div class="custom-tooltip">
                                <p><strong>Run ID:</strong> ${data[0].run_id}</p>
                                <p><strong>Run Date:</strong> ${data[0].run_date}</p>
                                <p><strong>Sample Name:</strong> ${data[0].sample_id}</p>
                                <p><strong>Sample Tags:</strong> ${data[0].sample_tags.join(" ,")}</p>
                                <p><strong>Sample Type:</strong> ${data[0].sample_type}</p>
                                <p><strong>Group:</strong> ${data[0].group}</p>
                                <p><strong>Taxon RPM:</strong> ${data[0].taxon_rpm.toFixed(2)}</p>
                                <p><strong>Host RPM:</strong> ${data[0].host_rpm.toFixed(2)}</p>
                            </div>
                        `;
                    } else {
                        return `
                            <div class="custom-tooltip">
                                <h3>Failed to retrieve data for this library</h3>
                            </div>
                        `
                    }
                }
            }
        };
        }

    $: {
        if (taxonHistoryData.length > 0) {
            historyPlotOptions = plotLogScatter ? getLogHistoryPlotOptions() : getHistoryPlotOptions();
        }
    }


</script>


<div class="py-12">          
    
    {#if loading}
        <div class="flex justify-center py-24">
            <ProgressRadial width="sm:w-12 md:w-24" stroke={20} meter="stroke-tertiary-500" track="stroke-tertiary-500/30" />
        </div>
    {:else}      
        {#if taxonHistoryData.length > 0}
            <div class="grid grid-cols-2 gap-4 flex">
                <div class="taxon-history-host-chart">
                    <ScatterChart 
                    data={taxonHistoryData} 
                    options={historyPlotOptions}
                    ></ScatterChart>
                </div>
                <button type="button" class="btn-icon btn-icon-sm variant-ghost" on:click={removeSelectedTaxon}>
                    <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" class="w-4 h-4">
                        <path stroke-linecap="round" stroke-linejoin="round" d="M6 18 18 6M6 6l12 12" />
                    </svg>                      
                </button>
            </div>
            
        {:else}
            <div class="flex justify-center">
                <p class="text-sm opacity-60 pt-24">
                    Taxon data not found
                </p>
            </div>
        {/if}
    {/if}

</div>


<style lang="postcss">
	/* Scoped globally but applied only within this component */
    :global(.taxon-history-host-chart .chart-grid-backdrop) {
        display: none;
    }
</style>