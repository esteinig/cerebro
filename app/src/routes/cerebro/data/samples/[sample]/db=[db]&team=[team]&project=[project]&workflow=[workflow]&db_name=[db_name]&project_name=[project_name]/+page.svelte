<script>
   
    
    //    $: serverFiltersActive = serverFilterConfig.domains.length > 0 ||
    //         serverFilterConfig.tags.length > 0 ||
    //         serverFilterConfig.kmer_min_reads > 0 ||
    //         serverFilterConfig.kmer_databases.length > 0 ||
    //         serverFilterConfig.alignment_min_reads > 0 ||
    //         serverFilterConfig.alignment_min_bases > 0 ||
    //         serverFilterConfig.alignment_min_regions > 0 ||
    //         serverFilterConfig.alignment_min_coverage > 0 ||
    //         serverFilterConfig.alignment_min_ref_length > 0 ||
    //         serverFilterConfig.assembly_min_contig_length > 0 ||
    //         serverFilterConfig.assembly_min_contig_identity > 0 ||
    //         serverFilterConfig.assembly_min_contig_coverage > 0
    
    //     $: clientFiltersActive = clientFilterConfig.domains.length > 0 ||
    //         clientFilterConfig.genera.length > 0 ||
    //         clientFilterConfig.species.length > 0 ||
    //         clientFilterConfig.modules.alignment ||
    //         clientFilterConfig.modules.kmer ||
    //         clientFilterConfig.modules.assembly ||
    //         clientFilterConfig.minimum.rpm > 0 ||
    //         clientFilterConfig.minimum.rpm_kmer > 0 ||
    //         clientFilterConfig.minimum.rpm_alignment > 0 ||
    //         clientFilterConfig.minimum.contigs > 0 ||
    //         clientFilterConfig.minimum.bases > 0
    
    
    //     let showServerSideFilters: boolean = false;
</script>


        <!-- {#if selectedView === "classification"}
            <p class="mb-1 mt-4">
                <span class="opacity-60">Taxonomy filters</span>
            </p>
            <div class="w-full border border-primary-500 rounded-md p-4" >
                
                <div class="grid grid-cols-3 sm:grid-cols-3 md:grid-cols-3 gap-x-4 align-center">
                    <div class="col-span-2">
                        <p class="opacity-60 text-xs pb-8 pr-4">
                            {#if showServerSideFilters}
                                Apply classification table filters by switching to client-side filters
                            {:else}
                                Apply fine-grained evidence filters by switching to server-side filters
                            {/if}
                            
                        </p>
                    </div>
                    <div class="col-span-1 text-xs">
                        <SlideToggle name="server-side-filters" bind:checked={showServerSideFilters} active="variant-filled-tertiary dark:variant-filled-tertiary" size="sm">Server-side</SlideToggle>
                    </div>
                </div>
                {#if showServerSideFilters}
                    <ServerFilterConfiguration bind:serverFilterConfig={serverFilterConfig}></ServerFilterConfiguration>
                    <div class="text-center py-8">
                        <button type="button" class="btn variant-outline-primary w-3/4" on:click={() => reloadTable()}>
                            <svg aria-hidden="true" fill="none" stroke="currentColor" class="w-5 h-5 mr-3" stroke-width="1.5" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
                                <path d="M16.023 9.348h4.992v-.001M2.985 19.644v-4.992m0 0h4.992m-4.993 0l3.181 3.183a8.25 8.25 0 0013.803-3.7M4.031 9.865a8.25 8.25 0 0113.803-3.7l3.181 3.182m0-4.991v4.99" stroke-linecap="round" stroke-linejoin="round"></path>
                            </svg>
                            Reload Table
                        </button>
                    </div>
                {:else}
                    <ClientFilterConfiguration bind:clientFilterConfig={clientFilterConfig}></ClientFilterConfiguration>
                {/if}
                
            </div>

            <p class="mb-1 mt-4">
                <span class="opacity-60">Taxonomy highlights</span>
            </p>
            <div class="w-full border border-primary-500 rounded-md p-4" >
                
                <TaxonomyHighlights title="Common contamination" bind:highlightConfig={contamHighlightConfig}></TaxonomyHighlights>
                <TaxonomyHighlights title="Syndromic pathogens" bind:highlightConfig={syndromeHighlightConfig}></TaxonomyHighlights>

            </div>

        {/if}
    </div>
    {#if selectedView === "qc"}
        <div class="col-span-2">
            <div class="w-full mb-5">
                <p class="mb-1"><span class="opacity-60">Quality Control</span></p>
                <div class="mb-4 border border-primary-500 p-8 rounded-md">
                    <QualityControl selectedWorkflowConfiguration={selectedWorkflowConfiguration} selectedModels={selectedModels} selectedQualityControlSummaries={selectedQualityControlSummaries} />
                </div>
            </div>
        </div>
    {:else if selectedView === "classification"}
        <div class="col-span-2">
            <div>
                <div class="w-full mb-5">
                    <div class="">
                        <span class="opacity-60">Classification overview</span>
                        {#if clientFiltersActive}
                            <span class="text-warning-500 opacity-60 text-xs ml-5">
                                Client-side table filters are active
                            </span>
                         {/if}
                        {#if serverFiltersActive}
                            <span class="text-warning-500 opacity-60 text-xs ml-5">
                                Server-side evidence filters are active
                            </span>
                        {/if}
                    </div>
                    <div class="mb-4 border border-primary-500 p-4 rounded-md">
                        <Classification selectedIdentifiers={selectedIdentifiers} selectedModels={selectedModels} clientFilterConfig={clientFilterConfig} serverFilterConfig={serverFilterConfig} taxonHighlightConfig={taxonHighlightConfig}></Classification>
                    </div>
                </div>
            </div>
        </div>
    {:else if selectedView === "candidates"}
        <div class="col-span-2">
            <div class="w-full mb-5">
                <p class="mb-1"><span class="opacity-60">Candidates</span></p>
                <div class="mb-4 border border-primary-500 p-4 rounded-md">
                    <Candidates selectedModels={selectedModels}></Candidates>
                </div>
            </div>
        </div>
    {:else if selectedView === "comments"}
        <div class="col-span-2">
            <div class="w-full mb-5">
                <p class="mb-1">
                    <span class="opacity-60">Sample comments</span>
                    <span class="opacity-60 text-xs ml-5">Comments are submitted for all libraries of <span class="font-semibold">{$page.params.sample}</span>. Other team members are able to view and submit comments.
                </p>
                <div class="mb-4 border border-primary-500 p-4 rounded-md">
                    <CommentBox selectedDatabaseId={$page.params.db} selectedProjectId={$page.params.project}></CommentBox>
                </div>
            </div>
        </div>
    {:else}
        <div class="col-span-2">
            <div class="w-full mb-5">
                <p class="mb-1"><span class="opacity-60">Sample reports</span>
                    <span class="opacity-60 text-xs ml-5"></p>
                <div class="mb-4 border border-primary-500 p-4 rounded-md">
                    <Reports selectedWorkflowConfiguration={selectedWorkflowConfiguration} selectedIdentifiers={selectedIdentifiers} selectedModels={selectedModels}></Reports>
                </div>
            </div>
        </div>
    {/if} -->


    