#!/usr/bin/env nextflow

/* 
vim: syntax=groovy
-*- mode: groovy;-*-

Cerebro: metagenomic and -transcriptomic diagnostics for clinical production environments
Copyright (C) 2023  Eike Steinig

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


/* 
===================
Workflow parameters
===================
*/

params.help                                          = false
params.monochrome                                    = false
params.outdir                                        = "cerebro_v${workflow.manifest.version}"

// File input from command-line

params.fastq_pe                                      = null
params.fastq_ont                                     = null

// ===================================
// Viral detection and assembly module
// ===================================

// Background depletion

params.virus_background_references                   = null                                                 // white-space separated string of background reference files to use in background depletion (minimap2)
params.virus_background_dbs                          = null                                                 // white-space separated string of database paths to use in background depletion (Kraken2)
params.virus_background_taxa                         = "Bacteria Archaea Eukaryota Holozoa Nucletmycea"     // white-space separated string of taxonomic identifiers or names to conduct sub-level depletion for (Scrubby) - must be below `Domain` level
params.virus_background_direct                       = "131567"                                             // white-space separated string of taxonomic identifiers or names for depletion of reads directly assigned to the taxon (Scrubby) [cellular organisms]

// Viral database and blacklist

params.virus_db_index                                = null
params.virus_db_fasta                                = null
params.virus_blacklist                               = null


params.subset_min_shared_hashes                     = 2 
params.subset_group_index                           = null     // null or int, recognized in process and switches to group selection [e.g 1 on seq index kraken fmt "taxid|1147|NC_014725.1"] -- picks the best genome by max shared hashes from taxid group
params.subset_group_sep                             = "|"      // format required for both bacterial and eukaryotic *|{taxid}|*

// Dummy parameters for domain-specific scanning stage in the alignment pipelines
// these are only used to configure processes with specific parameters for viral
// bacterial and eukaryotic alignment options

 


// Workflow imports

include { cerebro_production }              from './lib/workflows/subworkflows/production'
include { quality_control_illumina }        from './lib/workflows/subworkflows/quality_control'
include { quality_control_ont }             from './lib/workflows/subworkflows/quality_control'
include { pathogen_detection }              from './lib/workflows/pathogen_detection'
include { aneuploidy_detection_illumina }   from './lib/workflows/aneuploidy_detection'
include { culture_identification }          from './lib/workflows/culture_identification'
include { reference_alignment }             from './lib/workflows/validation/reference_alignment'

// Utility imports
include { PingServer; ProcessSamples; ProcessSamplesTaxonomy; QualityControlTable } from './lib/processes/cerebro'
include { RasusaReadsMultiple } from './lib/processes/rasusa' addParams(subdir: "subsample")
include { c; parse_file_params; init_msg; complete_msg; help_msg; get_paired_reads; get_single_reads; WriteConfig } from './lib/utils'

workflow {

    // Startup message and help
    params.help ? help_msg() : init_msg()

    if (params.test.dev.enabled){
        log.info "${c('red')}Development mode for utilities active, inputs and modules deactivated.${c('reset')}"
        exit 0
    }

    if (params.test.io.enabled) {
        log.info "${c('red')}Development test mode for inputs active, modules deactivated.${c('reset')}"
        def inputs = parse_file_params()
        exit 0
    }

    /*
    ==========================
    Input validation and reads
    ==========================
    */

    def inputs = parse_file_params()

    if (params.production.enabled) {

        if (params.production.api.upload.enabled) {
            PingServer(Channel.of("PING")) | collect                             // collect to await result before proceeding
        }
        
        if (params.ont.enabled) {
            data = get_single_reads(inputs.sample_sheet, true)
        } else {
            data = get_paired_reads(inputs.sample_sheet, true)
        }

        reads = data.pathogen                                                    // pathogen analysis uses all input reads
        reads_aneuploidy = data.aneuploidy                                       // host genome analysis if enabled uses subset of activated input reads
    } else {
        if (params.fastq_ont && params.ont.enabled) {
            reads = get_single_reads(null, false)                               
        } else if (params.fastq_pe) {
            reads = get_paired_reads(null, false)                               
            reads_aneuploidy = reads                                             // all input read files are also used for host genome analysis if enabled
        } else {
            error"${c('red')}Production mode is not active - you must provide either `--fastq_pe <GLOB>` or `--ont.enabled --fastq_ont <FILE-OR-GLOB>` for read inputs${c('reset')}"
        }
        
    }   

        
    if (params.subsample.enabled) {
        reads = RasusaReadsMultiple(reads, params.subsample.reads, 1..params.subsample.replicates)
    }

    /*
    =================
    Workflow modules
    =================
    */

    started = java.time.LocalDateTime.now()

    if (params.test.qc.enabled) {
        println "\n${c('red')}Quality control mode is active, other modules deactivated.${c('reset')}\n"

        // ===========================
        // Quality control subworkflow
        // ===========================

        
        if (params.ont.enabled) {
            quality_control_ont(
                reads, 
                inputs.ercc_fasta, 
                inputs.phage_fasta, 
                inputs.host_depletion_kraken2, 
                inputs.host_depletion_reference, 
                params.qc.host.depletion.enabled,
                params.qc.controls.phage.enabled
            )
        } else {
            quality_control_illumina(
                reads, 
                inputs.adapter_fasta, 
                inputs.ercc_fasta, 
                inputs.phage_fasta, 
                inputs.host_depletion_kraken2, 
                inputs.host_depletion_reference, 
                params.qc.deduplication.enabled,
                params.qc.deduplication.method,
                params.qc.host.depletion.enabled,
                params.qc.controls.phage.enabled
            )
        }

    } else {
        
        if (params.validation.enabled) {

            // ===============================
            // Validation experiment workflows
            // ===============================
            
            reference_alignment(inputs)

        } else {

            // =====================================
            // Host aneuploidy detection subworkflow
            // =====================================
            
            if (params.host.enabled && params.host.aneuploidy.enabled) {
                aneuploidy = aneuploidy_detection_illumina(reads_aneuploidy, inputs)
            } else {
                aneuploidy = Channel.empty()
            }

            // ==============================
            // Pathogen detection subworkflow
            // ==============================

            if (params.taxa.enabled) {

                pathogen_detection(reads, inputs, params.ont.enabled, params.taxa.enabled)

                WriteConfig(
                    inputs.sample_sheet, 
                    pathogen_detection.out.results | collect, 
                    started
                )

                if (params.production.enabled && params.production.api.upload.enabled) {
                    cerebro_production(
                        inputs.taxonomy_directory,
                        pathogen_detection.out.results,
                        inputs.sample_sheet,
                        WriteConfig.out.config
                    )
                } else {
                    if (params.process.enabled) {
                        if (params.process.taxa) {
                            samples = ProcessSamplesTaxonomy(pathogen_detection.out.results, inputs.taxonomy_directory)
                        } else {
                            samples = ProcessSamples(pathogen_detection.out.results)
                        }
                        samples | map { it -> it[1] } | QualityControlTable
                    }
                }
            }


            // ==============================
            // Pathogen detection subworkflow
            // ==============================

            if (params.taxa.enabled) {

            }

            // ==================================================================
            // Cultured isolate assembly for taxonomic identification subworkflow
            // ==================================================================

            if (params.culture.enabled) {
                ont_reads = get_single_reads(inputs.sample_sheet, false);
                pe_reads = get_paired_reads(inputs.sample_sheet, false);
                
                culture_identification(ont_reads, pe_reads, inputs);
            }
        }
    }
}   

workflow.onComplete {
    complete_msg()
}


workflow.onError {
    println "\n${c('red')}Pipeline execution stopped with the following message: ${workflow.errorMessage}${c('reset')}\n"
}
