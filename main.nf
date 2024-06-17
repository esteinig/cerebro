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



params.bacteria_mash_index                          = "db/bacteria.msh"
params.bacteria_fasta                               = "db/bacteria.fasta"

params.eukaryots_mash_index                         = "db/eukaryots_complete.msh"
params.eukaryots_fasta                              = "db/eukaryots_complete.fasta"

params.subset_min_shared_hashes                     = 2 
params.subset_group_index                           = null     // null or int, recognized in process and switches to group selection [e.g 1 on seq index kraken fmt "taxid|1147|NC_014725.1"] -- picks the best genome by max shared hashes from taxid group
params.subset_group_sep                             = "|"      // format required for both bacterial and eukaryotic *|{taxid}|*

// Dummy parameters for domain-specific scanning stage in the alignment pipelines
// these are only used to configure processes with specific parameters for viral
// bacterial and eukaryotic alignment options

params.vircov_scan_min_cov                           = null  
params.vircov_scan_min_len                           = null 
params.vircov_scan_min_mapq                          = null  // this applies to ungrouped alignment and filter reads that do not map uniquely to database references
params.vircov_scan_reads                             = null  
params.vircov_scan_coverage                          = null 
params.vircov_scan_regions                           = null
params.vircov_scan_regions_coverage                  = null
params.vircov_group_by                               = null
params.vircov_group_sep                              = null
params.vircov_group_select_by                        = null        
params.vircov_group_select_segment_field             = null
params.vircov_group_select_segment_field_nan         = null
params.vircov_remap_min_cov                          = null  
params.vircov_remap_min_len                          = null
params.vircov_remap_min_mapq                         = null
params.vircov_remap_reads                            = null
params.vircov_remap_coverage                         = null
params.vircov_remap_regions                          = null
params.vircov_remap_regions_coverage                 = null

// Viral consensus assembly (submodule)


// ======================
// K-mer profiling module
// ====================== 

params.kraken2_dbs                                   = null
params.kraken2_minimum_hit_groups                    = 3

params.bracken_read_length                           = 100
params.bracken_taxonomic_level                       = "S"
params.bracken_read_threshold                        = 3

// ==============================================
// Bacteria and Eukaryot subset alignment module
// ==============================================

// Bacterial and eukaryotic subset alignment

params.bacteria_min_cov                             = 0
params.bacteria_min_len                             = 50
params.bacteria_min_mapq                            = 60          
params.bacteria_min_regions                         = 0        
params.bacteria_min_reads                           = 0         

params.eukaryots_min_cov                            = 0
params.eukaryots_min_len                            = 50
params.eukaryots_min_mapq                           = 60        
params.eukaryots_min_regions                        = 0           
params.eukaryots_min_reads                          = 0       


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

    if (params.mode.dev.enabled){
        println "\n${c('red')}Development mode for utilities active, inputs and modules deactivated.${c('reset')}\n"
    } else {

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
                 println "\n${c('red')}Production mode is not active - you must provide either `--fastq_pe <GLOB>` or `--ont.enabled --fastq_ont <FILE-OR-GLOB>` for read inputs${c('reset')}\n"
            }
            
        }   

        if (params.mode.io.enabled){
            println "\n${c('red')}Development mode for inputs active, modules deactivated.${c('reset')}\n"
        } else {
            
            if (params.subsample.enabled) {
                reads = RasusaReadsMultiple(reads, params.subsample.reads, 1..params.subsample.replicates)
            }

            /*
            =================
            Workflow modules
            =================
            */

            started = java.time.LocalDateTime.now()

            if (params.mode.qc.enabled) {
                println "\n${c('red')}Quality control mode is active, other modules deactivated.${c('reset')}\n"

                // ===========================
                // Quality control subworkflow
                // ===========================

               
                if (params.ont.enabled) {
                    quality_control_ont(
                        reads, 
                        inputs.ercc_fasta, 
                        inputs.phage_fasta, 
                        inputs.host_depletion_dbs, 
                        inputs.host_depletion_references, 
                        params.qc.host.depletion.enabled,
                        params.qc.controls.phage.enabled
                    )
                } else {
                    quality_control_illumina(
                        reads, 
                        inputs.adapter_fasta, 
                        inputs.ercc_fasta, 
                        inputs.phage_fasta, 
                        inputs.host_depletion_dbs, 
                        inputs.host_depletion_references, 
                        params.qc.deduplication.enabled,
                        params.qc.deduplication.method,
                        params.qc.host.depletion.enabled,
                        params.qc.controls.phage.enabled
                    )
                }

            } else {
                
                if (params.validation.enabled) {

                    // =====================================
                    // Validation experiment workflows
                    // =====================================
                    
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

                        pathogen_detection(reads, inputs, params.ont.enabled)

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
    }
}

workflow.onComplete {
    complete_msg()
}


workflow.onError {
    println "\n${c('red')}Pipeline execution stopped with the following message: ${workflow.errorMessage}${c('reset')}\n"
}
