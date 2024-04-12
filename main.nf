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

===========
CEREBRO
===========

Version:    0.3.0
Inception:  2023-01-26
Maintainer: esteinig

*/

/*

Notes for this workflow:

    - Fastp deduplication is not deterministic due to low-bit hash algorithm and hash collisions, affects ~ 0.01% of reads: https://github.com/OpenGene/fastp
    
    - UMI-tools deduplication may also not be determistic despite specified seed value:  https://github.com/CGATOxford/UMI-tools/pull/550

    - UMI-tools with the current Neb UMI protocol against human reference is escalating RAM massively for some deep samples - 
      there may be interactions with the method (directional-adjacency), allowing multi-mapped read alignments, and generatign output stats (disabled)
      as described here: https://github.com/CGATOxford/UMI-tools/issues/561 where it seems that RAM usage escalates exponentially with increasing number 
      of unique UMIs per position (likely due to high read depth on human. This may be the case since the space of possible UMIs is 3^12 and we have
      high human coverage - this may still need to be tested properly.



*/

nextflow.enable.dsl=2

/* 
===================
Workflow parameters
===================
*/

params.help                                          = false
params.outdir                                        = "mgp_csf_v${workflow.manifest.version}"
params.databases                                     = 'db'
params.taxonomy                                      = 'db/ncbi_tax'

// File input

params.fastq                                         = null

params.fastq_illumina                                = null
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

// Viral scanning and remapping

params.vircov_scan_min_cov                           = 0  
params.vircov_scan_min_len                           = 50 
params.vircov_scan_min_mapq                          = 60                                                   // this applies to ungrouped alignment and filter reads that do not map uniquely to database references

params.vircov_scan_reads                             = 0  
params.vircov_scan_coverage                          = 0 
params.vircov_scan_regions                           = 0
params.vircov_scan_regions_coverage                  = 0

params.vircov_group_by                               = "taxid="
params.vircov_group_sep                              = ";"

params.vircov_group_select_by                        = "coverage"         
params.vircov_group_select_segment_field             = "segment="
params.vircov_group_select_segment_field_nan         = "segment=N/A"


params.vircov_remap_min_cov                          = 0  
params.vircov_remap_min_len                          = 50 
params.vircov_remap_min_mapq                         = 0  

params.vircov_remap_reads                            = 0
params.vircov_remap_coverage                         = 0
params.vircov_remap_regions                          = 0
params.vircov_remap_regions_coverage                 = 0

// Viral consensus assembly (submodule)

params.ivar_mpileup_args                             = "-d 50000"
params.ivar_min_qual                                 = 20
params.ivar_min_depth                                = 20
params.ivar_min_freq                                 = 0.75

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

params.bacteria_mash_index                          = "db/bacteria.msh"
params.bacteria_fasta                               = "db/bacteria.fasta"

params.eukaryots_mash_index                         = "db/eukaryots_complete.msh"
params.eukaryots_fasta                              = "db/eukaryots_complete.fasta"

params.subset_min_shared_hashes                     = 2 
params.subset_group_index                           = null     // null or int, recognized in process and switches to group selection [e.g 1 on seq index kraken fmt "taxid|1147|NC_014725.1"] -- picks the best genome by max shared hashes from taxid group
params.subset_group_sep                             = "|"      // format required for both bacterial and eukaryotic *|{taxid}|*

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


// ====================================================
// Metagenome and -transciptome assembly with alignment
// ====================================================

params.meta_blast_nt                                = null  // path of database directory without the 'nt' file/index specifier
params.meta_blast_nt_min_evalue                     = 0.000001
params.meta_blast_nt_min_identity                   = 95    // we only want high identity hits for identification, not looking for unusual things for now
params.meta_blast_nt_max_seqs                       = 100   // use the highest scoring hits ones for LCA computation

params.meta_diamond_nr                              = null  // path of database file created with diamond (including taxonomic information)
params.meta_diamond_nr_min_evalue                   = 0.000001
params.meta_diamond_nr_min_identity                 = 95    // we only want high identity hits for identification, not looking for unusual things for now
params.meta_diamond_nr_max_seqs                     = 100   // use the highest scoring hits ones for LCA computation

params.meta_diamond_nr_block_size                   = 2     // I think this improves runtime with higher memory (2(b+9b/c)MB where b = block size and c = index chunks)
params.meta_diamond_nr_index_chunks                 = 4

// ====================
// Other configurations
// ====================

params.monochrome                                   = false
params.minimap2_preset                              = "sr"       // Scrubby global settings for minimap2-preset 
params.covtobed_max_cov                             = 1000000    // Viral remapping if coverage really high (e.g. TWIST)



// Workflow imports

include { cerebro_production }              from './lib/workflows/subworkflows/production'
include { quality_control_illumina }        from './lib/workflows/subworkflows/quality_control'
include { pathogen_detection }              from './lib/workflows/pathogen_detection'
include { aneuploidy_detection_illumina }   from './lib/workflows/aneuploidy_detection'
include { culture_identification }          from './lib/workflows/culture_identification'

// Utility imports
include { PingServer} from './lib/processes/cerebro'
include { RasusaReadsMultiple } from './lib/processes/rasusa' addParams(subdir: "subsample")
include { c; parse_file_params; init_msg; complete_msg; help_msg; get_paired_reads; get_single_reads; WriteConfig } from './lib/utils'

workflow {

    // Startup message and help
    params.help ? help_msg() : init_msg()

    if (params.mode.dev.enabled){
        println "\n${c('red')}Development mode for utilities active, inputs and modules deactivated.${c('reset')}\n"
    } else {

        /*
        
        ================
        Input validation
        ================

        Necessary for input file checks and optional or paired input validation

        */
        
        def inputs = parse_file_params()

        if (params.production.enabled) {
            if (params.production.api.upload.enabled) {
                PingServer(Channel.of("PING")) | collect                             // collect to await result before proceeding
            }
            data = get_paired_reads(inputs.sample_sheet, true)                       // production requires sample sheet
            reads = data.pathogen                                                    // pathogen analysis uses all input reads
            reads_aneuploidy = data.aneuploidy                                       // host genome analysis if enabled uses subset of activated input reads
        } else {
            reads = get_paired_reads(inputs.sample_sheet, false)                     // pathogen analysis 
            reads_aneuploidy = reads                                                 // all input read files are also used for host genome analysis if enabled
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

                quality_control_illumina(
                    reads, 
                    inputs.adapter_fasta, 
                    inputs.ercc_fasta, 
                    inputs.phage_fasta, 
                    inputs.host_depletion_dbs, 
                    inputs.host_depletion_references, 
                    inputs.host_ercc_index, 
                    params.qc.deduplication.enabled,
                    params.qc.deduplication.method,
                    params.qc.host.depletion.enabled,
                    params.qc.controls.phage.enabled
                )  // uses only the pathogen analysis reads (host analysis is the same or a subset of pathogen analysis)

            } else {
                

                // =====================================
                // Host aneuploidy detection subworkflow
                // =====================================
                
                if (params.host.enabled && params.host.aneuploidy.enabled) {
                    aneuploidy = aneuploidy_detection(reads_aneuploidy, inputs)
                } else {
                    aneuploidy = Channel.empty()
                }

                // ==============================
                // Pathogen detection subworkflow
                // ==============================

                if (params.taxa.enabled) {
                    pathogen_detection(reads, inputs, false)

                    pathogen_detection.out.results | view

                    WriteConfig(
                        inputs.sample_sheet, 
                        pathogen_detection.out.results | collect, 
                        started
                    )

                    if (params.production.enabled) {
                        cerebro_production(
                            inputs.taxonomy,
                            pathogen_detection.out.results,
                            inputs.sample_sheet,
                            WriteConfig.out.config,
                        )
                    }

                }

                if ()
                
                // ===========================================
                // Cultured isolate identification subworkflow
                // ===========================================

                if (params.culture.enabled) {
                    ont_reads = get_single_reads(inputs.sample_sheet, false);
                    pe_reads = get_paired_reads(inputs.sample_sheet, false);
                    
                    culture_identification(ont_reads, pe_reads, inputs);
                    
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
