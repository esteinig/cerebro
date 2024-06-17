import groovy.json.JsonOutput
import java.nio.file.Files

/* 
=============
MSG FUNCTIONS
=============
*/


c_reset = params.monochrome ? '' : "\033[0m";
c_dim = params.monochrome ? '' : "\033[2m";
c_black = params.monochrome ? '' : "\033[0;30m";
c_green = params.monochrome ? '' : "\033[0;32m";
c_yellow = params.monochrome ? '' : "\033[0;33m";
c_blue = params.monochrome ? '' : "\033[0;34m";
c_purple = params.monochrome ? '' : "\033[0;35m";
c_cyan = params.monochrome ? '' : "\033[0;36m";
c_white = params.monochrome ? '' : "\033[0;37m";
c_red = params.monochrome ? '' : "\033[0;31m";

// 8-bit colors
c_indigo = params.monochrome ? '' : "\033[38;5;57m";
c_light_indigo = params.monochrome ? '' : "\033[38;5;63m";
c_light_blue = params.monochrome ? '' : "\033[38;5;33m";

// Get specific colors in other modules by 
// importing this function
def c(color) {
    if (color === "red"){
        return c_red
    } else {
        return c_reset
    }
}

def init_msg(){

    production_msg = params.production ? "${c_light_blue}Production mode active.${c_reset}" : ""

    log.info """
    ${c_indigo}=====================
    ${c_light_blue}M G P - C E R E B R O
    ${c_indigo}=====================${c_reset}

    ${workflow.manifest.description}

    Version:              ${c_light_indigo}v${workflow.manifest.version}${c_reset}
    Documentation:        ${c_light_blue}https://docs.meta-gp.org${c_reset}

    Workflow:             ${c_light_blue}${workflow.sessionId}${c_reset}
    Started:              ${c_light_indigo}${workflow.start}${c_reset}

    Profiles:             ${c_white}${workflow.profile}${c_reset}
    Working directory:    ${c_white}${workflow.workDir}${c_reset}

    ${production_msg}

    """.stripIndent()
}

def complete_msg(){

    log.info """

    ${c_indigo}=====================
    ${c_light_blue}M G P - C E R E B R O
    ${c_indigo}=====================${c_reset}

    Version:              ${c_light_indigo}v${workflow.manifest.version}${c_reset}
    Workflow:             ${c_light_blue}${workflow.sessionId}${c_reset}
    Completed:            ${c_light_indigo}${workflow.complete}${c_reset}

    ${c_indigo}=============================================================${c_reset}

    Please cite the following tools if used in the pipeline:

        - cerebro        0.7.0      https://github.com/esteinig/cerebro
        - umi-tools      1.1.4      https://github.com/CGATOxford/UMI-tools
        - calib          0.3.4      https://github.com/vpc-ccg/calib
        - covtobed       1.3.5      https://github.com/telatin/covtobed       
        - minimap2       2.24       https://github.com/lh3/minimap2          
        - samtools                  https://github.com/samtools/samtools      
        - kraken2        2.1.2      https://github.com/DerrickWood/kraken2    
        - fastp          0.23.2     https://github.com/OpenGene/fastp         
        - nextflow       22.10.4    https://github.com/nextflow-io/nextflow  
        - ivar           1.3.1      https://github.com/andersen-lab/ivar      
        - spades         3.15.5     https://github.com/ablab/spades           
        - strobealign    0.8.0      https://github.com/ksahlin/strobealign    
        - blast          2.13.0     https://github.com/ncbi                   
        - mash           2.3        https://github.com/marbl/Mash             
        - diamond        2.1.4      https://github.com/bbuchfink/diamond      
        - vircov         0.6.0      https://github.com/esteinig/vircov        
        - scrubby        0.3.0      https://github.com/esteinig/scrubby 
        - rasusa         0.7.1      https://github.com/mbhall88/rasusa
        - cnvkit         0.9.10     https://github.com/etal/cnvkit
        - nanoq          0.10.0     https://github.com/esteinig/nanoq

    Bibtex citations file can be found in the output directory for your convenience;
    please include citations of tools used in your workflow run when citing Cerebro
    in research publications.

    ${c_indigo}=============================================================${c_reset}

    ${c_light_blue}Cerebro${c_reset} depends on open-source libraries, including:
    
        - needletail                 https://github.com/onecodex/needletail
        - taxonomy                   https://github.com/onecodex/taxonomy 
        - niffler                    https://github.com/luizirber/niffler 
        - typst                      https://github.com/typst/typst

    ${c_indigo}=============================================================${c_reset}

    Cerebro is part of the Australian metagenomics 
    diagnostics consortium ${c_light_blue}META-GP${c_reset}. 

    For more information on accredited software for clinical  
    public health metagenomics, please see the documentation 
    repository at: 
    
    ${c_light_blue}https://docs.meta-gp.org${c_reset}
 
    ${c_indigo}=============================================================${c_reset}
    """.stripIndent()
}

// Help message

def help_msg(){

    log.info """
    ${c_indigo}=====================
    ${c_light_blue}M G P - C E R E B R O
    ${c_indigo}=====================${c_reset}

    Version:               ${c_light_indigo}v${workflow.mainfest.version}${c_reset}
    Documentation:        ${c_light_blue}https://docs.meta-gp.org${c_reset}

    """.stripIndent()
    System.exit(0)
}




/* 
================
CONFIG FUNCTIONS
================
*/

process WriteConfig {

    publishDir "$params.outdir", mode: "copy", pattern: "config.json"
    publishDir "$params.outdir", mode: "copy", pattern: "${workflow.sessionId}.csv", saveAs: { 'sample_sheet.csv' }

    input:
        path sample_sheet
        val  results
        val  started

    output:
        path 'config.json', emit: config
        path "${workflow.sessionId}.csv", emit: sample_sheet optional true 

    script:

        completed = java.time.LocalDateTime.now()

        config = [
            id: "$workflow.sessionId",
            name: "$workflow.runName",
            pipeline: "$workflow.manifest.name",
            version: "$workflow.manifest.version",
            started: "$started",
            completed: "$completed",
            params: params
        ]

        json = JsonOutput.toJson(config)
        json_pretty = JsonOutput.prettyPrint(json)

        if (sample_sheet) {
            """
            echo '${json_pretty}' > config.json
            cp $sample_sheet ${workflow.sessionId}.csv
            """
        } else {
            """
            echo '${json_pretty}' > config.json
            """
        }
        
}




/* 
===============
PARAM FUNCTIONS
===============
*/

def required_param_not_set_msg(param_name, param_value){
    println("\n${c_red}Required parameter `$param_name` not set ($param_value)${c_reset}\n")
    System.exit(1) 
}

/* Parse parameters that require file input and staging */
def parse_file_params(){
    
    // For files that we need staged into the processes, we use channels with collect 
    // and first methods to obtain `DataFlowVariables` - there may be a better way, 
    // but it's a little opaque to me what that might be.

    host_depletion_references = []; 
    host_depletion_dbs = [];
    if (params.qc.host.depletion.enabled) {
        if (params.qc.host.depletion.references && params.qc.host.depletion.databases) {
            host_depletion_references = Channel.fromPath(check_file_string(params.qc.host.depletion.references)).collect()
            host_depletion_dbs = Channel.fromPath(check_file_string(params.qc.host.depletion.databases)).collect()
        } else {
            println("\n${c_red}Host depletion is activated, but no references or database files provided (--qc.host.depletion.references | --qc.host.depletion.databases).${c_reset}\n")
            Thread.sleep(2000)
            System.exit(1) 
        }
    }


    virus_background_references = []; 
    virus_background_dbs = [];
    virus_db_index = []; 
    virus_db_fasta = [];

    if (params.taxa.enabled && params.taxa.alignment.enabled){
        if (params.virus_db_index && params.virus_db_fasta) {
            virus_db_index = Channel.fromPath(check_file(params.virus_db_index)).first()
            virus_db_fasta = Channel.fromPath(check_file(params.virus_db_fasta)).first() 
        } else {
            println("\n${c_red}Viral detection is activated, but no reference index or sequence file provides (--virus_db_index | --virus_db_fasta).${c_reset}\n")
            Thread.sleep(2000)
            System.exit(1) 
        }

        if (params.virus_background_references && params.virus_background_dbs) {
            virus_background_references = Channel.fromPath(check_file_string(params.virus_background_references)).collect()
            virus_background_dbs = Channel.fromPath(check_file_string(params.virus_background_dbs)).collect()
        } else {
            println("\n${c_red}Virus detection module is activated, but no background references or database files provided for depletion (--virus_background_dbs | --virus_background_references).${c_reset}\n")
            Thread.sleep(2000)
            System.exit(1) 
        }

    }
    

    ercc_fasta = [];
    if (params.qc.controls.ercc.enabled) {
        if (params.qc.controls.ercc.fasta) {
            ercc_fasta = Channel.fromPath(check_file(params.qc.controls.ercc.fasta)).first()
        } else {
            println("\n${c_red}ERCC control is activated, but no reference sequence file provided (--qc.controls.ercc.fasta).${c_reset}\n")
            Thread.sleep(2000)
            System.exit(1) 
        }
    } 

    phage_fasta = [];
    if (params.qc.controls.phage.enabled) {
        if (params.qc.controls.phage.fasta) {
            phage_fasta = Channel.fromPath(check_file(params.qc.controls.phage.fasta)).first()
        } else {
            println("\n${c_red}Phage control is activated, but no reference sequence file provided (--qc.controls.phage.fasta).${c_reset}\n")
            Thread.sleep(2000)
            System.exit(1) 
        }
    } 

    kraken2_dbs = [];
    if (params.taxa.enabled && params.taxa.kmer.enabled && params.taxa.kmer.kraken2.enabled){
        if (params.kraken2_dbs) {
            kraken2_dbs = Channel.fromPath(check_file_string(params.kraken2_dbs)).collect() 
        } else {
            println("\n${c_red}K-mer profiling is activated, but no reference databases provided (--kraken2_dbs).${c_reset}\n")
            Thread.sleep(2000)
            System.exit(1) 
        }
    }
    

    eukaryots_mash_index = [];
    eukaryots_fasta = [];
    if (params.taxa.enabled && params.taxa.alignment.enabled && params.taxa.alignment.eukaryots) {
        if (params.eukaryots_mash_index && params.eukaryots_fasta) {
            eukaryots_mash_index = Channel.fromPath(check_file(params.eukaryots_mash_index)).first()
            eukaryots_fasta = Channel.fromPath(check_file(params.eukaryots_fasta)).first()
            check_matching_reference_data(Channel.fromPath(params.eukaryots_mash_index), Channel.fromPath(params.eukaryots_fasta), "eukaryotic references")
        } else {
            println("\n${c_red}Eukaryotic alignment is activated but no references files were provided (--eukaryots_mash_index | --eukaryots_fasta).${c_reset}\n")
            Thread.sleep(2000)
            System.exit(1) 
        }
    }
    

    bacteria_mash_index = [];
    bacteria_fasta = [];
    if (params.taxa.enabled && params.taxa.alignment.enabled && params.taxa.alignment.bacteria) {
        if (params.bacteria_mash_index && params.bacteria_fasta) {
            bacteria_mash_index = Channel.fromPath(check_file(params.bacteria_mash_index)).first()
            bacteria_fasta = Channel.fromPath(check_file(params.bacteria_fasta)).first()
            check_matching_reference_data(Channel.fromPath(params.bacteria_mash_index), Channel.fromPath(params.bacteria_fasta), "bacterial references")
        } else {
            println("\n${c_red}Bacterial alignment is activated but no references files were provided (--bacteria_mash_index | --bacteria_fasta).${c_reset}\n")
            Thread.sleep(2000)
            System.exit(1) 
        }
    }

    meta_diamond_nr = [];
    meta_blast_nt = [];
    if (params.taxa.enabled && params.taxa.assembly.enabled && params.taxa.assembly.meta.enabled) {
        if (params.meta_blast_nt) {
            meta_blast_nt = Channel.fromPath(check_file(params.meta_blast_nt)).first()
        } else {
            println("\n${c_red}BLASTN for meta-assembly is activated, but no reference database path provided (--meta_blast_nt).${c_reset}\n")
            Thread.sleep(2000)
            System.exit(1) 
        }
        if (params.meta_diamond_nr) {
            meta_diamond_nr = Channel.fromPath(check_file(params.meta_diamond_nr)).first()
        } else {
            println("\n${c_red}DIAMOND for meta-assembly is activated, but no reference database path provided (--meta_diamond_nr).${c_reset}\n")
            Thread.sleep(2000)
            System.exit(1) 
        }
    }
    
    // Here we only need a file, as we create the channel in the sample sheet parsing function

    sample_sheet = [];
    if (params.production.enabled) {
        if (params.production.sample_sheet) {
            sample_sheet = file(check_file(params.production.sample_sheet))
        } else {
            println("\n${c_red}Production settings are activated, but no sample sheet was provided (--production.sample_sheet).${c_reset}\n")
            Thread.sleep(2000); // we need a small delay so the message is printed before we exit with error code
            System.exit(1) 
        }
    }

    // For production also check that the `nodes.dmp` and `names.dmp` files for taxonomy are available in the database directory

    taxonomy_directory = [];
    if (params.production.enabled || (params.process.enabled && params.process.taxa)) {
        if (params.taxonomy) {
            check_file("$params.taxonomy/nodes.dmp")
            check_file("$params.taxonomy/names.dmp")
            taxonomy_directory = Channel.fromPath(check_file(params.taxonomy)).first()
        } else {
            println("\n${c_red}Settings are activated that require a taxonomy (production or post-processing), but no taxonomy directory (containing `nodes.dmp` and `names.dmp` files) was provided (--database.taxonomy.directory).${c_reset}\n")
            Thread.sleep(2000); // we need a small delay so the message is printed before we exit with error code
            System.exit(1) 
        }
    }

    // Aneuploidy detection host analysis


    aneuploidy_reference_index = [];
    if (params.host.enabled && params.host.aneuploidy.enabled) {
        if (params.host.aneuploidy.reference_index) {
            aneuploidy_reference_index = Channel.fromPath(check_file(params.host.aneuploidy.reference_index)).first()
        } else {
            println("\n${c_red}Aneuploidy detection is activated, but no host reference index for alignment was provided (--host.aneuploidy.reference_index).${c_reset}\n")
            Thread.sleep(2000); // we need a small delay so the message is printed before we exit with error code
            System.exit(1) 
        }
    }
    
    aneuploidy_controls = [];
    if (params.host.enabled && params.host.aneuploidy.enabled && params.host.aneuploidy.cnvkit.enabled) {
        if (params.host.aneuploidy.cnvkit.normal_control) {
            aneuploidy_controls = Channel.fromPath(check_file_string(params.host.aneuploidy.cnvkit.normal_control)).collect()
        } else {
            println("\n${c_red}Aneuploidy detection is activated, but no normal control alignments were provided (--host.aneuploidy.normal_control).${c_reset}\n")
            Thread.sleep(2000); // we need a small delay so the message is printed before we exit with error code
            System.exit(1) 
        }
    }

    aneuploidy_reference_fasta = [];
    if (params.host.enabled && params.host.aneuploidy.enabled && params.host.aneuploidy.cnvkit.enabled) {
        if (params.host.aneuploidy.cnvkit.reference_fasta) {
            aneuploidy_reference_fasta = Channel.fromPath(check_file(params.host.aneuploidy.cnvkit.reference_fasta)).first()
        } else {
            println("\n${c_red}Aneuploidy detection is activated, but no host reference sequence (.fasta) was provided (--host.aneuploidy.reference_fasta).${c_reset}\n")
            Thread.sleep(2000); // we need a small delay so the message is printed before we exit with error code
            System.exit(1) 
        }
    }

    // Entirely optional parameters

    virus_blacklist = params.virus_blacklist ? Channel.fromPath(check_file(params.virus_blacklist)).first() : []        
    adapter_fasta   = params.qc.reads.fastp.adapter_fasta ? Channel.fromPath(check_file(params.qc.reads.fastp.adapter_fasta)).first() : []    

    return [
        host_depletion_references: host_depletion_references,
        host_depletion_dbs: host_depletion_dbs, 
        virus_background_references: virus_background_references, 
        virus_background_dbs: virus_background_dbs, 
        virus_db_index: virus_db_index, 
        virus_db_fasta: virus_db_fasta, 
        virus_blacklist: virus_blacklist, 
        adapter_fasta: adapter_fasta, 
        ercc_fasta: ercc_fasta, 
        phage_fasta: phage_fasta,
        kraken2_dbs: kraken2_dbs, 
        bacteria_mash_index: bacteria_mash_index, 
        bacteria_fasta: bacteria_fasta, 
        eukaryots_mash_index: eukaryots_mash_index, 
        eukaryots_fasta: eukaryots_fasta, 
        meta_blast_nt: meta_blast_nt,
        meta_diamond_nr: meta_diamond_nr, 
        sample_sheet: sample_sheet,
        taxonomy_directory: taxonomy_directory,
        aneuploidy_reference_index: aneuploidy_reference_index,
        aneuploidy_reference_fasta: aneuploidy_reference_fasta,
        aneuploidy_controls: aneuploidy_controls
    ]
}


/* 
================
HELPER FUNCTIONS
================
*/


// Helper function to check that matching bacterial index (e.g. Mash or other) 
// and reference sequence files for database subsetting exist.
def check_matching_reference_data(index, fasta, msg) {
    index.concat(fasta) 
    | map { tuple(it.baseName, it) } 
    | groupTuple 
    | map { it[1] } 
    | filter { it.size() == 2 } 
    | ifEmpty{ exit 1, "\n${c_red}Could not detect matching index and fasta files ($msg must match by file base name).${c_reset}\n" }
}

// Helper function to check if file exists and 
// return the Nextflow file(...) for staging
def check_file(file_path) {
    def fp = new File(file_path)
    if (!fp.exists()){
        println("\n${c_red}File path not found: ${file_path} ${c_reset}\n")
        Thread.sleep(2000);
        System.exit(1)
    }
    return file_path
}

// Helper function to check a white-space separated 
// string of paths, see if file exists and return 
// it as an array of file(...) for staging
def check_file_string(file_path_string) {
    def file_path_strings = file_path_string.tokenize(',');
    def file_paths = file_path_strings.each {  
        def fp = new File(it)
        if (!fp.exists()){
            println("\n${c_red}File path not found: ${fp} ${c_reset}\n")
            Thread.sleep(2000); 
            System.exit(1)
        }
        it
    }
    return file_paths
}

/* 
===========
FILE INPUTS
===========
*/

def get_paired_reads(sample_sheet, production) {

    if (params.production.sample_sheet){
        reads = from_sample_sheet_illumina(sample_sheet, production)
    } else if (params.fastq_pe) {
        reads = channel.fromFilePairs(params.fastq_pe, flat: true, checkIfExists: true)
    } else {
        println "\n${c('red')}Either `--fastq_pe` or `--production.sample_sheet` have to be specified for paired read input${c('reset')}\n"
        Thread.sleep(2000); 
        System.exit(1)
    }
    return reads
}

def get_single_reads(sample_sheet, production) {

    if (params.production.sample_sheet){
        reads = from_sample_sheet_ont(sample_sheet, production)
    } else if (params.fastq_ont) {
        reads = channel.fromPath(params.fastq_ont, checkIfExists: true) | map { tuple(it.getSimpleName(), it) } 
    } else {
        println "\n${c('red')}Either `--fastq_ont` or `--production.sample_sheet` have to be specified for nanopore read input${c('reset')}\n"
        Thread.sleep(2000); 
        System.exit(1)
    }
    return reads
}

def from_sample_sheet_illumina(file, production){

    def row_number = 2 // with header

    fastq_files = channel.fromPath("$file") | splitCsv(header:true, strip:true) | map { row -> 

        // Check that required columns are present
        if (row.sample_id === null || row.forward_path === null || row.reverse_path === null) {
            println "\n${c('red')}Sample sheet did not contain required columns (sample_id, forward_path, reverse_path)${c('reset')}\n"
            Thread.sleep(2000); 
            System.exit(1)
        }

        // Check that required values for this row are set
        if (row.sample_id.isEmpty() || row.forward_path.isEmpty() || row.reverse_path.isEmpty() ) {
            println "\n${c('red')}Sample sheet did not contain required values (sample_id, forward_path, reverse_path) in row: ${row_number}${c('reset')}\n"
            Thread.sleep(2000); 
            System.exit(1)
        }

        forward = new File(row.forward_path)
        reverse = new File(row.reverse_path)

        if (!forward.exists()){
            println("Forward read file does not exist: ${forward}")
            return 
        }
        if (!reverse.exists()){
            println("Reverse read file does not exist: ${reverse}")
            return 
        }

        if (production) {
             // Check that required columns are present
            if (row.run_id === null || row.run_date === null) {
                println "\n${c('red')}Production mode is activated, but the sample sheet did not contain required columns (run_date, run_id)${c('reset')}\n"
                Thread.sleep(2000); 
                System.exit(1)
            }
            // Check that run date and identifier are set for production models
            if (row.run_id.isEmpty() || row.run_date.isEmpty()) {
                println "\n${c('red')}Production mode is activated, but the sample sheet did not contain a sequence run identifier or date for sample: ${row.sample_id}${c('reset')}\n"
                Thread.sleep(2000); 
                System.exit(1)
            }
            // Check if the aneuploidy column is present
            if (row.aneuploidy === null || row.aneuploidy.isEmpty()) {
                println "\n${c('red')}Production mode is activated, but the sample sheet did not specify the aneuploidy detection/consent column for sample: ${row.sample_id}${c('reset')}\n"
                Thread.sleep(2000); 
                System.exit(1)
            }
        }

        row_number++

        return tuple(row.sample_id, row.forward_path, row.reverse_path, row.aneuploidy.toBoolean())

    }

    fastq_files | ifEmpty { exit 1, "\n${c_red}Could not find read files specified in sample sheet.${c_reset}\n" }

    return [
        aneuploidy: fastq_files | filter { it[3] } | map { tuple(it[0], it[1], it[2]) },
        pathogen: fastq_files | map { tuple(it[0], it[1], it[2]) }
    ]

}


def from_sample_sheet_ont(file, production){

    def row_number = 2 // with header

    fastq_files = channel.fromPath("$file") | splitCsv(header:true, strip:true) | map { row -> 

        // Check that required columns are present
        if (row.sample_id === null || row.fastq == null) {
            println "\n${c('red')}Sample sheet did not contain required columns (sample_id, forward_path, reverse_path)${c('reset')}\n"
            Thread.sleep(2000); 
            System.exit(1)
        }

        // Check that required values for this row are set
        if (row.sample_id.isEmpty() || row.fastq.isEmpty()) {
            println "\n${c('red')}Sample sheet did not contain required values (sample_id, forward_path, reverse_path) in row: ${row_number}${c('reset')}\n"
            Thread.sleep(2000); 
            System.exit(1)
        }

        fastq = new File(row.fastq)

        if (!fastq.exists()){
            println("Fastq read file does not exist: ${forward}")
            return 
        }

        if (production) {
             // Check that required columns are present
            if (row.run_id === null || row.run_date === null) {
                println "\n${c('red')}Production mode is activated, but the sample sheet did not contain required columns (run_date, run_id)${c('reset')}\n"
                Thread.sleep(2000); 
                System.exit(1)
            }
            // Check that run date and identifier are set for production models
            if (row.run_id.isEmpty() || row.run_date.isEmpty()) {
                println "\n${c('red')}Production mode is activated, but the sample sheet did not contain a sequence run identifier or date for sample: ${row.sample_id}${c('reset')}\n"
                Thread.sleep(2000); 
                System.exit(1)
            }
            // Check if the aneuploidy column is present
            if (row.aneuploidy === null || row.aneuploidy.isEmpty()) {
                println "\n${c('red')}Production mode is activated, but the sample sheet did not specify the aneuploidy detection/consent column for sample: ${row.sample_id}${c('reset')}\n"
                Thread.sleep(2000); 
                System.exit(1)
            }
        }

        row_number++

        return tuple(row.sample_id, row.fastq, row.aneuploidy.toBoolean())

    }

    fastq_files | ifEmpty { exit 1, "\n${c_red}Could not find read files specified in sample sheet.${c_reset}\n" }

    return [
        aneuploidy: fastq_files | filter { it[3] } | map { tuple(it[0], it[1]) },
        pathogen: fastq_files | map { tuple(it[0], it[1]) }
    ]

}

def from_sample_sheet_reference_alignment_illumina(file){

    def row_number = 2 // with header

    fastq_files = channel.fromPath("$file") | splitCsv(header:true, strip:true) | map { row -> 

        // Check that required columns are present
        if (row.sample_id === null || row.forward_path === null || row.reverse_path === null || row.reference_path === null) {
            println "\n${c('red')}Sample sheet did not contain required columns (sample_id, forward_path, reverse_path, reference_path)${c('reset')}\n"
            Thread.sleep(2000); 
            System.exit(1)
        }

        // Check that required values for this row are set
        if (row.sample_id.isEmpty() || row.forward_path.isEmpty() || row.reverse_path.isEmpty() || row.reference_path.isEmpty()) {
            println "\n${c('red')}Sample sheet did not contain required values (sample_id, forward_path, reverse_path, reference_path) in row: ${row_number}${c('reset')}\n"
            Thread.sleep(2000); 
            System.exit(1)
        }

        forward = new File(row.forward_path)
        reverse = new File(row.reverse_path)
        reference = new File(row.reference_path)

        if (!forward.exists()){
            println("Forward read file does not exist: ${forward}")
            return 
        }
        if (!reverse.exists()){
            println("Reverse read file does not exist: ${reverse}")
            return 
        }
        if (!reference.exists()){
            println("Reference genome file does not exist: ${reference}")
            return 
        }

        row_number++

        return tuple(row.sample_id, row.forward_path, row.reverse_path, row.reference_path)

    }

    fastq_files | ifEmpty { exit 1, "\n${c_red}Could not find read files specified in sample sheet.${c_reset}\n" }

    return fastq_files

}


def from_sample_sheet_reference_alignment_ont(file){

    def row_number = 2 // with header

    fastq_files = channel.fromPath("$file") | splitCsv(header:true, strip:true) | map { row -> 

        // Check that required columns are present
        if (row.sample_id === null || row.fastq === null || row.reference_path === null) {
            println "\n${c('red')}Sample sheet did not contain required columns (sample_id, fastq, reference_path)${c('reset')}\n"
            Thread.sleep(2000); 
            System.exit(1)
        }

        // Check that required values for this row are set
        if (row.sample_id.isEmpty() || row.fastq.isEmpty() || row.reference_path.isEmpty()) {
            println "\n${c('red')}Sample sheet did not contain required values (sample_id, fastq, reference_path) in row: ${row_number}${c('reset')}\n"
            Thread.sleep(2000); 
            System.exit(1)
        }

        fastq = new File(row.fastq)

        if (!fastq.exists()){
            println("Forward read file does not exist: ${forward}")
            return 
        }
        if (!reference.exists()){
            println("Reference genome file does not exist: ${reference}")
            return 
        }

        row_number++

        return tuple(row.sample_id, row.fastq, row.reference_path)

    }

    fastq_files | ifEmpty { exit 1, "\n${c_red}Could not find read files specified in sample sheet.${c_reset}\n" }

    return fastq_files

}