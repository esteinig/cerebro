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

        - cerebro        1.0.0      https://github.com/esteinig/cerebro      
        - minimap2       2.24       https://github.com/lh3/minimap2         
        - bowtie2        2.24       https://github.com/        
        - samtools                  https://github.com/samtools/samtools      
        - kraken2        2.1.2      https://github.com/DerrickWood/kraken2    
        - bracken        3.0.0      https://github.com/
        - fastp          0.23.2     https://github.com/OpenGene/fastp         
        - nextflow       24.04      https://github.com/nextflow-io/nextflow  
        - ivar           1.3.1      https://github.com/andersen-lab/ivar    
        - megahit        2.1.3      https://github.com/ksahlin/megahit     
        - spades         4.0.0      https://github.com/ablab/spades           
        - strobealign    0.13.0     https://github.com/ksahlin/strobealign    
        - blast          2.13.0     https://github.com/ncbi                   
        - mash           2.3        https://github.com/marbl/Mash             
        - diamond        2.1.4      https://github.com/bbuchfink/diamond      
        - vircov         1.0.0      https://github.com/esteinig/vircov        
        - scrubby        1.0.0      https://github.com/esteinig/scrubby 
        - rasusa         2.0.0      https://github.com/mbhall88/rasusa
        - cnvkit         0.9.10     https://github.com/etal/cnvkit
        - nanoq          0.10.0     https://github.com/esteinig/nanoq

    Bibtex citations file can be found in the output directory for your convenience;
    please include citations of tools used in your workflow when citing Cerebro in
    research publications.

    ${c_indigo}=============================================================${c_reset}

    ${c_light_blue}Cerebro${c_reset} depends on open-source libraries, including:
    
        - needletail                 https://github.com/onecodex/needletail
        - taxonomy                   https://github.com/onecodex/taxonomy 
        - niffler                    https://github.com/luizirber/niffler 
        - typst                      https://github.com/typst/typst

    ${c_indigo}=============================================================${c_reset}

    Cerebro is part of the Australian metagenomics diagnostics consortium ${c_light_blue}META-GP${c_reset}. 

    For more information on accredited software for clinical public health metagenomics, 
    please see the documentation: 
    
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
    error "Required parameter `$param_name` not set ($param_value)"
}

/* Parse parameters that require file input and staging */
def parse_file_params(){
    
    // Here we only need a file, as we create the channel in the sample sheet parsing function

    sample_sheet = [];
    if (params.production.enabled) {
        if (params.production.sample_sheet) {
            sample_sheet = file(check_file(params.production.sample_sheet))
        } else {
            error "Production settings are activated, but no sample sheet was provided (--production.sample_sheet)."
        }
    }

    // For production also check that the `nodes.dmp` and `names.dmp` files for taxonomy are available in the database directory

    taxonomy_directory = [];
    if (params.process.enabled && params.process.taxa) {
        if (params.database.taxonomy) {
            check_file("$params.database.taxonomy/nodes.dmp")
            check_file("$params.database.taxonomy/names.dmp")
            taxonomy_directory = Channel.fromPath(check_file(params.taxonomy)).first()
        } else {
            error "Settings are activated that require a taxonomy (production or post-processing), but no taxonomy directory (NCBI-style) was provided (--database.taxonomy.directory)"
        }
    }

    // Quality control inputs

    host_depletion_reference = []; 
    host_depletion_kraken2 = [];
    if (params.qc.host.depletion.enabled) {
        if (params.qc.host.depletion.reference && params.qc.host.depletion.kraken2) {
            host_depletion_reference = Channel.fromPath(check_file_string(params.qc.host.depletion.reference)).collect()
            host_depletion_kraken2 = Channel.fromPath(check_file_string(params.qc.host.depletion.kraken2)).collect()
        } else {
            error "Host depletion is activated, but no FASTA reference/alignment index or database files were provided (--qc.host.depletion.fasta | --qc.host.depletion.kraken2)"
        }
    }

    background_depletion_reference = []; 
    background_depletion_kraken2 = [];
    if (params.qc.background.depletion.enabled) {
        if (params.qc.background.depletion.reference && params.qc.background.depletion.kraken2) {
            background_depletion_reference = Channel.fromPath(check_file_string(params.qc.background.depletion.reference)).collect()
            background_depletion_kraken2 = Channel.fromPath(check_file_string(params.qc.background.depletion.kraken2)).collect()
        } else {
            error "Background depletion is activated, but no FASTA reference/alignment index or Kraken2 database files were provided (--qc.background.depletion.fasta | --qc.background.depletion.kraken2)"
        }
    }

    ercc_fasta = [];
    if (params.qc.controls.ercc.enabled) {
        if (params.qc.controls.ercc.fasta) {
            ercc_fasta = Channel.fromPath(check_file(params.qc.controls.ercc.fasta)).first()
        } else {
            error "ERCC control is activated, but no FASTA reference/alignment index was provided (--qc.controls.ercc.fasta)."
        }
    } 

    phage_fasta = [];
    if (params.qc.controls.phage.enabled) {
        if (params.qc.controls.phage.fasta) {
            phage_fasta = Channel.fromPath(check_file(params.qc.controls.phage.fasta)).first()
        } else {
            error "Phage control is activated, but no FASTA reference/alignment index was provided (--qc.controls.phage.fasta)."
        }
    } 

    // Taxonomic classification

    alignment_index = []; 
    alignment_fasta = [];
    if (params.taxa.enabled && params.taxa.alignment.enabled){
        if (params.databases.alignment.index && params.databases.alignment.fasta) {
            alignment_index = Channel.fromPath(check_file(params.databases.alignment.index)).first()
            alignment_fasta = Channel.fromPath(check_file(params.databases.alignment.fasta)).first() 
        } else {
            error "Alignment pathogen detection is activated, but no reference index or sequence file were provided (--databases.alignment.index | --databases.alignment.fasta)."
        }
    }

    kraken2_dbs = [];
    if (params.taxa.enabled && params.taxa.kmer.enabled && params.taxa.kmer.kraken2.enabled){
        if (params.databases.kmer.kraken2) {
            kraken2_dbs = Channel.fromPath(check_file_string(params.databases.kmer.kraken2)).collect() 
        } else {
            error "K-mer profiling with Kraken2 is activated, but no reference databases were provided (--databases.kmer.kraken2)."
        }
    }

    metabuli_dbs = [];
    if (params.taxa.enabled && params.taxa.kmer.enabled && params.taxa.kmer.metabuli.enabled){
        if (params.databases.kmer.metabuli) {
            metabuli_dbs = Channel.fromPath(check_file_string(params.databases.kmer.metabuli)).collect() 
        } else {
            error "K-mer profiling with Metabuli is activated, but no reference databases were provided (--databases.kmer.metabuli)."
        }
    }

    meta_diamond_nr = [];
    meta_blast_nt = [];
    if (params.taxa.enabled && params.taxa.assembly.enabled && params.taxa.assembly.meta.enabled) {
        if (params.databases.assembly.blast) {
            meta_blast_nt = Channel.fromPath(check_file(params.databases.assembly.blast)).first()
        } else {
            error "BLAST nucleotide classification for metagenome assembly is activated, but no reference database was provided (--databases.assembly.blast)."
        }
        if (params.databases.assembly.diamond) {
            meta_diamond_nr = Channel.fromPath(check_file(params.databases.assembly.diamond)).first()
        } else {
            error "DIAMOND protein classification for metagenome assembly is activated, but no reference database was provided (--databases.assembly.diamond)."
        }
    }
    
    // Panviral background depletion and alignment analysis


    panviral_background_depletion_reference = []; 
    panviral_background_depletion_kraken2 = [];
    if (params.panviral.enabled) {
        if (params.databases.panviral.background.reference && params.databases.panviral.background.kraken2) {
            panviral_background_depletion_reference = Channel.fromPath(check_file_string(params.databases.panviral.background.reference)).collect()
            panviral_background_depletion_kraken2 = Channel.fromPath(check_file_string(params.databases.panviral.background.kraken2)).collect()
        } else {
            error "Panviral background depletion is activated, but no FASTA reference/alignment index or Kraken2 database files were provided (--databases.panviral.background.reference & --databases.panviral.background.kraken2)"
        }
    }

    panviral_alignment_index = []; 
    panviral_alignment_fasta = [];
    if (params.panviral.enabled){
        if (params.databases.alignment.index && params.databases.alignment.fasta) {
            panviral_alignment_index = Channel.fromPath(check_file(params.databases.panviral.index)).first()
            panviral_alignment_fasta = Channel.fromPath(check_file(params.databases.panviral.fasta)).first() 
        } else {
            error "Alignment pathogen detection is activated, but no reference index or sequence file were provided (--databases.panviral.index | --databases.panviral.fasta)."
        }
    }


    panviral_blacklist = params.databases.panviral.blacklist ? Channel.fromPath(check_file(params.databases.panviral.blacklist)).first() : []   

    // Aneuploidy detection host analysis


    aneuploidy_reference_index = [];
    if (params.host.enabled && params.host.aneuploidy.enabled) {
        if (params.host.aneuploidy.reference_index) {
            aneuploidy_reference_index = Channel.fromPath(check_file(params.host.aneuploidy.reference_index)).first()
        } else {
            error "Aneuploidy detection is activated, but no host reference index was provided (--host.aneuploidy.reference_index)."
        }
    }
    
    aneuploidy_controls = [];
    if (params.host.enabled && params.host.aneuploidy.enabled && params.host.aneuploidy.cnvkit.enabled) {
        if (params.host.aneuploidy.cnvkit.normal_control) {
            aneuploidy_controls = Channel.fromPath(check_file_string(params.host.aneuploidy.cnvkit.normal_control)).collect()
        } else {
            error "Aneuploidy detection is activated, but no normal control alignments were provided (--host.aneuploidy.normal_control)."
        }
    }

    aneuploidy_reference_fasta = [];
    if (params.host.enabled && params.host.aneuploidy.enabled && params.host.aneuploidy.cnvkit.enabled) {
        if (params.host.aneuploidy.cnvkit.reference_fasta) {
            aneuploidy_reference_fasta = Channel.fromPath(check_file(params.host.aneuploidy.cnvkit.reference_fasta)).first()
        } else {
            error "Aneuploidy detection is activated, but no host reference sequence (.fasta) was provided (--host.aneuploidy.reference_fasta)."
        }
    }

    // Entirely optional parameters

    virus_blacklist = params.virus_blacklist ? Channel.fromPath(check_file(params.virus_blacklist)).first() : []        
    adapter_fasta   = params.qc.reads.fastp.adapter_fasta ? Channel.fromPath(check_file(params.qc.reads.fastp.adapter_fasta)).first() : []    

    return [
        sample_sheet: sample_sheet,
        taxonomy_directory: taxonomy_directory, 
        adapter_fasta: adapter_fasta, 
        ercc_fasta: ercc_fasta, 
        phage_fasta: phage_fasta,
        host_depletion_reference: host_depletion_reference,
        host_depletion_kraken2: host_depletion_kraken2, 
        background_depletion_reference: background_depletion_reference,
        background_depletion_kraken2: background_depletion_kraken2,
        alignment_index: alignment_index,
        alignment_fasta: alignment_fasta,
        kraken2_dbs: kraken2_dbs, 
        metabuli_dbs: metabuli_dbs, 
        meta_blast_nt: meta_blast_nt,
        meta_diamond_nr: meta_diamond_nr, 
        aneuploidy_reference_index: aneuploidy_reference_index,
        aneuploidy_reference_fasta: aneuploidy_reference_fasta,
        aneuploidy_controls: aneuploidy_controls,
        panviral_blacklist: panviral_blacklist
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
    | ifEmpty{ exit 1, "Could not detect matching index and fasta files ($msg must match by file base name)." }
}

// Helper function to check if file exists and 
// return the Nextflow file(...) for staging
def check_file(file_path) {
    def fp = new File(file_path)
    if (!fp.exists()){
        error "File path not found: ${file_path} "
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
            error "File path not found: ${fp} "
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
        error "Either `--fastq_pe` or `--production.sample_sheet` have to be specified for paired read input"
    }
    return reads
}

def get_single_reads(sample_sheet, production) {

    if (params.production.sample_sheet){
        reads = from_sample_sheet_ont(sample_sheet, production)
    } else if (params.fastq_ont) {
        reads = channel.fromPath(params.fastq_ont, checkIfExists: true) | map { tuple(it.getSimpleName(), it) } 
    } else {
        error "Either `--fastq_ont` or `--production.sample_sheet` have to be specified for nanopore read input"
    }
    return reads
}

def from_sample_sheet_illumina(file, production){

    def row_number = 2 // with header

    fastq_files = channel.fromPath("$file") | splitCsv(header:true, strip:true) | map { row -> 

        // Check that required columns are present
        if (row.sample_id === null || row.forward_path === null || row.reverse_path === null) {
            error "Sample sheet did not contain required columns (sample_id, forward_path, reverse_path)"
        }

        // Check that required values for this row are set
        if (row.sample_id.isEmpty() || row.forward_path.isEmpty() || row.reverse_path.isEmpty() ) {
            error "Sample sheet did not contain required values (sample_id, forward_path, reverse_path) in row: ${row_number}"
        }

        forward = new File(row.forward_path)
        reverse = new File(row.reverse_path)

        if (!forward.exists()){
            error "Forward read file does not exist: ${forward}"
        }
        if (!reverse.exists()){
            error "Reverse read file does not exist: ${reverse}"
        }

        if (production) {
             // Check that required columns are present
            if (row.run_id === null || row.run_date === null) {
                error "Production mode is activated, but the sample sheet did not contain required columns (run_date, run_id)"
            }
            // Check that run date and identifier are set for production models
            if (row.run_id.isEmpty() || row.run_date.isEmpty()) {
                error "Production mode is activated, but the sample sheet did not contain a sequence run identifier or date for sample: ${row.sample_id}"
            }
            // Check if the aneuploidy column is present
            if (row.aneuploidy === null || row.aneuploidy.isEmpty()) {
                error "Production mode is activated, but the sample sheet did not specify the aneuploidy detection/consent column for sample: ${row.sample_id}"
            }
        }

        row_number++

        return tuple(row.sample_id, row.forward_path, row.reverse_path, row.aneuploidy.toBoolean())

    }

    fastq_files | ifEmpty { exit 1, "Could not find read files specified in sample sheet." }

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
            error "Sample sheet did not contain required columns (sample_id, forward_path, reverse_path)"
        }

        // Check that required values for this row are set
        if (row.sample_id.isEmpty() || row.fastq.isEmpty()) {
            error "Sample sheet did not contain required values (sample_id, forward_path, reverse_path) in row: ${row_number}"
        }

        fastq = new File(row.fastq)

        if (!fastq.exists()){
            error "Fastq read file does not exist: ${forward}"
        }

        if (production) {
             // Check that required columns are present
            if (row.run_id === null || row.run_date === null) {
                error "Production mode is activated, but the sample sheet did not contain required columns (run_date, run_id)"
            }
            // Check that run date and identifier are set for production models
            if (row.run_id.isEmpty() || row.run_date.isEmpty()) {
                error "Production mode is activated, but the sample sheet did not contain a sequence run identifier or date for sample: ${row.sample_id}"
            }
            // Check if the aneuploidy column is present
            if (row.aneuploidy === null || row.aneuploidy.isEmpty()) {
                error "Production mode is activated, but the sample sheet did not specify the aneuploidy detection/consent column for sample: ${row.sample_id}"
            }
        }

        row_number++

        return tuple(row.sample_id, row.fastq, row.aneuploidy.toBoolean())

    }

    fastq_files | ifEmpty { exit 1, "Could not find read files specified in sample sheet." }

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
            error "Sample sheet did not contain required columns (sample_id, forward_path, reverse_path, reference_path)"
        }

        // Check that required values for this row are set
        if (row.sample_id.isEmpty() || row.forward_path.isEmpty() || row.reverse_path.isEmpty() || row.reference_path.isEmpty()) {
            error "Sample sheet did not contain required values (sample_id, forward_path, reverse_path, reference_path) in row: ${row_number}"
        }

        forward = new File(row.forward_path)
        reverse = new File(row.reverse_path)
        reference = new File(row.reference_path)

        if (!forward.exists()){
            error "Forward read file does not exist: ${forward}"
        }
        if (!reverse.exists()){
            error "Reverse read file does not exist: ${reverse}"
        }
        if (!reference.exists()){
            error "Reference genome file does not exist: ${reference}"
        }

        row_number++

        return tuple(row.sample_id, row.forward_path, row.reverse_path, row.reference_path)

    }

    fastq_files | ifEmpty { exit 1, "Could not find read files specified in sample sheet." }

    return fastq_files

}


def from_sample_sheet_reference_alignment_ont(file){

    def row_number = 2 // with header

    fastq_files = channel.fromPath("$file") | splitCsv(header:true, strip:true) | map { row -> 

        // Check that required columns are present
        if (row.sample_id === null || row.fastq === null || row.reference_path === null) {
            error "Sample sheet did not contain required columns (sample_id, fastq, reference_path)"
        }

        // Check that required values for this row are set
        if (row.sample_id.isEmpty() || row.fastq.isEmpty() || row.reference_path.isEmpty()) {
            error "Sample sheet did not contain required values (sample_id, fastq, reference_path) in row: ${row_number}"
        }

        fastq = new File(row.fastq)

        if (!fastq.exists()){
            error "Fastq read file does not exist: ${forward}"
        }
        if (!reference.exists()){
            error "Reference genome file does not exist: ${reference}"
        }

        row_number++

        return tuple(row.sample_id, row.fastq, row.reference_path)

    }

    fastq_files | ifEmpty { exit 1, "Could not find read files specified in sample sheet." }

    return fastq_files

}