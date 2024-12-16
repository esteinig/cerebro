import groovy.json.JsonOutput

/* Production */

def getProductionConfig() {

    return [
        apiUrl: params.cerebroProduction.authConfig.apiUrl ? params.cerebroProduction.authConfig.apiUrl : error("URL for Cerebro API was not provided (-> params.cerebroProduction.authConfig.apiUrl)"),
        authToken: params.cerebroProduction.authConfig.tokenEnvironmentVariable ? System.getenv(params.cerebroProduction.authConfig.tokenEnvironmentVariable): error("No token environment variable was provided for production configuration (-> params.cerebroProduction.authConfig.tokenEnvironmentVariable)"),
     ]
}

/* Panviral Enrichment */

def getPanviralEnrichmentDatabases() {

    return [
        panviralEnrichment: getPanviralEnrichmentVirusDatabases(),
        qualityControl: getQualityControlDatabases(),
    ]
    
}

def getPanviralEnrichmentVirusDatabases() {

    return [
        virusDatabase: getAlignmentReferenceIndex(
            params.panviralEnrichment.virusReference, 
            params.panviralEnrichment.virusIndex, 
            params.panviralEnrichment.virusAligner,
            "panviral enrichment :: virus"
        ),
        taxonomy: params.cerebroProduction.enabled ? getPanviralEnrichmentTaxonomy() : Channel.empty() // necessary for production only
    ]
}

def getPanviralEnrichmentTaxonomy() {

    return getFilePath(
        params.panviralEnrichment.taxonomy, 
        "panviral enrichment :: taxonomy"
    )
}

/* Pathogen Detection */

def getPathogenDetectionDatabases() {
    
    return [
        qualityControl: getQualityControlDatabases(),
        taxonomicProfile: getTaxonomicProfileDatabases(),
        metagenomeAssembly: getMetagenomeAssemblyDatabases(),
    ]
}


def getQualityControlDatabases() {

    return [
        hostDepletion:       params.qualityControl.hostDepletion       ? getQualityHostDatabase(params.qualityControl)          :  Channel.empty(),
        internalControls:    params.qualityControl.internalControls    ? getQualityInternalControls(params.qualityControl)      :  Channel.empty(),
        backgroundDepletion: params.qualityControl.backgroundDepletion ? getQualityBackgroundDatabase(params.qualityControl)    :  Channel.empty(),
    ]
}

/* Pathogen metagenome assembly profiling databases */

def getMetagenomeAssemblyDatabases() {

    def magParams = params.pathogenDetection.metagenomeAssembly;

    return [
        ncbiDatabase:               magParams.ncbiDatabase    ?  getPathogenAssemblyNcbiDatabase(magParams)                 :  Channel.empty(),
        contigProfile:              magParams.contigProfile   ?  getPathogenAssemblyContigProfileDatabase(magParams)       :  Channel.empty(),
    ]
}


def getPathogenAssemblyNcbiDatabase(magParams) {

    return getFilePath(
        magParams.ncbiDatabaseIndex, 
        "pathogen detection :: assembly :: ncbiDatabase"
    )
}


def getPathogenAssemblyContigProfileDatabase(magParams) {

    return getFilePath(
        magParams.contigProfileIndex, 
        "pathogen detection :: assembly :: profileDatabase"
    )
}


def getPathogenAssemblySylphDatabase(magParams) {

    return getFilePath(
        magParams.contigProfileIndex, 
        "pathogen detection :: assembly :: sylph"
    )
}
def getPathogenAssemblySylphDatabaseMetadata(magParams) {

    return getFilePath(
        magParams.contigProfileMetadata,
        "pathogen detection :: assembly :: sylph meta"
    )
}

/* Pathogen taxonomic profiling databases */


def getTaxonomicProfileDatabases() {

    def profileParams = params.pathogenDetection.taxonomicProfile;

    return [
        vircovDatabase:     profileParams.alignment                                                                                                         ?  getPathogenProfileVircovDatabase(profileParams)        :  Channel.empty(),
        krakenDatabase:     profileParams.classifier && profileParams.classifierMethod.contains("kraken2")                                                  ?  getPathogenProfileKrakenDatabase(profileParams)        :  Channel.empty(),
        metabuliDatabase:   profileParams.classifier && profileParams.classifierMethod.contains("metabuli")                                                 ?  getPathogenProfileMetabuliDatabase(profileParams)      :  Channel.empty(),
        sylphDatabase:      profileParams.profiler && profileParams.profilerMethod.contains("sylph")                                                        ?  getPathogenProfileSylphDatabase(profileParams)         :  Channel.empty(),
        sylphMetadata:      profileParams.profiler && profileParams.profilerMethod.contains("sylph")                                                        ?  getPathogenProfileSylphDatabaseMetadata(profileParams) :  Channel.empty(),
        kmcpDatabase:       profileParams.profiler && (profileParams.classifierMethod.contains("kmcp") || profileParams.profilerMethod.contains("kmcp"))    ?  getPathogenProfileKmcpDatabase(profileParams)          :  Channel.empty(),
        ganonDatabase:      profileParams.profiler && (profileParams.classifierMethod.contains("ganon") || profileParams.profilerMethod.contains("ganon"))  ?  getPathogenProfileGanonDatabase(profileParams)         :  Channel.empty(),
        taxonomy:           profileParams.taxonomy                                                                                                          ?  getPathogenProfileTaxonomy(profileParams)              :  Channel.empty(),
    ]
}

def getPathogenProfileTaxonomy(profileParams) {

    return getFilePath(
        profileParams.taxonomy, 
        "pathogen detection :: tax profile :: taxonomy"
    )
}

def getPathogenProfileVircovDatabase(profileParams) {

    return getAlignmentReferenceIndex(
        profileParams.alignmentReference, 
        profileParams.alignmentIndex, 
        profileParams.alignmentMethod,
        "pathofen detection :: tax profile :: vircov"
    )
}
def getPathogenProfileKrakenDatabase(profileParams) {

    return getFilePath(
        profileParams.krakenIndex, 
        "pathogen detection :: tax profile :: kraken"
    )
}
def getPathogenProfileMetabuliDatabase(profileParams) {

    return getFilePath(
        profileParams.metabuliIndex, 
        "pathogen detection :: tax profile :: metabuli"
    )
}
def getPathogenProfileSylphDatabase(profileParams) {

    return getFilePath(
        profileParams.sylphIndex, 
        "pathogen detection :: tax profile :: sylph"
    )
}
def getPathogenProfileSylphDatabaseMetadata(profileParams) {

    return getFilePath(
        profileParams.sylphMetadata,
        "pathogen detection :: tax profile :: sylph meta"
    )
}
def getPathogenProfileKmcpDatabase(profileParams) {

    return getFilePath(
        profileParams.kmcpIndex, 
        "pathogen detection :: tax profile :: kmcp"
    )
}
def getPathogenProfileGanonDatabase(profileParams) {

    return getFilePath(
        profileParams.ganonIndex, 
        "pathogen detection :: tax profile :: ganon"
    )
}

/* Quality control databases */

def getQualityHostDatabase(qualityControlParams) {

    return getAlignmentIndex(
        qualityControlParams.hostDepletionIndex, 
        qualityControlParams.hostDepletionAligner,
        "pathogen detection :: quality control :: host"
    )
}
def getQualityBackgroundDatabase(qualityControlParams) {

    return getAlignmentReferenceIndex(
        qualityControlParams.backgroundDepletionReference, 
        qualityControlParams.backgroundDepletionIndex, 
        qualityControlParams.backgroundDepletionAligner,
        "pathogen detection :: quality control :: background"
    )
}
def getQualityInternalControls(qualityControlParams) {

    return getAlignmentReferenceIndex(
        qualityControlParams.internalControlsReference, 
        qualityControlParams.internalControlsIndex, 
        qualityControlParams.internalControlsAligner,
        "pathogen detection :: quality control :: internal controls"
    )
}

/**
 * Retrieves the Bowtie2 or Bowtie21 index files based on the provided index prefix.
 *
 * This function first checks for the presence of Bowtie2 index files with extensions
 * such as `.1.bt2`, `.2.bt2`, etc. If any of these files are missing, it then checks
 * for an alternative set of Bowtie21 index files with extensions like `.1.bt21`, `.2.bt21`, etc.
 * 
 * The function ensures that either all `.bt2` or all `.bt21` files exist. If neither
 * set is fully available, an exception is thrown.
 * 
 * @param index The prefix of the Bowtie2 index files (e.g., the base name without the extension).
 * @return A list of absolute paths to the existing Bowtie2 or Bowtie21 index files.
 * @throws RuntimeException if neither the complete set of Bowtie2 nor Bowtie21 index files are found.
 */
def getBowtie2IndexFiles(String index) {
    // List of expected Bowtie2 index file extensions
    def bowtie2Extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    
    // List of expected Bowtie21 index file extensions as an alternative
    def bowtie21Extensions = ['.1.bt21', '.2.bt21', '.3.bt21', '.4.bt21', '.rev.1.bt21', '.rev.2.bt21']

    // Helper closure to check the existence of files with given extensions
    def checkFiles = { extensions ->
        extensions.collect { ext -> 
            def file = new File("${index}${ext}") // Construct the file path
            if (!file.exists()) {
                return null // Return null if the file does not exist
            }
            return file.absolutePath // Return the absolute path of the file
        }
    }

    // Attempt to retrieve the Bowtie2 files
    def bowtie2Files = checkFiles(bowtie2Extensions)

    // If any Bowtie2 file is missing, attempt to retrieve the Bowtie21 files
    if (bowtie2Files.any { it == null }) {
        bowtie2Files = checkFiles(bowtie21Extensions)
        
        // If any Bowtie21 file is missing, throw an exception
        if (bowtie2Files.any { it == null }) {
            error "Neither Bowtie2 nor Bowtie21 index files exist for: ${index}"
        }
    }

    // Return the list of valid file paths
    return bowtie2Files
}


def getAlignmentIndex(String index, String aligner, String description) {

    if (index == null) {
        error "Index path for '${description}' cannot be null"
    }

    // Check if the aligner is Bowtie2 and retrieve the appropriate index files
    if (aligner == "bowtie2") {
        try {
            indexPaths = getBowtie2IndexFiles(index)
        } catch (RuntimeException e) {
            error "Error retrieving Bowtie2 index files for '${description}' database: ${e.message}"
        }
    } else {
        // For other aligners, simply use the index as a single file path
        indexPaths = [new File(index).absolutePath]
        if (!new File(index).exists()) {
            error "Index file for '${description}' database does not exist: ${index}"
        }
    }
    return indexPaths
    
}


def getFilePath(String file, String description) {

    if (file == null) {
        error "File path for '${description}' cannot be null"
    }
    filePath = new File(file).absolutePath
    if (!new File(file).exists()) {
        error "File for '${description}' database does not exist: ${file}"
    }
    return filePath
    
}

def getClassifierIndex(String index, String description) {

    if (index == null) {
        error "Index path for '${description}' cannot be null"
    }

    // For other aligners, simply use the index as a single file path
    indexPaths = [new File(index).absolutePath]
    if (!new File(index).exists()) {
        error "Index file for '${description}' database does not exist: ${index}"
    }
    return indexPaths
    
}

/* Wrapper alignment infex/reference database functions */

def getAlignmentReferenceIndex(String reference, String index, String aligner, String description) {

    indexPaths = getAlignmentIndex(index, aligner, description)

    referenceFile = new File(reference)
    if (!referenceFile.exists()) {
        error "Reference sequence file for '${description}' database does not exist: ${referenceFile}"
    }
    
    return [indexPaths, referenceFile.absolutePath]

}


def getClassifierReferenceIndex(String reference, String index, String description) {

    indexPaths = getClassifierIndex(index, description)

    referenceFile = new File(reference)
    if (!referenceFile.exists()) {
        error "Reference sequence file for '${description}' database does not exist: ${referenceFile}"
    }
    
    return [indexPaths, referenceFile.absolutePath]

}

/* Reads and sample sheet inputs */

def getReads(String fastqPaired, String fastqNanopore, String sampleSheet, Boolean sampleSheetProduction) {
    
    def pairedReads = fastqPaired ? fastqPaired.split(/\s+/) as List<String> : fastqPaired
    def nanoporeReads = fastqNanopore ? fastqNanopore.split(/\s+/) as List<String> : fastqNanopore

    if (sampleSheet){
        return readSampleSheet(sampleSheet, false)
    } else if (fastqPaired) {
        return channel.fromFilePairs(pairedReads, flat: true, checkIfExists: true)
    } else if (fastqNanopore) {
        return channel.fromPath(nanoporeReads, checkIfExists: true) | map { tuple(it.getSimpleName(), it) } 
    } else {
        error "Either one or both of '--fastqPaired' and '--fastqNanopore' or '--sampleSheet' have to be specified for read input"
    }
}


def readSampleSheet(file, production){

    def row_number = 2 // with header

    fastqFiles = channel.fromPath("$file") | splitCsv(header:true, strip:true) | map { row -> 

        // Check that required columns are present
        if (row.sample_id === null || row.forward_path === null || row.reverse_path === null) {
            error "Sample sheet did not contain required columns (sample_id, forward_path, reverse_path)"
        }

        // Check that required values for this row are set
        if (row.sample_id.isEmpty() || row.forward_path.isEmpty() ) {
            error "Sample sheet did not contain required values (sample_id, forward_path) in row: ${row_number}"
        }

        forward = new File(row.forward_path)
        if (!forward.exists()){
            error "Forward read file does not exist: ${forward}"
        }

        if (!row.reverse_path.isEmpty()) {
            reverse = new File(row.reverse_path)
            if (!reverse.exists()){
                error "Reverse read file does not exist: ${reverse}"
            }
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

    fastqFiles | ifEmpty { exit 1, "Could not find read files specified in sample sheet." }

    return fastqFiles | map { tuple(it[0], it[1], it[2]) }

}

/* Formatted messages for users */

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
c_indigo = params.monochrome ? '' : "\033[38;5;57m";
c_light_indigo = params.monochrome ? '' : "\033[38;5;63m";
c_light_blue = params.monochrome ? '' : "\033[38;5;33m";

def initMessage(){


    log.info """
    ${c_indigo}=====================
    ${c_light_blue}C E R E B R O
    ${c_indigo}=====================${c_reset}

    ${workflow.manifest.description}

    Version:              ${c_light_indigo}v${workflow.manifest.version}${c_reset}
    Documentation:        ${c_light_blue}https://docs.meta-gp.org${c_reset}

    Workflow:             ${c_light_blue}${workflow.sessionId}${c_reset}
    Started:              ${c_light_indigo}${workflow.start}${c_reset}

    Profiles:             ${c_white}${workflow.profile}${c_reset}
    Working directory:    ${c_white}${workflow.workDir}${c_reset}

    """.stripIndent()
}

def completionMessage(){

    log.info """

    ${c_indigo}=====================
    ${c_light_blue}C E R E B R O
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
    
    ${c_light_blue}https://cerebro.meta-gp.org${c_reset}
 
    ${c_indigo}=============================================================${c_reset}
    """.stripIndent()
}

def helpMessage(){

    log.info """
    ${c_indigo}=====================
    ${c_light_blue}C E R E B R O
    ${c_indigo}=====================${c_reset}

    Version:              ${c_light_indigo}v${workflow.mainfest.version}${c_reset}
    Documentation:        ${c_light_blue}https://cerebro.meta-gp.org${c_reset}

    """.stripIndent()
    System.exit(0)
}


process PipelineConfig {

    publishDir "$params.outputDirectory", mode: "copy", pattern: "config.json"

    input:
        val samples
        val cerebroWorkflow
        val started

    output:
        path 'config.json', emit: config

    script:

        completed = java.time.LocalDateTime.now()

        config = [
            id: "$workflow.sessionId",
            name: "$workflow.runName",
            pipeline: "$workflow.manifest.name",
            version: "$workflow.manifest.version",
            started: "$started",
            completed: "$completed",
            workflow: "$cerebroWorkflow",
            params: params
        ]

        json = JsonOutput.toJson(config)
        json_pretty = JsonOutput.prettyPrint(json)

        """
        echo '${json_pretty}' > config.json
        """
        
}