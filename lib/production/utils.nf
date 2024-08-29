
process StageInputFiles {

    input:
    path(stageJson)

    output:
    tuple env(PIPELINE), env(SAMPLE_ID), path("*.gz")

    script:

    """
    SAMPLE_ID=\$(cerebro-fs stage --json $stageJson --outdir . --pipeline pipeline.txt)
    PIPELINE=\$(cat pipeline.txt)
    """
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
            throw new RuntimeException("Neither Bowtie2 nor Bowtie21 index files exist for: ${index}")
        }
    }

    // Return the list of valid file paths
    return bowtie2Files
}

def getPanviralEnrichmentDatabases() {
    return [
        host: getPanviralEnrichmentHostDatabase(),
        virus: getPanviralEnrichmentVirusDatabase(),
        control:  getPanviralEnrichmentControlDatabase()
    ]
}

def getPanviralEnrichmentHostDatabase() {

    hostIndex = params.panviralEnrichment.hostIndex

    // Check if the host aligner is Bowtie2 and retrieve the appropriate index files
    if (params.panviralEnrichment.hostAligner == "bowtie2") {
        try {
            hostIndexPaths = getBowtie2IndexFiles(hostIndex)
        } catch (RuntimeException e) {
            throw new RuntimeException("Error retrieving Bowtie2 index files for host database: ${e.message}")
        }
    } else {
        // For other aligners, simply use the hostIndex as a single file path
        hostIndexPaths = [new File(hostIndex).absolutePath]
        if (!new File(hostIndex).exists()) {
            throw new RuntimeException("Host index file does not exist: ${hostIndex}")
        }
    }
    return hostIndexPaths
}


def getPanviralEnrichmentVirusDatabase() {

    virusReference = params.panviralEnrichment.virusReference
    virusIndex = params.panviralEnrichment.virusIndex

    // Check if the host aligner is Bowtie2 and retrieve the appropriate index files
    if (params.panviralEnrichment.virusAligner == "bowtie2") {
        try {
            virusIndexPaths = getBowtie2IndexFiles(virusIndex)
        } catch (RuntimeException e) {
            throw new RuntimeException("Error retrieving Bowtie2 index files for virus database: ${e.message}")
        }
    } else {
        // For other aligners, simply use the hostIndex as a single file path
        virusIndexFile = new File(virusIndex)
        virusIndexPaths = [virusIndexFile.absolutePath]
        if (!virusIndexFile.exists()) {
            throw new RuntimeException("Virus index file does not exist: ${virusIndex}")
        }
    }
    
    virusReferenceFile = new File(virusReference)
    if (!virusReferenceFile.exists()) {
        throw new RuntimeException("Virus reference file does not exist: ${virusReference}")
    }

    return [virusIndexPaths, virusReferenceFile.absolutePath]
}



def getPanviralEnrichmentControlDatabase() {

    controlReference = params.panviralEnrichment.controlReference
    controlIndex = params.panviralEnrichment.controlIndex

    // Check if the host aligner is Bowtie2 and retrieve the appropriate index files
    if (params.panviralEnrichment.controlAligner == "bowtie2") {
        try {
            controlIndexPaths = getBowtie2IndexFiles(controlIndex)
        } catch (RuntimeException e) {
            throw new RuntimeException("Error retrieving Bowtie2 index files for control database: ${e.message}")
        }
    } else {
        // For other aligners, simply use the hostIndex as a single file path
        controlIndexFile = new File(controlIndex)
        controlIndexPaths = [controlIndexFile.absolutePath]
        if (!controlIndexFile.exists()) {
            throw new RuntimeException("Virus index file does not exist: ${controlIndex}")
        }
    }

    controlReferenceFile = new File(controlReference)
    if (!controlReferenceFile.exists()) {
        throw new RuntimeException("Virus reference file does not exist: ${controlReference}")
    }

    return [controlIndexPaths, controlReferenceFile.absolutePath]
}