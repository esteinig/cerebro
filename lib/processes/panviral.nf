

process VirusRecovery {
    
    label "panviralVirusRecovery"
    tag { sampleID }

    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.viruses.tsv"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
    val(aligner)

    output:
    tuple (val(sampleID), path(forward), path(reverse), emit: reads)
    tuple (val(sampleID), path("${sampleID}.viruses.tsv"), emit: results)

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    vircov run -i $forward -i $reverse -o ${sampleID}.viruses.tsv --index $alignmentIndex --reference vircov__reference --scan-threads $task.cpus --remap-threads $params.resources.threads.vircovRemap --parallel $params.resources.threads.vircovParallel --workdir data/
    """
    
}


process ProcessOutput {
    
    label "cerebro"
    tag { sampleID }

    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.qc.json"

    input:
    tuple val(sampleID), path(result_files)

    output:
    tuple val(sampleID), path("${sampleID}.qc.json")

    script:

    """
    cerebro-pipe process panviral --id ${sampleID} --qc ${sampleID}.qc.json 
    """
    
}
