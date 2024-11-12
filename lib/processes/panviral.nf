

process VirusRecovery {
    
    label "panviralVirusRecovery"
    tag { sampleID }

    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.tsv"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
    val(aligner)
    val(vircovArgs)

    output:
    tuple (val(sampleID), path(forward), path(reverse), emit: reads)
    tuple (val(sampleID), path("${sampleID}.tsv"), emit: results)

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    vircov run -i $forward -i $reverse -o ${sampleID}.tsv --index $alignmentIndex --reference vircov__reference --scan-threads $task.cpus --remap-threads $params.resources.threads.vircovRemap --parallel $params.resources.threads.vircovParallel --workdir data/ $vircovArgs
    """
    
}


process ProcessOutput {
    
    label "panviralVirusTools"
    tag { sampleID }

    publishDir "$params.outputDirectory/results", mode: "copy", pattern: "viruses.vircov.tsv"

    input:
    tuple val(sampleID), path(result_files)

    output:
    tuple val(sampleID), path("viruses.vircov.tsv")

    script:

    """
    vircov tools concat-output --input *.tsv --output viruses.vircov.tsv --file-id
    """
    
}
