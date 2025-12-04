

process VirusRecovery {
    
    label "panviralVirusRecovery"
    tag { sampleID }

    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.vircov.tsv"
    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}/*"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), path(reference, name: 'vircov__reference')  // index and reference can be the same
    val(aligner)
    val(vircovArgs)

    output:
    tuple (val(sampleID), path(forward), path(reverse), emit: reads)
    tuple (val(sampleID), path("${sampleID}.vircov.tsv"), emit: results)
    tuple (val(sampleID), path("${sampleID}/*"), emit: output)

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    vircov run -i $forward -i $reverse -o ${sampleID}.vircov.tsv --aligner $aligner --index $alignmentIndex --reference vircov__reference --scan-threads $task.cpus --remap-threads $params.resources.threads.vircovRemap --remap-parallel $params.resources.threads.vircovParallel --workdir ${sampleID} --keep $vircovArgs
    """
    
}


process PanviralTable {
    
    label "panviralVirusTools"

    publishDir "$params.outputDirectory/results", mode: "copy", pattern: "panviral.tsv"

    input:
    path(vircov_outputs)

    output:
    path("panviral.tsv")

    script:

    """
    vircov tools concat-output --input *.tsv --output panviral.tsv --file-id   
    """
    
}

process ProcessOutput {
    
    tag { sampleID }
    label "cerebro"

    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.qc.json"
    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.pv.json"

    publishDir "$params.outputDirectory/results/samples/$sampleID", mode: "copy", pattern: "${sampleID}.qc.json"
    publishDir "$params.outputDirectory/results/samples/$sampleID", mode: "copy", pattern: "${sampleID}.pv.json"

    input:
    tuple val(sampleID), path(result_files)

    output:
    tuple (val(sampleID), path("${sampleID}.qc.json"), path("${sampleID}.pv.json"), emit: results)
    val(sampleID), emit: samples

    script:

    """
    cerebro-pipeline process panviral --id ${sampleID} --qc ${sampleID}.qc.json --panviral ${sampleID}.pv.json --paired-end --qc-fail-ok
    """
    
}