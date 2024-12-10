

process VirusRecovery {
    
    label "panviralVirusRecovery"
    tag { sampleID }

    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.vircov.tsv"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
    val(aligner)
    val(vircovArgs)

    output:
    tuple (val(sampleID), path(forward), path(reverse), emit: reads)
    tuple (val(sampleID), path("${sampleID}.vircov.tsv"), emit: results)

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    vircov run -i $forward -i $reverse -o ${sampleID}.vircov.tsv --aligner $aligner --index $alignmentIndex --reference vircov__reference --scan-threads $task.cpus --remap-threads $params.resources.threads.vircovRemap --remap-parallel $params.resources.threads.vircovParallel --workdir data/ $vircovArgs
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

    input:
    tuple val(sampleID), path(result_files)

    output:
    tuple (val(sampleID), path("${sampleID}.qc.json"), path("${sampleID}.pv.json"), emit: results)

    script:

    """
    cerebro-pipe process panviral --id ${sampleID} --qc ${sampleID}.qc.json --panviral ${sampleID}.pv.json --paired-end
    """
    
}


process UploadOutput {
    
    tag { sampleID }
    label "cerebro"

    publishDir "$params.outputDirectory/panviral/models", mode: "copy", pattern: "${sampleID}.cerebro.json"

    input:
    tuple val(sampleID), path(quality), path(panviral)
    path(taxonomyDirectory)
    val(teamName)
    val(databaseName)
    val(projectName)

    output:
    path("${sampleID}.cerebro.json")

    script:

    """
    cerebro-client --team $teamName --token ${System.getenv(params.cerebroProduction.apiToken)} upload-panviral --quality $quality --panviral $panviral --taxonomy $taxonomyDirectory --team-name $teamName --db-name $databaseName --project-name $projectName --run-id TESTRUN --model-dir .
    """
    
}
