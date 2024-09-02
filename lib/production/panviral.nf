/* ==================================
 * PANVIRAL ENRICHMENT FOR PRODUCTION
 * ==================================
*/



process ReadQuality {

    label "panviralReadQuality"
    tag { sampleID }


    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.reads.json"

    input:
    tuple val(sampleID), path(forward), path(reverse)

    output:
    tuple (val(sampleID), path("${sampleID}__reads__R1.fq.gz"), path("${sampleID}__reads__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), path("${sampleID}.reads.json"), emit: results)
   

    script:
    
    adapter_sequence = params.panviralEnrichment.adapterForward && params.panviralEnrichment.adapterReverse ? "--adapter_sequence=$params.panviralEnrichment.adapterForward --adapter_sequence_r2=$params.panviralEnrichment.adapterReverse" : ""
    
    """
    fastp -i $forward -I $reverse -o ${sampleID}__reads__R1.fq.gz -O ${sampleID}__reads__R2.fq.gz --thread $task.cpus --json ${sampleID}.reads.json --length_required 50 --cut_tail --cut_tail_mean_quality 20 --low_complexity_filter --complexity_threshold 30 --detect_adapter_for_pe --trim_poly_g --poly_g_min_len 10 $adapter_sequence
    """

}


process HostDepletion {

    label "panviralHostDepletion"
    tag { sampleID }


    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.host.json"
    
    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(index)
    val(aligner)

    output:
    tuple (val(sampleID), path("${sampleID}__host__R1.fq.gz"), path("${sampleID}__host__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), path("${sampleID}.host.json"), emit: results)
    
    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    scrubby reads -i $forward -i $reverse --index $alignmentIndex --aligner $aligner --threads $task.cpus -o ${sampleID}__host__R1.fq.gz -o ${sampleID}__host__R2.fq.gz --json ${sampleID}.host.json
    """
}


process InternalControls {
    
    label "panviralInternalControls"
    tag { sampleID }

    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.controls.tsv"
    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.controls.json"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
    val(aligner)

    output:
    tuple (val(sampleID), path("${sampleID}__controls__R1.fq.gz"), path("${sampleID}__controls__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), path("${sampleID}.controls.tsv"), emit: results)
    path("${sampleID}.controls.json")

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    vircov coverage -i $forward -i $reverse -o ${sampleID}.controls.tsv --aligner $aligner --index $alignmentIndex --reference vircov__reference --threads $task.cpus --workdir data/ --zero --read-id reads.txt
    scrubby alignment -i $forward -i $reverse -a reads.txt -o ${sampleID}__controls__R1.fq.gz -o ${sampleID}__controls__R2.fq.gz --json ${sampleID}.controls.json
    """
    
}

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
    vircov run -i $forward -i $reverse -o ${sampleID}.viruses.tsv --index $alignmentIndex --reference vircov__reference --scan-threads $task.cpus --remap-threads $params.panviralEnrichment.remapThreads --parallel $params.panviralEnrichment.remapParallel --workdir data/
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

workflow PanviralEnrichment {

    take:
        reads
        hostDatabase
        virusDatabase
        controlDatabase
    main:

        ReadQuality(reads)

        HostDepletion(
            ReadQuality.out.reads, 
            hostDatabase, 
            params.panviralEnrichment.hostAligner
        )

        InternalControls(
            HostDepletion.out.reads,
            controlDatabase,
            params.panviralEnrichment.controlAligner
        )

        VirusRecovery(
            InternalControls.out.reads, 
            virusDatabase, 
            params.panviralEnrichment.virusAligner
        )

        results = ReadQuality.out.results.mix(
            HostDepletion.out.results, 
            InternalControls.out.results,
            VirusRecovery.out.results
        )

        results | groupTuple | ProcessOutput

}