/* ==================================
 * PANVIRAL ENRICHMENT FOR PRODUCTION
 * ==================================
*/



process Fastp {

    label "fastp"
    tag { sampleID }


    publishDir "$params.outdir/panviral_enrichment/$sampleID/$stageID", mode: "copy", pattern: "${sampleID}.fastp.json"
    publishDir "$params.outdir/panviral_enrichment/$sampleID/$stageID", mode: "copy", pattern: "${stageJson}"

    input:
    tuple val(sampleID), val(stageID), path(stageJson), path(forward), path(reverse)

    output:
    tuple (val(sampleID), val(stageID), path(stageJson), path("${sampleID}__qc__R1.fq.gz"), path("${sampleID}__qc__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), val(stageID), path("${sampleID}.fastp.json"), emit: results)
   

    script:
    
    adapter_sequence = params.panviralEnrichment.adapterForward && params.panviralEnrichment.adapterReverse ? "--adapter_sequence=$params.panviralEnrichment.adapterForward --adapter_sequence_r2=$params.panviralEnrichment.adapterReverse" : ""
    

    """
    fastp -i $forward -I $reverse -o ${sampleID}__qc__R1.fq.gz -O ${sampleID}__qc__R2.fq.gz --thread $task.cpus --json ${sampleID}.fastp.json --length_required 50 --cut_tail --cut_tail_mean_quality 20 --low_complexity_filter --complexity_threshold 30 --detect_adapter_for_pe --trim_poly_g --poly_g_min_len 10 $adapter_sequence
    """

}

process Scrubby {

    label "scrubby"
    tag { sampleID }

    publishDir "$params.outdir/panviral_enrichment/$sampleID/$stageID", mode: "copy", pattern: "${sampleID}.scrubby.json"
    
    input:
    tuple val(sampleID), val(stageID), path(stageJson), path(forward), path(reverse)
    path(index)
    val(aligner)

    output:
    tuple (val(sampleID), val(stageID), path(stageJson), path("${sampleID}__host__R1.fq.gz"), path("${sampleID}__host__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), val(stageID), path("${sampleID}.scrubby.json"), emit: results)
    
    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    scrubby reads -i $forward -i $reverse -o ${sampleID}__host__R1.fq.gz -o ${sampleID}__host__R2.fq.gz --aligner $aligner --index $alignmentIndex --threads $task.cpus --json ${sampleID}.scrubby.json
    """
}



process Vircov {
    
    label "vircov"
    tag { sampleID }

    publishDir "$params.outdir/panviral_enrichment/$sampleID", mode: "copy", pattern: "${sampleID}.vircov.tsv"

    input:
    tuple val(sampleID), val(stageID), path(stageJson), path(forward), path(reverse)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
    val(aligner)

    output:
    tuple (val(sampleID), val(stageID), path(stageJson), path(forward), path(reverse), emit: reads)
    tuple (val(sampleID), val(stageID), path("${sampleID}.vircov.tsv"), emit: results)

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    vircov run -i $forward -i $reverse -o ${sampleID}.vircov.tsv --index $alignmentIndex --reference vircov__reference --scan-threads $task.cpus --remap-threads $params.panviralEnrichment.remapThreads --parallel $params.panviralEnrichment.remapParallel --workdir data/
    """
    
}

workflow PanviralEnrichment {

    take:
        reads
        hostIndex
        virusIndex
    main:

        readQualityControl = Fastp(reads)

        hostDepletion = Scrubby(
            readQualityControl.reads, 
            hostIndex, 
            params.panviralEnrichment.hostAligner
        )

        virusRecovery = Vircov(
            hostDepletion.reads, 
            virusIndex, 
            params.panviralEnrichment.virusAligner
        )

        results = readQualityControl.results.mix(
            hostDepletion.results, 
            virusRecovery.results
        )

    emit:
        results = results

}