/* ==================================
 * PANVIRAL ENRICHMENT FOR PRODUCTION
 * ==================================
*/



process QualityControl {

    label "fastp"
    tag { sampleID }


    publishDir "$params.outputDirectory/panviral_enrichment/$sampleID", mode: "copy", pattern: "${sampleID}.fastp.json"

    input:
    tuple val(sampleID), path(forward), path(reverse)

    output:
    tuple (val(sampleID), path("${sampleID}__qc__R1.fq.gz"), path("${sampleID}__qc__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), path("${sampleID}.fastp.json"), emit: results)
   

    script:
    
    adapter_sequence = params.panviralEnrichment.adapterForward && params.panviralEnrichment.adapterReverse ? "--adapter_sequence=$params.panviralEnrichment.adapterForward --adapter_sequence_r2=$params.panviralEnrichment.adapterReverse" : ""
    

    """
    fastp -i $forward -I $reverse -o ${sampleID}__qc__R1.fq.gz -O ${sampleID}__qc__R2.fq.gz --thread $task.cpus --json ${sampleID}.fastp.json --length_required 50 --cut_tail --cut_tail_mean_quality 20 --low_complexity_filter --complexity_threshold 30 --detect_adapter_for_pe --trim_poly_g --poly_g_min_len 10 $adapter_sequence
    """

}


process HostDepletion {

    label "scrubby"
    tag { sampleID }

    publishDir "$params.outputDirectory/panviral_enrichment/$sampleID", mode: "copy", pattern: "${sampleID}.scrubby.json"
    
    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(index)
    val(aligner)

    output:
    tuple (val(sampleID), path("${sampleID}__host__R1.fq.gz"), path("${sampleID}__host__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), path("${sampleID}.scrubby.json"), emit: results)
    
    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    scrubby reads -i $forward -i $reverse -o ${sampleID}__host__R1.fq.gz -o ${sampleID}__host__R2.fq.gz --aligner $aligner --index $alignmentIndex --threads $task.cpus --json ${sampleID}.scrubby.json
    """
}


process InternalControls {
    
    label "vircov_scrubby"
    tag { sampleID }

    publishDir "$params.outputDirectory/panviral_enrichment/$sampleID", mode: "copy", pattern: "${sampleID}.controls.tsv"
    publishDir "$params.outputDirectory/panviral_enrichment/$sampleID", mode: "copy", pattern: "${sampleID}.controls.json"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
    val(aligner)

    output:
    tuple (val(sampleID), path("${sampleID}__controls__R1.fq.gz"), path("${sampleID}__controls__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), path("${sampleID}.controls.tsv"), path("${sampleID}.controls.json"), emit: results)

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    vircov coverage -i $forward -i $reverse -o ${sampleID}.controls.tsv --aligner $aligner --index $alignmentIndex --reference vircov__reference --threads $task.cpus --workdir data/ --zero --read-id reads.txt
    scrubby alignment -i $forward -i $reverse -a reads.txt -o ${sampleID}__controls__R1.fq.gz -o ${sampleID}__controls__R2.fq.gz --json ${sampleID}.controls.json
    """
    
}

process VirusRecovery {
    
    label "vircov"
    tag { sampleID }

    publishDir "$params.outputDirectory/panviral_enrichment/$sampleID", mode: "copy", pattern: "${sampleID}.vircov.tsv"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
    val(aligner)

    output:
    tuple (val(sampleID), path(forward), path(reverse), emit: reads)
    tuple (val(sampleID), path("${sampleID}.vircov.tsv"), emit: results)

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
        hostDB
        virusDB
        controlDB
    main:

        QualityControl(reads)

        HostDepletion(
            QualityControl.out.reads, 
            hostDB, 
            params.panviralEnrichment.hostAligner
        )

        InternalControls(
            HostDepletion.out.reads,
            controlDB,
            params.panviralEnrichment.controlAligner
        )

        VirusRecovery(
            InternalControls.out.reads, 
            virusDB, 
            params.panviralEnrichment.virusAligner
        )

        QualityControl.out.results.mix(
            HostDepletion.out.results, 
            InternalControls.out.results,
            VirusRecovery.out.results
        )
}