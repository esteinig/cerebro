
process InputScan {

    tag { sampleID }
    label "qualityReadScan"


    publishDir "$params.outputDirectory/quality/$sampleID", mode: "copy", pattern: "${sampleID}.input.json"

    input:
    tuple val(sampleID), path(forward), path(reverse)

    output:
    tuple (val(sampleID), path("${sampleID}.input.json"), emit: results)

    script:
    
    """
    cerebro-pipe tools scan-reads -i $forward -i $reverse --json ${sampleID}.input.json
    """

}

process OutputScan {

    tag { sampleID }
    label "qualityReadScan"


    publishDir "$params.outputDirectory/quality/$sampleID", mode: "copy", pattern: "${sampleID}.output.json"

    input:
    tuple val(sampleID), path(forward), path(reverse)

    output:
    tuple (val(sampleID), path("${sampleID}.output.json"), emit: results)

    script:
    
    """
    cerebro-pipe tools scan-reads -i $forward -i $reverse --json ${sampleID}.output.json
    """

}

process SyntheticControls { 
    
    tag { sampleID }
    label "qualitySyntheticControls"

    publishDir "$params.outputDirectory/quality/$sampleID", mode: "copy", pattern: "${sampleID}.synthetic.tsv"
    publishDir "$params.outputDirectory/quality/$sampleID", mode: "copy", pattern: "${sampleID}.synthetic.json"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
    val(aligner)

    output:
    tuple (val(sampleID), path("${sampleID}__synthetic__R1.fq.gz"), path("${sampleID}__synthetic__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), path("${sampleID}.synthetic.tsv"), emit: results)
    path("${sampleID}.synthetic.json")

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    vircov coverage -i $forward -i $reverse -o ${sampleID}.synthetic.tsv --aligner $aligner --index $alignmentIndex --reference vircov__reference --threads $task.cpus --workdir data/ --zero --read-id reads.txt
    scrubby alignment -i $forward -i $reverse -a reads.txt -o ${sampleID}__synthetic__R1.fq.gz -o ${sampleID}__synthetic__R2.fq.gz --json ${sampleID}.synthetic.json
    """
    
}


process ReadQuality {

    tag { sampleID }
    label "qualityReadQuality"


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


process Deduplication {

    tag { sampleID }
    label "qualityDeduplication"


    publishDir "$params.outputDirectory/panviral/$sampleID", mode: "copy", pattern: "${sampleID}.dedup.json"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    val(head)
    val(deterministic)

    output:
    tuple (val(sampleID), path("${sampleID}__dedup__R1.fq.gz"), path("${sampleID}__dedup__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), path("${sampleID}.dedup.json"), emit: results)
   

    script:
    
    """
    cerebro-pipe tools umi-dedup-naive -i $forward -i $reverse -o ${sampleID}__dedup__R1.fq.gz -o ${sampleID}__dedup__R2.fq.gz --head $head --reads $forward --json ${sampleID}.dedup.json
    """

}

process HostDepletion {

    tag { sampleID }
    label "qualityHostDepletion"


    publishDir "$params.outputDirectory/quality/$sampleID", mode: "copy", pattern: "${sampleID}.host.json"
    
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
    
    tag { sampleID }
    label "qualityInternalControls"

    publishDir "$params.outputDirectory/quality/$sampleID", mode: "copy", pattern: "${sampleID}.controls.tsv"
    publishDir "$params.outputDirectory/quality/$sampleID", mode: "copy", pattern: "${sampleID}.controls.json"

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



process BackgroundDepletion {

    tag { sampleID }
    label "qualityBackgroundDepletion"

    publishDir "$params.outputDirectory/quality/$sampleID", mode: "copy", pattern: "${sampleID}.background.tsv"
    publishDir "$params.outputDirectory/quality/$sampleID", mode: "copy", pattern: "${sampleID}.background.json"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
    val(aligner)

    output:
    tuple (val(sampleID), path("${sampleID}__background__R1.fq.gz"), path("${sampleID}__background__R2.fq.gz"), emit: reads)
    tuple (val(sampleID), path("${sampleID}.background.tsv"), emit: results)
    path("${sampleID}.background.json")

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    """
    vircov coverage -i $forward -i $reverse -o ${sampleID}.background.tsv --aligner $aligner --index $alignmentIndex --reference vircov__reference --threads $task.cpus --workdir data/ --zero --read-id reads.txt
    scrubby alignment -i $forward -i $reverse -a reads.txt -o ${sampleID}__background__R1.fq.gz -o ${sampleID}__background__R2.fq.gz --json ${sampleID}.controls.json
    """
    
}


process ProcessOutput {
    
    tag { sampleID }
    label "cerebro"

    publishDir "$params.outputDirectory/quality/$sampleID", mode: "copy", pattern: "${sampleID}.qc.json"

    input:
    tuple val(sampleID), path(result_files)

    output:
    path("${sampleID}.qc.json")

    script:

    qc_bg_param = params.cerebroConfig.qualityControlBackgroundOnly ? "--qc-background" : ""

    """
    cerebro-pipe process quality --id ${sampleID} --output ${sampleID}.qc.json $qc_bg_param
    """
    
}


process QualityControlTables {
    
    tag { sampleID }
    label "cerebro"

    publishDir "$params.outputDirectory/quality", mode: "copy", pattern: "*.tsv"

    input:
    path(result_files)

    output:
    tuple path("qc_reads.tsv"), path("qc_bg.tsv"), path("qc_ctrl.tsv")

    script:

    """
    cerebro-pipe tables quality --json *.qc.json --reads qc_reads.tsv --background qc_bg.tsv --controls qc_ctrl.tsv
    """
    
}