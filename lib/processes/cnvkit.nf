process CnvKitAneuploidy {

    tag { "$id : $idx_name" }
    label "cnvkit_aneuploidy"

    publishDir "$params.outdir/workflow/$params.subdir/${id}", mode: "copy", pattern: "${id}-scatter.png"
    publishDir "$params.outdir/workflow/$params.subdir/${id}", mode: "copy", pattern: "${id}.cns"
    publishDir "$params.outdir/workflow/$params.subdir/${id}", mode: "copy", pattern: "${id}.cnr"

    input:
    tuple val(id), val(idx_name), path(sample_bam)
    path(control_bam)
    path(fasta)

    output:
    tuple (val(id), path("results"), emit: results)
    path("${id}-scatter.png")
    path("${id}.cns")
    path("${id}.cnr")
    


    script:

    drop_low_coverage = params.host.aneuploidy.cnvkit.drop_low_coverage ? "--drop-low-coverage" : ""
    target_average_bin_size = params.host.aneuploidy.cnvkit.target_avg_size ? "--target-avg-size $params.host.aneuploidy.cnvkit.target_avg_size" : ""

    """
    cnvkit.py batch $sample_bam -n $control_bam -m wgs -f $fasta -p $task.cpus \
        --output-reference ${id}_${idx_name}.cnn --output-dir results/ --diagram --scatter $drop_low_coverage $target_average_bin_size
    mv results/${id}-scatter.png .
    mv results/${id}.cns .
    mv results/${id}.cnr .
    """

}