process RasusaReadsMultiple {

    tag { "$id : $reads : $replicate" }
    label "rasusa"

    publishDir "$params.outdir/workflow/$params.subdir/fastq", mode: "symlink", pattern: "${id}_${reads}_${replicate}_*.fq.gz"

    input:
    tuple val(id), path(forward), path(reverse)
    each reads
    each replicate

    output:
    tuple val("${id}_${reads}_${replicate}"), path("${id}_${reads}_${replicate}_1.fq.gz"), path("${id}_${reads}_${replicate}_2.fq.gz")

    script:

    """
    rasusa -i $forward -i $reverse --num $reads -o ${id}_${reads}_${replicate}_1.fq.gz -o ${id}_${reads}_${replicate}_2.fq.gz
    """

}