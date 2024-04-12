process MetaSpades {

    label "spades_meta"
    tag { id }

    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.fasta"

    input:
    tuple val(id), path(forward), path(reverse)

    output:
    tuple(val(id), path("${id}.fasta"), emit: contigs)

    script:

    """
    spades.py --meta -t $task.cpus -1 $forward -2 $reverse -o ${id}_meta_spades -k $params.meta.spades.k
    mv ${id}_meta_spades/contigs.fasta ${id}.fasta
    """

}