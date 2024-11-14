process Dragonflye {

    label "dragonflye"
    
    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.contigs.fa", saveAs: { "assembly__dragonflye__fasta" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "symlink", pattern: "${id}_qc.fq.gz"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.json"

    input:
    tuple val(id), path(reads)

    output:
    tuple (val(id), path("${id}_qc.fq.gz"), emit: reads)
    tuple (val(id), path("qc__nanoq__reads"), emit: results)
    path("${id}.json")
    
    publishDir "$params.outdir/culture_id/ont/assembly", mode: "copy", pattern: "*.contigs.fa"
    publishDir "$params.outdir/culture_id/ont/assembly", mode: "copy", pattern: "*.assembly.tsv"

    input:
    tuple val(id), path(reads), path(reference)

    output:
    tuple (val(id), path("${id}.contigs.fa"), emit: assembly)
    tuple (val(id), path("${id}.assembly.tsv"), emit: results)

    script:

    """
    dragonflye --reads $reads --outdir ${id}_assembly --cpus $task.cpus $params.dragonflye.args
    mv ${id}_assembly/contigs.fa ${id}.contigs.fa
    assembly-scan ${id}.contigs.fa > ${id}.assembly.tsv
    """

}