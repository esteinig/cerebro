process Metabuli {

    tag { id }
    label "metabuli"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}_${idx_name}.metabuli.report", saveAs: { "kmer__metabuli__${idx_name}" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.metabuli.tsv"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}.metabuli.report"

    input:
    tuple val(id), path(forward), path(reverse)
    each path(metabuli_db)

    output:
    tuple(val(id), path("kmer__metabuli__${idx_name}"), emit: results)
    tuple path("${id}_${idx_name}.metabuli.report"), path("${id}_${idx_name}.metabuli")

    script:

    mem = task.memory.split()[0]
    idx_name = metabuli_db.baseName

    min_score = params.taxa.metabuli.precision.enabled && params.taxa.metabuli.precision.min_score ? "--min-score $params.taxa.metabuli.precision.min_score" : ""
    min_sp_score = params.taxa.metabuli.precision.enabled && params.taxa.metabuli.precision.min_sp_score ? "--min-sp-score $params.taxa.metabuli.precision.min_sp_score" : ""

    """
    metabuli classify $min_score $min_sp_score --max-ram $mem --threads $task.cpus $forward $reverse $metabuli_db classified/ $id
    
    cp classified/${id}_classifications.tsv ${id}_${idx_name}.metabuli.tsv
    rm classified/${id}_classifications.tsv

    cp classified/${id}_report.tsv ${id}_${idx_name}.metabuli.report
    """

}

process MetabuliOnt {

    tag { id }
    label "metabuli"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}_${idx_name}.metabuli.report", saveAs: { "kmer__metabuli__${idx_name}" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.metabuli.tsv"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}.metabuli.report"

    input:
    tuple val(id), path(reads)
    each path(metabuli_db)

    output:
    tuple(val(id), path("kmer__metabuli__${idx_name}"), emit: results)
    tuple path("${id}_${idx_name}.metabuli.report"), path("${id}_${idx_name}.metabuli")

    script:

    mem = task.memory.split()[0]
    idx_name = metabuli_db.baseName

    min_score = params.taxa.metabuli.precision.enabled && params.taxa.metabuli.precision.min_score ? "--min-score $params.taxa.metabuli.precision.min_score" : ""
    min_sp_score = params.taxa.metabuli.precision.enabled && params.taxa.metabuli.precision.min_sp_score ? "--min-sp-score $params.taxa.metabuli.precision.min_sp_score" : ""

    """
    metabuli classify $min_score $min_sp_score --max-ram $mem --threads $task.cpus $reads $metabuli_db classified/ $id
    
    cp classified/${id}_classifications.tsv ${id}_${idx_name}.metabuli.tsv
    rm classified/${id}_classifications.tsv

    cp classified/${id}_report.tsv ${id}_${idx_name}.metabuli.report
    """

}