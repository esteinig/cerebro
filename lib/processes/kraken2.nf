process Kraken2Uniq {

    tag { id }
    label "kraken2uniq"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}_${idx_name}.kraken2uniq.report", saveAs: { "kmer__kraken2uniq__${idx_name}" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.kraken2uniq"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}.kraken2uniq.report"

    input:
    tuple val(id), path(forward), path(reverse)
    each path(kraken2_db)

    output:
    tuple(val(id), path("kmer__kraken2uniq__${idx_name}"), emit: results)
    tuple path("${id}_${idx_name}.kraken2uniq.report"), path("${id}_${idx_name}.kraken2uniq")

    script:
    idx_name = kraken2_db.baseName

    """
    kraken2 --db $kraken2_db --minimum-hit-groups $params.kraken2_minimum_hit_groups --report-minimizer-data --threads $task.cpus --output ${id}_${idx_name}.kraken2uniq --report ${id}_${idx_name}.kraken2uniq.report --paired $forward $reverse
    cp ${id}_${idx_name}.kraken2uniq.report kmer__kraken2uniq__${idx_name}
    """

}

process Kraken2Bracken {

    tag { id }
    label "kraken2bracken"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}_${idx_name}.kraken2bracken.report", saveAs: { "kmer__kraken2bracken__${idx_name}" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.kraken2"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}.kraken2bracken.report"

    input:
    tuple val(id), path(forward), path(reverse)
    each path(kraken2_db)

    output:
    tuple(val(id), path("kmer__kraken2bracken__${idx_name}"), emit: results)
    tuple path("${id}_${idx_name}.kraken2bracken.report"), path("${id}_${idx_name}.kraken2bracken")

    script:
    idx_name = kraken2_db.baseName

    """
    kraken2 --db $kraken2_db --confidence $params.kraken2_confidence --threads $task.cpus --output ${id}_${idx_name}.kraken2 --report ${id}_${idx_name}.kraken2.report --paired $forward $reverse
    bracken -d $kraken2_db -i ${id}_${idx_name}.kraken2.report -r $params.bracken_read_length -l $params.bracken_taxonomic_level -t $params.bracken_read_threshold -o {id}_${idx_name}.kraken2bracken.report -w sample.breport
    cp ${id}_${idx_name}.kraken2bracken.report kmer__kraken2bracken__${idx_name}
    """

}

process Kraken2UniqOnt {

    tag { id }
    label "kraken2uniq"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}_${idx_name}.kraken2uniq.report", saveAs: { "kmer__kraken2uniq__${idx_name}" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.kraken2uniq"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}.kraken2uniq.report"

    input:
    tuple val(id), path(reads)
    each path(kraken2_db)

    output:
    tuple(val(id), path("kmer__kraken2uniq__${idx_name}"), emit: results)
    tuple path("${id}_${idx_name}.kraken2uniq.report"), path("${id}_${idx_name}.kraken2uniq")

    script:
    idx_name = kraken2_db.baseName

    """
    kraken2 --db $kraken2_db --minimum-hit-groups $params.kraken2_minimum_hit_groups --report-minimizer-data --threads $task.cpus --output ${id}_${idx_name}.kraken2uniq --report ${id}_${idx_name}.kraken2uniq.report $reads
    cp ${id}_${idx_name}.kraken2uniq.report kmer__kraken2uniq__${idx_name}
    """

}

process Kraken2BrackenOnt {

    tag { id }
    label "kraken2bracken"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}_${idx_name}.kraken2bracken.report", saveAs: { "kmer__kraken2bracken__${idx_name}" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "symlink", pattern: "${id}_${idx_name}.kraken2"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}.kraken2bracken.report"

    input:
    tuple val(id), path(reads)
    each path(kraken2_db)

    output:
    tuple(val(id), path("kmer__kraken2bracken__${idx_name}"), emit: results)
    tuple path("${id}_${idx_name}.kraken2bracken.report"), path("${id}_${idx_name}.kraken2bracken")

    script:
    idx_name = kraken2_db.baseName

    """
    kraken2 --db $kraken2_db --confidence $params.kraken2_confidence --threads $task.cpus --output ${id}_${idx_name}.kraken2 --report ${id}_${idx_name}.kraken2.report $reads
    bracken -d $kraken2_db -i ${id}_${idx_name}.kraken2.report -r $params.bracken_read_length -l $params.bracken_taxonomic_level -t $params.bracken_read_threshold -o {id}_${idx_name}.kraken2bracken.report -w sample.breport
    cp ${id}_${idx_name}.kraken2bracken.report kmer__kraken2bracken__${idx_name}
    """

}