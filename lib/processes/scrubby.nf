// Input 1 and Input 2 are used for R1 and R2 so that sequential depletions (e.g. ERCC -> Host) that would otherwise have the same name are not overwritten

process ScrubbyReadsKrakenDepletion {

    label "scrubby_reads"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.json", saveAs: { "$params.result_file" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.json"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.log"
    publishDir "$params.outdir/workflow/$params.subdir/fastq", mode: "symlink", pattern: "${id}_depleted*.fastq"

    input:
    tuple val(id), path("input_1"), path("input_2")
    path "dbs/*"

    output:
    tuple val(id), path("${id}_depleted_1.fastq"), path("${id}_depleted_2.fastq")
    tuple path("${id}.json"), path("${id}.log")

    script:

    """
    scrubby scrub-reads -i input_1 input_2 -o ${id}_depleted_1.fastq ${id}_depleted_2.fastq --kraken-db dbs/* --kraken-taxa $params.kraken_taxa --kraken-taxa-direct $params.kraken_taxa_direct --kraken-threads $task.cpus --json ${id}.json 2> ${id}.log
    """
}

process ScrubbyReadsKrakenMinimapDepletion {

    label "scrubby_reads"
    tag { "$id" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.json", saveAs: { "$params.result_file" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.json"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.log"
    publishDir "$params.outdir/workflow/$params.subdir/fastq", mode: "symlink", pattern: "${id}_depleted*.fq.gz"

    input:
    tuple val(id), path("input_1"), path("input_2")
    path "dbs/*"
    path "refs/*"

    output:
    tuple (val(id), path("${id}_depleted_1.fq.gz"), path("${id}_depleted_2.fq.gz"), emit: reads)
    tuple (val(id), path("$params.result_file"), emit: results)
    tuple path("${id}.log"), path("${id}.json")

    script:

    """
    scrubby scrub-reads -i input_1 input_2 -o ${id}_depleted_1.fq.gz ${id}_depleted_2.fq.gz --kraken-db dbs/* --kraken-taxa $params.kraken_taxa --kraken-taxa-direct $params.kraken_taxa_direct --kraken-threads $task.cpus --minimap2-index refs/* --min-cov $params.deplete_min_cov --min-len $params.deplete_min_len --min-mapq $params.deplete_min_mapq --minimap2-threads $task.cpus --json ${id}.json 2> ${id}.log
    cp ${id}.json $params.result_file
    """
}


process ScrubbyReadsKrakenMinimapDepletionOnt {

    label "scrubby_reads"
    tag { "$id" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.json", saveAs: { "$params.result_file" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.json"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.log"
    publishDir "$params.outdir/workflow/$params.subdir/fastq", mode: "symlink", pattern: "${id}_depleted.fq.gz"

    input:
    tuple val(id), path("input_1")
    path "dbs/*"
    path "refs/*"

    output:
    tuple (val(id), path("${id}_depleted.fq.gz"), emit: reads)
    tuple (val(id), path("$params.result_file"), emit: results)
    tuple path("${id}.log"), path("${id}.json")

    script:

    """
    scrubby scrub-reads -i input_1 -o ${id}_depleted.fq.gz --kraken-db dbs/* --kraken-taxa $params.kraken_taxa --kraken-taxa-direct $params.kraken_taxa_direct --kraken-threads $task.cpus --minimap2-index refs/* --min-cov $params.deplete_min_cov --min-len $params.deplete_min_len --min-mapq $params.deplete_min_mapq --minimap2-threads $task.cpus --json ${id}.json 2> ${id}.log
    cp ${id}.json $params.result_file
    """
}

process ScrubbyAlignmentDepletion {

    // Auto detection based on SAM/BAM/CRAM/PAF extension

    label "scrubby_single"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.json", saveAs: { "$params.result_file" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.json"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.log"
    publishDir "$params.outdir/workflow/$params.subdir/fastq", mode: "symlink", pattern: "${id}_depleted*.fq.gz"

    input:
    tuple val(id), path("input_1"), path("input_2")
    tuple val(id), val(idx_name), path(alignment)

    output:
    tuple (val(id), path("${id}_depleted_1.fq.gz"), path("${id}_depleted_2.fq.gz"), emit: reads)
    tuple (val(id), path("$params.result_file"), emit: results)
    tuple path("${id}.log"), path("${id}.json")

    script:

    """
    scrubby scrub-alignment -i input_1 input_2 -o ${id}_depleted_1.fq.gz ${id}_depleted_2.fq.gz --alignment $alignment --alignment-name ${idx_name} --min-cov $params.deplete_min_cov --min-len $params.deplete_min_len --min-mapq $params.deplete_min_mapq --json ${id}.json 2> ${id}.log
    cp ${id}.json $params.result_file
    """
}


process ScrubbyAlignmentDepletionOnt {

    // Auto detection based on SAM/BAM/CRAM/PAF extension

    label "scrubby_single"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.json", saveAs: { "$params.result_file" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.json"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.log"
    publishDir "$params.outdir/workflow/$params.subdir/fastq", mode: "symlink", pattern: "${id}_depleted.fq.gz"

    input:
    tuple val(id), path("input_1")
    tuple val(id), val(idx_name), path(alignment)

    output:
    tuple (val(id), path("${id}_depleted.fq.gz"), emit: reads)
    tuple (val(id), path("$params.result_file"), emit: results)
    tuple path("${id}.log"), path("${id}.json")

    script:

    """
    scrubby scrub-alignment -i input_1 -o ${id}_depleted.fq.gz --alignment $alignment --alignment-name ${idx_name} --min-cov $params.deplete_min_cov --min-len $params.deplete_min_len --min-mapq $params.deplete_min_mapq --json ${id}.json 2> ${id}.log
    cp ${id}.json $params.result_file
    """
}

process ScrubbyExtractVircovReads {

    // TXT file of read identifiers is supported by Scrubby

    label "scrubby_single"
    tag { "$id : $db_name : $idx_name" }

    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: { "${id}_${idx_name}.json" }
    publishDir "$params.outdir/workflow/$params.subdir/fastq", mode: "symlink", pattern: "${id}_${idx_name}_extracted*.fq.gz"

    input:
    tuple val(id), path("input_1"), path("input_2")
    tuple val(id), val(idx_name), path(read_path)
    tuple val(id), val(db_name), path(vircov_path), path(ref_fasta), path(ref_alignment)

    output:
    tuple(val(id), val(db_name), val(idx_name), path(ref_fasta), path(ref_alignment), path(vircov_path), path("${id}_${idx_name}_extracted_1.fq.gz"), path("${id}_${idx_name}_extracted_2.fq.gz"), emit: aligned)
    path("${id}_${idx_name}.json")

    script:

    """
    scrubby scrub-alignment --alignment $read_path --alignment-name ${idx_name} -i input_1 input_2 -o ${id}_${idx_name}_extracted_1.fq.gz ${id}_${idx_name}_extracted_2.fq.gz --extract --json ${id}_${idx_name}.json
    """
    
}

process ScrubbyExtractVircovReadsOnt {

    // TXT file of read identifiers is supported by Scrubby

    label "scrubby_single"
    tag { "$id : $db_name : $idx_name" }

    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: { "${id}_${idx_name}.json" }
    publishDir "$params.outdir/workflow/$params.subdir/fastq", mode: "symlink", pattern: "${id}_${idx_name}_extracted.fq.gz"

    input:
    tuple val(id), path("input_1")
    tuple val(id), val(idx_name), path(read_path)
    tuple val(id), val(db_name), path(vircov_path), path(ref_fasta), path(ref_alignment)

    output:
    tuple(val(id), val(db_name), val(idx_name), path(ref_fasta), path(ref_alignment), path(vircov_path), path("${id}_${idx_name}_extracted.fq.gz"), emit: aligned)
    path("${id}_${idx_name}.json")

    script:

    """
    scrubby scrub-alignment --alignment $read_path --alignment-name ${idx_name} -i input_1 -o ${id}_${idx_name}_extracted.fq.gz --extract --json ${id}_${idx_name}.json
    """
    
}