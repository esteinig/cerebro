process VircovReferenceSelection {

    tag { "$id : $db_name" }
    label "vircov"

    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_scan_${db_name}.select.tsv"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_scan_${db_name}.grouped.tsv"

    input:
    tuple val(id), path(forward), path(reverse)
    tuple val(id_idx), val(db_name), path(alignment)
    path(alignment_fasta)
    path(blacklist)

    output:
    tuple (val(id), path(forward), path(reverse), emit: reads)
    tuple (val(id), path(forward), path(reverse), val(db_name), path("references/*.fasta"), emit: references) optional true 
    tuple (val(id), val(db_name), path("align__vircov__${db_name}__scan"), emit: results)
    path("${id}_scan_${db_name}.select.tsv")

    script:
    
    segment_field = params.vircov_group_select_segment_field ? "--segment-field '$params.vircov_group_select_segment_field'" : ""
    segment_field_nan = params.vircov_group_select_segment_field_nan ? "--segment-field-nan '$params.vircov_group_select_segment_field_nan'" : ""
    exclude = blacklist ? "--exclude $params.virus_blacklist" : ""

    """
    vircov coverage --alignment $alignment --fasta $alignment_fasta --min-len $params.vircov_scan_min_len --min-cov $params.vircov_scan_min_cov --min-mapq $params.vircov_scan_min_mapq --reads $params.vircov_scan_reads --coverage $params.vircov_scan_coverage --regions $params.vircov_scan_regions --regions-coverage $params.vircov_scan_regions_coverage --group-by "$params.vircov_group_by" --group-sep "$params.vircov_group_sep" --group-select-by "$params.vircov_group_select_by" --group-select-split references --group-select-order --group-select-data "${id}_scan_${db_name}.grouped.tsv" $segment_field $segment_field_nan $exclude > ${id}_scan_${db_name}.select.tsv
    cp ${id}_scan_${db_name}.select.tsv align__vircov__${db_name}__scan
    """
    
}

process VircovReferenceSelectionOnt {

    tag { "$id : $db_name" }
    label "vircov"

    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_scan_${db_name}.select.tsv"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_scan_${db_name}.grouped.tsv"

    input:
    tuple val(id), path(reads)
    tuple val(id_idx), val(db_name), path(alignment)
    path(alignment_fasta)
    path(blacklist)

    output:
    tuple (val(id), path(reads), emit: reads)
    tuple (val(id), path(reads), val(db_name), path("references/*.fasta"), emit: references) optional true 
    tuple (val(id), val(db_name), path("align__vircov__${db_name}__scan"), emit: results)
    path("${id}_scan_${db_name}.select.tsv")

    script:
    
    segment_field = params.vircov_group_select_segment_field ? "--segment-field '$params.vircov_group_select_segment_field'" : ""
    segment_field_nan = params.vircov_group_select_segment_field_nan ? "--segment-field-nan '$params.vircov_group_select_segment_field_nan'" : ""
    exclude = blacklist ? "--exclude $params.virus_blacklist" : ""

    """
    vircov coverage --alignment $alignment --fasta $alignment_fasta --min-len $params.vircov_scan_min_len --min-cov $params.vircov_scan_min_cov --min-mapq $params.vircov_scan_min_mapq --reads $params.vircov_scan_reads --coverage $params.vircov_scan_coverage --regions $params.vircov_scan_regions --regions-coverage $params.vircov_scan_regions_coverage --group-by "$params.vircov_group_by" --group-sep "$params.vircov_group_sep" --group-select-by "$params.vircov_group_select_by" --group-select-split references --group-select-order --group-select-data "${id}_scan_${db_name}.grouped.tsv" $segment_field $segment_field_nan $exclude > ${id}_scan_${db_name}.select.tsv
    cp ${id}_scan_${db_name}.select.tsv align__vircov__${db_name}__scan
    """
    
}


process VircovRealign {

    tag { "$id : $idx_name : $db_name" }
    label "vircov"

    publishDir "$params.outdir/workflow/$params.subdir/reads", mode: "copy", pattern: "${id}_${idx_name}.txt"
    publishDir "$params.outdir/workflow/$params.subdir/references/${id}", mode: "copy", pattern: "*.fasta"

    input:
    tuple val(id), path(forward), path(reverse)
    tuple val(id_idx), val(db_name), val(idx_name), path(fasta), path(alignment)

    output:
    tuple (val(id), path(forward), path(reverse), emit: reads)
    tuple (val(id), val(idx_name), path("${id}_${idx_name}.txt"), emit: read_ids)
    tuple (val(id), val(db_name), path("${id}_${idx_name}.tsv"), path(fasta), path(alignment), emit: aligned)

    script:
    
    """
    vircov coverage --alignment $alignment --fasta $fasta --min-len $params.vircov_remap_min_len --min-cov $params.vircov_remap_min_cov --min-mapq $params.vircov_remap_min_mapq --reads $params.vircov_remap_reads --coverage $params.vircov_remap_coverage --regions $params.vircov_remap_regions --regions-coverage $params.vircov_remap_regions_coverage --read-ids ${id}_${idx_name}.txt -v > ${id}_${idx_name}.tsv
    """

}


process VircovRealignOnt {

    tag { "$id : $idx_name : $db_name" }
    label "vircov"

    publishDir "$params.outdir/workflow/$params.subdir/reads", mode: "copy", pattern: "${id}_${idx_name}.txt"
    publishDir "$params.outdir/workflow/$params.subdir/references/${id}", mode: "copy", pattern: "*.fasta"

    input:
    tuple val(id), path(reads)
    tuple val(id_idx), val(db_name), val(idx_name), path(fasta), path(alignment)

    output:
    tuple (val(id), path(reads), emit: reads)
    tuple (val(id), val(idx_name), path("${id}_${idx_name}.txt"), emit: read_ids)
    tuple (val(id), val(db_name), path("${id}_${idx_name}.tsv"), path(fasta), path(alignment), emit: aligned)

    script:
    
    """
    vircov coverage --alignment $alignment --fasta $fasta --min-len $params.vircov_remap_min_len --min-cov $params.vircov_remap_min_cov --min-mapq $params.vircov_remap_min_mapq --reads $params.vircov_remap_reads --coverage $params.vircov_remap_coverage --regions $params.vircov_remap_regions --regions-coverage $params.vircov_remap_regions_coverage --read-ids ${id}_${idx_name}.txt -v > ${id}_${idx_name}.tsv
    """

}

process ConcatenateVircov {
    
    // Unlimited number of Vircov output files to concatenate using `find`

    tag { "$id" }
    label "vircov"

    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_remap_${db_name}.tsv"
    
    input:
    tuple val(id), val(db_name), path(vircov_files), path(fasta_files), path(bam_files)

    output:
    tuple (val(id), val(db_name), path("align__vircov__${db_name}__remap"), emit: results)
    path("${id}_remap_${db_name}.tsv")
    
    """
    vircov=\$(find . -name "${id}*.tsv")
    cat \$vircov > ${id}_remap_${db_name}.tsv
    cp ${id}_remap_${db_name}.tsv align__vircov__${db_name}__remap
    """

}

process VircovZero {
    
    // Simple assessment of aligned reads with zero counts, needs reference sequences

    tag { "$id : $idx_name : $alignment" }
    label "vircov"

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.tsv", saveAs: { "$params.result_file" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), path(forward), path(reverse)
    tuple val(id_idx), val(idx_name), path(alignment)
    path(fasta)

    output:
    tuple (val(id), path(forward), path(reverse), emit: reads)
    tuple (val(id), path("$params.result_file"), emit: results)
    path("${id}.tsv")c
    
    script:
    
    """
    vircov coverage --alignment $alignment --fasta $fasta --zero -v > ${id}.tsv
    cp ${id}.tsv $params.result_file
    """
    
}

process VircovAlignReferenceZero {
    
    // Simple assessment of aligned reads with zero counts, needs reference sequences

    tag { "$id : $idx_name : $alignment" }
    label "vircov"

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.tsv", saveAs: { "$params.result_file" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), path(forward), path(reverse)
    tuple val(id_idx), val(idx_name), path(alignment), path(fasta)

    output:
    tuple (val(id), path(forward), path(reverse), emit: reads)
    tuple (val(id), path("$params.result_file"), emit: results)
    path("${id}.tsv")
    
    script:
    
    """
    vircov coverage --alignment $alignment --fasta $fasta --zero -v > ${id}.tsv
    cp ${id}.tsv $params.result_file
    """
    
}



process VircovZeroOnt {
    
    // Simple assessment of aligned reads with zero counts, needs reference sequences

    tag { "$id : $idx_name : $alignment" }
    label "vircov"

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.tsv", saveAs: { "$params.result_file" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.tsv"

    input:
    tuple val(id), path(reads)
    tuple val(id_idx), val(idx_name), path(alignment)
    path(fasta)

    output:
    tuple (val(id), path(reads), emit: reads)
    tuple (val(id), path("$params.result_file"), emit: results)
    path("${id}.tsv")
    
    script:
    
    """
    vircov coverage --alignment $alignment --fasta $fasta --zero -v > ${id}.tsv
    cp ${id}.tsv $params.result_file
    """
    
}


process VircovSubsetAlign {

    tag { "$id : $idx_name" }
    label "vircov"

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}_${idx_name}_subset.tsv", saveAs: { "align__vircov__${idx_name}__subset" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}_subset.tsv"

    input:
    tuple val(id), path(forward), path(reverse)
    tuple val(id_idx), val(idx_name), path(alignment), path(fasta)

    output:
    tuple(val(id), path(forward), path(reverse), emit: reads)
    tuple(val(id), path("align__vircov__${idx_name}__subset"), emit: results)
    path("${id}_${idx_name}_subset.tsv")

    script:
    
    """
    vircov coverage --alignment $alignment --fasta $fasta --min-len $params.vircov_min_len --min-cov $params.vircov_min_cov --min-mapq $params.vircov_min_mapq --reads $params.vircov_min_reads --regions $params.vircov_min_regions -v > ${id}_${idx_name}_subset.tsv
    cp ${id}_${idx_name}_subset.tsv align__vircov__${idx_name}__subset
    """

}



process VircovSubsetAlignOnt {

    tag { "$id : $idx_name" }
    label "vircov"

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}_${idx_name}_subset.tsv", saveAs: { "align__vircov__${idx_name}__subset" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}_subset.tsv"

    input:
    tuple val(id), path(reads)
    tuple val(id_idx), val(idx_name), path(alignment), path(fasta)

    output:
    tuple(val(id), path(reads), emit: reads)
    tuple(val(id), path("align__vircov__${idx_name}__subset"), emit: results)
    path("${id}_${idx_name}_subset.tsv")

    script:
    
    """
    vircov coverage --alignment $alignment --fasta $fasta --min-len $params.vircov_min_len --min-cov $params.vircov_min_cov --min-mapq $params.vircov_min_mapq --reads $params.vircov_min_reads --regions $params.vircov_min_regions -v > ${id}_${idx_name}_subset.tsv
    cp ${id}_${idx_name}_subset.tsv align__vircov__${idx_name}__subset
    """

}
