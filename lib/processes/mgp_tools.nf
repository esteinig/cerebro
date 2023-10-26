process Anonymize {

    label "mgp_tools"
    tag { id }

    publishDir "$params.outdir/workflow/$params.subdir/${id}", mode: "symlink", pattern: "*.fq.gz"

    input:
    tuple val(id), path(forward), path(reverse)
    
    output:
    tuple val(id), path("${uuid}_R1.fq.gz"), path("${uuid}_R2.fq.gz")

    script: 

    uuid = UUID.randomUUID().toString().substring(0, 8)

    """
    cerebro tools utils anonymize -i $forward -i $reverse -o ${uuid}_R1.fq.gz -o ${uuid}_R2.fq.gz --fake-illumina-header
    """

}

process MashDatabaseSubset {

    label "mgp_tools"
    tag { id }

    publishDir "$params.outdir/workflow/$params.subdir/subsets/${id}", mode: "symlink", pattern: "${idx_name}_subset.fasta"

    input:
    tuple val(id), val(idx_name), path(forward), path(reverse), path(mash_screen)  
    path(fasta)

    output:
    tuple val(id), val(idx_name), path(forward), path(reverse), path("${idx_name}_subset.fasta")

    script: 

    if (params.subset_group_index) {
        """
        cerebro tools utils subset --fasta $fasta --mash $mash_screen --min-identity 0 --min-shared-hashes $params.min_shared_hashes --group-index $params.subset_group_index --group-sep "$params.subset_group_sep" --output ${idx_name}_subset.fasta 
        """
    } else {
        """
        cerebro tools utils subset --fasta $fasta --mash $mash_screen --min-identity 0 --min-shared-hashes $params.min_shared_hashes --output ${idx_name}_subset.fasta
        """
    }

}