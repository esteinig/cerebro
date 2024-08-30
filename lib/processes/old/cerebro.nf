// Always uses environment variable for authentication token

process PingServer {

    tag { "$params.production.api.url/status" }
    label "cerebro"
    
    input:
    val(results)  // Wait on result files before pinging server in production pipeline

    output:
    val("SUCCESS")

    script:

    """
    cerebro --token-env $params.production.api.token --api-url $params.production.api.url client ping-api
    """

}

process ProcessSamplesTaxonomy {

    publishDir "$params.outdir/cerebro", mode: "copy", pattern: "${id}.json"

    label "cerebro"
    
    input:
    tuple val(id), path(results)
    path(taxonomy)

    output:
    tuple(val(id), path("${id}.json"), emit: cerebro)

    script:

    """
    cerebro workflow process --input . --sample-id ${id} --taxonomy $taxonomy --output ${id}.json 
    """
}

process ProcessSamples {

    publishDir "$params.outdir/cerebro", mode: "copy", pattern: "${id}.json"

    label "cerebro"
    
    input:
    tuple val(id), path(results)

    output:
    tuple(val(id), path("${id}.json"), emit: cerebro)

    script:

    """
    cerebro workflow process --input . --sample-id ${id} --output ${id}.json 
    """
    
}



process QualityControlTable {

    publishDir "$params.outdir", mode: "copy", pattern: "qc.tsv"

    label "cerebro"
    
    input:
    path(json_files)
    
    output:
    path("qc.tsv")

    script:

    ercc_mass = params.process.ercc_mass ? "--ercc-mass $params.process.ercc_mass" : ""

    """
    cerebro workflow quality --input *.json --output qc.tsv --header $ercc_mass
    """
  

}

process UploadSample {

    publishDir "$params.outdir/cerebro", mode: "copy", pattern: "${id}.json"

    tag { "$params.cerebro_url/cerebro" }
    label "cerebro"
    
    input:
    tuple val(id), path(sample_json)
    path(sample_sheet)
    path(config_json)

    output:
    val("SUCCESS")

    script:

    """
    cerebro --token-env  $params.production.api.token --api-url $params.production.api.url client upload --input $sample_json --sample-sheet $sample_sheet --workflow-config $config_json --team-name $params.production.api.upload.team_name --project-name $params.production.api.upload.project_name
    """

}


process AlignmentTable {

    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.tsv"
    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.tsv", saveAs: { "align__vircov_remap__${db_name}" }

    tag { "$id : $db_name" }
    label "cerebro"
    
    input:
    tuple val(id), val(db_name), path(files_1), path(files_2), path(files_2)

    output:
    tuple (val(id), path("${id}.tsv"), emit: results)

    script:

    // Consensus files may not exist - accounts for empty file input

    """

    if [ ! -f "align__vircov__${db_name}__scan" ]; then
        touch align__vircov__${db_name}__scan
    fi

    if [ ! -f "align__vircov__${db_name}__remap" ]; then
        touch align__vircov__${db_name}__remap
    fi
        
    fasta_count=\$(find . -name "*.fasta" | wc -l)

    if [ \$fasta_count -gt 0 ]; then
        cov_count=\$(find . -name "*.txt" | wc -l)

        cat *.consensus.fasta > ${id}.consensus.fasta

        if [ \$cov_count -gt 0 ]; then
            cat *.txt > ${id}.coverage.tsv
            cerebro workflow tools scan-remap --id ${id} --db ${db_name} --scan align__vircov__${db_name}__scan --remap align__vircov__${db_name}__remap --consensus ${id}.consensus.fasta --coverage ${id}.coverage.tsv --output ${id}.tsv -H
        else
            cerebro workflow tools scan-remap --id ${id} --db ${db_name} --scan align__vircov__${db_name}__scan --remap align__vircov__${db_name}__remap --consensus ${id}.consensus.fasta --output ${id}.tsv -H
        fi

    else
        cov_count=\$(find . -name "*.txt" | wc -l)

        if [ \$cov_count -gt 0 ]; then
            cat *.txt > ${id}.coverage.tsv
            cerebro workflow tools scan-remap --id ${id} --db ${db_name} --scan align__vircov__${db_name}__scan --remap align__vircov__${db_name}__remap --coverage ${id}.coverage.tsv --output ${id}.tsv -H
        else
            cerebro workflow tools scan-remap --id ${id} --db ${db_name} --scan align__vircov__${db_name}__scan --remap align__vircov__${db_name}__remap --output ${id}.tsv -H
        fi

    fi    
    cp ${id}.tsv align__vircov_remap__${db_name}
    """

}

    
process MashDatabaseSubset {

    label "cerebro"
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
        cerebro workflow tools subset-fasta --fasta $fasta --mash $mash_screen --min-identity 0 --min-shared-hashes $params.taxa.alignment.subset.min_shared_hashes --group-index $params.taxa.alignment.subset.group_index --group-sep "$params.taxa.alignment.subset.group_sep" --output ${idx_name}_subset.fasta 
        """
    } else {
        """
        cerebro workflow tools subset-fasta --fasta $fasta --mash $mash_screen --min-identity 0  --min-shared-hashes $params.taxa.alignment.subset.min_shared_hashes --output ${idx_name}_subset.fasta
        """
    }

}

process MashDatabaseSubsetOnt {

    label "cerebro"
    tag { id }

    publishDir "$params.outdir/workflow/$params.subdir/subsets/${id}", mode: "symlink", pattern: "${idx_name}_subset.fasta"

    input:
    tuple val(id), val(idx_name), path(mash_screen)  
    path(fasta)

    output:
    tuple val(id), val(idx_name), path("${idx_name}_subset.fasta")

    script: 

    if (params.subset_group_index) {
        """
        cerebro workflow tools subset-fasta --fasta $fasta --mash $mash_screen --min-identity 0 --min-shared-hashes $params.taxa.alignment.subset.min_shared_hashes --group-index $params.taxa.alignment.subset.group_index --group-sep "$params.taxa.alignment.subset.group_sep" --output ${idx_name}_subset.fasta 
        """
    } else {
        """
        cerebro workflow tools subset-fasta --fasta $fasta --mash $mash_screen --min-identity 0 --min-shared-hashes $params.taxa.alignment.subset.min_shared_hashes --output ${idx_name}_subset.fasta
        """
    }
}