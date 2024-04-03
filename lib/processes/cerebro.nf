// Always uses environment variable for authentication token

process PingServer {

    tag { "$params.cerebro_api_url/status" }
    label "cerebro"
    
    input:
    val(results)  // Wait on result files before pinging server in production pipeline

    output:
    val("SUCCESS")

    script:

    """
    cerebro --token-env $params.cerebro_token_env --api-url $params.cerebro_api_url client ping-api
    """

}

process ParseSample {

    publishDir "$params.outdir/cerebro", mode: "copy", pattern: "${id}.json"

    label "cerebro"
    
    input:
    tuple val(id), path(results)
    path(taxonomy)

    output:
    tuple(val(id), path("${id}.json"), emit: cerebro)

    script:

    """
    cerebro workflow parse-sample --input . --sample-id ${id} --taxonomy $taxonomy --output ${id}.json 
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
    cerebro --token-env $params.cerebro_token_env --api-url $params.cerebro_api_url client upload --input $sample_json --sample-sheet $sample_sheet --workflow-config $config_json --team-name $params.cerebro_team_name --project-name $params.cerebro_project_name
    """

}


process ProcessVirusAlignment {

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

    