process DiamondNR {

    label "diamond_nr"
    tag { id }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}_ncbi_nr.tsv", saveAs: { "assembly__diamond__ncbi_nr" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_ncbi_nr.tsv"

    input:
    tuple val(id), path(contigs)
    path(database)

    output:
    tuple(val(id), path("assembly__diamond__ncbi_nr"), emit: results)
    path("${id}_ncbi_nr.tsv")

    script:

    // TODO: Empty file guard from assembly

    """
    diamond blastx -d $database -p $task.cpus -q $contigs -f 6 qseqid qlen qstart qend sseqid slen sstart send length nident pident evalue bitscore staxids sscinames stitle \
    --max-target-seqs $params.meta_diamond_nr_max_seqs --index-chunks $params.meta_diamond_nr_index_chunks --block-size $params.meta_diamond_nr_block_size --evalue $params.meta_diamond_nr_min_evalue --id $params.meta_diamond_nr_min_identity -o ${id}_ncbi_nr.tsv
    cp ${id}_ncbi_nr.tsv assembly__diamond__ncbi_nr 
    """ 

}