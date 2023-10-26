process BlastNT {

    label "blast_nt"
    tag { id }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}_ncbi_nt.tsv", saveAs: { "assembly__blastn__ncbi_nt" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_ncbi_nt.tsv"

    input:
    tuple val(id), path(contigs)
    path(database)

    output:
    tuple(val(id), path("assembly__blastn__ncbi_nt"), emit: results)
    path("${id}_ncbi_nt.tsv")

    script:

    // Set the execution environment variable BLASTDB to database path to enable the taxonomic assignments

    """
    BLASTDB=$database blastn -num_threads $task.cpus -query $contigs -perc_identity $params.meta_blast_nt_min_identity -evalue $params.meta_blast_nt_min_evalue -db ${database}/nt \
    -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length nident pident evalue bitscore staxid ssciname stitle' -max_target_seqs $params.meta_blast_nt_max_seqs > ${id}_ncbi_nt.tsv
    cp ${id}_ncbi_nt.tsv assembly__blastn__ncbi_nt
    """

}