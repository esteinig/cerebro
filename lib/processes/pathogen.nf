


process Kraken2 {

    tag { sampleID }
    label "pathogenProfileKraken2"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kraken2.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kraken2.report"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(krakenDatabase)
    val(krakenConfidence)

    output:

    tuple (val(sampleID), path(krakenDatabase), path("${sampleID}.kraken2.report"), emit: bracken)
    tuple (val(sampleID), path("${sampleID}.kraken2.report"), path("${sampleID}.kraken2.tsv"), emit: results)

    script:

    """
    kraken2 --db $krakenDatabase --confidence $krakenConfidence --threads $task.cpus --output ${sampleID}.kraken2.tsv --report ${sampleID}.kraken2.report --paired $forward $reverse
    """

}


process Bracken {

    tag { sampleID }
    label "pathogenProfileBracken"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.bracken.report"

    input:
    tuple val(sampleID), path(krakenDatabase), path(krakenReport)
    val(brackenReadLength)
    val(brackenRank)
    val(brackenReads)

    output:
    tuple (val(sampleID), path("${sampleID}.bracken.report"), emit: results)

    script:

    """
    bracken -d $kraken2_db -i $krakenReport -r $brackenReadLength -l $brackenTaxRank -t $brackenMinReads -o ${sampleID}.bracken.report
    """

}


process Sylph {

    tag { sampleID }
    label "pathogenProfileSylph"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.sylph.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.sylphmpa"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(sylphDatabase)
    path(sylphDatabaseMetadata)

    output:
    tuple (path("${sampleID}.sylph.tsv"), emit: results)

    script:

    """
    sylph profile $sylphDatabase -1 $forward -2 $reverse -c 100 --min-number-kmers 20 -t $task.cpu > ${sampleID}.sylph.tsv
    sylph_to_taxprof.py -m $sylphDatabaseMetadata -s ${sampleID}.sylph.tsv -o ${sampleID}.sylphmpa
    """

}


process Metabuli {

    tag { sampleID }
    label "pathogenProfileMetabuli"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metabuli.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metabuli.report"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(metabuliDatabase)

    output:
    tuple (path("${sampleID}.metabuli.report"), path("${sampleID}.metabuli.tsv"), emit: results)

    script:

    memoryLimit = task.memory.split()[0]

    """
    metabuli classify --max-ram $memoryLimit --threads $task.cpus $forward $reverse $metabuliDatabase classified/ $sampleID
    cp classified/${sampleID}_classifications.tsv ${sampleID}.metabuli.tsv
    rm classified/${sampleID}_classifications.tsv
    cp classified/${sampleID}_report.tsv ${sampleID}.metabuli.report
    """

}


process Kmcp {

    tag { sampleID }
    label "pathogenProfileKmcp"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(kmcpDatabase)
    val(kmcpMode)


    script:

    """
    kmcp search -d $kmcpDatabase -1 $forward -2 $reverse -o ${sampleId}.reads.tsv.gz --threads ${task.cpus}
    kmcp profile --taxid-map $kmcpDatabase/taxid.map --taxdump $kmcpDatabase/taxdump --mode $kmcpMode
    """

}