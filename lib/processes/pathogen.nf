process Kraken2 {

    tag { sampleID }
    label "pathogenProfileKraken2"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kraken2.reads.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kraken2.reads.report"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(krakenDatabase)
    val(krakenConfidence)
    val(krakenMemoryMapping)

    output:

    tuple (val(sampleID), path(krakenDatabase), path("${sampleID}.kraken2.reads.report"), emit: bracken)
    tuple (val(sampleID), path("${sampleID}.kraken2.reads.report"), path("${sampleID}.kraken2.reads.tsv"), emit: results)

    script:

    memoryMapping = krakenMemoryMapping ? "--memory-mapping" : ""

    """
    kraken2 --db $krakenDatabase --confidence $krakenConfidence --threads $task.cpus --output ${sampleID}.kraken2.reads.tsv --report ${sampleID}.kraken2.reads.report $memoryMapping --paired $forward $reverse 

    if [ ! -f "${sampleID}.kraken2.reads.tsv" ]; then
        touch "${sampleID}.kraken2.reads.tsv"
    fi

    """

}

process Kraken2Nanopore {

    tag { sampleID }
    label "pathogenProfileKraken2"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kraken2.reads.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kraken2.reads.report"

    input:
    tuple val(sampleID), path(reads)
    path(krakenDatabase)
    val(krakenConfidence)
    val(krakenMemoryMapping)

    output:

    tuple (val(sampleID), path(krakenDatabase), path("${sampleID}.kraken2.reads.report"), emit: bracken)
    tuple (val(sampleID), path("${sampleID}.kraken2.reads.report"), path("${sampleID}.kraken2.reads.tsv"), emit: results)

    script:

    memoryMapping = krakenMemoryMapping ? "--memory-mapping" : ""

    """
    kraken2 --db $krakenDatabase --confidence $krakenConfidence --threads $task.cpus --output ${sampleID}.kraken2.reads.tsv --report ${sampleID}.kraken2.reads.report $memoryMapping $reads
    
    if [ ! -f "${sampleID}.kraken2.reads.tsv" ]; then
        touch "${sampleID}.kraken2.reads.tsv"
    fi
    """

}


process Bracken {

    tag { sampleID }
    label "pathogenProfileBracken"

    errorStrategy 'ignore'  // empty files e.g. BLANK or those with < brackenMinReads at brackenRank

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.bracken.abundance.report"

    input:
    tuple val(sampleID), path(krakenDatabase), path(krakenReport)
    val(brackenReadLength)
    val(brackenRank)
    val(brackenMinReads)

    output:
    tuple (val(sampleID), path("${sampleID}.bracken.abundance.report"), emit: results)

    script:

    """
    bracken -d $krakenDatabase -i $krakenReport -r $brackenReadLength -l $brackenRank -t $brackenMinReads -o ${sampleID}.bracken.abundance.report
    """
}


process Sylph {

    tag { sampleID }
    label "pathogenProfileSylph"

    errorStrategy 'ignore'  // empty files e.g. BLANK

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.sylph.abundance.report"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(sylphDatabase)
    path(sylphMetadata)
    val(sylphMinNumberKmers)
    val(sylphQueryCompression)

    output:
    tuple (val(sampleID), path("${sampleID}.sylph.abundance.report"), emit: results)

    script:

    """
    sylph profile $sylphDatabase -1 $forward -2 $reverse -c $sylphQueryCompression --min-number-kmers $sylphMinNumberKmers -t $task.cpus > ${sampleID}.sylph.tsv
    python $baseDir/lib/scripts/sylph_to_taxprof.py -m $sylphMetadata -s ${sampleID}.sylph.tsv -o "" 
    
    if [ -f "${forward}.sylphmpa" ]; then
        tail -n +2 ${forward}.sylphmpa > ${sampleID}.sylph.abundance.report
    else
        touch "${sampleID}.sylph.abundance.report"
    fi

    """

}


process SylphNanopore {

    tag { sampleID }
    label "pathogenProfileSylph"

    errorStrategy 'ignore'  // empty files e.g. BLANK

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.sylph.abundance.report"

    input:
    tuple val(sampleID), path(reads)
    path(sylphDatabase)
    path(sylphMetadata)
    val(sylphMinNumberKmers)
    val(sylphQueryCompression)

    output:
    tuple (val(sampleID), path("${sampleID}.sylph.abundance.report"), emit: results)

    script:

    """
    sylph profile $sylphDatabase $reads -c $sylphQueryCompression --min-number-kmers $sylphMinNumberKmers -t $task.cpus > ${sampleID}.sylph.tsv
    python $baseDir/lib/scripts/sylph_to_taxprof.py -m $sylphMetadata -s ${sampleID}.sylph.tsv -o "" 
    
    if [ -f "${reads}.sylphmpa" ]; then
        tail -n +2 ${reads}.sylphmpa > ${sampleID}.sylph.abundance.report
    else
        touch "${sampleID}.sylph.abundance.report"
    fi

    """

}


process Metabuli {

    tag { sampleID }
    label "pathogenProfileMetabuli"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metabuli.reads.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metabuli.reads.report"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(metabuliDatabase)

    output:
    tuple (val(sampleID), path("${sampleID}.metabuli.reads.report"), path("${sampleID}.metabuli.reads.tsv"), emit: results)

    script:

    memoryLimit = "${task.memory}".split()[0]

    """
    if [[ $forward == *.gz ]]; then
        line_count=\$(zcat "$forward" | wc -l)
    else
        line_count=\$(wc -l < "$forward")
    fi

    if [[ \$line_count -eq 0 ]]; then
        touch ${sampleID}.metabuli.reads.tsv
        touch ${sampleID}.metabuli.reads.report
    else
        metabuli classify --max-ram $memoryLimit --threads $task.cpus $forward $reverse $metabuliDatabase classified/ $sampleID
        cp classified/${sampleID}_classifications.tsv ${sampleID}.metabuli.reads.tsv
        cp classified/${sampleID}_report.tsv ${sampleID}.metabuli.reads.report
    fi
    """
    

}


process MetabuliNanopore {

    tag { sampleID }
    label "pathogenProfileMetabuli"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metabuli.reads.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metabuli.reads.report"

    input:
    tuple val(sampleID), path(reads)
    path(metabuliDatabase)

    output:
    tuple (val(sampleID), path("${sampleID}.metabuli.reads.report"), path("${sampleID}.metabuli.reads.tsv"), emit: results)

    script:

    memoryLimit = "${task.memory}".split()[0]

    """
    metabuli classify --max-ram $memoryLimit --threads $task.cpus --seq-mode 3 $reads $metabuliDatabase classified/ $sampleID

    cp classified/${sampleID}_classifications.tsv ${sampleID}.metabuli.reads.tsv
    cp classified/${sampleID}_report.tsv ${sampleID}.metabuli.reads.report
    """

}



process Kmcp { 

    tag { sampleID }
    label "pathogenProfileKmcp"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kmcp.reads.report"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kmcp.abundance.report"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kmcp.reads.tsv.gz"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(kmcpDatabase)
    val(kmcpMode)
    val(kmcpLevel)
    val(kmcpMinQueryCoverage)

    output:
    tuple (val(sampleID), path("${sampleID}.kmcp.reads.tsv.gz"), path("${sampleID}.kmcp.reads.report"), path("${sampleID}.kmcp.abundance.report"), emit: results)


    script:

    """
    kmcp search -d $kmcpDatabase -1 $forward -2 $reverse -o ${sampleID}.reads.tsv.gz --threads $task.cpus
    kmcp profile --level $kmcpLevel --taxid-map $kmcpDatabase/taxid.map --taxdump $kmcpDatabase/taxonomy --mode $kmcpMode -t $kmcpMinQueryCoverage -o ${sampleID}.k.report -B ${sampleID} -C ${sampleID} ${sampleID}.reads.tsv.gz
    

    if [ -f "${sampleID}.binning.gz" ]; then
        mv ${sampleID}.binning.gz ${sampleID}.kmcp.reads.tsv.gz
    else
        touch "${sampleID}.kmcp.reads.tsv.gz"
    fi

    if [ -f "${sampleID}.profile" ]; then
        tail -n +6 ${sampleID}.profile > ${sampleID}.kmcp.abundance.report
    else
        touch "${sampleID}.kmcp.abundance.report"
    fi
    
    if [ -f "${sampleID}.k.report" ]; then
        mv ${sampleID}.k.report ${sampleID}.kmcp.reads.report
    else
        touch "${sampleID}.kmcp.reads.report"
    fi

    """

}


process KmcpNanopore { 

    tag { sampleID }
    label "pathogenProfileKmcp"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kmcp.reads.report"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kmcp.abundance.report"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kmcp.reads.tsv.gz"

    input:
    tuple val(sampleID), path(reads)
    path(kmcpDatabase)
    val(kmcpMode)
    val(kmcpLevel)
    val(kmcpMinQueryCoverage)

    output:
    tuple (val(sampleID), path("${sampleID}.kmcp.reads.tsv.gz"), path("${sampleID}.kmcp.reads.report"), path("${sampleID}.kmcp.abundance.report"), emit: results)


    script:

    """
    seqkit sliding -s 100 -W 300 -o ${sampleID}.reads.short.fastq $reads
    kmcp search -d $kmcpDatabase ${sampleID}.reads.short.fastq -o ${sampleID}.reads.tsv.gz --threads $task.cpus
    kmcp profile --level $kmcpLevel --taxid-map $kmcpDatabase/taxid.map --taxdump $kmcpDatabase/taxonomy --mode $kmcpMode -t $kmcpMinQueryCoverage -o ${sampleID}.k.report -B ${sampleID} -C ${sampleID} ${sampleID}.reads.tsv.gz
    mv ${sampleID}.binning.gz ${sampleID}.kmcp.reads.tsv.gz
    tail -n +6 ${sampleID}.profile > ${sampleID}.kmcp.abundance.report
    mv ${sampleID}.k.report ${sampleID}.kmcp.reads.report
    """

}

process GanonReads { 

    tag { sampleID }
    label "pathogenProfileGanon"

    errorStrategy 'ignore'  // empty files e.g. BLANK

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.ganon.reads.report"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.ganon.reads.tsv"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(ganonDatabase)
    val(ganonDatabasePrefix)
    val(ganonMultipleMatches)

    output:
    tuple (val(sampleID), path("${sampleID}.ganon.reads.report"), path("${sampleID}.ganon.reads.tsv"), emit: results)


    script:

    // Sequence abundance configuration for report (--binning)

    """
    ganon classify --db-prefix $ganonDatabase/$ganonDatabasePrefix --paired-reads $forward $reverse --output-prefix $sampleID --threads $task.cpus --binning --multiple-matches $ganonMultipleMatches --output-one

    mv ${sampleID}.tre ${sampleID}.ganon.reads.report
    mv ${sampleID}.one ${sampleID}.ganon.reads.tsv
    """

}

process RemapKmerReads {


    tag { sampleID }
    label "KmerRemap"


    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple val(sampleID), path(report_file), path(reads_file)
    path(remapKmerDatabase)
    path(taxonomy)


    output:
    tuple (val(sampleID), path(forward), path(reverse), emit: reads)
    tuple (val(sampleID), path("${sampleID}.vircov.tsv"), emit: results)

    script:

    """
    cerebro-pipeline pathogen remap-kmer-reads --reads-file $reads_file --report-file $report_file --r1 $forward --r2 $reverse --taxonomy $taxonomy --rank genus --db $remapKmerDatabase
    """

}


process GanonReadsNanopore { 

    tag { sampleID }
    label "pathogenProfileGanon"

    errorStrategy 'ignore'  // empty files e.g. BLANK

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.ganon.reads.report"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.ganon.reads.tsv"

    input:
    tuple val(sampleID), path(reads)
    path(ganonDatabase)
    val(ganonDatabasePrefix)
    val(ganonMultipleMatches)

    output:
    tuple (val(sampleID), path("${sampleID}.ganon.reads.report"), path("${sampleID}.ganon.reads.tsv"), emit: results)


    script:

    // Sequence abundance configuration for report (--binning)

    """
    ganon classify --db-prefix $ganonDatabase/$ganonDatabasePrefix --single-reads $reads --output-prefix $sampleID --threads $task.cpus --binning --multiple-matches $ganonMultipleMatches --output-one

    mv ${sampleID}.tre ${sampleID}.ganon.reads.report
    mv ${sampleID}.one ${sampleID}.ganon.reads.tsv
    """

}

process GanonProfile { 

    tag { sampleID }
    label "pathogenProfileGanon"

    errorStrategy 'ignore'  // empty files e.g. BLANK

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.ganon.abundance.report"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(ganonDatabase)
    val(ganonDatabasePrefix)
    val(ganonMultipleMatches)

    output:
    tuple (val(sampleID), path("${sampleID}.ganon.abundance.report"), emit: results)


    script:

    // Taxonomic abundance configuration for report (default)

    """
    ganon classify --db-prefix $ganonDatabase/$ganonDatabasePrefix --paired-reads $forward $reverse --output-prefix $sampleID --threads $task.cpus --multiple-matches $ganonMultipleMatches

    mv ${sampleID}.tre ${sampleID}.ganon.abundance.report
    """

}

process GanonProfileNanopore { 

    tag { sampleID }
    label "pathogenProfileGanon"

    errorStrategy 'ignore'  // empty files e.g. BLANK

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.ganon.abundance.report"

    input:
    tuple val(sampleID), path(reads)
    path(ganonDatabase)
    val(ganonDatabasePrefix)
    val(ganonMultipleMatches)

    output:
    tuple (val(sampleID), path("${sampleID}.ganon.abundance.report"), emit: results)


    script:

    // Taxonomic abundance configuration for report (default)

    """
    ganon classify --db-prefix $ganonDatabase/$ganonDatabasePrefix --single-reads $reads --output-prefix $sampleID --threads $task.cpus --multiple-matches $ganonMultipleMatches

    mv ${sampleID}.tre ${sampleID}.ganon.abundance.report
    """

}

process MetaSpades {

    tag { sampleID }
    label "pathogenAssemblyMetaspades"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metaspades.fasta"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    val(kmerList)
    val(minContigLength)
    val(args)

    output:
    tuple(val(sampleID), val("metaspades"), path("${sampleID}.metaspades.fasta"), emit: contigs)
    tuple(val(sampleID), path(forward), path(reverse), emit: reads)

    script:

    """
    spades.py --meta -t $task.cpus -1 $forward -2 $reverse -o assembly/ -k $kmerList $args
    mv assembly/contigs.fasta ${sampleID}.metaspades.fasta
    """

}


process MetaSpadesNanopore {

    tag { sampleID }
    label "pathogenAssemblyMetaspades"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metaspades.fasta"

    input:
    tuple val(sampleID), path(reads)
    val(kmerList)
    val(minContigLength)
    val(args)

    output:
    tuple(val(sampleID), val("metaspades"), path("${sampleID}.metaspades.fasta"), emit: contigs)
    tuple(val(sampleID), path(reads), emit: reads)

    script:

    """
    spades.py --meta -t $task.cpus -s $reads -o assembly/ -k $kmerList $args
    mv assembly/contigs.fasta ${sampleID}.metaspades.fasta
    """

}

process Megahit {

    tag { sampleID }
    label "pathogenAssemblyMegahit"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.megahit.fasta"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    val(kmerList)
    val(minContigLength)
    val(args)

    output:
    tuple(val(sampleID), val("megahit"), path("${sampleID}.megahit.fasta"), emit: contigs)
    tuple(val(sampleID), path(forward), path(reverse), emit: reads)

    script:

    """
    megahit -t $task.cpus -1 $forward -2 $reverse -o assembly --k-list $kmerList --min-contig-len $minContigLength $args
    mv assembly/final.contigs.fa ${sampleID}.megahit.fasta
    """
}


process MegahitNanopore {

    tag { sampleID }
    label "pathogenAssemblyMegahit"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.megahit.fasta"

    input:
    tuple val(sampleID), path(reads)
    val(kmerList)
    val(minContigLength)
    val(args)

    output:
    tuple(val(sampleID), val("megahit"), path("${sampleID}.megahit.fasta"), emit: contigs)
    tuple(val(sampleID), path(reads), emit: reads)

    script:

    """
    megahit -t $task.cpus -r $reads -o assembly --k-list $kmerList --min-contig-len $minContigLength $args
    mv assembly/final.contigs.fa ${sampleID}.megahit.fasta
    """
    
}


process MetaMdbgNanopore {

    tag { sampleID }
    label "pathogenAssemblyMetaMdbg"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metamdbg.fasta"

    input:
    tuple val(sampleID), path(reads)
    val(minReadOverlap)

    output:
    tuple(val(sampleID), val("metamdbg"), path("${sampleID}.metamdbg.fasta"), emit: contigs)
    tuple(val(sampleID), path(reads), emit: reads)

    script:

    """
    metaMDBG --out-dir assembly/ --in-ont $reads --threads $task.cpus
    mv assembly/contigs.fasta ${sampleID}.metamdbg.fasta
    """

}

process ContigCoverage {
    
    tag { sampleID }
    label "pathogenAssemblyContigCoverage"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "symlink", pattern: "${sampleID}.${assembler}.bam*"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    tuple val(sampleID2), path(forward), path(reverse)

    output:
    tuple(val(sampleID), val(assembler), path(contigs), emit: contigs)
    tuple(val(sampleID), path("${sampleID}.${assembler}.bam"), path("${sampleID}.${assembler}.bam.bai"), emit: coverage)

    script:

    """
    minimap2 -ax sr -t $task.cpus $contigs $forward $reverse | samtools view -@ $task.cpus -hbF 4 - | samtools sort -@ $task.cpus - > ${sampleID}.${assembler}.bam
    samtools index ${sampleID}.${assembler}.bam
    """
}


process ContigCoverageNanopore {
    
    tag { sampleID }
    label "pathogenAssemblyContigCoverage"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    tuple val(sampleID2), path(reads)

    output:
    tuple(val(sampleID), val(assembler), path(contigs), emit: contigs)
    tuple(val(sampleID), path("${sampleID}.contigs.bam"), path("${sampleID}.contigs.bam.bai"), emit: coverage)

    script:

    """
    minimap2 -ax map-ont -t $task.cpus $contigs $reads | samtools view -@ $task.cpus -hbF 4 - | samtools sort -@ $task.cpus - > ${sampleID}.contigs.bam
    samtools index ${sampleID}.contigs.bam
    """
}


process Concoct {
    
    tag { sampleID }
    label "pathogenAssemblyConcoct"

    errorStrategy 'ignore'  // empty files e.g. BLANK

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "concoct_bins"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    tuple val(sampleID2), path(contigBam), path(contigBai)
    val(chunkSize)
    val(readLength)
    val(minBinSize)
    val(minContigLength)

    output:
    tuple(val(sampleID), path("concoct_bins"), emit: bins)

    script:

    """
    cut_up_fasta.py $contigs -c $chunkSize -o 0 --merge_last -b contigs_${chunkSize}.bed > contigs_${chunkSize}.fasta
    concoct_coverage_table.py contigs_${chunkSize}.bed $contigBam > coverage_table.tsv
    concoct --read_length $readLength --length_threshold $minContigLength --threads $task.cpus --composition_file contigs_${chunkSize}.fasta --coverage_file coverage_table.tsv -b concoct/
    merge_cutup_clustering.py concoct/clustering_gt${minContigLength}.csv > concoct/clustering_merged.csv

    mkdir -p concoct_bins
    extract_fasta_bins.py $contigs concoct/clustering_merged.csv --output_path concoct_bins
    """
}

process Metabat2 {
    
    tag { sampleID }
    label "pathogenAssemblyMetabat2"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "metabat_bins"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    tuple val(sampleID2), path(contigBam), path(contigBai)
    val(minBinSize)
    val(minContigLength)


    script:

    """
    runMetaBat.sh --numThreads $task.cpus --minContig $minContigLength --minClsSize $minBinSize $contigs $contigBam
    """

}

process SemiBin2 {
    
    tag { sampleID }
    label "pathogenAssemblySemiBin2"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "symlink", pattern: "semibin2"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    tuple val(sampleID2), path(contigBam), path(contigBai)
    val(minContigLength)

    script:

    """
    SemiBin2 single_easy_bin -i $contigs -b $contigBam -o semibin2 --environment global --threads $task.cpus --min-len $minContigLength
    """

}

process Vamb {
    
    tag { sampleID }
    label "pathogenAssemblySemiBin2"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "symlink", pattern: "semibin2"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    tuple val(sampleID2), path(contigBam), path(contigBai)

    script:

    """
    SemiBin2 single_easy_bin -i $contigs -b $contigBam -o semibin2 --environment global --threads $task.cpus
    """

}

process Vircov {
    
    label "pathogenProfileVircov"
    tag { sampleID }

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.vircov.tsv"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), path(reference, name: 'vircov__reference')  // index and reference can be the same
    val(aligner)
    val(secondary)
    val(remapThreads)
    val(remapParallel)
    val(vircovArgs)

    output:
    tuple (val(sampleID), path(forward), path(reverse), emit: reads)
    tuple (val(sampleID), path("${sampleID}.vircov.tsv"), emit: results)

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    secondaryFlag = secondary ? "--secondary" : "" 

    """
    vircov run -i $forward -i $reverse -o ${sampleID}.vircov.tsv --aligner $aligner --index $alignmentIndex --reference vircov__reference --workdir data/ \
    --scan-threads $task.cpus --remap-threads $remapThreads --remap-parallel $remapParallel $secondaryFlag $vircovArgs
    """
    
}

process VircovNanopore {
    
    label "pathogenProfileVircov"
    tag { sampleID }

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.alignment.tsv"

    input:
    tuple val(sampleID), path(reads)
    tuple path(index), path(reference, name: 'vircov__reference')  // index and reference can be the same
    val(aligner)
    val(secondary)
    val(remapThreads)
    val(remapParallel)

    output:
    tuple (val(sampleID), path(reads), emit: reads)
    tuple (val(sampleID), path("${sampleID}.alignment.tsv"), emit: results)

    script:

    alignmentIndex = index[0]

    secondaryFlag = secondary ? "--secondary" : "" 

    """
    vircov run -i $reads -o ${sampleID}.alignment.tsv --aligner minimap2 --preset map-ont --index $alignmentIndex --reference vircov__reference --workdir data/ \
    --scan-threads $task.cpus --remap-threads $remapThreads --remap-parallel $remapParallel $secondaryFlag
    """
    
}

process BlastContigs {

    label "pathogenAssemblyBlast"
    tag { sampleID }

    publishDir  "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.blast.assembly.tsv"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    path(database)
    val(databasePrefix)
    val(minPercentIdentity)
    val(minEvalue)
    val(maxTargetSeqs)


    output:
    tuple(val(sampleID), path("${sampleID}.blast.assembly.tsv"), emit: results)

    script:
    
    // Set the execution environment variable BLASTDB to database path to enable the taxonomic assignments
    
    """
    BLASTDB=$database blastn -num_threads $task.cpus -query $contigs -perc_identity $minPercentIdentity -evalue $minEvalue -max_target_seqs $maxTargetSeqs -db ${database}/${databasePrefix} -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length nident pident evalue bitscore staxid ssciname stitle' > ${sampleID}.blast.assembly.tsv
    """

}

process BlastContigsBitscoreStream {

    label "pathogenAssemblyBlast"
    tag { sampleID }

    publishDir  "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.blast.assembly.tsv"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    path(database)
    val(databasePrefix)
    val(minPercentIdentity)
    val(minEvalue)
    val(maxTargetSeqs)


    output:
    tuple(val(sampleID), path("${sampleID}.blast.bitscore.tsv"), emit: results)

    script:

    """
    BLASTDB="$database" blastn \
        -num_threads "$task.cpus" \
        -query "$contigs" \
        -perc_identity "$minPercentIdentity" \
        -evalue "$minEvalue" \
        -max_hsps 1 \
        -db "${database}/${databasePrefix}" \
        -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length nident pident evalue bitscore staxid ssciname stitle' \
        | awk 'BEGIN{FS=OFS="\t"} {
            k=\$1                          # qseqid
            bs=\$13                        # bitscore
            if(!(k in max) || bs>max[k]) {max[k]=bs; line[k]=\$0}
        }
        END{for(k in line) print line[k]}' \
    > "${sampleID}.blast.assembly.tsv"
    """
}

process ProcessOutputIllumina {
    
    tag { sampleID }
    label "cerebroPipeline"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.pd.json"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.qc.json"
    
    publishDir "$params.outputDirectory/results/samples/$sampleID", mode: "copy", pattern: "*"

    input:
    tuple val(sampleID), path(result_files)

    output:
    tuple (path("${sampleID}.qc.json"), path("${sampleID}.pd.json"), emit: results)
    val(sampleID), emit: samples

    script:

    """
    cerebro-pipeline process pathogen --id ${sampleID} --qc ${sampleID}.qc.json --pathogen ${sampleID}.pd.json --paired-end --qc-fail-ok
    """  
}

process ProcessOutputNanopore {
    
    tag { sampleID }
    label "cerebro"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.pd.json"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.qc.json"

    publishDir "$params.outputDirectory/results/samples/$sampleID", mode: "copy", pattern: "${sampleID}.pd.json"
    publishDir "$params.outputDirectory/results/samples/$sampleID", mode: "copy", pattern: "${sampleID}.qc.json"

    input:
    tuple val(sampleID), path(result_files)

    output:
    tuple (path("${sampleID}.qc.json"), path("${sampleID}.pd.json"), emit: results)

    script:

    """
    cerebro-pipeline process pathogen --id ${sampleID} --qc ${sampleID}.qc.json --pathogen ${sampleID}.pd.json --qc-fail-ok
    """
}

process PathogenDetectionTable {
    
    label "cerebro"

    publishDir "$params.outputDirectory/results", mode: "copy", pattern: "species.tsv"

    input:
    path(result_files)
    path(taxonomy_directory)

    output:
    path("species.tsv")

    script:

    """
    echo '{"ranks": ["Species"]}' > filter.json
    cerebro-pipeline table pathogen-detection --json *.pd.json --output species.tsv --taxonomy $taxonomy_directory --filter-json filter.json
    """
}
