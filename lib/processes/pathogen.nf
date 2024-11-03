


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

process Kraken2Nanopore {

    tag { sampleID }
    label "pathogenProfileKraken2"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kraken2.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kraken2.report"

    input:
    tuple val(sampleID), path(reads)
    path(krakenDatabase)
    val(krakenConfidence)

    output:

    tuple (val(sampleID), path(krakenDatabase), path("${sampleID}.kraken2.report"), emit: bracken)
    tuple (val(sampleID), path("${sampleID}.kraken2.report"), path("${sampleID}.kraken2.tsv"), emit: results)

    script:

    """
    kraken2 --db $krakenDatabase --confidence $krakenConfidence --threads $task.cpus --output ${sampleID}.kraken2.tsv --report ${sampleID}.kraken2.report $reads
    
    if [ ! -f "${sampleID}.kraken2.tsv" ]; then
        touch "${sampleID}.kraken2.tsv"
    fi
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
    val(brackenMinReads)

    output:
    tuple (val(sampleID), path("${sampleID}.bracken.report"), emit: results)

    script:

    if (krakenReport.size() == 0) {
        """
        touch ${sampleID}.bracken.report
        """
    } else {
        """
        bracken -d $krakenDatabase -i $krakenReport -r $brackenReadLength -l $brackenRank -t $brackenMinReads -o ${sampleID}.bracken.report
        """
    }
}


process Sylph {

    tag { sampleID }
    label "pathogenProfileSylph"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.sylph.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.sylph.report"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(sylphDatabase)
    path(sylphMetadata)
    path(sylphMinNumberKmers)
    path(sylphQueryCompression)

    output:
    tuple (val(sampleID), path("${sampleID}.sylph.tsv"), path("${sampleID}.sylph.report"), emit: results)

    script:

    """
    sylph profile $sylphDatabase -1 $forward -2 $reverse -c $sylphQueryCompression --min-number-kmers $sylphMinNumberKmers -t $task.cpus > ${sampleID}.sylph.tsv
    python $baseDir/lib/scripts/sylph_to_taxprof.py -m $sylphMetadata -s ${sampleID}.sylph.tsv -o "" 
    
    if [ -f "${reads}.sylphmpa" ]; then
        mv ${reads}.sylphmpa ${sampleID}.sylph.report
    else
        touch "${sampleID}.sylph.report"
    fi

    """

}


process SylphNanopore {

    tag { sampleID }
    label "pathogenProfileSylph"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.sylph.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.sylph.report"

    input:
    tuple val(sampleID), path(reads)
    path(sylphDatabase)
    path(sylphMetadata)
    path(sylphMinNumberKmers)
    path(sylphQueryCompression)

    output:
    tuple (val(sampleID), path("${sampleID}.sylph.tsv"), path("${sampleID}.sylph.report"), emit: results)

    script:

    """
    sylph profile $sylphDatabase $reads -c $sylphQueryCompression --min-number-kmers $sylphMinNumberKmers -t $task.cpus > ${sampleID}.sylph.tsv
    python $baseDir/lib/scripts/sylph_to_taxprof.py -m $sylphMetadata -s ${sampleID}.sylph.tsv -o "" 
    
    if [ -f "${reads}.sylphmpa" ]; then
        mv ${reads}.sylphmpa ${sampleID}.sylph.report
    else
        touch "${sampleID}.sylph.report"
    fi

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
    tuple (val(sampleID), path("${sampleID}.metabuli.report"), path("${sampleID}.metabuli.tsv"), emit: results)

    script:

    memoryLimit = "${task.memory}".split()[0]

    """
    metabuli classify --max-ram $memoryLimit --threads $task.cpus $forward $reverse $metabuliDatabase classified/ $sampleID

    cp classified/${sampleID}_classifications.tsv ${sampleID}.metabuli.tsv
    rm classified/${sampleID}_classifications.tsv
    cp classified/${sampleID}_report.tsv ${sampleID}.metabuli.report
    """

}


process MetabuliNanopore {

    tag { sampleID }
    label "pathogenProfileMetabuli"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metabuli.tsv"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.metabuli.report"

    input:
    tuple val(sampleID), path(reads)
    path(metabuliDatabase)

    output:
    tuple (val(sampleID), path("${sampleID}.metabuli.report"), path("${sampleID}.metabuli.tsv"), emit: results)

    script:

    memoryLimit = "${task.memory}".split()[0]

    """
    metabuli classify --max-ram $memoryLimit --threads $task.cpus --seq-mode 3 $reads $metabuliDatabase classified/ $sampleID

    cp classified/${sampleID}_classifications.tsv ${sampleID}.metabuli.tsv
    rm classified/${sampleID}_classifications.tsv
    cp classified/${sampleID}_report.tsv ${sampleID}.metabuli.report
    """

}



process Kmcp { 

    tag { sampleID }
    label "pathogenProfileKmcp"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kmcp.abundance.report"
    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.kmcp.reads.report"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(kmcpDatabase)
    val(kmcpMode)
    val(kmcpLevel)

    output:
    tuple (val(sampleID), path("${sampleID}.reads.kmcp.report"), path("${sampleID}.abundance.kmcp.report"), emit: results)


    script:

    """
    kmcp search -d $kmcpDatabase -1 $forward -2 $reverse -o ${sampleID}.reads.tsv.gz --threads $task.cpus
    kmcp profile --level $kmcpLevel --taxid-map $kmcpDatabase/taxid.map --taxdump $kmcpDatabase/taxonomy --mode $kmcpMode -B ${sampleID} -C ${sampleID}.abundance.kmcp.report ${sampleID}.reads.tsv.gz
    mv ${sampleID}.binning.gz ${sampleID}.reads.kmcp.report
    """

}

process GanonReads { 

    tag { sampleID }
    label "pathogenProfileGanon"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.reads.ganon.report"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(ganonDatabase)
    val(ganonMultipleMatches)

    output:
    tuple (val(sampleID), path("${sampleID}.reads.ganon.report"), emit: results)


    script:

    // Sequence abundance configuration for report (--binning

    """
    ganon classify --db-prefix --paired-reads $forward $reverse --output-prefix $sampleID --threads $task.cpus --binning --multiple-matches $ganonMultipleMatches
    mv ${sampleID}.rep ${sampleID}.reads.ganon.report
    """

}

process GanonProfile { 

    tag { sampleID }
    label "pathogenProfileGanon"

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.rep", saveAs: "${sampleID}.abundance.ganon.report"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    path(ganonDatabase)
    val(ganonMultipleMatches)

    output:
    tuple (val(sampleID), path("${sampleID}.abundance.ganon.report"), emit: results)


    script:

    // Taxonomic abundance configuration for report (default)

    """
    ganon classify --db-prefix --paired-reads $forward $reverse --output-prefix $sampleID --threads $task.cpus --multiple-matches $ganonMultipleMatches
    mv ${sampleID}.rep ${sampleID}.abundance.ganon.report
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
    tuple(val(sampleID), path(forward), path(reverse), emit: reads)

    script:

    """
    metaMDBG --out-dir assembly/ --in-ont $reads --threads $task.cpus
    mv assembly/contigs.fasta ${sampleID}.metamdbg.fasta
    """

}

process ContigCoverage {
    
    tag { sampleID }
    label "pathogenAssemblyContigCoverage"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    tuple val(sampleID), path(forward), path(reverse)

    output:
    tuple(val(sampleID), val(assembler), path(contigs), emit: contigs)
    tuple(val(sampleID), path("${sampleID}.contigs.bam"), path("${sampleID}.contigs.bam.bai"), emit: coverage)

    """
    minimap2 -ax sr -t $task.cpus $contigs $forward $reverse | samtools view -@ $task.cpus -hbF 12 - | samtools sort -@ $task.cpus - > ${sampleID}.contigs.bam
    samtools index ${sampleID}.contigs.bam
    """
}


process ContigCoverageNanopore {
    
    tag { sampleID }
    label "pathogenAssemblyContigCoverage"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    tuple val(sampleID), path(reads)

    output:
    tuple(val(sampleID), val(assembler), path(contigs), emit: contigs)
    tuple(val(sampleID), path("${sampleID}.contigs.bam"), path("${sampleID}.contigs.bam.bai"), emit: coverage)

    """
    minimap2 -ax map-ont -t $task.cpus $contigs $reads| samtools view -@ $task.cpus -hbF 12 - | samtools sort -@ $task.cpus - > ${sampleID}.contigs.bam
    samtools index ${sampleID}.contigs.bam
    """
}


process Concoct {
    
    tag { sampleID }
    label "pathogenAssemblyConcoct"


    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "concoct_bins"

    input:
    tuple val(sampleID), val(assembler), path(contigs)
    tuple val(sampleID), path(contigBam), path(contigBai)
    val(chunkSize)
    val(readLength)
    val(minBinSize)
    val(minContigLength)

    output:
    tuple(val(sampleID), path("concoct_bins"), emit: bins)

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
    tuple val(sampleID), path(contigBam), path(contigBai)
    val(minBinSize)
    val(minContigLength)

    output:


    """
    runMetaBat.sh --numThreads $task.cpus --minContig $minContigLength --minClsSize $minBinSize $contigs $contigBam
    """

}

process Vircov {
    
    label "pathogenProfileVircov"
    tag { sampleID }

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.alignment.tsv"

    input:
    tuple val(sampleID), path(forward), path(reverse)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
    val(aligner)
    val(secondary)
    val(remapThreads)
    val(remapParallel)

    output:
    tuple (val(sampleID), path(forward), path(reverse), emit: reads)
    tuple (val(sampleID), path("${sampleID}.alignment.tsv"), emit: results)

    script:

    indexName = index[0].getSimpleName()
    alignmentIndex = aligner == "bowtie2" ? indexName : index[0]

    secondaryFlag = secondary ? "--secondary" : "" 

    """
    vircov run -i $forward -i $reverse -o ${sampleID}.alignment.tsv --aligner $aligner --index $alignmentIndex --reference vircov__reference --workdir data/ \
    --scan-threads $task.cpus --remap-threads $remapThreads --parallel $remapParallel $secondaryFlag
    """
    
}


process VircovNanopore {
    
    label "pathogenProfileVircov"
    tag { sampleID }

    publishDir "$params.outputDirectory/pathogen/$sampleID", mode: "copy", pattern: "${sampleID}.alignment.tsv"

    input:
    tuple val(sampleID), path(reads)
    tuple path(index), (path(reference) , stageAs: 'vircov__reference')  // index and reference can be the same
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
    vircov run -i $forward -i $reverse -o ${sampleID}.alignment.tsv --aligner minimap2 --preset map-ont --index $alignmentIndex --reference vircov__reference --workdir data/ \
    --scan-threads $task.cpus --remap-threads $remapThreads --parallel $remapParallel $secondaryFlag
    """
    
}