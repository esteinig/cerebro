/* Ensure params.referenceFasta and params.normalControlBam exists then execute with:
    
    nextflow run main.nf -p mamba --outdir test --pairedReads "*_{R1,R2}.fastq.gz"

*/ 

process MinimapAneuploidy {

    tag { "$id : $indexName" }
    label "minimap2"

    publishDir "$params.outdir/${id}", mode: "symlink", pattern: "${id}.bam"
    publishDir "$params.outdir/${id}", mode: "copy", pattern: "${id}_coverage.txt"
    publishDir "$params.outdir/${id}", mode: "copy", pattern: "${id}_markdup.json"

    input:
    tuple val(id), path(forward), path(reverse)
    path(index)

    output:
    tuple (val(id), path(forward), path(reverse), emit: reads)
    tuple (val(id), val(indexName), path("${id}.bam"), emit: alignment)
    tuple (val(id), val(indexName), path("${id}_coverage.txt"), emit: coverage)
    tuple (val(id), val(indexName), path("${id}_markdup.json"), optional: true, emit: deduplicate)

    script:

    indexName = index.getSimpleName()

    if (params.deduplicate){
        
        // Optical duplicate setting `-d 100` for HiSeq style platforms, use `-d 2500` for NovaSeq style platforms

        """
        minimap2 -t $task.cpus --sam-hit-only -ax sr $index $forward $reverse | samtools view -Su - | samtools collate -@ $task.cpus -O -u - | samtools fixmate -@ $task.cpus -m -u - -  | samtools sort -@ $task.cpus -u - | samtools markdup -@ $task.cpus -f ${id}_markdup.json --json -S -d $params.markdupDistance --mode s --include-fails -r - ${id}.bam
        samtools coverage ${id}.bam --no-header > ${id}_coverage.txt
        """
    } else {
        """
        minimap2 -t $task.cpus --sam-hit-only -ax sr $index $forward $reverse | samtools view -Su - | samtools sort -@ $task.cpus - -o ${id}.bam
        samtools coverage ${id}.bam- --no-header > ${id}_coverage.txt
        """
    }
    
}

process CnvKitAneuploidy {

    tag { "$id : $indexName" }
    label "cnvkit"

    publishDir "$params.outdir/${id}", mode: "copy", pattern: "${id}.png"
    publishDir "$params.outdir/${id}", mode: "copy", pattern: "${id}.cns"
    publishDir "$params.outdir/${id}", mode: "copy", pattern: "${id}.cnr"
    publishDir "$params.outdir/${id}", mode: "symlink", pattern: "${id}.cnn"

    input:
    tuple val(id), val(indexName), path(bam)
    path(controlBam)
    path(fasta)

    output:
    path("${id}.png")
    path("${id}.cns")
    path("${id}.cnr")
    path("${id}.cnn")
    
    script:

    dropLowCoverage = params.dropLowCoverage ? "--drop-low-coverage" : ""
    targetAverageBinSize = params.targetAverageBinSize ? "--target-avg-size $params.targetAverageBinSize" : ""

    """
    cnvkit.py batch $bam -n $controlBam -m wgs -f $fasta -p $task.cpus --output-reference ${id}.cnn --output-dir results/ --diagram --scatter $dropLowCoverage $targetAverageBinSize
    
    mv results/${id}-scatter.png ${id}.png
    mv results/${id}.cns .
    mv results/${id}.cnr .
    """

}

workflow {

        pairedReads = channel.fromFilePairs(
            params.pairedReads, 
            flat: true, 
            checkIfExists: true
        )

        referenceFasta = channel.fromPath(params.referenceFasta, checkIfExists: true).first()
        normalControlBam = channel.fromPath(params.normalControlBam, checkIfExists: true).first()

        pairedReads | view

        MinimapAneuploidy(
            pairedReads, 
            referenceFasta
        )

        CnvKitAneuploidy(
            MinimapAneuploidy.out.alignment, 
            normalControlBam, 
            referenceFasta, 
        )

}

