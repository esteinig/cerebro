
params.outdir                               = "aneuploidy"
params.pairedReads                          = "*_{R1,R2}.fastq.gz"
params.deduplicate                          = true
params.host.aneuploidy.deduplicate          = true                           // samtools markdup removal on host alignment bam
params.markdupDistance                      = 100                            // http://www.htslib.org/algorithms/duplicate.html
params.referenceFasta                       = "db/chm13v2.fa"                
params.normalControl                        = "db/HG007.CHM13v2.5x.bam"      // control alignment reference must be referenceFasta

params.dropLowCoverage  = false
params.targetAverageBinSize = null

process MinimapAneuploidy {

    tag { "$id : $idx_name" }
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
    tuple (val(id), val(indexName), path("coverage.txt"), emit: coverage)
    tuple (val(id), val(indexName), path("markdup.json"), optional: true, emit: deduplicate)

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

    tag { "$id : $idx_name" }
    label "cnvkit"

    publishDir "$params.outdir/${id}", mode: "copy", pattern: "${id}-scatter.png"
    publishDir "$params.outdir/${id}", mode: "copy", pattern: "${id}.cns"
    publishDir "$params.outdir/${id}", mode: "copy", pattern: "${id}.cnr"
    publishDir "$params.outdir/${id}", mode: "symlink", pattern: "${id}_${indexName}.cnn"

    input:
    tuple val(id), val(indexName), path(bam)
    path(controlBam)
    path(fasta)

    output:
    path("${id}-scatter.png")
    path("${id}.cns")
    path("${id}.cnr")
    path("${id}_${indexName}.cnn")
    
    script:

    dropLowCoverage = params.dropLowCoverage ? "--drop-low-coverage" : ""
    targetAverageBinSize = params.targetAverageBinSize ? "--target-avg-size $params.targetAverageBinSize" : ""

    """
    cnvkit.py batch $bam -n $controlBam -m wgs -f $fasta -p $task.cpus --output-reference ${id}_${indexName}.cnn --output-dir results/ --diagram --scatter $dropLowCoverage $targetAverageBinSize
    
    mv results/${id}-scatter.png .
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

        MinimapAneuploidy(
            pairedReads, 
            params.referenceFasta
        )

        CnvKitAneuploidy(
            MinimapAneuploidy.out.alignment, 
            params.normalControl, 
            params.referenceFasta
        )

}

