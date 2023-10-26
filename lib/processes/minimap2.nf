process MinimapAlignPAF {

    tag { "$id : $idx_name" }
    label "minimap2"

    publishDir "$params.outdir/workflow/$params.subdir", mode: "symlink", pattern: "${id}.paf"

    input:
    tuple val(id), path(forward), path(reverse)
    path(index)

    output:
    tuple val(id), path(forward), path(reverse)
    tuple val(id), val(idx_name), path("${id}.paf")

    script:

    idx_name = index.baseName

    """
    minimap2 -t $task.cpus -c -x sr ${index} $forward $reverse > ${id}.paf
    """

}

process MinimapRealignBAM {

    tag { "$id : $idx_name" }
    label "minimap2_realign"

    publishDir "$params.outdir/workflow/$params.subdir/", mode: "symlink", pattern: "${id}_${idx_name}.bam"
    publishDir "$params.outdir/workflow/$params.subdir/coverage", mode: "copy", pattern: "${id}_${idx_name}.bed"
    publishDir "$params.outdir/workflow/$params.subdir/coverage", mode: "copy", pattern: "${id}_${idx_name}.txt"

    input:
    tuple val(id), path(forward), path(reverse), val(db_name), path(fasta)

    output:
    tuple(val(id), path(forward), path(reverse), emit: reads)
    tuple(val(id), val(db_name), val(idx_name), path(fasta), path("${id}_${idx_name}.bam"), emit: aligned)
    tuple(val(id), val(db_name), path("${id}_${idx_name}.txt"), emit: coverage)
    path("${id}_${idx_name}.bed")
    
    script:

    idx_name = fasta.baseName

    """
    minimap2 -t $task.cpus --sam-hit-only -ax sr ${fasta} $forward $reverse | samtools view -@ $task.cpus -Sb - | samtools sort -@ $task.cpus - -o ${id}_${idx_name}.bam
    covtobed --max-cov $params.covtobed_max_cov ${id}_${idx_name}.bam > ${id}_${idx_name}.bed
    samtools coverage ${id}_${idx_name}.bam --no-header > coverage.txt
    grep "^>" $fasta | cut -c2- | cut -d' ' -f2- > descr.txt
    paste coverage.txt descr.txt > ${id}_${idx_name}.txt
    """

}


process MinimapIndexSubset {

    tag { "$idx_name" }
    label "minimap2"

    publishDir "$params.outdir/workflow/$params.subdir/subsets/${id}", mode: "symlink", pattern: "${idx_name}.mmi"

    input:
    tuple val(id), val(idx_name), path(forward), path(reverse), path(fasta)

    output:
    tuple val(id), path(forward), path(reverse)
    tuple val(id), val(idx_name), path("${idx_name}.mmi"), path(fasta)

    script:

    """
    minimap2 -t $task.cpus -x sr -d ${idx_name}.mmi $fasta
    """

}

process MinimapAlignSubsetPAF {

    tag { "$id : $idx_name" }
    label "minimap2"

    publishDir "$params.outdir/workflow/$params.subdir/subsets/${id}/", mode: "symlink", pattern: "${idx_name}.paf"

    input:
    tuple val(id), path(forward), path(reverse)
    tuple val(id), val(idx_name), path(index), path(fasta)

    output:
    tuple val(id), path(forward), path(reverse)
    tuple val(id), val(idx_name), path("${idx_name}.paf"), path(fasta)

    script:

    """
    minimap2 -t $task.cpus -c -x sr ${index} $forward $reverse > ${idx_name}.paf
    """

}


process MinimapAneuploidy {

    tag { "$id : $idx_name" }
    label "minimap2_aneuploidy"

    publishDir "$params.outdir/workflow/$params.subdir/${id}", mode: "symlink", pattern: "${id}.bam"
    publishDir "$params.outdir/workflow/$params.subdir/${id}", mode: "copy", pattern: "coverage.txt"
    publishDir "$params.outdir/workflow/$params.subdir/${id}", mode: "copy", pattern: "markdup.json"

    input:
    tuple val(id), path(forward), path(reverse)
    path(index)

    output:
    tuple (val(id), path(forward), path(reverse), emit: reads)
    tuple (val(id), val(idx_name), path("${id}.bam"), emit: alignment)
    tuple (val(id), val(idx_name), path("coverage.txt"), emit: coverage)
    tuple (val(id), val(idx_name), path("markdup.json"), optional: true, emit: deduplicate)

    script:

    idx_name = index.getSimpleName()

    if (params.host.aneuploidy.deduplicate){
        
        // Optical duplicate setting `-d 100` for HiSeq style platforms, use `-d 2500` for NovaSeq style platforms

        """
        minimap2 -t $task.cpus --sam-hit-only -ax sr $index $forward $reverse | samtools view -Su - | \
            samtools collate -@ $task.cpus -O -u - | samtools fixmate -@ $task.cpus -m -u - -  | samtools sort -@ $task.cpus -u - | \
            samtools markdup -@ $task.cpus -f markdup.json --json -S -d $params.host.aneuploidy.markdup_distance --mode s --include-fails -r - ${id}.bam
        samtools coverage ${id}.bam --no-header > coverage.txt
        """
    } else {
        """
        minimap2 -t $task.cpus --sam-hit-only -ax sr $index $forward $reverse | samtools view -Su - | samtools sort -@ $task.cpus - -o ${id}.bam
        samtools coverage ${id}.bam- --no-header > coverage.txt
        """
    }
    

}