process MashScreenWinner {

    label "mash"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/workflow/$params.subdir/${id}", mode: "copy", pattern: "${idx_name}.txt"

    input:
    tuple val(id), path(forward), path(reverse)
    path(index)

    output:
    tuple val(id), val(idx_name), path(forward), path(reverse), path("${idx_name}.txt")

    script:

    idx_name = index.baseName

    """
    mash screen -w -p $task.cpus $index $forward $reverse > ${idx_name}.txt
    """
    

}
