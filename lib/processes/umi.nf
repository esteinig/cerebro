// Calib was designed for sparse UMI space utilisation: https://github.com/vpc-ccg/calib/issues/47 - we therefore set edit-distance to 0

process CalibDeduplication {

    label "umi_calib"
    tag { "$id" }

    input:
    tuple val(id), path(forward), path(reverse)

    output:
    tuple (val(id), path("${id}_calib_dedup_1.fq.gz"), path("${id}_calib_dedup_2.fq.gz"), emit: reads)

    // Parameters must be selected for specific read lengths (150 bp reads) and UMI lengths (fixed to 12 for protocol here). Input must be gzipped.
    
    """
    cerebro tools umi prepare-calib --input $forward --output ${id}_calib_1.fq.gz 
    calib -f ${id}_calib_1.fq.gz -r $reverse -o ${id}. -c $task.cpus -l1 12 -l2 0 -e 0 -t 2 -k 8 -m 7 --gzip-input 
    cerebro tools umi dedup-calib --input $forward $reverse --output ${id}_calib_dedup_1.fq.gz ${id}_calib_dedup_2.fq.gz --clusters ${id}.cluster
    rm ${id}_calib_1.fq.gz ${id}.cluster
    """
    
}

// Uses as much memory as the size of the uncompressed read input due to non-hashed 
// deduplication to avoid potential hash collisions and ensure deterministic
// outputs (necessary for accreditation) - note that Fastp deduplication is not
// deterministic and causes hash collisions in 0.01% of reads (see repository)
// which makes the deduplication process unusable for our workflow!

process NaiveDeduplication {

    label "umi_cerebro"
    tag { "$id" }

    input:
    tuple val(id), path(forward), path(reverse)
    val(umi)

    output:
    tuple (val(id), path("${id}_dedup_1.fq.gz"), path("${id}_dedup_2.fq.gz"), emit: reads)
    
    script:

    no_umi_opt = umi ? "" : "--no-umi"

    """
    cerebro tools umi dedup-naive --input $forward $reverse --reads $forward --output ${id}_dedup_1.fq.gz ${id}_dedup_2.fq.gz $no_umi_opt
    """
   
}