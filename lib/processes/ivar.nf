process IvarConsensus {

    label "ivar"
    tag { "$id : $idx_name" }

    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}.variants.tsv"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}_${idx_name}.consensus.fasta"

    input:
    tuple val(id), val(db_name), val(idx_name), path(reference), path(bam), path(vircov), path(fwd_aligned), path(rev_aligned)
    
    output:
    tuple (val(id), val(db_name), path("${id}_${idx_name}.consensus.fasta"), emit: consensus) optional true
    path("${id}_${idx_name}.variants.tsv") optional true

    script:

    // If no reads were extracted, do not continue with consensus assembly.

    // For segmented genomes with multiple sequences in the reference file, we must
    // split the sequences into their own files first, realign the extracted reads to 
    // each reference sequence and construct the consensus for each before merging again
    // - if there is only a single sequence, the re-mapping step is ommitted and the 
    // input alignment used directly

    // This is all somewhat inefficient, but at least automates this process without 
    // splitting the segments in previous steps. Note that the VirecovRemapping step
    // due to filtering of the alignment can leave unpaired reads, which are excluded
    // in the realignment loop for the consensus sequences.

    // Note: > 2 as the top level directory (.) is counted by `find`

    // At the end we replace the otherwise not very informative header (by default constructed 
    // from the input sequence file) in the consensus file with the original reference alignment
    // header so that identifier and annotation are restored (e.g. to trace segments)

    if (fwd_aligned.size() > 0 && rev_aligned.size() > 0) {

        """
        cerebro tools split --input $reference --outdir sequences/ 
        count=\$(find sequences/ | wc -l)
        
        if [ \$count -gt 2 ]
        then
            for seq in sequences/*; do
                seq_name=\$(basename \${seq%%.*})

                minimap2 -t $task.cpus --sam-hit-only -ax sr \${seq} $fwd_aligned $rev_aligned | samtools view -@ $task.cpus -Sb - | samtools sort -@ $task.cpus - -o tmp.bam
                                
                samtools mpileup $params.ivar_mpileup_args -A -B -Q 0 tmp.bam | ivar consensus -p \${seq_name}.consensus \
                    -q $params.ivar_min_qual \
                    -t $params.ivar_min_freq \
                    -m $params.ivar_min_depth \
                    -n N
                
                samtools mpileup $params.ivar_mpileup_args -A -B -Q 0 tmp.bam | ivar variants -p \${seq_name}.variants \
                    -q $params.ivar_min_qual \
                    -t $params.ivar_min_freq \
                    -m $params.ivar_min_depth \
                    -r \$seq 

                rm tmp.bam

                HEADER=\$(head -1 \$seq)
                sed -i "1s~.*~\$HEADER~" \${seq_name}.consensus.fa

            done
            cat *.consensus.fa > ${id}_${idx_name}.consensus.fasta
            cat *.variants.tsv > ${id}_${idx_name}.variants.tsv
        else 
            samtools mpileup $params.ivar_mpileup_args -A -B -Q 0 $bam | ivar consensus -p ${id}_${idx_name}.consensus \
                    -q $params.ivar_min_qual \
                    -t $params.ivar_min_freq \
                    -m $params.ivar_min_depth \
                    -n N
                
            samtools mpileup $params.ivar_mpileup_args -A -B -Q 0 $bam | ivar variants -p ${id}_${idx_name}.variants \
                -q $params.ivar_min_qual \
                -t $params.ivar_min_freq \
                -m $params.ivar_min_depth \
                -r $reference 


            HEADER=\$(head -1 $reference)
            sed -i "1s~.*~\$HEADER~" ${id}_${idx_name}.consensus.fa

            mv ${id}_${idx_name}.consensus.fa ${id}_${idx_name}.consensus.fasta
        fi
        """
    } else {
        """
        echo "No reads found in aligned input read files"
        """
    }

}

