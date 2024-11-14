
include { MetaSpades } from '../../processes/spades' addParams(
    subdir: "assembly/spades_meta/"
)
include { BlastNT } from '../../processes/blast' addParams(
    subdir: "assembly/spades_meta/blast"
)
include { DiamondNR } from '../../processes/diamond' addParams(
    subdir: "assembly/spades_meta/diamond"
)

workflow metagenome_assembly {
    take:
        reads                                                                                     
        blast_db
        diamond_db
    main:
        MetaSpades(reads)
        BlastNT(MetaSpades.out.contigs, blast_db)
        DiamondNR(MetaSpades.out.contigs, diamond_db)
    emit: 
        reads = reads
        contigs = MetaSpades.out.contigs
        results = BlastNT.out.results.mix(DiamondNR.out.results)
}