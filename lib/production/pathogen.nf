/* ==================================
 * PATHOGEN DETECTION FOR PRODUCTION
 * ==================================
*/


include { QualityControl } from "./quality";
include { Kraken2; Bracken; Metabuli; Sylph; Kmcp } from "../processes/pathogen";


workflow PathogenDetection {

    take:
        reads
        qualityControlDatabases
        taxonomicProfileDatabases
    main:

        /* Read and background control module */

        QualityControl(
            reads, 
            qualityControlDatabases
        )

        /* Taxonomic classification and profiling module */

        TaxonomicProfile(
            QualityControl.out.reads, 
            taxonomicProfileDatabases
        )

        /* Metagenome assembly and taxonomic profiling module */




}


workflow TaxonomicProfile {
    take:
        reads
        databases
    main:
        profileParams = params.pathogenDetection.taxonomicProfile

        if (profileParams.classifier.contains("kraken2")) {
            Kraken2(
                reads,
                databases.krakenDatabase,
                profileParams.krakenConfidence
            )

            if (profileParams.classifierBrackenProfile) {
                Bracken(
                    Kraken2.out.bracken,
                    profileParams.brackenReadLength,
                    profileParams.brackenRank,
                    profileParams.brackenMinReads
                )
            }
        }
        
        if (profileParams.classifier.contains("metabuli")) {
            Metabuli(
                reads,
                databases.metabuliDatabase
            )
        }

        if (profileParams.classifier.contains("sylph")) {
            Sylph(
                reads,
                databases.sylphDatabase,
                databases.sylphMetadata
            )
        }

        if (profileParams.classifier.contains("kmcp")) {
            Kmcp(
                reads,
                databases.kmcpDatabase,
                profileParams.kmcpMode
            )
        }
}

workflow MetagenomeAssembly {
    take:
        reads
        databases
    main:
        magParams = params.pathogenDetection.metagenomeAssembly

        if (magParams.assembler.contains("metaspades")) {
            contigs = MetaSpades(
                reads,
                magParams.metaspadesKmer
            ).contigs
        } else if (magParams.assembler.contains("megahit") {
            contigs = Megahit(
                reads,
                magParams.megahitKmer
            ).contigs
        } else {
            contigs = Channel.empty()
        }

        

        
}
