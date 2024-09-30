/* ==================================
 * PATHOGEN DETECTION FOR PRODUCTION
 * ==================================
*/


include { QualityControl } from "./quality";
include { Kraken2; Bracken; Metabuli; Sylph; Kmcp } from "../processes/pathogen";
include { MetaSpades; Megahit } from "../processes/pathogen";
include { ContigCoverage as MetaSpadesCoverage; ContigCoverage as MegahitCoverage } from "../processes/pathogen";
include { Concoct as MetaSpadesConcoct; Concoct as MegahitConcoct } from "../processes/pathogen";
include { Metabat2 as MetaSpadesMetabat2; Metabat2 as MegahitMetabat2 } from "../processes/pathogen";


workflow PathogenDetection {

    take:
        reads
        qualityControlDatabases
        taxonomicProfileDatabases
        metagenomeAssemblyDatabases
    main:

        /* Read and background controls module */

        QualityControl(
            reads, 
            qualityControlDatabases
        )

        /* Taxonomic read classification and profiling module */

        TaxonomicProfile(
            QualityControl.out.reads, 
            taxonomicProfileDatabases
        )

        /* Metagenome assembly and taxonomic profiling module */

        MetagenomeAssembly(
            QualityControl.out.reads,
            metagenomeAssemblyDatabases
        )



}


workflow TaxonomicProfile {
    take:
        reads
        databases
    main:
        profileParams = params.pathogenDetection.taxonomicProfile

        if (profileParams.classifierMethod.contains("kraken2")) {
            Kraken2(
                reads,
                databases.krakenDatabase,
                profileParams.krakenConfidence
            )
            if (profileParams.profilerMethod.contains("bracken")) {
                Bracken(
                    Kraken2.out.bracken,
                    profileParams.brackenReadLength,
                    profileParams.brackenRank,
                    profileParams.brackenMinReads
                )
            }
        }
        
        if (profileParams.classifierMethod.contains("metabuli")) {
            Metabuli(
                reads,
                databases.metabuliDatabase
            )
        }

        if (profileParams.profilerMethod.contains("sylph")) {
            Sylph(
                reads,
                databases.sylphDatabase,
                databases.sylphMetadata
            )
        }

        if (profileParams.profilerMethod.contains("kmcp")) {
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

        if (magParams.assemblyMethod.contains("metaspades")) {
            
            metaspadesAssemblyCoverage = MetaSpades(
                reads,
                magParams.assemblyKmerList,
                magParams.assemblyMinContigLength,
                magParams.assemblyArgs
            ) | MetaSpadesCoverage
            
            if (magParams.binningMethod.contains("concoct")) {
                MetaSpadesConcoct(
                    metaspadesAssemblyCoverage,
                    magParams.binningChunkSize,
                    magParams.binningReadLength,
                    magParams.binningMinBinSize,
                    magParams.binningMinContigLength
                )
            }
            if (magParams.binningMethod.contains("metabat2")) {
                MetaSpadesMetabat2(
                    metaspadesAssemblyCoverage,
                    magParams.binningChunkSize,
                    magParams.binningReadLength,
                    magParams.binningMinBinSize,
                    magParams.binningMinContigLength
                )
            }
        }

        if (magParams.assemblyMethod.contains("megahit")) {
            
            megahitAssemblyCoverage = Megahit(
                reads,
                magParams.assemblyKmerList,
                magParams.assemblyMinContigLength,
                magParams.assemblyArgs
            ) | MegahitCoverage
            
            if (magParams.binningMethod.contains("concoct")) {
                MegahitConcoct(
                    megahitAssemblyCoverage,
                    magParams.binningChunkSize,
                    magParams.binningReadLength,
                    magParams.binningMinBinSize,
                    magParams.binningMinContigLength
                )
            }

            if (magParams.binningMethod.contains("metabat2")) {
                MegahitMetabat2(
                    megahitAssemblyCoverage,
                    magParams.binningChunkSize,
                    magParams.binningReadLength,
                    magParams.binningMinBinSize,
                    magParams.binningMinContigLength,
                )
            }
        }

        


        

        
}
