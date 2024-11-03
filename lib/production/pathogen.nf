/* ==================================
 * PATHOGEN DETECTION FOR PRODUCTION
 * ==================================
*/


include { QualityControl; QualityControlNanopore } from "./quality";
include { Vircov; VircovNanopore } from "../processes/pathogen";
include { Kraken2; Bracken; Metabuli; Sylph; Kmcp; GanonReads; GanonProfile } from "../processes/pathogen";

include { Kraken2Nanopore; Bracken as BrackenNanopore; MetabuliNanopore; SylphNanopore } from "../processes/pathogen";

include { MetaSpades; Megahit; MetaSpadesNanopore; MegahitNanopore } from "../processes/pathogen";
include { ContigCoverage as MetaSpadesCoverage; ContigCoverage as MegahitCoverage } from "../processes/pathogen";
include { ContigCoverageNanopore as MetaSpadesCoverageNanopore; ContigCoverageNanopore as MegahitCoverageNanopore } from "../processes/pathogen";
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

        if (params.pathogenDetection.taxonomicProfile.enabled) {
            TaxonomicProfile(
                QualityControl.out.reads, 
                taxonomicProfileDatabases
            )
        }
        
        /* Metagenome assembly and taxonomic profiling module */

        if (params.pathogenDetection.metagenomeAssembly.enabled) {
            MetagenomeAssembly(
                QualityControl.out.reads,
                metagenomeAssemblyDatabases
            )
        }


}

workflow PathogenDetectionNanopore {

    take:
        reads
        qualityControlDatabases
        taxonomicProfileDatabases
        metagenomeAssemblyDatabases
    main:

        /* Read and background controls module */

        QualityControlNanopore(
            reads, 
            qualityControlDatabases
        )

        /* Taxonomic read classification and profiling module */

        if (params.pathogenDetection.taxonomicProfile.enabled) {
            TaxonomicProfileNanopore(
                QualityControlNanopore.out.reads, 
                taxonomicProfileDatabases
            )
        }
        
        /* Metagenome assembly and taxonomic profiling module */

        if (params.pathogenDetection.metagenomeAssembly.enabled) {
            MetagenomeAssemblyNanopore(
                QualityControlNanopore.out.reads,
                metagenomeAssemblyDatabases
            )
        }


}


workflow TaxonomicProfile {
    take:
        reads
        databases
    main:
        profileParams = params.pathogenDetection.taxonomicProfile

        if (profileParams.alignment) {
            Vircov(
                reads,
                databases.vircovDatabase,
                profileParams.alignmentMethod,
                profileParams.alignmentSecondary,
                params.resources.threads.vircovRemap,
                params.resources.threads.vircovParallel
            )
        }

        if (profileParams.classifier && profileParams.classifierMethod.contains("kraken2")) {
            Kraken2(
                reads,
                databases.krakenDatabase,
                profileParams.krakenConfidence
            )
            if (profileParams.profiler && profileParams.profilerMethod.contains("bracken")) {
                Bracken(
                    Kraken2.out.bracken,
                    profileParams.brackenReadLength,
                    profileParams.brackenRank,
                    profileParams.brackenMinReads
                )
            }
        }
        
        if (profileParams.classifier && profileParams.classifierMethod.contains("metabuli")) {
            Metabuli(
                reads,
                databases.metabuliDatabase
            )
        }

        if (profileParams.classifier && profileParams.classifierMethod.contains("ganon")) {
            GanonReads(
                reads,
                databases.ganonDatabase,
                profileParams.ganonMultipleMatches
            )
        }

        if (profileParams.profiler && profileParams.profilerMethod.contains("kmcp")) || 
           (profileParams.classifier && profileParams.classifierMethod.contains("kmcp")) {
            Kmcp(
                reads,
                databases.kmcpDatabase,
                profileParams.kmcpMode,
                profileParams.kmcpLevel
            )
        }

        if (profileParams.profiler && profileParams.profilerMethod.contains("ganon")) {
            GanonProfile(
                reads,
                databases.ganonDatabase,
                profileParams.ganonMultipleMatches
            )
        }


        if (profileParams.profiler && profileParams.profilerMethod.contains("sylph")) {
            Sylph(
                reads,
                databases.sylphDatabase,
                databases.sylphMetadata,
                profileParams.sylphMinNumberKmers,
                profileParams.sylphQueryCompression
            )
        }
}



workflow TaxonomicProfileNanopore {
    take:
        reads
        databases
    main:
        profileParams = params.pathogenDetection.taxonomicProfile

        if (profileParams.alignment) {
            VircovNanopore(
                reads,
                databases.vircovDatabase,
                profileParams.alignmentMethod,
                profileParams.alignmentSecondary,
                params.resources.threads.vircovRemap,
                params.resources.threads.vircovParallel
            )
        }

        if (profileParams.classifier && profileParams.classifierMethod.contains("kraken2")) {
            Kraken2Nanopore(
                reads,
                databases.krakenDatabase,
                profileParams.krakenConfidence
            )
            if (profileParams.profiler && profileParams.profilerMethod.contains("bracken")) {
                BrackenNanopore(
                    Kraken2Nanopore.out.bracken,
                    profileParams.brackenReadLength,
                    profileParams.brackenRank,
                    profileParams.brackenMinReads
                )
            }
        }
        
        if (profileParams.classifier && profileParams.classifierMethod.contains("metabuli")) {
            MetabuliNanopore(
                reads,
                databases.metabuliDatabase
            )
        }

        if (profileParams.profiler && profileParams.profilerMethod.contains("sylph")) {
            SylphNanopore(
                reads,
                databases.sylphDatabase,
                databases.sylphMetadata,
                profileParams.sylphMinNumberKmers,
                profileParams.sylphQueryCompression
            )
        }
}



workflow MetagenomeAssemblyNanopore {
    take:
        reads
        databases
    main:

        magParams = params.pathogenDetection.metagenomeAssembly

        if (magParams.assemblyMethod.contains("metaspades")) {
            
            metaspadesAssemblyCoverage = MetaSpadesNanopore(
                reads,
                magParams.assemblyKmerList,
                magParams.assemblyMinContigLength,
                magParams.assemblyArgs
            ) | MetaSpadesCoverageNanopore
            
            if (magParams.binningMethod.contains("concoct")) {
                MetaSpadesConcoct(
                    metaspadesAssemblyCoverage,
                    magParams.binningChunkSize,
                    magParams.binningReadLength,
                    magParams.binningMinBinSize,
                    0
                )
            }
            if (magParams.binningMethod.contains("metabat2")) {
                MetaSpadesMetabat2(
                    metaspadesAssemblyCoverage,
                    magParams.binningMinBinSize,
                    magParams.binningMinContigLength
                )
            }
        }

        if (magParams.assemblyMethod.contains("megahit")) {
            
            megahitAssemblyCoverage = MegahitNanopore(
                reads,
                magParams.assemblyKmerList,
                magParams.assemblyMinContigLength,
                magParams.assemblyArgs
            ) | MegahitCoverageNanopore
            
            if (magParams.binningMethod.contains("concoct")) {
                MegahitConcoct(
                    megahitAssemblyCoverage,
                    magParams.binningChunkSize,
                    magParams.binningReadLength,
                    magParams.binningMinBinSize,
                    0
                )
            }

            if (magParams.binningMethod.contains("metabat2")) {
                MegahitMetabat2(
                    megahitAssemblyCoverage,
                    magParams.binningMinBinSize,
                    magParams.binningMinContigLength,
                )
            }
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
                    0
                )
            }
            if (magParams.binningMethod.contains("metabat2")) {
                MetaSpadesMetabat2(
                    metaspadesAssemblyCoverage,
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
                    0
                )
            }

            if (magParams.binningMethod.contains("metabat2")) {
                MegahitMetabat2(
                    megahitAssemblyCoverage,
                    magParams.binningMinBinSize,
                    magParams.binningMinContigLength,
                )
            }
        }

        


        

        
}
