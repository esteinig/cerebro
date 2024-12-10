/* ==================================
 * PATHOGEN DETECTION FOR PRODUCTION
 * ==================================
*/


include { QualityControl; QualityControlNanopore } from "./quality";
include { Vircov; VircovNanopore } from "../processes/pathogen";
include { Kraken2; Bracken; Metabuli; Sylph; Kmcp; GanonReads; GanonProfile; ProcessOutputIllumina; PathogenDetectionTable } from "../processes/pathogen";

include { Kraken2Nanopore; Bracken as BrackenNanopore; MetabuliNanopore; SylphNanopore; KmcpNanopore; GanonReadsNanopore; GanonProfileNanopore; ProcessOutputNanopore; PathogenDetectionTable as PathogenDetectionTableNanopore  } from "../processes/pathogen";

include { MetaSpades; Megahit; MetaSpadesNanopore; MegahitNanopore } from "../processes/pathogen";
include { ContigCoverage as MetaSpadesCoverage; ContigCoverage as MegahitCoverage } from "../processes/pathogen";
include { ContigCoverageNanopore as MetaSpadesCoverageNanopore; ContigCoverageNanopore as MegahitCoverageNanopore } from "../processes/pathogen";
include { Concoct as MetaSpadesConcoct; Concoct as MegahitConcoct } from "../processes/pathogen";
include { Metabat2 as MetaSpadesMetabat2; Metabat2 as MegahitMetabat2 } from "../processes/pathogen";
include { SemiBin2 as MegahitSemiBin2; SemiBin2 as MetaSpadesSemiBin2 } from "../processes/pathogen";
include { BlastContigs as MegahitBlastNcbi; BlastContigs as MetaSpadesBlastNcbi } from "../processes/pathogen";

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
                taxonomicProfileDatabases,
                QualityControl.out.results
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
                taxonomicProfileDatabases,
                QualityControlNanopore.out.results
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
        qualityControlResults
    main:
        profileParams = params.pathogenDetection.taxonomicProfile

        if (profileParams.alignment) {
            Vircov(
                reads,
                databases.vircovDatabase,
                profileParams.alignmentMethod,
                profileParams.alignmentSecondary,
                params.resources.threads.vircovRemap,
                params.resources.threads.vircovParallel,
                profileParams.vircovArgs
            )
        }

        if (profileParams.classifier && profileParams.classifierMethod.contains("kraken2")) {
            Kraken2(
                reads,
                databases.krakenDatabase,
                profileParams.krakenConfidence,
                profileParams.krakenMemoryMapping
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
                profileParams.ganonDatabasePrefix,
                profileParams.ganonMultipleMatches
            )
        }

        if ((profileParams.profiler && profileParams.profilerMethod.contains("kmcp")) || (profileParams.classifier && profileParams.classifierMethod.contains("kmcp"))) {
            Kmcp(
                reads,
                databases.kmcpDatabase,
                profileParams.kmcpMode,
                profileParams.kmcpLevel,
                profileParams.kmcpMinQueryCoverage
            )
        }

        if (profileParams.profiler && profileParams.profilerMethod.contains("ganon")) {
            GanonProfile(
                reads,
                databases.ganonDatabase,
                profileParams.ganonDatabasePrefix,
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


        vircovResults = profileParams.alignment  ? Vircov.out.results : Channel.empty()

        results = vircovResults.mix(
            (profileParams.profiler && profileParams.profilerMethod.contains("kmcp")) || (profileParams.classifier && profileParams.classifierMethod.contains("kmcp")) ? Kmcp.out.results : Channel.empty(),
            (profileParams.classifier && profileParams.classifierMethod.contains("kraken2")) && (profileParams.profiler && profileParams.profilerMethod.contains("bracken")) ? Bracken.out.results : Channel.empty(),
            (profileParams.classifier && profileParams.classifierMethod.contains("metabuli")) ? Metabuli.out.results : Channel.empty(),
            (profileParams.classifier && profileParams.classifierMethod.contains("ganon")) ? GanonReads.out.results : Channel.empty(),
            (profileParams.classifier && profileParams.classifierMethod.contains("kraken2")) ? Kraken2.out.results : Channel.empty(),
            (profileParams.profiler && profileParams.profilerMethod.contains("ganon")) ? GanonProfile.out.results : Channel.empty(),
            (profileParams.profiler && profileParams.profilerMethod.contains("sylph")) ? Sylph.out.results : Channel.empty()
        )

        json = results.mix(qualityControlResults) | groupTuple | map { d -> [d[0], d[1..-1].flatten()] } | ProcessOutputIllumina
        tables = PathogenDetectionTable(json | collect, databases.taxonomy)
    
    emit:
        reads   = reads
        results = results
        tables  = tables

}



workflow TaxonomicProfileNanopore {
    take:
        reads
        databases
        qualityControlResults
    main:
        profileParams = params.pathogenDetection.taxonomicProfile

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
                profileParams.krakenConfidence,
                profileParams.krakenMemoryMapping
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

        if (profileParams.classifier && profileParams.classifierMethod.contains("ganon")) {
            GanonReadsNanopore(
                reads,
                databases.ganonDatabase,
                profileParams.ganonDatabasePrefix,
                profileParams.ganonMultipleMatches
            )
        }

        if ((profileParams.profiler && profileParams.profilerMethod.contains("kmcp")) || (profileParams.classifier && profileParams.classifierMethod.contains("kmcp"))) {
            KmcpNanopore(
                reads,
                databases.kmcpDatabase,
                profileParams.kmcpMode,
                profileParams.kmcpLevel,
                profileParams.kmcpMinQueryCoverage
            )
        }

        if (profileParams.profiler && profileParams.profilerMethod.contains("ganon")) {
            GanonProfileNanopore(
                reads,
                databases.ganonDatabase,
                profileParams.ganonDatabasePrefix,
                profileParams.ganonMultipleMatches
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


        vircovResults = profileParams.alignment  ? VircovNanopore.out.results : Channel.empty()

        results = vircovResults.mix(
            (profileParams.profiler && profileParams.profilerMethod.contains("kmcp")) || (profileParams.classifier && profileParams.classifierMethod.contains("kmcp")) ? KmcpNanopore.out.results : Channel.empty(),
            (profileParams.classifier && profileParams.classifierMethod.contains("kraken2")) && (profileParams.profiler && profileParams.profilerMethod.contains("bracken")) ? BrackenNanopore.out.results : Channel.empty(),
            (profileParams.classifier && profileParams.classifierMethod.contains("metabuli")) ? MetabuliNanopore.out.results : Channel.empty(),
            (profileParams.classifier && profileParams.classifierMethod.contains("ganon")) ? GanonReadsNanopore.out.results : Channel.empty(),
            (profileParams.classifier && profileParams.classifierMethod.contains("kraken2")) ? Kraken2Nanopore.out.results : Channel.empty(),
            (profileParams.profiler && profileParams.profilerMethod.contains("ganon")) ? GanonProfileNanopore.out.results : Channel.empty(),
            (profileParams.profiler && profileParams.profilerMethod.contains("sylph")) ? SylphNanopore.out.results : Channel.empty()
        )


        // process results to json and get tables
        json = results.mix(qualityControlResults) | groupTuple | map { d -> [d[0], d[1..-1].flatten()] } | ProcessOutputNanopore
        tables = PathogenDetectionTableNanopore(json | collect, databases.taxonomy)
    
    emit:
        reads   = reads
        results = results
        tables  = tables
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
            if (magParams.binningMethod.contains("semibin2")) {
                MetaSpadesSemiBin2(
                    metaspadesAssemblyCoverage,
                    magParams.binningMinContigLength
                )
            }

            if (magParams.ncbiDatabase && magParams.ncbiDatabaseMethod.contains("blast")) {
                MetaSpadesBlastNcbi(
                    MetaSpades.out.contigs,
                    databases.ncbiDatabase
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

            if (magParams.binningMethod.contains("semibin2")) {
                MegahitSemiBin2(
                    megahitAssemblyCoverage,
                    magParams.binningMinContigLength
                )
            }

            if (magParams.ncbiDatabase && magParams.ncbiDatabaseMethod.contains("blast")) {

                MegahitBlastNcbi(
                    Megahit.out.contigs,
                    databases.ncbiDatabase
                )
            }
        }

        


        

        
}
