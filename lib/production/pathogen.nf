/* ==================================
 * PATHOGEN DETECTION FOR PRODUCTION
 * ==================================
*/


include { QualityControl; QualityControlNanopore } from "./quality";
include { Vircov; VircovNanopore } from "../processes/pathogen";
include { Kraken2; Bracken; Metabuli; Sylph; Kmcp; GanonReads; GanonProfile; Megahit; MetaSpades; ProcessOutputIllumina; PathogenDetectionTable } from "../processes/pathogen";
include { ContigCoverage as MetaSpadesCoverage; ContigCoverage as MegahitCoverage } from "../processes/pathogen";
include { Concoct as MetaSpadesConcoct; Concoct as MegahitConcoct } from "../processes/pathogen";
include { Metabat2 as MetaSpadesMetabat2; Metabat2 as MegahitMetabat2 } from "../processes/pathogen";
include { SemiBin2 as MegahitSemiBin2; SemiBin2 as MetaSpadesSemiBin2 } from "../processes/pathogen";
include { BlastContigs as MegahitBlast; BlastContigs as MetaSpadesBlast } from "../processes/pathogen";

include { PipelineConfig } from "./utils";

workflow PathogenDetection {

    take:
        reads
        qualityControlDatabases
        taxonomicProfileDatabases
        metagenomeAssemblyDatabases
        productionConfig
        stagedFileData
    main:

        cerebroWorkflow = "pathogen-detection"
        workflowStarted = java.time.LocalDateTime.now()

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


        results = QualityControl.out.results.mix(
            params.pathogenDetection.taxonomicProfile.enabled ? TaxonomicProfile.out.results : Channel.empty(),
            params.pathogenDetection.metagenomeAssembly.enabled ? MetagenomeAssembly.out.results : Channel.empty()
        )


        results | groupTuple | map { d -> [d[0], d[1..-1].flatten()] } | ProcessOutputIllumina
        
        PathogenDetectionTable(
            ProcessOutputIllumina.out.results | collect, 
            taxonomicProfileDatabases.taxonomy
        )

        PipelineConfig(
            ProcessOutputIllumina.out.samples | collect,
            cerebroWorkflow, 
            workflowStarted
        )

        if (params.cerebroProduction.enabled) {
            
            UploadOutput(
                stagedFileData.mix(ProcessOutputIllumina.out.results) | groupTuple | map { d -> d.flatten() }, 
                taxonomicProfileDatabases.taxonomy, 
                PipelineConfig.out.config,
                productionConfig.apiUrl,
                productionConfig.authToken
            )

        }

}

workflow TaxonomicProfile {
    take:
        reads
        databases
    main:
        profileParams = params.pathogenDetection.taxonomicProfile
        
        log.info "TaxonomicProfile parameters: ${profileParams}"

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



        // if (profileParams.profiler && profileParams.profilerMethod.contains("ganon")) {
        //     GanonProfile(
        //         reads,
        //         databases.ganonDatabase,
        //         profileParams.ganonDatabasePrefix,
        //         profileParams.ganonMultipleMatches
        //     )
        // }


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
            // (profileParams.profiler && profileParams.profilerMethod.contains("ganon")) ? GanonProfile.out.results : Channel.empty(),
            (profileParams.profiler && profileParams.profilerMethod.contains("sylph")) ? Sylph.out.results : Channel.empty()
        )

    emit:
        results = results

}


workflow MetagenomeAssembly {
    take:
        reads
        databases
    main:
    
        magParams = params.pathogenDetection.metagenomeAssembly

        filteredReads = reads.filter { read_id, fq1, fq2 ->
            !params.skipAssembly.any { skip -> read_id.contains(skip) }
        }

        if (magParams.assemblyMethod.contains("metaspades")) {
            
            metaspadesAssemblyCoverage = MetaSpades(
                filteredReads,
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

            
            if (magParams.contigProfile && magParams.contigProfileMethod.contains("blast")) {
                MetaSpadesBlast(
                    MetaSpades.out.contigs,
                    databases.contigProfile,
                    magParams.contigProfileBlastPrefix,
                    magParams.contigProfileBlastMinPercentIdentity,
                    magParams.contigProfileBlastMinEvalue,
                    magParams.contigProfileBlastMaxTargetSeqs,
                )
            }
        }

        if (magParams.assemblyMethod.contains("megahit")) {
            
            megahitAssemblyCoverage = Megahit(
                filteredReads,
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

            if (magParams.contigProfile && magParams.contigProfileMethod.contains("blast")) {
                MegahitBlast(
                    Megahit.out.contigs,
                    databases.contigProfile,
                    magParams.contigProfileBlastPrefix,
                    magParams.contigProfileBlastMinPercentIdentity,
                    magParams.contigProfileBlastMinEvalue,
                    magParams.contigProfileBlastMaxTargetSeqs,
                )
            }
        }

        results = Channel.empty().mix(
            (magParams.assemblyMethod.contains("megahit") && magParams.contigProfile && magParams.contigProfileMethod.contains("blast")) ? MegahitBlast.out.results : Channel.empty(),
            (magParams.assemblyMethod.contains("metaspades") && magParams.contigProfile && magParams.contigProfileMethod.contains("blast")) ? MetaSpadesBlast.out.results : Channel.empty()
        )

    emit:
        results = results
}