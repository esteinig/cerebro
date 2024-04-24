
/* 
===========================
BACTERIAL SUBSET ALIGNMENT
===========================
*/

include { MashScreenWinner as MashScreenBacteria } from '../../processes/mash' addParams(
    subdir: "alignment/bacteria/mash_subset"
)
include { MashDatabaseSubset as MashDatabaseSubsetBacteria } from '../../processes/cerebro' addParams(
    subdir: "alignment/bacteria/mash_subset",
    min_shared_hashes: params.subset_min_shared_hashes
)
include { MinimapIndexSubset as MinimapIndexSubsetBacteria } from '../../processes/minimap2' addParams(
    subdir: "alignment/bacteria/mash_subset"
)
include { MinimapAlignSubsetPAF as MinimapAlignSubsetBacteria } from '../../processes/minimap2' addParams(
    subdir: "alignment/bacteria/mash_subset"
)
include { VircovSubsetAlign as VircovCoverageSubsetBacteria } from '../../processes/vircov' addParams(
    subdir: "alignment/bacteria/mash_subset",
    vircov_min_reads: params.bacteria_min_reads,
    vircov_min_regions: params.bacteria_min_regions,
    vircov_min_cov: params.bacteria_min_cov,
    vircov_min_len: params.bacteria_min_len,
    vircov_min_mapq: params.bacteria_min_mapq
)

workflow bacteria_subset_alignment {
    take:
        reads                                                                                     
        bacteria_mash_index                                                                    
        bacteria_fasta
    main:
        MashScreenBacteria(reads, bacteria_mash_index)                                 
        MashDatabaseSubsetBacteria(MashScreenBacteria.out, bacteria_fasta)       
        MinimapIndexSubsetBacteria(MashDatabaseSubsetBacteria.out)                     
        MinimapAlignSubsetBacteria(MinimapIndexSubsetBacteria.out)
        VircovCoverageSubsetBacteria(MinimapAlignSubsetBacteria.out)         
    emit:
        reads = reads
        results = VircovCoverageSubsetBacteria.out.results
}

/* 
===========================
BACTERIAL SUBSET ALIGNMENT
===========================
*/

include { MashScreenWinner as MashScreenEukaryots } from '../../processes/mash' addParams(
    subdir: "alignment/eukaryots/mash_subset"
)
include { MashDatabaseSubset as MashDatabaseSubsetEukaryots } from '../../processes/cerebro' addParams(
    subdir: "alignment/eukaryots/mash_subset",
    min_shared_hashes: params.subset_min_shared_hashes
)
include { MinimapIndexSubset as MinimapIndexSubsetEukaryots } from '../../processes/minimap2' addParams(
    subdir: "alignment/eukaryots/mash_subset"
)
include { MinimapAlignSubsetPAF as MinimapAlignSubsetEukaryots } from '../../processes/minimap2' addParams(
    subdir: "alignment/eukaryots/mash_subset"
)
include { VircovSubsetAlign as VircovCoverageSubsetEukaryots } from '../../processes/vircov' addParams(
    subdir: "alignment/eukaryots/mash_subset",
    vircov_min_reads: params.eukaryots_min_reads,
    vircov_min_regions: params.eukaryots_min_regions,
    vircov_min_cov: params.eukaryots_min_cov,
    vircov_min_len: params.eukaryots_min_len,
    vircov_min_mapq: params.eukaryots_min_mapq
)

workflow eukaryots_subset_alignment {
    take:
        reads                                                                                     
        eukaryots_mash_index                                                                    
        eukaryots_fasta
    main:
        MashScreenEukaryots(reads, eukaryots_mash_index)                                 
        MashDatabaseSubsetEukaryots(MashScreenEukaryots.out, eukaryots_fasta)       
        MinimapIndexSubsetEukaryots(MashDatabaseSubsetEukaryots.out)                     
        MinimapAlignSubsetEukaryots(MinimapIndexSubsetEukaryots.out)
        VircovCoverageSubsetEukaryots(MinimapAlignSubsetEukaryots.out)      
    emit:
        reads = reads  
        results = VircovCoverageSubsetEukaryots.out.results  
}
