/* 
===================
PRODUCTION WORKFLOW
===================
*/

include { PingServer; ProcessSamples; UploadSample } from '../../processes/cerebro'

workflow cerebro {
    take:
        taxonomy_directory
        result_files
        sample_sheet
        config_file
    main:
        if (params.cerebro_upload){
            PingServer(result_files) | collect  // await results before proceeding
        }

        ProcessSamples(result_files, taxonomy_directory)
        
        if (params.cerebro_upload){
            UploadSample(ProcessSamples.out.cerebro, sample_sheet, config_file)
        }
}

workflow cerebro_production {
    take:
        taxonomy_directory
        result_files
        sample_sheet
        config_file
    main:
        cerebro(taxonomy_directory, result_files, sample_sheet, config_file)
}