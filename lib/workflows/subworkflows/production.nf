/* 
===================
PRODUCTION WORKFLOW
===================
*/

include { PingServer; ParseSample; UploadSample } from '../../processes/cerebro'

workflow cerebro {
    take:
        taxonomy_directory
        result_files
        sample_sheet
        config_file
    main:
        if (params.cerebro_upload){
            PingServer(result_files) | collect  // await result before proceeding
        }
        ParseSample(result_files, taxonomy_directory)
        if (params.cerebro_upload){
            UploadSample(ParseSample.out.cerebro, sample_sheet, config_file)
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