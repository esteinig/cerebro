

workflow bacterial_species_id {
    take:
        ont_reads
        pe_reads
    main:
        if (ont_reads && pe_reads == null) {
            assemblies = NanoqScan(ont_reads) | Nanoq | Dragonflye 
        } else if (ont_reads && pe_reads) {
            ont_reads = NanoqScan(ont_reads) | Nanoq
            pe_reads = 
        }

        annotations = Bakta(assemblies) | PlotAnnotations
        gtdb_species = GtdbSpecies(assemblies) | PlotTaxonomy
        
    
    emit:
         
}