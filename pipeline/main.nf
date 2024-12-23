nextflow.enable.dsl = 2

include { HIFIASM } from './modules/hifiasm'
include { PURGE_DUPS } from './modules/purge_dups'
include { HAP_HIC } from './modules/hap_hic'
include { RAGTAG } from './modules/ragtag'

workflow {
    input_ch = Channel.fromPath(params.samples)
        .splitCsv(header:true)
        .map { row -> tuple(
            row.sample,
            row.ploidy,
            file(row.hifi),
            file(row.hic_r1),
            file(row.hic_r2)
        )}
        
    reference = file(params.reference)
    
    HIFIASM(input_ch)
    PURGE_DUPS(HIFIASM.out.gfa)
    HAP_HIC(PURGE_DUPS.out.purged_assembly
        .join(input_ch.map { tuple(it[0], it[3], it[4]) }))
    RAGTAG(HAP_HIC.out.scaffolds, reference)
}