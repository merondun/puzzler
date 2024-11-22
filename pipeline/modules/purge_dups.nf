process PURGE_DUPS {
    tag "${sample}"
    publishDir "${params.wd}/${sample}/HiCHiFi/02_PurgeDups", mode: 'copy'

    label 'purge_dups'
    
    input:
        tuple val(sample), path(gfa), path(hifi)
    
    output:
        tuple val(sample), path("purged.fa"), emit: purged_assembly
    
    script:
        """
        gfatools gfa2fa ${gfa} > ${sample}.hic.p_ctg.fa

        split_fa ${sample}.hic.p_ctg.fa > asm.split
        minimap2 -t ${task.cpus} -xasm5 -DP asm.split asm.split | pigz -c > asm.split.self.paf.gz
        purge_dups -M1000 -E1000 asm.split.self.paf.gz > dups.bed
        get_seqs dups.bed ${sample}.hic.p_ctg.fa
        """
}
