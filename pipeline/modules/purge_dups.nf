process PURGE_DUPS {
    publishDir "${params.wd}/${sample}/HiCHiFi/02_PurgeDups", mode: 'copy'
    
    input:
        tuple val(sample), path(gfa), path(hifi)
    
    output:
        tuple val(sample), path("purged.fa"), emit: purged_assembly
        path "purged/*", emit: busco_results
    
    script:
        """
        gfatools gfa2fa ${gfa} > ${sample}.hic.p_ctg.fa

        split_fa ${sample}.hic.p_ctg.fa > asm.split
        minimap2 -t 64 -xasm5 -DP asm.split asm.split | pigz -c > asm.split.self.paf.gz
        purge_dups -M1000 -E1000 asm.split.self.paf.gz > dups.bed
        get_seqs dups.bed ${sample}.hic.p_ctg.fa
        
        busco -i purged.fa -l embryophyta_odb10 -m genome -c 64 -f --out_path . --out purged
        """
}
