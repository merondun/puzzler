process HAP_HIC {
    publishDir "${params.wd}/${sample}/HiCHiFi/03_HapHiC", mode: 'copy'
    
    input:
        tuple val(sample), path(purged_assembly), path(hic_r1), path(hic_r2)
    
    output:
        tuple val(sample), path("04.build/scaffolds.fa"), emit: scaffolds
    
    script:
        """
        bwa index ${purged_assembly}
        bwa mem -5SP -t 64 ${purged_assembly} ${hic_r1} ${hic_r2} | \\
            samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o HiC.bam
        
        /project/coffea_pangenome/Software/Merondun/HapHiC/utils/filter_bam \\
            HiC.bam 1 --nm 3 --threads 64 | \\
            samtools view - -b -@ 64 -o HiC.filtered.bam
        
        haphic pipeline ${purged_assembly} HiC.filtered.bam 28 \\
            --RE "GATC,GANTC" --correct_nrounds 2 \\
            --threads 64 --processes 64 --max_inflation 7.0
        """
}

