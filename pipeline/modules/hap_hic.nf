process HAP_HIC {
    tag "${sample}"
    publishDir "${params.wd}/${sample}/HiCHiFi/03_HapHiC", mode: 'copy'

    label 'hap_hic'
    
    input:
        tuple val(sample), path(purged_assembly), path(hic_r1), path(hic_r2)
    
    output:
        tuple val(sample), path("04.build/scaffolds.fa"), emit: scaffolds
    
    script:
        """
        bwa index ${purged_assembly}
        bwa mem -5SP -t ${task.cpus} ${purged_assembly} ${hic_r1} ${hic_r2} | \\
            samblaster | samtools view - -@ ${task.cpus} -S -h -b -F 3340 -o HiC.bam
        
        /project/coffea_pangenome/Software/Merondun/HapHiC/utils/filter_bam \\
            HiC.bam 1 --nm 3 --threads ${task.cpus} | \\
            samtools view - -b -@ ${task.cpus} -o HiC.filtered.bam
        
        haphic pipeline ${purged_assembly} HiC.filtered.bam 28 \\
            --RE "GATC,GANTC" --correct_nrounds 2 \\
            --threads ${task.cpus} --processes ${task.cpus} --max_inflation 7.0
        """
}

