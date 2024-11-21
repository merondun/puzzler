process HIFIASM {
    publishDir "${params.wd}/${sample}/HiCHiFi", mode: 'copy'

    input:
        tuple val(sample), val(ploidy), path(hifi), path(hic_r1), path(hic_r2)
    
    output:
        tuple val(sample), path("${sample}.hic.p_ctg.gfa"), emit: gfa
    
    script:
        """
        hifiasm --n-hap ${ploidy} -s 0.9 -t 64 -o ${sample} \\
            --h1 ${hic_r1} --h2 ${hic_r2} ${hifi}
        """
}
