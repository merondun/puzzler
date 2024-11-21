process RAGTAG {
    publishDir "${params.wd}/Primary_Assemblies", mode: 'copy'

    input:
        tuple val(sample), path(scaffolds)
        path reference
    
    output:
        tuple val(sample), path("${sample}.HifiasmHifiHiC-PurgeDups-HapHiC-RagTag.fa"), emit: final_assembly
        path "${sample}_final/*", emit: busco_results
    
    script:
        """
        ragtag.py scaffold ${reference} ${scaffolds} -o ragtag_output --aligner minimap2 -t 64
        sed 's/_RagTag//g' ragtag_output/ragtag.scaffold.fasta > \\
            ${sample}.HifiasmHifiHiC-PurgeDups-HapHiC-RagTag.fa
        
        busco -i ${sample}.HifiasmHifiHiC-PurgeDups-HapHiC-RagTag.fa \\
            -l embryophyta_odb10 -m genome -c 64 -f \\
            --out_path . --out ${sample}_final
        """
}
