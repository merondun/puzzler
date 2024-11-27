process RAGTAG {
    errorStrategy 'ignore'
    afterScript 'echo "Sample: $sample, Exit: $exit_status" >> failed_samples.log'
    tag "${sample}"

    publishDir path: "${params.wd}/Primary_Assemblies", mode: 'copy',
        pattern: "*.HifiasmHifiHiC-PurgeDups-HapHiC-RagTag.fa"
    publishDir path: "${params.wd}/${sample}/HiCHiFi/04_RagTag", mode: 'copy',
        pattern: "ragtag_output"

    label 'ragtag'

    input:
        tuple val(sample), path(scaffolds)
        path reference
    
    output:
        tuple val(sample), path("${sample}.HifiasmHifiHiC-PurgeDups-HapHiC-RagTag.fa"), optional: true, emit: final_assembly
        path "ragtag_output", optional: true
    
    script:
        """
        ragtag.py scaffold ${reference} ${scaffolds} -o ragtag_output --aligner minimap2 -t ${task.cpus}
        sed 's/_RagTag//g' ragtag_output/ragtag.scaffold.fasta > \\
            ${sample}.HifiasmHifiHiC-PurgeDups-HapHiC-RagTag.fa
        """
}
