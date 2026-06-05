process bowtie2 {
    tag "bowtie2 on ${sample_name}"
    label "process_high"

    publishDir "${launchDir}/analysis/bowtie2"

    module 'bowtie2/2.3.4.1'
    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.mapped_unmapped.bam"), emit: bowtie2_mapped_unmapped_bam
    tuple val(sample_name), path("${sample_name}.both_unmapped.bam"), val(is_SE), emit: bowtie2_bam_both_unmapped_bam
    path("${sample_name}.mapped_unmapped.stats"), emit: samtools_stats

    script:
    def index = params.bowtie2.(params.genome)
    def input_reads = is_SE ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def unmapped_filter = is_SE ? "-f 4 -F 256" : "-f 12 -F 256"
    """
    bowtie2 -p ${task.cpus} -x ${index} \
    ${input_reads} \
    | samtools sort -O BAM -o ${sample_name}.mapped_unmapped.bam

    samtools index ${sample_name}.mapped_unmapped.bam
    samtools stats ${sample_name}.mapped_unmapped.bam > ${sample_name}.mapped_unmapped.stats

    samtools view -b ${unmapped_filter} \
    ${sample_name}.mapped_unmapped.bam > ${sample_name}.both_unmapped.bam
    """
}
