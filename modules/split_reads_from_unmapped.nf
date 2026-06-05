process split_reads_from_unmapped {
    tag "split reads - ${sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/split_reads"

    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(bam_file), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.host_remove.R{1,2}.fastq.gz"), val(is_SE), emit: split_reads

    script:
    if (is_SE)
    """
        samtools sort -n ${bam_file} -o ${sample_name}.sorted.bam
        samtools fastq ${sample_name}.sorted.bam \
            -0 ${sample_name}.host_remove.R1.fastq.gz -n
        touch ${sample_name}.host_remove.R2.fastq.gz
    """
    else
    """
        samtools sort -n ${bam_file} -o ${sample_name}.sorted.bam
        samtools fastq ${sample_name}.sorted.bam \
            -1 ${sample_name}.host_remove.R1.fastq.gz \
            -2 ${sample_name}.host_remove.R2.fastq.gz \
            -0 /dev/null -s /dev/null -n
    """
}
