process bwa_mem {
    tag "${sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/bwa_mem"

    module "bwa/0.7.17"
    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.bam"), emit: bam_file

    // -M: mark shorter split hits as secondary
    script:
    def index = params.bwa_index
    def input_reads = is_SE ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    """
    bwa mem -M \
    -t ${task.cpus} \
    ${index} \
    ${input_reads} \
    | samtools sort \
    -O "BAM" \
    -o ${sample_name}.bam -
    """
}
