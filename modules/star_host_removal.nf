process star_host_removal {
    tag "star on ${sample_name}"
    label "process_high"

    publishDir "${launchDir}/analysis/star"

    module "STAR/2.7.10a"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.host_remove.R{1,2}.fastq.gz"), val(is_SE), emit: split_reads
    path("${sample_name}.Log.final.out"), emit: star_log

    script:
    def index = params.star_index.(params.genome)
    if (is_SE)
    """
        STAR --runThreadN ${task.cpus} \
            --genomeDir ${index} \
            --readFilesIn ${reads[0]} \
            --readFilesCommand zcat \
            --outSAMtype None \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${sample_name}.

        gzip -c ${sample_name}.Unmapped.out.mate1 > ${sample_name}.host_remove.R1.fastq.gz
        touch ${sample_name}.host_remove.R2.fastq.gz
    """
    else
    """
        STAR --runThreadN ${task.cpus} \
            --genomeDir ${index} \
            --readFilesIn ${reads[0]} ${reads[1]} \
            --readFilesCommand zcat \
            --outSAMtype None \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${sample_name}.

        gzip -c ${sample_name}.Unmapped.out.mate1 > ${sample_name}.host_remove.R1.fastq.gz
        gzip -c ${sample_name}.Unmapped.out.mate2 > ${sample_name}.host_remove.R2.fastq.gz
    """
}
