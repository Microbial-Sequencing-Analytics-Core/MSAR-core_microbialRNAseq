nextflow.enable.dsl=2

include { fastp                    } from './modules/fastp'
include { fastqc                   } from './modules/fastqc'
include { multiqc                  } from './modules/multiqc'
include { fastq_screen             } from './modules/fastq_screen'
include { bowtie2                  } from './modules/bowtie2'
include { split_reads_from_unmapped } from './modules/split_reads_from_unmapped'
include { star_host_removal        } from './modules/star_host_removal'
include { bwa_mem                  } from './modules/bwa_mem'
include { samtools                 } from './modules/samtools'
include { filt_bam                 } from './modules/filt_bam'
include { get_mapped_reads         } from './modules/get_mapped_reads'
include { kraken2                  } from './modules/kraken2'
include { merge_idxstats           } from './modules/merge_idxstats'
include { kraken_biom              } from './modules/kraken_biom'


ch_reads = Channel
    .fromPath(params.samplesheet, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
        def is_singleEnd = !row.fq2
        def meta  = [sample_name: row.sample, single_end: is_singleEnd]
        def reads = is_singleEnd
            ? [ file(row.fq1, checkIfExists: true) ]
            : [ file(row.fq1, checkIfExists: true), file(row.fq2, checkIfExists: true) ]
        [meta, reads]
    }


workflow {
    println "Reference FASTA: ${params.bwa_index}"

    fastqc(ch_reads)
    fastp(ch_reads)
    fastq_screen(fastp.out.trim_reads)

    if (params.host_removal == "RNA") {
        star_host_removal(fastp.out.trim_reads)
        ch_for_bwa = star_host_removal.out.split_reads
    } else if (params.host_removal == "DNA") {
        bowtie2(fastp.out.trim_reads)
        split_reads_from_unmapped(bowtie2.out.bowtie2_bam_both_unmapped_bam)
        ch_for_bwa = split_reads_from_unmapped.out.split_reads
    } else {
        ch_for_bwa = fastp.out.trim_reads
    }

    bwa_mem(ch_for_bwa)
    samtools(bwa_mem.out.bam_file)
    merge_idxstats(samtools.out.idxstats_out.collect())
    multiqc(
        fastq_screen.out.fastq_screen_out
            .mix(fastqc.out.zips, fastp.out.fastp_json, samtools.out.stats_out)
            .collect()
    )

    if (params.run_kraken) {
        filt_bam(bwa_mem.out.bam_file)
        get_mapped_reads(filt_bam.out.filt_bam)
        kraken2(get_mapped_reads.out.mapped_reads)
        kraken_biom(kraken2.out.kraken2_report.collect())
    }
}

workflow.onComplete {
    println(workflow.success ? "\nDone!" : "Oops .. something went wrong")
}
