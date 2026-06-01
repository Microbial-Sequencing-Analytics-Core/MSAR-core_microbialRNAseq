nextflow.enable.dsl=2

process fastp {
    tag "${meta.sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/fastp/", mode : "copy"

    module "fastp/1.0.1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta.sample_name), path("${meta.sample_name}_trimmed_R{1,2}.fastq.gz"), val(meta.single_end), emit: trim_reads
    path("${meta.sample_name}.fastp.json"), emit: fastp_json

    // --trim_front2 1
    script:
    def adapter = "/mnt/beegfs/kimj32/reference/adapters.fa"
    //--adapter_fasta $adapter \
    if(!meta.single_end) {
    """
        fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        --thread ${task.cpus} \
        -o ${meta.sample_name}_trimmed_R1.fastq.gz \
        -O ${meta.sample_name}_trimmed_R2.fastq.gz \
        --json ${meta.sample_name}.fastp.json
    """
    } else {
    """
        fastp \
        -i ${reads} \
        --thread ${task.cpus} \
        -o ${meta.sample_name}_trimmed_R1.fastq.gz \
        --json ${meta.sample_name}.fastp.json
        """
    }
}


process fastqc {
    tag "Fastqc on ${meta.sample_name}"
    label "process_low"

    publishDir "${launchDir}/analysis/fastqc/", mode: "copy"

    module 'FastQC/0.11.9'

    input:
    tuple val(meta), path(reads)

    output:
    path ("*.zip"), emit: zips
    path ("*.html"), emit: htmls

    script:
    """
        fastqc --threads ${task.cpus} ${reads}
    """
}


process multiqc {
    tag "Multiqc on this project"
    label "process_dual"

    publishDir "${launchDir}/analysis/multiqc/", mode : "copy"

    module "python/3.11.1"

    input:
    path(files)

    output:
    path("*.html"), emit: multiqc_output

    script:
    config_yaml = "/home/kimj32/config_defaults.yaml"
    """
        multiqc ${files} --filename "multiqc_report.html" --config ${config_yaml}
    """
}


process fastq_screen {
    tag "Fastq-screen on ${sample_name}"
    label "process_low"

    publishDir "${launchDir}/analysis/fastq_screen", mode : "copy"

    module 'FastQScreen/0.14.1'
    module 'bowtie2/2.3.4.1'

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    path("*.html")
    path("*.txt"), emit: fastq_screen_out

    // threads option is already defined in fastq_screeN_conf
    script:
    conf = "/mnt/beegfs/kimj32/polymerase/polymeraseDependencies/FastQ_Screen_Genomes/fastq_screen.conf"
    """
        fastq_screen --aligner bowtie2 \
        --conf ${conf} \
        ${reads[0]} \
        --outdir ./  \
        --threads ${task.cpus}
    """
}


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
    // https://www.metagenomics.wiki/tools/short-read/remove-host-sequences#h.p2wn17ntls79
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


process split_reads_from_unmapped {
    tag "split reads - ${sample_name}"
    label "process_medium"

    publishDir "${launchDir}/analysis/split_reads"

    module "samtools/1.16.1"

    input:
    tuple val(sample_name), path(bam_file), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.host_remove.R{1,2}.fastq.gz"), val(is_SE),  emit: split_reads

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


process metaphlan {
    tag "metaphlan on ${sample_name}"
    label "memory_medium"

    publishDir "${launchDir}/analysis/metaphlan", mode : "copy"

    module "MetaPhlAn/4.0"
    // singularity

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    path("${sample_name}.bowtie2.bz2")
    tuple val(sample_name), path("${sample_name}.profiled_metagenome.txt"), emit: profiled_metagenome

    script:
    def bowtie2db = "/mnt/beegfs/kimj32/reference/metaphlan4/metaphlan_databases/"
    // def bowtie2db = "/mnt/beegfs/root/MetaPhlAn/" // this is from the HPC - directory is not writable
    def input_reads = is_SE ? "${reads[0]}" : "${reads[0]},${reads[1]}"
    """
        metaphlan ${input_reads} \
        --nproc ${task.cpus} \
        --input_type fastq \
        -x mpa_vOct22_CHOCOPhlAnSGB_202212 \
        --bowtie2db  ${bowtie2db} \
        -t rel_ab_w_read_stats \
        -o ${sample_name}.profiled_metagenome.txt \
        --bowtie2out ${sample_name}.bowtie2.bz2
    """
}


process concat_fq {
    // See https://forum.biobakery.org/t/humann3-paired-end-reads/862
    tag "concat ${sample_name}"
    label "process_dual"

    publishDir "${launchDir}/analysis/concat"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    tuple val(sample_name), path("${sample_name}.concat.fastq.gz"), emit: concat_reads

    script:
    if (is_SE)
    """
        cp ${reads[0]} ${sample_name}.concat.fastq.gz
    """
    else
    """
        cat ${reads[0]} ${reads[1]} > ${sample_name}.concat.fastq.gz
    """
}


process megahit {
    tag "${sample_name}"
    label "memory_medium"

    publishDir (
        path: "${launchDir}/analysis/megahit",
        pattern: "**",
        saveAS: { filename -> filename.endsWith("final.contigs.fa") ? "${sample_name}/${sample_name}.contigs.fa" : filename}
    )

    module "megahit/1.2.9"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    path("${sample_name}", type: "dir")
    tuple val(sample_name), path("${sample_name}/final.contigs.fa"), emit: final_contigs_fa
    path("${sample_name}.*.log"), emit: "log"

    script:
    def input_reads = is_SE ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
        megahit ${input_reads} \
            -o ${sample_name} \
            -t ${task.cpus} \
            --k-min 27 \
            --k-max 47 \
            1> ${sample_name}.stdout.log \
            2> ${sample_name}.stderr.log
    """
}


process humann {
    tag "humann on ${sample_name}"
    label "memory_high"

    publishDir "${launchDir}/analysis/humann", mode : "copy"

    module "HUMAnN/3.0"
    // singularity

    input:
    tuple val(sample_name), path(reads)

    output:
    path("*")

    // --nucleotide-database : /mnt/beegfs/kimj32/reference/Humann/chocophlan
    // --protein-database : /mnt/beegfs/kimj32/reference/Humann/uniref
    // /cm/shared/apps/MetaPhlAn/4.0/envs/mpa/bin/metaphlan
    script:
    def nucl_db = "/mnt/beegfs/kimj32/reference/Humann/chocophlan"
    def prot_db = "/mnt/beegfs/kimj32/reference/Humann/uniref"
    def bowtie2_db = "/cm/shared/apps/HUMAnN/3.0/lib/python3.9/site-packages/metaphlan/metaphlan_databases"
    // def cpus = task.cpus - 2
    """
        humann --input ${reads} \
        --threads ${task.cpus} \
        --memory-use maximum \
        -o ${sample_name}_humann_out \
        --metaphlan-options "--bowtie2db ${bowtie2_db} --offline"
    """
}


process kraken2{
    tag "kraken2 on ${sample_name}"
    label "process_high"

    publishDir "${launchDir}/analysis/kraken2"

    module "kraken/2.1.2"

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    path("*_kraken.txt")
    path("*_report.txt"), emit: kraken2_report

    script:
    def kraken2_db = "/mnt/beegfs/kimj32/reference/KRAKEN_DB"
    def paired_flag = is_SE ? "" : "--paired"
    def input_reads = is_SE ? "${reads[0]}" : "${reads}"
    """
        kraken2 --db ${kraken2_db} \
        --threads ${task.cpus} \
        --use-names \
        --output ${sample_name}_kraken.txt \
        --report ${sample_name}_report.txt \
        ${paired_flag} \
        ${input_reads}
     """
}


process kraken_biom {
    tag "kraken_biom"
    label "process_low"

    publishDir "${launchDir}/analysis/kraken_biom", mode: "copy"

    input:
    path(reports)

    output:
    path("biom_table.tsv"), emit: biom_table

    script:
    def kraken_biom = "~/beegfs/python_env/kraken2/bin/kraken-biom"
    """
    ${kraken_biom} ${reports} --fmt tsv -o biom_table.tsv
    """
}


process count_reads {
    tag "${sample_name} - ${stage}"
    label "process_low"

    publishDir "${launchDir}/analysis/read_counts", mode: "copy"

    input:
    tuple val(sample_name), path(reads), val(stage)

    output:
    path("${sample_name}.${stage}.count.csv"), emit: count_csv

    script:
    """
    read_count=\$(( \$(zcat ${reads[0]} | wc -l) / 4 ))
    reads_million=\$(echo \$read_count | awk '{printf "%.2f", \$1 / 1000000}')
    echo "sample,stage,reads,reads_million" > ${sample_name}.${stage}.count.csv
    echo "${sample_name},${stage},\$read_count,\$reads_million" >> ${sample_name}.${stage}.count.csv
    """
}


process merge_read_counts {
    label "process_low"

    publishDir "${launchDir}/analysis/read_counts", mode: "copy"

    input:
    path(csv_files)  // all per-sample count CSVs collected into one list

    output:
    path("read_counts.csv")

    script:
    """
    # Pivot long-format per-sample CSVs into one wide-format table (one row per sample)
    awk -F',' '
        FNR==1 { next }
        \$2=="raw"          { raw[\$1]=\$3 }
        \$2=="post_fastp"   { post_fastp[\$1]=\$3 }
        \$2=="host_removed" { host_removed[\$1]=\$3 }
        END {
            print "sample,raw,post_fastp,host_removed"
            for (s in raw) print s "," raw[s] "," post_fastp[s] "," host_removed[s]
        }
    ' ${csv_files} > read_counts.csv
    """
}


process make_taxa_counts {
    tag "${sample_name}"
    label "process_low"

    publishDir "${launchDir}/analysis/taxa_counts", mode: "copy"

    input:
    tuple val(sample_name), path(profiled_metagenome)

    output:
    path("${sample_name}.species.txt")
    path("${sample_name}.genus.txt")

    script:
    """
    grep 's__' ${profiled_metagenome} | grep -v 't__' | cut -f1,5 | sed 's/.*|//g' > ${sample_name}.species.txt
    grep 'g__' ${profiled_metagenome} | grep -v 't__' | grep -v 's__' | cut -f1,5 | sed 's/.*|//g' > ${sample_name}.genus.txt
    """
}


// https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4.1
process StrainPhlAn{
    tag "StrainPhlAn on ${sample_name}"
    label "memory_medium"

    publishDir "${launchDir}/analysis/StrainPhlAn"

    module "MetaPhlAn/4.1.1"
    module 'bowtie2/2.3.4.1'
    // singularity

    input:
    tuple val(sample_name), path(reads), val(is_SE)

    output:
    // path("*")
    tuple val(sample_name), path("${sample_name}.sam.bz2"), emit: StrainPhlAn_sam_bz2

    script:
    def bowtie2db = "/mnt/beegfs/kimj32/reference/metaphlan4/metaphlan_databases/"
    // def bowtie2db = "/mnt/beegfs/root/MetaPhlAn/" // this is from the HPC - directory is not writable
    """
        metaphlan ${reads[0]},${reads[1]} \
        --nproc ${task.cpus} \
        --input_type fastq \
        -s ${sample_name}.sam.bz2 \
        -o ${sample_name}.profiled.tsv \
        --bowtie2out ${sample_name}.bowtie2.bz2
    """
}


process concensus_markers {
    tag "concensus_markers on ${sample_name}"
    label "memory_medium"

    publishDir "${launchDir}/analysis/concensus_markers", mode : "copy"

    module "MetaPhlAn/4.0"
    module 'bowtie2/2.3.4.1'

    input:
    tuple val(sample_name), path(StrainPhlAn_sam_bz2)

    output:
    path("consensus_markers")

    script:
    """
        sample2markers.py -i ${StrainPhlAn_sam_bz2} -o consensus_markers -n ${task.cpus}
    """


}

// See https://bioinformatics.stackexchange.com/questions/20227/how-does-one-account-for-both-single-end-and-paired-end-reads-as-input-in-a-next
ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)

// adapted from https://bioinformatics.stackexchange.com/questions/20227/how-does-one-account-for-both-single-end-and-paired-end-reads-as-input-in-a-next
ch_reads = ch_samplesheet.splitCsv(header:true).map {

    // This is the read1 and read2 entry
    r1 = it['fq1']
    r2 = it['fq2']

    // Detect wiether single-end or paired-end
    is_singleEnd = r2.toString() =='' ? true : false

    // The "meta" map, which is a Nextflow/Groovy map with id (the sample name) and a single_end logical entry
    meta = [sample_name: it['sample'], single_end: is_singleEnd]

    // We return a nested map, the first entry is the meta map, the second one is the read(s)
    r2.toString()=='' ? [meta, [r1]] : [meta, [r1, r2]]

}

workflow {
    fastqc(ch_reads)
    fastp(ch_reads)
    fastq_screen(fastp.out.trim_reads)

    // host removal: bowtie2 for DNA, STAR for RNA
    if (params.host_removal == "RNA") {
        star_host_removal(fastp.out.trim_reads)
        ch_host_removed    = star_host_removal.out.split_reads
        ch_host_removal_qc = star_host_removal.out.star_log
    } else {
        bowtie2(fastp.out.trim_reads)
        split_reads_from_unmapped(bowtie2.out.bowtie2_bam_both_unmapped_bam)
        ch_host_removed    = split_reads_from_unmapped.out.split_reads
        ch_host_removal_qc = bowtie2.out.samtools_stats
    }

    // Count reads at raw, post-fastp, and host-removed stages
    ch_raw_counts     = ch_reads.map { meta, reads -> [meta.sample_name, reads, "raw"] }
    ch_fastp_counts   = fastp.out.trim_reads.map { sn, reads, se -> [sn, reads, "post_fastp"] }
    ch_hostrem_counts = ch_host_removed.map { sn, reads, se -> [sn, reads, "host_removed"] }
    count_reads( ch_raw_counts.mix(ch_fastp_counts, ch_hostrem_counts) )
    merge_read_counts( count_reads.out.count_csv.collect() )

    if (params.run_metaphlan) {
        metaphlan(ch_host_removed)
        make_taxa_counts(metaphlan.out.profiled_metagenome)
    }
    if (params.run_megahit) {
        megahit(ch_host_removed)
        multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips, fastp.out.fastp_json, ch_host_removal_qc, megahit.out.log).collect() )
    }
    if (params.run_humann) {
        concat_fq(ch_host_removed)
        humann(concat_fq.out.concat_reads)
    }
    if (params.run_kraken2) {
        kraken2(ch_host_removed)
        kraken_biom(kraken2.out.kraken2_report.collect())
        multiqc( fastq_screen.out.fastq_screen_out.mix(fastqc.out.zips, fastp.out.fastp_json, ch_host_removal_qc).collect() )
    }
    if (params.run_StrainPhlAn) {
        StrainPhlAn(ch_host_removed)
    }
}

workflow.onComplete {
    println ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}