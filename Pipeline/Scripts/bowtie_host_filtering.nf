nextflow.enable.dsl = 2

params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory!"
}


/************************** 
-----------CALL------------
**************************/

// $ nextflow run ../scratch/rki_bowtie_host_filter/Pipeline/Scripts/rki_bowtie_host_filtering.nf \
// -profile rki_slurm,rki_mamba \
// -c Scripts/rki_profile.config \
// -with-report logs/test.$(date +%T).report.html \


/************************** 
DEFINE VARIABLES
**************************/

projectDir = "/home/rykalinav/scratch/rki_bowtie_host_filter/Pipeline"

// Parameters for host db
params.index_prefix = "GRCh38_noalt_as"


/************************** 
---------WORKFLOW----------
**************************/
ch_fastq = Channel.fromFilePairs("${projectDir}/RawData/*_R{1,2}*.fastq.gz", checkIfExists: true)
ch_indeces = Channel.fromPath("${projectDir}/References/GRCh38_noalt_as/*.bt2", checkIfExists: true)
        


workflow {
    ch_mapped_reads = MAP(params.index_prefix, ch_fastq, ch_indeces.collect())
    ch_unmapped_reads = FILTER(ch_mapped_reads)
    ch_sorted_unmapped_reads = SORT(ch_unmapped_reads)
    ch_host_removed_reads = SPLIT(ch_sorted_unmapped_reads)
    ch_marked_header_reads = MARK_HEADER(ch_host_removed_reads)
}

/************************** 
PROCESSES
**************************/

// bowtie2
process MAP {
    label "bowtie_samtools"
    conda "${projectDir}/Environments/bowtie2_v2.5.2.yml"
    publishDir "${params.outdir}/01_all_mapped_reads", mode: "copy", overwrite: true

    input:
        val (prefix)
        tuple val(id), path(reads)
        path (indeces)

    output:
        tuple val(id), path("${id}_mapped_and_unmapped.bam")

    /* 
    bowtie2 mapping against host sequence database, 
    keeping both aligned and unaligned reads (paired-end reads)
    */

    script:
        """
        bowtie2 \
          -p ${task.cpus} \
          -x ${prefix} \
          -1 ${reads[0]} \
          -2 ${reads[1]} |\
          
        samtools view \
          -b \
          -S \
          - \
          -o ${id}_mapped_and_unmapped.bam
        """    
}

// samtools
process FILTER {
    conda "${projectDir}/Environments/bowtie_samtools.yml"
    publishDir "${params.outdir}/02_all_unmapped_reads", mode: "copy", overwrite: true

    input:
        tuple val(id), path(mapped_reads)

    output:
        tuple val(id), path("${id}_both_reads_unmapped.bam"),  emit: all_unmapped

    /* 
    Get unmapped pairs (both reads R1 and R2 unmapped) 
    Extract only (-f 12) alignments with both reads unmapped: <read unmapped><mate unmapped>
    Do not (-F 256) extract alignments which are: <not primary alignment>
    */

    script:
        """
        samtools view \
          -b \
          -f 12 \
          -F 256 \
          ${mapped_reads} \
          -o ${id}_both_reads_unmapped.bam
        """    
}

// samtools
process SORT {
    label "samtools_sort"
    conda "${projectDir}/Environments/bowtie2_v2.5.2.yml"
    publishDir "${params.outdir}/03_sorted_unmapped_reads", mode: "copy", overwrite: true

    input:
        tuple val(id), path(unmapped_reads)

    output:
        tuple val(id), path("${id}_unmapped_sorted.bam")

    /* 
    Sort bam file by read name (-n) to have paired reads next to each other
    */

    script:
        """
        samtools sort \
          -n \
          -@ 2 \
          ${unmapped_reads} \
          -o ${id}_unmapped_sorted.bam
        """    
}

// samtools
process SPLIT {
    label "samtools_fastq"
    conda "${projectDir}/Environments/bowtie2_v2.5.2.yml"
    publishDir "${params.outdir}/04_host_removed_reads", mode: "copy", overwrite: true

    input:
        tuple val(id), path(sorted_reads)

    output:
        tuple val(id), path("${id}_host_removed_R{1,2}.fastq.gz")

    /* 
    Split paired-end reads into separated fastq files id_{R1,R2}.fastq.gz
    */

    script:
        """
        samtools fastq \
           -@ 8 \
           -n \
           ${sorted_reads} \
           -1 ${id}_host_removed_R1.fastq.gz \
           -2 ${id}_host_removed_R2.fastq.gz        
        """    
}

// awk
process MARK_HEADER {
    publishDir "${params.outdir}/05_header_marked_reads", mode: "copy", overwrite: true

    input:
        tuple val(id), path(host_removed_reads)

    output:
        tuple val(id), path("${id}_R{1,2}.fastq.gz")

    /* 
    Add 1 and 2 to the end of fastq read header correspondingly (needed for shiver)
    */

    script:
        """
        zcat ${host_removed_reads[0]} |\
        awk '{if (NR%4 == 1) {print \$1 " 1"} else print}' |\
        gzip -c > ${id}_R1.fastq.gz
        rm ${host_removed_reads[0]}

        zcat ${host_removed_reads[1]} |\
        awk '{if (NR%4 == 1) {print \$1 " 2"} else print}' |\
        gzip -c > ${id}_R2.fastq.gz
        rm ${host_removed_reads[1]}   
       """    
}