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
// -with-report logs/test.$(date +%T).report.html \
// # use this parameter for an empty test run
// -stub

/************************** 
DEFINE VARIABLES
**************************/

projectDir = "/home/rykalinav/scratch/rki_bowtie_host_filter/Pipeline"

// Parameters for host db
params.index_prefix = "GRCh38_noalt_as"


/************************** 
---------WORKFLOW----------
**************************/
ch_fastq = Channel
        .fromFilePairs("${projectDir}/RawData/*_R{1,2}*.fastq.gz", checkIfExists: true)
        .ifEmpty { exit 1, "Fastq files are not found: " }

ch_indeces = Channel
        .fromPath("${projectDir}/References/GRCh38_noalt_as/*.bt2")
        .ifEmpty { exit 1, "Bowtie2 index directory not found: " }


workflow {
    ch_mapped_reads = MAP(params.index_prefix, ch_fastq, ch_indeces.collect())
    ch_unmapped_reads = FILTER(ch_mapped_reads)
}

/************************** 
PROCESSES
**************************/

// bowtie2
process MAP {
    label "bowtie_samtools"
    conda "${projectDir}/Environments/bowtie_samtools.yml"
    publishDir "${params.outdir}/01_all_mapped_reads", mode: "copy", overwrite: true

    input:
        val (prefix)
        tuple val(id), path(reads)
        path (indeces)

    output:
        tuple val(id), path("${id}_mapped_and_unmapped.bam"),  emit: all_reads_bam

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
          -2 ${reads[1]} | samtools view \
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
        tuple val(id), path("${id}_both_reads_unmapped.bam"),  emit: all_reads_bam

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

