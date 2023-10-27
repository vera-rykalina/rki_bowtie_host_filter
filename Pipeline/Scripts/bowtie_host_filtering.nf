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
params.host_db = "${rojectDir}/References/GRCh38_noalt_as/"


// Get software versions
process SoftwareVersions {
    publishDir "${params.outdir}/01_software_versions", mode: 'copy'

    output:
        path "software_versions.txt"

    script:
    """
    echo "software\tversion\tbuild\tchannel" > tempfile
    
    conda list | \
    grep "bowtie2\\|samtools\\|" >> tempfile

    # replace blanks with tab for easier processing downstream
    tr -s '[:blank:]' '\t' < tempfile > software_versions.txt
    """
}


/************************** 
---------WORKFLOW----------
**************************/
ch_software_versions = SoftwareVersions()
//ch_infiles = channel.fromFilePairs("${projectDir}/RawData/*_R{1,2}*.fastq.gz")

workflow {

    //ch_classified = CLASSIFY(ch_infiles, params.krakendb)
    //ch_extracted = EXTRACT(ch_classified.classified_fastq, ch_classified.kraken2_output)
    //ch_merged_ags = ch_classified.unclassified_fastq.combine(ch_extracted.filtered, by:0)
    //ch_merged_compressed = MERGE(ch_merged_ags)
    

}

/************************** 
PROCESSES
**************************/

// bowtie2
process MAP {

    label "bowtie"
    conda "${projectDir}/Environments/bowtie2.yml"
    publishDir "${params.outdir}/01_mapped_reads/${id}", mode: "copy", overwrite: true

    input:
        tuple val(id), path(reads)
  
    output:
        tuple val(id), path("${id}_classified_R*.fastq"),     emit: classified_fastq
        tuple val(id), path("${id}_unclassified_R*.fastq"),   emit: unclassified_fastq
        tuple val(id), path("${id}_kraken2_out.txt"),         emit: kraken2_output
        tuple val(id), path("${id}_kraken2_report.txt"),      emit: kraken2_report

    script:
        """
        bowtie2 \
        -p 8 \
        -x ${host_db} \
        -1 ${read[0]} \
        -2 ${read[1]} \
        -S ${id}_mapped_and_unmapped.sam"
        
        """
     
}

// krakenTools
process EXTRACT {
    label "krakentools"
    conda "${projectDir}/Environments/krakentools.yml"
    publishDir "${params.outdir}/02_homo_filtered_reads", mode: "copy", overwrite: true

    input:
        tuple val(id), path(reads)
        tuple val(id), path(kraken2_output)
    
    output:
        tuple val(id), path("${id}_filtered_R*.fastq"), emit: filtered
    
    script:
        """
            extract_kraken_reads.py \
                -k ${kraken2_output} \
                --taxid 9606 \
                --exclude \
                -s1 ${reads[0]} \
                -s2 ${reads[1]} \
                -o ${id}_filtered_R1.fastq \
                -o2 ${id}_filtered_R2.fastq \
                --fastq-output   
        """
}

// merge unclassified and filtered reads
process MERGE {
    publishDir "${params.outdir}/03_merged_reads", failOnError: true, mode: "copy", overwrite: true

    input:
        tuple val(id), path(unclassified), path(filtered)

    output:
        tuple val("${id}"), path("${id}_merged_R*.fastq.gz"), emit: merged_fastq_gz

    script:
        """
        gzip -c ${unclassified[0]} ${filtered[0]} > ${id}_merged_R1.fastq.gz
        gzip -c ${unclassified[1]} ${filtered[1]} > ${id}_merged_R2.fastq.gz
        """
}
