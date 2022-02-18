include { initOptions; saveFiles } from '../../../lib/nf/functions'
options = initOptions(params.containsKey('options') ? params.options : [:], 'gather_samples')

process GATHER_SAMPLES {
    tag "${meta.id}"
    label "process_low"

    // TODO: publishDir needs to be changed to accomdate params.publish_dir_mode and params.force
    publishDir "${params.outdir}/${meta.id}", mode: "copy", overwrite: false 
       // saveAs: { filename -> saveFiles(filename:filename, opts:options) }
    
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("fastqs/${meta.id}*.fastq.gz"), path("extra/*.gz"), path("${meta.id}-genome-size.txt"), emit: raw_fastq, optional: true
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions
    path "*-{error,merged}.txt", optional: true

    shell:
    fcov = params.coverage.toInteger() == 0 ? 150 : Math.round(params.coverage.toInteger() * 1.5)
    r1 = reads[0]
    r2 = reads[1]

    '''
    
    mkdir -p fastqs
    mkdir -p extra

    cp -L !{r1} fastqs/!{meta.id}_R1.fastq.gz
    cp -L !{r2} fastqs/!{meta.id}_R2.fastq.gz
    touch extra/empty.fna.gz

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        art: $(echo $(art_illumina --help 2>&1) | sed 's/^.*Version //;s/ .*$//')
        fastq-dl: $(echo $(fastq-dl --version 2>&1) | sed 's/fastq-dl //')
        fastq-scan: $(echo $(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        mash: $(echo $(mash --version 2>&1))
        ncbi-genome-download: $(echo $(ncbi-genome-download --version 2>&1))
        pigz: $(echo $(pigz --version 2>&1) | sed 's/pigz //')
    END_VERSIONS
    '''

    
}