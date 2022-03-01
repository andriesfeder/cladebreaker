//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]

import java.nio.file.Files
import java.nio.file.Paths

def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample.replace("_T1","")
    meta.single_end   = row.single_end.toBoolean()

    if (params.force && file(params.outdir).exists()) {
            file("${params.outdir}/${meta.id}/fastqs").deleteDir()
        }

    String sample = row.sample.toString().replace("_T1","")
    String outdir = "${params.outdir}/${sample}"
    final File dir = new File("${outdir}/fastqs")
    dir.mkdirs()

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {

        Path source = Paths.get(row.fastq_1)
        Path target = Paths.get("${params.outdir}/${sample}/fastqs/${sample}.fastq.gz")

        Files.copy(source, target)

        array = [ meta, [ file(target) ], file(outdir) ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }

        Path source1 = Paths.get(row.fastq_1)
        Path source2 = Paths.get(row.fastq_2)
        Path target1 = Paths.get("${params.outdir}/${sample}/fastqs/${sample}_R1.fastq.gz")
        Path target2 = Paths.get("${params.outdir}/${sample}/fastqs/${sample}_R2.fastq.gz")

        Files.copy(source1, target1)
        Files.copy(source2, target2)

        array = [ meta, [ file(target1), file(target2) ], file(outdir) ]
    }
    return array
}
