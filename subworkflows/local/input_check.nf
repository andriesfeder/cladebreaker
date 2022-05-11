//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    force_imp()
    valid_sheet = SAMPLESHEET_CHECK ( samplesheet )
    valid_sheet
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set{ reads }
    valid_sheet
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_assembly_channels(it) }
        .set{ assemblies }
    // TODO: remove these
    // reads.view()
    // assemblies.view()
    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    assemblies                                // channel: [ val(meta), [ assemblies] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]

import java.nio.file.Files
import java.nio.file.Paths

def force_imp() {
    if (params.force && file(params.outdir).exists()) {
            file("${params.outdir}").deleteDir()
        }
}

def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample.replace("_T1","")
    meta.single_end   = row.single_end.toBoolean()
    meta.assembly     = row.assembly.toBoolean()

    //if (params.force && file(params.outdir).exists()) {
    //        file("${params.outdir}/${meta.id}/input_data").deleteDir()
    //    }

    String sample = row.sample.toString().replace("_T1","")
    String outdir = "${params.outdir}/${sample}"

    final File dir = new File("${outdir}/input_data")
    dir.mkdirs()

    def array = []
    def array_assembly = []
    if (!file(row.file_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.file_1}"
    }
    if (meta.single_end) {

        Path source = Paths.get(row.file_1)
        Path target = Paths.get("${params.outdir}/${sample}/input_data/${sample}.fastq.gz")

        Files.copy(source, target)

        array = [ meta, [ file(target) ]]
    } else if (!meta.single_end && !meta.assembly){
        if (!file(row.file_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.file_2}"
        }

        Path source1 = Paths.get(row.file_1)
        Path source2 = Paths.get(row.file_2)
        Path target1 = Paths.get("${params.outdir}/${sample}/input_data/${sample}_R1.fastq.gz")
        Path target2 = Paths.get("${params.outdir}/${sample}/input_data/${sample}_R2.fastq.gz")

        Files.copy(source1, target1)
        Files.copy(source2, target2)

        array = [ meta, [ file(target1), file(target2) ]]
    } else {

        return null
    }

    return array
}

def create_assembly_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample.replace("_T1","")
    meta.single_end   = row.single_end.toBoolean()
    meta.assembly     = row.assembly.toBoolean()

    // if (params.force && file(params.outdir).exists()) {
    //        file("${params.outdir}/${meta.id}/input_data").deleteDir()
    //    }

    String sample = row.sample.toString().replace("_T1","")
    String outdir = "${params.outdir}/${sample}"

    final File dir = new File("${outdir}/input_data")
    dir.mkdirs()

    def array = []
    if (!file(row.file_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> assembly file does not exist!\n${row.file_1}"
    }
    if (meta.assembly) {

        Path source = Paths.get(row.file_1)
        Path target = Paths.get("${params.outdir}/${sample}/input_data/${sample}.fna")

        Files.copy(source,target)

        array = [meta,  file(target) ]
    } else {

        return null
    }

    return array
}
