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
    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    assemblies                                // channel: [ val(meta), [ assemblies] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

def force_imp() {
    if (params.force && file(params.outdir).exists()) {
        file("${params.outdir}").deleteDir()
    }
}

def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample.replace("_T1","")
    meta.single_end = row.single_end.toBoolean()
    meta.assembly   = row.assembly.toBoolean()

    def sample = row.sample.toString().replace("_T1","")
    def outdir = "${params.outdir}/${sample}"

    def dir = new File("${outdir}/input_data")
    dir.mkdirs()

    def array = []
    if (!file(row.file_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.file_1}"
    }
    if (meta.single_end) {
        def source = java.nio.file.Paths.get(row.file_1)
        def target = java.nio.file.Paths.get("${params.outdir}/${sample}/input_data/${sample}.fastq.gz")
        java.nio.file.Files.copy(source, target)
        array = [ meta, [ file(target) ]]
    } else if (!meta.single_end && !meta.assembly) {
        if (!file(row.file_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.file_2}"
        }
        def source1 = java.nio.file.Paths.get(row.file_1)
        def source2 = java.nio.file.Paths.get(row.file_2)
        def target1 = java.nio.file.Paths.get("${params.outdir}/${sample}/input_data/${sample}_R1.fastq.gz")
        def target2 = java.nio.file.Paths.get("${params.outdir}/${sample}/input_data/${sample}_R2.fastq.gz")
        java.nio.file.Files.copy(source1, target1)
        java.nio.file.Files.copy(source2, target2)
        array = [ meta, [ file(target1), file(target2) ]]
    } else {
        return null
    }

    return array
}

def create_assembly_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample.replace("_T1","")
    meta.single_end = row.single_end.toBoolean()
    meta.assembly   = row.assembly.toBoolean()

    def sample = row.sample.toString().replace("_T1","")
    def outdir = "${params.outdir}/${sample}"

    def dir = new File("${outdir}/input_data")
    dir.mkdirs()

    def array = []
    if (!file(row.file_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> assembly file does not exist!\n${row.file_1}"
    }
    if (meta.assembly) {
        def source = java.nio.file.Paths.get(row.file_1)
        def target = java.nio.file.Paths.get("${params.outdir}/${sample}/input_data/${sample}.fna")
        java.nio.file.Files.copy(source, target)
        array = [ meta, file(target) ]
    } else {
        return null
    }

    return array
}
