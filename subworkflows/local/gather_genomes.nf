//
// Check input samplesheet and get read channels
//

include { PROKKA             } from '../../modules/nf-core/modules/prokka/main'
include { BAKTA              } from '../../modules/local/bakta/main'
include { NCBIGENOMEDOWNLOAD } from '../../modules/nf-core/modules/ncbigenomedownload/main'

workflow GATHER_GENOMES {

    take:
    ncbi_genomes // path(accessions)

    main:
    ch_versions = Channel.empty()
    ch_paths = Channel.empty()
    ch_paths = ncbi_genomes.splitText()

    ch_paths
        .map { create_gca_channels( it ) }
        .unique()
        .set { gca }

    NCBIGENOMEDOWNLOAD (
        gca
    )
    ch_versions = ch_versions.mix(NCBIGENOMEDOWNLOAD.out.versions.first())
    annotation_gff = Channel.empty()
    if( params.ref == null ){
        ncbi_fna = NCBIGENOMEDOWNLOAD.out.fna.transpose()

        if ( params.annotator == 'bakta' ) {
            BAKTA (
                ncbi_fna,
                file(params.bakta_db)
            )
            ch_versions    = ch_versions.mix(BAKTA.out.versions.first())
            annotation_gff = annotation_gff.mix(BAKTA.out.gff)
        } else {
            prokka_ncbi = ncbi_fna
                .combine(Channel.fromPath( params.proteins ))
                .combine(Channel.fromPath( params.prodigal_tf ))
            PROKKA (
                prokka_ncbi
            )
            ch_versions    = ch_versions.mix(PROKKA.out.versions.first())
            annotation_gff = annotation_gff.mix(PROKKA.out.gff)
        }
    }

    emit:
    ncbi           = NCBIGENOMEDOWNLOAD.out.fna
    annotation_gff                      // channel: [ val(meta), path(gff) ]
    versions       = ch_versions        // channel: [ versions.yml ]
}

def create_gca_channels( String gca ) {
    gca = gca.trim()
    def meta = [:]
    meta.id         = gca
    meta.single_end = false
    meta.assembly   = true

    def filePath = java.nio.file.Paths.get("${workflow.workDir}/tmp/${gca}.txt")
    def fr = new java.io.FileWriter("${workflow.workDir}/tmp/${gca}.txt")
    fr.write(gca)
    fr.close()

    def array = [ meta, file(filePath) ]
    return array
}

def gather_gff(String t) {
    def meta = [:]
    meta.id = "ALL_SAMPLES"

    new File("${workflow.workDir}/tmp/gff/").mkdirs()
    def filePath = java.nio.file.Paths.get("${workflow.workDir}/tmp/gff/")

    def dir = new File("${params.outdir}")
    def listing = dir.listFiles()
    listing.each { sample ->
        if (sample.isDirectory()) {
            print(sample)
        }
        def annotation_path = new File("${sample}/annotation/*.gff")
        if (annotation_path.exists()) {
            def path = "${workflow.workDir}/tmp/gff/"
            java.nio.file.Files.copy(annotation_path.toPath(), (new File(path + annotation_path.getName())).toPath())
        }
    }
    def array = [meta, file(filePath)]
    return array
}
