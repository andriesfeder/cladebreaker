//
// Check input samplesheet and get read channels
//

include { PROKKA                      } from '../../modules/nf-core/modules/prokka/main'
include { NCBIGENOMEDOWNLOAD          } from '../../modules/nf-core/modules/ncbigenomedownload/main'

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
    prokka_gff = Channel.empty()
    if( params.ref == null ){
        prokka_ncbi =  Channel.empty()
        prokka_ncbi = prokka_ncbi.mix(NCBIGENOMEDOWNLOAD.out.fna)
        prokka_ncbi = prokka_ncbi.transpose()
        prokka_ncbi = prokka_ncbi.combine(Channel.fromPath( params.proteins )).combine(Channel.fromPath( params.prodigal_tf ))

        PROKKA (
            prokka_ncbi
        )
        ch_versions = ch_versions.mix(PROKKA.out.versions.first())
        prokka_gff = prokka_emit.mix(PROKKA.out.gff)
    }

    emit:
    ncbi = NCBIGENOMEDOWNLOAD.out.fna
    prokka_gff           // channel: [ val(meta), [ assemblies] ]
    versions = ch_versions             // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]

import java.nio.file.Files
import java.nio.file.Paths
import java.io.FileWriter
import java.nio.channels.FileChannel


def create_gca_channels( String gca ) {
    gca = gca.trim()
    def meta = [:]
    meta.id           = gca
    meta.single_end   = false
    meta.assembly     = true

    Path filePath = Paths.get("${workflow.workDir}/tmp/${gca}.txt")
    FileWriter fr = new FileWriter("${workflow.workDir}/tmp/${gca}.txt")
    fr.write(gca)
    fr.close()

    array = [ meta, file(filePath) ]

    return array
}
def gather_gff(String t ) {
    def meta = [:]
    meta.id           = "ALL_SAMPLES"

    new File("${workflow.workDir}/tmp/gff/").mkdirs()
    Path filePath = Paths.get("${workflow.workDir}/tmp/gff/")

    File dir = new File("${params.outdir}")
    File[] listing = dir.listFiles()
    for (File sample : listing) {
        if (sample.isDirectory()){
            print(sample)
        }
        File annotation_path = new File("${sample}/annotation/*.gff")
        if (annotation_path.exists()) {
            String path = "${workflow.workDir}/tmp/gff/"
            Files.copy(annotation_path.toPath(), (new File(path + annotation_path.getName())).toPath())
        }
    }
    array = [meta, file(filePath)]

    return array

}


