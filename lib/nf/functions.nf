def saveFiles(Map args) {
    /* Modeled after nf-core/modules saveFiles function */
    def final_output = null
    def filename = ""
    def found_ignore = false
    def logs_subdir = args.containsKey('logs_subdir') ? args.logs_subdir : args.opts.logs_subdir
    def process_name = args.opts.process_name
    def publish_to_base = args.opts.publish_to_base.getClass() == Boolean ? args.opts.publish_to_base : false
    def publish_to_base_list = args.opts.publish_to_base.getClass() == ArrayList ? args.opts.publish_to_base : []
    if (args.filename) {
        if (args.filename.startsWith('.command')) {
            // Its a Nextflow process file, rename to "nf-<PROCESS_NAME>.*"
            ext = args.filename.replace(".command.", "")
            final_output = "logs/${process_name}/${logs_subdir}/nf-${process_name}.${ext}"
        } else if (args.filename.endsWith('.stderr.txt') || args.filename.endsWith('.stdout.txt') || args.filename.endsWith('.log')  || args.filename.endsWith('.err') || args.filename.equals('versions.yml')) {
            // Its a version file or  program specific log files
            final_output = "logs/${process_name}/${logs_subdir}/${args.filename}"
        } else {
            // Its a program output
            filename = args.filename
            if (filename.startsWith("results/")) {
                filename = filename.replace("results/","")
            }

            // *-error.txt should be at the base dir and 'blastdb' should go in blast folder
            if (filename.endsWith("-error.txt") || filename.endsWith("-genome-size.txt") || publish_to_base == true) {
                final_output = filename
            } else if (filename.startsWith("blastdb/")) {
                final_output = "blast/${filename}"
            } else if (filename.startsWith("total_contigs_")) {
                final_output = null
            } else if (params.publish_dir.containsKey(process_name)) {
                final_output = "${params.publish_dir[process_name]}/${filename}"
                if (final_output.startsWith("/")) {
                    final_output = filename
                }
            } else {
                if (args.opts.is_module) {
                    final_output = filename
                } else {
                    final_output = "${process_name}/${filename}"
                }
            }

            // Exclude files that should be ignored
            args.opts.ignore.each {
                if (filename.endsWith("${it}")) {
                    final_output = null
                }
            }

            // Publish specific files to base
            publish_to_base_list.each {
                if (filename.endsWith("${it}")) {
                    final_output = filename
                }
            }
        }
        return final_output == null ? null : final_output.replace("//", "/")
    }
}

def initOptions(Map args, String process_name) {
    /* Function to initialise default values and to generate a Groovy Map of available options for nf-core modules */
    def Map options = [:]
    options.args            = args.args ?: ''
    options.ignore          = args.ignore ?: []
    options.is_module       = args.is_module ?: false
    options.logs_subdir     = args.logs_subdir ?: ''
    options.process_name    = args.process_name ?: process_name
    options.publish_to_base = args.publish_to_base ?: false
    options.suffix          = args.suffix ?: ''

    return options
}