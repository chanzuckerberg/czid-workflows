version 1.1


workflow benchmark {
    input {
        File? taxon_counts
        File? contig_summary
        File? taxadb 
        String workflow_type 
    }

    if(workflow_type == "short-read-mngs") { 
        ## Coercing inputs to non-optional values. Will hopefully cause failures early 
        call short_read_mngs { 
            input: 
                taxon_counts = select_first([taxon_counts]),
                contig_summary = select_first([contig_summary]),
                taxadb = select_first([taxadb])
        }

    }
}
task short_read_mngs {
    input {
        File taxon_counts
        File contig_summary
        File taxadb
    }

    command <<<
        set -euxo pipefail
        echo "Hello world"
    >>>
}