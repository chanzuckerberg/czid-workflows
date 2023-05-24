version 1.1

import './short-read-mngs-benchmark.wdl' as sr

workflow benchmark {
    input {
        # mngs inputs
        File? taxon_counts_run_1
        File? contig_summary_run_1
        File? contig_fasta_run_1

        Array[File]? step_counts_run_1 

        File? taxon_counts_run_2
        File? contig_summary_run_2
        File? contig_fasta_run_2

        Array[File]? step_counts_run_2

        File? ground_truth

        String workflow_type 
        String docker_image_id
    }

    if(workflow_type == "short-read-mngs") { 
        ## Coercing inputs to non-optional values. Will hopefully cause failures early 
        call sr.short_read_mngs_benchmark { 
            input: 
                taxon_counts_run_1 = select_first([taxon_counts_run_1]),
                contig_summary_run_1 = select_first([contig_summary_run_1]),
                contig_fasta_run_1 = select_first([contig_fasta_run_1]),
                taxon_counts_run_2,
                contig_summary_run_2,
                contig_fasta_run_2,
                step_counts_run_1,
                step_counts_run_2,
                ground_truth = ground_truth,
                docker_image_id 
        }
    }
}
