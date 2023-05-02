version 1.1


workflow short_read_mngs_benchmark  {
    input { 
        File taxon_counts_run_1
        File contig_summary_run_1
        File contig_fasta_run_1

        File? taxon_counts_run_2
        File? contig_summary_run_2
        File? contig_fasta_run_2
        
        # ground_truth is a tsv with taxid, read count, 
        # proportion, order, and scientific name
        # compared to both run_1 and run_2
        File? ground_truth 

        String docker_image_id
    }

    if (defined(ground_truth)) { 
        call make_taxadb {
            input: 
                docker_image_id
        }
        call truth_benchmark as truth_benchmark_nt {
            input: 
                taxon_counts = taxon_counts_run_1,
                contig_summary = contig_summary_run_1,
                contig_fasta = contig_fasta_run_1,
                ground_truth = select_first([ground_truth]),
                taxadb_sqlite = make_taxadb.taxadb_sqlite,
                dbtype = "NT",
                docker_image_id
        }
        call truth_benchmark as truth_benchmark_nr {
            input: 
                taxon_counts = taxon_counts_run_1,
                contig_summary = contig_summary_run_1,
                contig_fasta = contig_fasta_run_1,
                ground_truth = select_first([ground_truth]),
                taxadb_sqlite = make_taxadb.taxadb_sqlite,
                dbtype = "NR",
                docker_image_id
        }
    }

    output {
        File? truth_nt = truth_benchmark_nt.ground_truth_output
        File? truth_nr = truth_benchmark_nr.ground_truth_output
    }
}

task make_taxadb { 
    input {
        String docker_image_id
    }
    command <<<
    taxadb download -o taxadb --type taxa
    taxadb create -i taxadb --dbname taxadb.sqlite || true
    >>>

    output {
        File? taxadb_sqlite = "taxadb.sqlite"
    }
    runtime {
        docker: docker_image_id
    }
}

task truth_benchmark { 
    input{
        File taxon_counts
        File contig_summary
        File contig_fasta 
        File ground_truth
        File? taxadb_sqlite
        String dbtype
        String docker_image_id 
    }
    command <<<
        set -euxo pipefail
        python3 <<CODE

        from benchmark_helpers import harvest, metrics
        from taxadb.taxid import TaxID
        import json

        taxadb = TaxID(dbtype="sqlite", dbname="~{taxadb_sqlite}")
        truth = harvest.read_truth_file("~{ground_truth}")

        tc = harvest.harvest_sample_taxon_counts(
            "~{taxon_counts}", 
            "~{contig_summary}", 
            "~{contig_fasta}", 
            "~{dbtype}", 
            taxadb
            )

        # TODO: figure out how to determine if data is paired or not
        total_reads = sum([i["reads_dedup"] for i in tc.values()])*2 #TODO: review
        relative_abundance_nt = {k: v["reads_dedup"]/total_reads for k, v in tc.items()}
        

        
        ground_truth = {
            "aupr": metrics.truth_aupr(relative_abundance_nt, truth),
            "l2_norm": metrics.truth_l2_norm(relative_abundance_nt, truth)
        }
        with open("ground_truth_~{dbtype}.json", "w") as f:
            json.dump(ground_truth, f)
        CODE
    >>>
    output { 
        File ground_truth_output = "ground_truth_~{dbtype}.json"
    }
    runtime {
        docker: docker_image_id
    }
}
