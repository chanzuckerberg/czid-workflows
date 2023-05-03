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
    call make_taxadb {
        input: 
            docker_image_id
    }
    
    call preprocess_taxa as preprocess_taxa_nt{
        input: 
            taxon_counts = taxon_counts_run_1,
            contig_summary = contig_summary_run_1,
            contig_fasta = contig_fasta_run_1,
            taxadb_sqlite = make_taxadb.taxadb_sqlite,
            dbtype = "NT",
            docker_image_id
    }
    call preprocess_taxa as preprocess_taxa_nr{
        input: 
            taxon_counts = taxon_counts_run_1,
            contig_summary = contig_summary_run_1,
            contig_fasta = contig_fasta_run_1,
            taxadb_sqlite = make_taxadb.taxadb_sqlite,
            dbtype = "NR",
            docker_image_id
    }

    if (defined(ground_truth)) { 

        call truth_benchmark as truth_benchmark_nt {
            input: 
                preprocessed_taxa = preprocess_taxa_nt.preprocessed_taxa,
                ground_truth = select_first([ground_truth]),
                dbtype = "NT",
                docker_image_id
        }
        call truth_benchmark as truth_benchmark_nr {
            input: 
                preprocessed_taxa = preprocess_taxa_nr.preprocessed_taxa,
                ground_truth = select_first([ground_truth]),
                dbtype = "NR",
                docker_image_id
        }
    }
    call notebook as test_notebook { 
            input: 
                preprocessed_taxa_nt = preprocess_taxa_nt.preprocessed_taxa,
                preprocessed_taxa_nr = preprocess_taxa_nr.preprocessed_taxa,
                docker_image_id
        }
    if (defined(taxon_counts_run_2) || defined(contig_summary_run_2) || defined(contig_fasta_run_2)) {
        call notebook { 
            input: 
                preprocessed_taxa_nt = preprocess_taxa_nt.preprocessed_taxa,
                preprocessed_taxa_nr = preprocess_taxa_nr.preprocessed_taxa,
                docker_image_id
        }

    }

    output {
        File preprocessed_nt = preprocess_taxa_nt.preprocessed_taxa
        File preprocessed_nr = preprocess_taxa_nr.preprocessed_taxa
        File benchmark_notebook = test_notebook.benchmark_notebook
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

task preprocess_taxa { 
    input{
        File taxon_counts
        File contig_summary
        File contig_fasta 
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

        tc = harvest.harvest_sample_taxon_counts(
            "~{taxon_counts}", 
            "~{contig_summary}", 
            "~{contig_fasta}", 
            "~{dbtype}", 
            taxadb
            )

        with open("preprocessed_taxa.json", "w") as f:
            json.dump(tc, f)
        
        CODE
    >>>
    output { 
        File preprocessed_taxa = "preprocessed_taxa.json"
    }
    runtime {
        docker: docker_image_id
    }
}

task notebook { 
    input{
        File preprocessed_taxa_nt
        File preprocessed_taxa_nr 
        String docker_image_id
    }
    command <<<
    jq -n '[(input | {NT: .}), (input | {NR: .})] | reduce .[] as $item ({}; . * $item)' \
        "~{preprocessed_taxa_nt}" "~{preprocessed_taxa_nr}" > combined_taxa.json
    
    export HARVEST_DATA="combined_taxa.json"
    cp /home/jovyan/notebooks/short-read-mngs-benchmarks.ipynb short-read-mngs-benchmarks.ipynb
    jupyter nbconvert --to notebook --execute --inplace short-read-mngs-benchmarks.ipynb

    >>>
    output { 
        File combined = "combined_taxa.json"
        File benchmark_notebook = "short-read-mngs-benchmarks.ipynb"
    }
    runtime { 
        docker: docker_image_id 
    }
}

task truth_benchmark { 
    input{
        File preprocessed_taxa
        File ground_truth
        String dbtype
        String docker_image_id 
    }
    command <<<
        set -euxo pipefail
        python3 <<CODE

        from benchmark_helpers import harvest, metrics
        from taxadb.taxid import TaxID
        import json

        with open("~{preprocessed_taxa}", "r") as f:
            tc = json.load(f) 

        truth = harvest.read_truth_file("~{ground_truth}")

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
