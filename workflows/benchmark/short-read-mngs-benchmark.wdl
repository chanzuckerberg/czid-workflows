version 1.1
import "./mngs-benchmark-tasks.wdl" as mbt


workflow short_read_mngs_benchmark  {
    # TODO: fix coloration potentially
    # TODO: output html file  
    input { 
        File taxon_counts_run_1
        File contig_summary_run_1
        File contig_fasta_run_1

        Array[File]? step_counts_run_1 

        File? taxon_counts_run_2
        File? contig_summary_run_2
        File? contig_fasta_run_2

        Array[File]? step_counts_run_2
        
        # ground_truth is a tsv with taxid, read count, 
        # proportion, order, and scientific name
        # compared to both run_1 and run_2
        File? ground_truth 

        String docker_image_id
    }
    call mbt.make_taxadb as make_taxadb {
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

    ## Compare submitted file to a ground truth tsv file
    ## returns the AUPR and L2 norm 
    ## TODO: actualy do something with the L2_norm

    if (defined(ground_truth)) { 

        call mbt.truth_benchmark as truth_benchmark_nt {
            input: 
                preprocessed_taxa = preprocess_taxa_nt.preprocessed_taxa,
                ground_truth = select_first([ground_truth]),
                dbtype = "NT",
                docker_image_id
        }
        call mbt.truth_benchmark as truth_benchmark_nr {
            input: 
                preprocessed_taxa = preprocess_taxa_nr.preprocessed_taxa,
                ground_truth = select_first([ground_truth]),
                dbtype = "NR",
                docker_image_id
        }
    }

    ## Read in step count files 

    if (defined(step_counts_run_1) && defined(step_counts_run_2)) {
        call mbt.read_step_counts as read_step_counts_run_1 {
            input:
                sc = select_first([step_counts_run_1]),
                docker_image_id,
        }
        
        call mbt.read_step_counts as read_step_counts_run_2 {
            input:
                sc = select_first([step_counts_run_2]),
                docker_image_id,
        }
        call mbt.merge_step_counts {
            input:
                step_counts = select_all([read_step_counts_run_1.step_counts, read_step_counts_run_2.step_counts]),
                docker_image_id,
        }
    }
  
    if (defined(taxon_counts_run_2) || defined(contig_summary_run_2) || defined(contig_fasta_run_2)) {
        call preprocess_taxa as preprocess_taxa_ref_nt {
            input:
                taxon_counts = select_first([taxon_counts_run_2]),
                contig_summary = select_first([contig_summary_run_2]),
                contig_fasta = select_first([contig_fasta_run_2]),
                taxadb_sqlite = make_taxadb.taxadb_sqlite,
                dbtype = "NT",
                docker_image_id 
        }
        
        call preprocess_taxa as preprocess_taxa_ref_nr {
            input:
                taxon_counts = select_first([taxon_counts_run_2]),
                contig_summary = select_first([contig_summary_run_2]),
                contig_fasta = select_first([contig_fasta_run_2]),
                taxadb_sqlite = make_taxadb.taxadb_sqlite,
                dbtype = "NR",
                docker_image_id 
        }

    }
    call notebook as test_notebook { 
        input: 
            preprocessed_taxa_nt = preprocess_taxa_nt.preprocessed_taxa,
            preprocessed_taxa_nr = preprocess_taxa_nr.preprocessed_taxa,
            preprocessed_taxa_ref_nt = preprocess_taxa_ref_nt.preprocessed_taxa,
            preprocessed_taxa_ref_nr = preprocess_taxa_ref_nr.preprocessed_taxa,
            ground_truth_nt = truth_benchmark_nt.ground_truth_output,
            ground_truth_nr = truth_benchmark_nr.ground_truth_output,
            step_counts = merge_step_counts.step_count_tsv,
            docker_image_id
    }

    output {
        File preprocessed_nt = preprocess_taxa_nt.preprocessed_taxa
        File preprocessed_nr = preprocess_taxa_nr.preprocessed_taxa
        File benchmark_notebook = test_notebook.benchmark_notebook
        File? step_counts_run_1_json = read_step_counts_run_1.step_counts
        File? step_counts_run_2_json = read_step_counts_run_2.step_counts
        File? step_count_tsv = merge_step_counts.step_count_tsv
        File? truth_nt = truth_benchmark_nt.ground_truth_output
        File? truth_nr = truth_benchmark_nr.ground_truth_output
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
        File? preprocessed_taxa_ref_nt
        File? preprocessed_taxa_ref_nr
        File? ground_truth_nt
        File? ground_truth_nr
        File? step_counts
        String docker_image_id
    }
    String jq_cmd = '[(input | {NT: .}), (input | {NR: .})] | reduce .[] as $item ({}; . * $item)'

    command <<<
    set -euxo pipefail
    jq -n '~{jq_cmd}' \
        "~{preprocessed_taxa_nt}" "~{preprocessed_taxa_nr}" > combined_taxa.json

    export HARVEST_DATA="combined_taxa.json"
    
    if [[ "~{defined(preprocessed_taxa_ref_nt)}" == "true" &&  "~{defined(preprocessed_taxa_ref_nr)}" == "true" ]]; then
        jq -n '~{jq_cmd}' \
            "~{preprocessed_taxa_ref_nt}" "~{preprocessed_taxa_ref_nr}" > ref_data.json
        export REF_LIB="ref_data.json";
    fi

    if [[ "~{defined(step_counts)}" == "true" ]]; then
        export STEP_COUNTS="~{step_counts}"
    fi

    if [[ "~{defined(ground_truth_nt)}" == "true" && "~{defined(ground_truth_nr)}" == "true" ]]; then 
        jq -n '~{jq_cmd}' \
            "~{ground_truth_nt}" "~{ground_truth_nr}" > ground_truth.json
        export GROUND_TRUTH="ground_truth.json";
    fi
    
    # TODO: handle empty NR or NT or both output
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

