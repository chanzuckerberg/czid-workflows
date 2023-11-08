version 1.1

task merge_step_counts {
    input {
        Array[File] step_counts
        String docker_image_id
    }
    command <<<
        set -euxo pipefail
        python3 <<CODE
        import pandas as pd
        import json 

        df_run_1 = pd.read_json("~{step_counts[0]}", orient="values", typ="series")
        df_run_1.name = "run_1_step_counts"

        df_run_2 = pd.read_json("~{step_counts[1]}", orient="values", typ="series")
        df_run_2.name = "ref_step_counts"

        df = pd.merge(df_run_1, df_run_2, left_index=True, right_index=True, how='outer').sort_values(by="run_1_step_counts", ascending=False)
        df.to_csv("step_count.tsv", sep="\t")
        CODE
    >>>
    output {
        File step_count_tsv = "step_count.tsv"
    }
    runtime {
        docker: docker_image_id
    }
}

task read_step_counts {
    input {
        Array[File] sc
        String docker_image_id
    }
    command <<<
        jq -s 'add | walk(if type == "string" then tonumber else . end)' ~{sep=' ' sc} > step_counts.json
    >>>
    output {
        File step_counts = "step_counts.json"
    }
    runtime {
        docker: docker_image_id
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

        total_reads = sum([i["reads_dedup"] for i in tc.values()])
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
