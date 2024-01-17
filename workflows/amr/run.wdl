version 1.1

import "../short-read-mngs/host_filter.wdl" as host_filter

struct RawSample {
    Array[File]+ raw_reads
}

struct FilteredSample {
    Array[File]+ subsampled_reads
    Array[File]+ non_host_reads
    File clusters
    File cluster_sizes
    File contigs
}


workflow amr {
    input {
        RawSample? raw_sample
        FilteredSample? filtered_sample
        String docker_image_id
        String sample_name
        String host_filtering_docker_image_id = "czid-short-read-mngs" # default local value
        File card_json = "s3://czid-public-references/card/2023-05-22/card.json"
        File card_ontology = "s3://czid-public-references/amr_v2/ontology/initial/ontology.json"
        File kmer_db = "s3://czid-public-references/card/2023-05-22/61_kmer_db.json"
        File amr_kmer_db = "s3://czid-public-references/card/2023-05-22/all_amr_61mers.txt"
        File wildcard_data = "s3://czid-public-references/card/2023-05-22/wildcard_database_v4.0.0.fasta"
        File wildcard_index = "s3://czid-public-references/card/2023-05-22/index-for-model-sequences.txt"
        Int min_contig_length = 100
        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }

    if (defined(raw_sample)) {
        RawSample raw_sample_in = select_first([raw_sample])
        Array[File]+ raw_reads = select_all(raw_sample_in.raw_reads)

        call host_filter.czid_host_filter as host_filter_stage { 
            input:
            fastqs_0 = raw_reads[0],
            fastqs_1 = if length(raw_reads) > 1 then raw_reads[1] else None,
            docker_image_id = host_filtering_docker_image_id,
            s3_wd_uri = s3_wd_uri
        }

        Array[File]+ host_filtered_reads = if (defined(host_filter_stage.hisat2_human_filtered1_fastq))
            then select_all([host_filter_stage.hisat2_human_filtered1_fastq, host_filter_stage.hisat2_human_filtered2_fastq])
            else select_all([host_filter_stage.hisat2_host_filtered1_fastq, host_filter_stage.hisat2_host_filtered2_fastq])

        call RunRedup as RunRedupRaw {
            input:
            host_filtered_reads = host_filtered_reads,
            subsampled_reads = select_all([host_filter_stage.subsampled_out_subsampled_1_fa, host_filter_stage.subsampled_out_subsampled_2_fa]),
            clusters = host_filter_stage.czid_dedup_out_duplicate_clusters_csv,
            cluster_sizes = host_filter_stage.czid_dedup_out_duplicate_cluster_sizes_tsv,
            docker_image_id = docker_image_id,
        }
        call RunSpades {
            input:
            reduplicated_reads = RunRedupRaw.redups_fa,
            min_contig_length = min_contig_length,
            docker_image_id = host_filtering_docker_image_id,
        }
    }

    if (defined(filtered_sample)) {
        FilteredSample filtered_sample_in = select_first([filtered_sample])
        File filtered_sample_contigs = filtered_sample_in.contigs

        call RunRedup as RunRedupFiltered {
            input:
            host_filtered_reads = filtered_sample_in.non_host_reads,
            subsampled_reads = filtered_sample_in.subsampled_reads,
            clusters = filtered_sample_in.clusters,
            cluster_sizes = filtered_sample_in.cluster_sizes,
            docker_image_id = docker_image_id,
        }
    }

    Array[File]+ non_host_reads_fa = select_first([RunRedupRaw.redups_fa, RunRedupFiltered.redups_fa])
    File contigs_fa = select_first([RunSpades.contigs, filtered_sample_contigs])

    call RunRgiBwtKma {
        input:
        non_host_reads_fa = non_host_reads_fa,
        card_json = card_json, 
        docker_image_id = docker_image_id
    }
    call RunRgiMain {
        input:
        contigs_fa = contigs_fa,
        card_json = card_json, 
        docker_image_id = docker_image_id
    }
    call MakeGeneCoverage {
        input:
        main_amr_results = RunRgiMain.main_amr_results,
        main_output_json = RunRgiMain.output_json,
        docker_image_id = docker_image_id
    }
    call RunRgiKmerBwt { 
        input:
        output_sorted_length_100 = RunRgiBwtKma.output_sorted_length_100,
        card_json = card_json,
        kmer_db = kmer_db,
        amr_kmer_db = amr_kmer_db,
        wildcard_data = wildcard_data,
        wildcard_index = wildcard_index,
        docker_image_id = docker_image_id
    }
    call RunRgiKmerMain { 
        input:
        main_output_json = RunRgiMain.output_json,
        card_json = card_json,
        kmer_db = kmer_db,
        amr_kmer_db = amr_kmer_db,
        wildcard_data = wildcard_data,
        wildcard_index = wildcard_index,
        docker_image_id = docker_image_id
    }
    call RunResultsPerSample { 
        input: 
        main_output = RunRgiMain.main_amr_results,
        main_species_output = RunRgiKmerMain.species_calling,
        kma_output = RunRgiBwtKma.kma_amr_results,
        kma_species_output = RunRgiKmerBwt.kma_species_calling,
        gene_coverage = MakeGeneCoverage.output_gene_coverage,
        card_ontology = card_ontology,
        docker_image_id = docker_image_id,
        sample_name = sample_name
    }

    call tsvToSam { 
        input: 
        contigs_fa = contigs_fa,
        final_summary = RunResultsPerSample.final_summary,
        docker_image_id = docker_image_id
    }

    call ZipOutputs {
        input:
        contigs_fa = contigs_fa,
        non_host_reads_fa = non_host_reads_fa,
        main_reports = select_all(
            [
                RunResultsPerSample.final_summary,
                RunResultsPerSample.synthesized_report,
            ]
        ),
        raw_reports = select_all(
            [
                RunRgiKmerBwt.kma_species_calling,
                RunRgiKmerBwt.sr_species_allele,
                RunRgiKmerMain.species_calling,
                RunRgiMain.main_amr_results,
                RunRgiBwtKma.kma_amr_results,
                RunRgiBwtKma.gene_mapping_data,
            ]
        ),
        intermediate_files = select_all(
            [
                RunRgiBwtKma.artifacts_mapping_stats,
                RunRgiBwtKma.overall_mapping_stats,
                RunRgiBwtKma.reference_mapping_stats,
                RunRgiBwtKma.output_sorted_length_100,
                RunRgiBwtKma.output_sorted_length_100_bai,
                RunRgiMain.output_json,
                RunRgiBwtKma.kma_amr_results_json,
            ]
        ),
        docker_image_id = docker_image_id
    }

}

task RunRedup {
    input {
        Array[File]+ host_filtered_reads
        Array[File]+ subsampled_reads
        File clusters
        File cluster_sizes
        String docker_image_id
    }
    command <<<
        set -euxo pipefail
        # exit if no duplicate reads
        if [[ ! $(grep -v "^1\t" "~{cluster_sizes}") ]]; then
            counter=1
            for fasta in ~{sep=' ' subsampled_reads}; do
                cp $fasta redups_$counter.fa
                ((counter++))
            done
            exit 0
        fi

        grep -v "^1\t" "~{cluster_sizes}" | cut -f2 > duplicated-reads.txt
        grep -h ">" ~{sep=' ' subsampled_reads} | sed "s/^>//" > passed_filters.txt

        python3 <<CODE
        pair_values = []
        passed_filters = set()
        duplicates = set()
        with open("passed_filters.txt", "r") as pf:
            passed_filters.update(pf.read().splitlines())
        with open("duplicated-reads.txt", "r") as dr:
            duplicates.update(dr.read().splitlines())
        with open("~{clusters}", "r") as clusters:
            for line in clusters:
                key, value = line.strip().split(",")
                if key in duplicates and key in passed_filters and key != value:
                    pair_values.append(value)

        with open("duplicated-pairs.txt", "w+") as f:
            for value in pair_values:
                f.write(value)
                f.write("\n")
        CODE

        # Extract the duplicate reads from the host filtered reads file in fasta format,
        # then concatenate those reads to the end of the corresponding subsampled reads file
        export HOST_FILTERED_READS_FILES=(~{sep=' ' host_filtered_reads})
        export SUBSAMPLED_READS_FILES=(~{sep=' ' subsampled_reads})

        for index in ${!HOST_FILTERED_READS_FILES[@]}; do
            seqkit grep -f duplicated-pairs.txt ${HOST_FILTERED_READS_FILES[$index]} | seqkit fq2fa -o duplicate_reads_$index.fasta
            cat ${SUBSAMPLED_READS_FILES[$index]} duplicate_reads_$index.fasta > redups_$index.fasta
        done
    >>>
    output {
        Array[File]+ redups_fa = glob("redups*.fasta")
    }
    runtime {
        docker: docker_image_id
    }
}

task RunResultsPerSample {
    input {
        File main_output
        File main_species_output
        File kma_output
        File kma_species_output
        File gene_coverage
        File card_ontology
        String docker_image_id
        String sample_name
    }
    command <<< 
        set -euxo pipefail

        python3 <<CODE
        import json
        import pandas as pd
        def clean_aro(df, column_name):
            """modifies dataframe inplace to clean the ARO string"""
            df[column_name] = df[column_name].map(lambda x: x.strip())


        def append_suffix_to_colname(df, suffix):
            """append suffix to column name in place"""
            df.rename(lambda x: x + suffix, axis="columns", inplace=True)


        def this_or_that(df, this_colname, that_colname):
            """returns this if not nan, else returns that"""
            thisorthat = df.apply(
                lambda x: x[this_colname] if not pd.isna(x[this_colname]) else x[that_colname],
                axis=1,
            )
            if thisorthat.empty:
                return ""
            
            return thisorthat


        main_output = pd.read_csv("~{main_output}", sep="\t")
        clean_aro(main_output, "Best_Hit_ARO")
        

        main_output.sort_values(by='Best_Hit_ARO', inplace=True)

        main_species_output = pd.read_csv("~{main_species_output}", sep="\t")
        clean_aro(main_species_output, "Best_Hit_ARO")

        main_species_output.sort_values(by = 'Best_Hit_ARO', inplace=True)

        kma_output = pd.read_csv("~{kma_output}", sep="\t")
        clean_aro(kma_output, "ARO Term")

        kma_output.sort_values(by = 'ARO Term', inplace=True)

        kma_species_output = pd.read_csv("~{kma_species_output}", sep="\t")
        clean_aro(kma_species_output, "ARO term")
        kma_species_output.sort_values(by='ARO term', inplace=True)

        # rename main-amr and main-species columns to account for duplciate column names
        append_suffix_to_colname(main_output, "_contig_amr")
        append_suffix_to_colname(main_species_output, "_contig_sp")

        # vv get rid of weird additional space in some of the contig amr fields to make matching work downstream
        main_output["Contig_contig_amr"] = main_output["Contig_contig_amr"].map(
            lambda x: x.strip()
        )
        # merge the data where Best_Hit_ARO and Contig name must match
        merge_a = main_output.merge(main_species_output, left_on = ['Best_Hit_ARO_contig_amr', 'Contig_contig_amr'],
                                                                    right_on = ['Best_Hit_ARO_contig_sp', 'Contig_contig_sp'], 
                                                                    how='outer',
                                                                    suffixes = [None, None])
        merge_a["ARO_contig"] = this_or_that(
            merge_a, "Best_Hit_ARO_contig_amr", "Best_Hit_ARO_contig_sp"
        )
        merge_a.to_csv("merge_a.tsv", index=None, sep="\t")

        merge_a['ARO_contig_amr'] = merge_a['ARO_contig_amr'].astype(str)
        
        append_suffix_to_colname(kma_output, "_kma_amr")  # ARO Term
        append_suffix_to_colname(kma_species_output, "_kma_sp")  # ARO term

        merge_b = kma_output.merge(kma_species_output, left_on = 'ARO Term_kma_amr', right_on = 'ARO term_kma_sp',
                                   how = 'outer', suffixes = [None, None])
        
        merge_b['ARO Accession_kma_amr'] = merge_b['ARO Accession_kma_amr'].astype(str)

        
        # remove kma results from variant models (because these results are inaccurate)
        merge_b = merge_b[merge_b['Reference Model Type_kma_amr'] != "protein variant model"] # remove protein variant model
        merge_b = merge_b[merge_b['Reference Model Type_kma_amr'] != "rRNA gene variant model"] # remove rRNA variant model
        merge_b["ARO_kma"] = this_or_that(
            merge_b, "ARO Term_kma_amr", "ARO term_kma_sp"
        )
        merge_b.to_csv("merge_b.tsv", index=None, sep="\t")

        # final merge of MAIN and KMA combined results
        merge_x = merge_a.merge(merge_b, left_on = 'ARO_contig_amr',
                                right_on = 'ARO Accession_kma_amr', how='outer',
                                suffixes = [None, None])

        merge_x["ARO_accession"] = this_or_that(merge_x, "ARO_contig_amr", "ARO Accession_kma_amr")
        merge_x["ARO_overall"] = this_or_that(merge_x, "ARO_contig", "ARO_kma")
        merge_x["Gene_Family_overall"] = this_or_that(
            merge_x, "AMR Gene Family_contig_amr", "AMR Gene Family_kma_amr"
        )
        merge_x["Drug_Class_overall"] = this_or_that(
            merge_x, "Drug Class_contig_amr", "Drug Class_kma_amr"
        )
        merge_x["model_overall"] = this_or_that(
            merge_x, "Model_type_contig_amr", "Reference Model Type_kma_amr"
        )
        merge_x["Resistance_Mechanism_overall"] = this_or_that(
            merge_x, "Resistance Mechanism_contig_amr", "Resistance Mechanism_kma_amr"
        )
        merge_x.to_csv("comprehensive_AMR_metrics.tsv", index=None, sep='\t')

        df = pd.read_csv("comprehensive_AMR_metrics.tsv", sep='\t')
        df["ARO_contig_amr"] = df["ARO_contig_amr"].apply(lambda x: f'{x:.0f}' if x==x else x)  # Convert non-nan to int float
        
        big_table = df[['ARO_overall', 'Gene_Family_overall', 'Drug_Class_overall', 'Resistance_Mechanism_overall', 'model_overall', 'All Mapped Reads_kma_amr', 'Percent Coverage_kma_amr','Depth_kma_amr', 'CARD*kmer Prediction_kma_sp', 'Cut_Off_contig_amr', 'Percentage Length of Reference Sequence_contig_amr', 'Best_Identities_contig_amr', 'CARD*kmer Prediction_contig_sp']]
        big_table.sort_values(by='Gene_Family_overall', inplace=True)

        big_table.to_csv("bigtable_report.tsv", sep='\t', index=None)

        if big_table.empty: # if no outputs, simply return
            open("primary_AMR_report.tsv", "a").close()
            exit(0)

        gene_coverage = pd.read_csv("~{gene_coverage}", sep='\t')
        def remove_na(input_set):
            set_list = list(input_set)
            return([i for i in set_list if i == i])

        ontology = json.load(open("~{card_ontology}"))
        def get_high_level_classes(drug_class_string_list):
            drug_classes = set()
            for drug_class_string in drug_class_string_list:
                drug_classes |= set(drug_class_string.split('; '))
            high_level_classes = set()
            for drug_class in drug_classes:
                if drug_class in ontology:
                    high_level_classes |= set(ontology[drug_class].get('highLevelDrugClasses', []))
            return high_level_classes

        this_list = list(set(df['ARO_overall']))
        result_df = {}
        for index in this_list:#[0:1]:
            sub_df = df[df['ARO_overall'] == index]
            result = {}
            gf = remove_na(set(sub_df['AMR Gene Family_contig_amr']).union(set(sub_df['AMR Gene Family_kma_amr'])))
            result['sample_name'] = "~{sample_name}"
            result['gene_family'] = ';'.join(gf) if len(gf) > 0 else None

            accession = list(set(sub_df['ARO_accession'].dropna().apply(lambda x: f'ARO:{int(x)}')))
            result['aro_accession'] = ';'.join(accession) if len(accession) > 0 else None
        
            dc = remove_na(set(sub_df['Drug Class_contig_amr']).union(set(sub_df['Drug Class_kma_amr'])))
            result['drug_class'] = ';'.join(dc) if len(dc) > 0 else None

            hldc = get_high_level_classes(dc)
            # Join classes with a semicolon and a space since we have a python list instead of a nice pandas union
            result['high_level_drug_class'] = '; '.join(hldc) if len(hldc) > 0 else None

            rm = remove_na(set(sub_df['Resistance Mechanism_contig_amr']).union(set(sub_df['Resistance Mechanism_kma_amr'])))
            result['resistance_mechanism'] = ';'.join(rm) if len(rm) > 0 else None

            
            sub_df.loc[(sub_df['Cut_Off_contig_amr'] == 'Strict') & (sub_df['Nudged_contig_amr'] == True), 'Cut_Off_contig_amr'] = "Nudged"
            co = remove_na(set(sub_df['Cut_Off_contig_amr']))

            result['cutoff'] = ';'.join(co) if len(co) > 0 else None
            
            mt = remove_na(set(sub_df['Model_type_contig_amr']).union(set(sub_df['Reference Model Type_kma_amr'])))
            result['model_type'] = ';'.join([' '.join(i.split(' ')[0:-1]) for i in mt]) if len(mt) > 0 else None
            
            rcb = remove_na(set(sub_df['Percent Coverage_kma_amr']))
            result['read_coverage_breadth'] = max(rcb) if len(rcb) > 0 else None

            gene_id = ";".join(remove_na(set(sub_df["ARO_contig_amr"].unique())))
            if gene_id:
                contig_coverage = gene_coverage[gene_coverage["ID"].astype('str') == gene_id]["gene_coverage_perc"].iloc[0]
            else:
                contig_coverage = None
            result["contig_coverage_breadth"] = contig_coverage
            
            cd = remove_na(set(sub_df['Depth_kma_amr']))
            result['read_coverage_depth'] = max(cd) if len(cd) > 0 else None
            
            result['num_contigs'] = len(remove_na(set(sub_df['Contig_contig_amr'])))
            
            nr = remove_na(set(sub_df['All Mapped Reads_kma_amr']))
            result['num_reads'] = max(nr) if len(nr) > 0 else None

            read_gi = remove_na(set(sub_df["Reference Sequence_kma_amr"]))
            result["read_gene_id"] = ";".join(read_gi) if len(read_gi) > 0 else None
            
            pid = remove_na(set(sub_df['Best_Identities_contig_amr']))
            result['contig_percent_id'] = sum(pid) / len(pid) if len(pid) > 0 else None
            
            contig_species = ' '.join(remove_na(set(sub_df['CARD*kmer Prediction_contig_sp'])))
            result['contig_species'] = contig_species
            
            read_species = ' '.join(remove_na(set(sub_df['CARD*kmer Prediction_kma_sp'])))
            result['read_species'] = read_species
            
            
            
            sp = remove_na(set(sub_df['CARD*kmer Prediction_contig_sp']).union(set(sub_df['CARD*kmer Prediction_kma_sp'])))
            final_species = {} 
            for s in sp:
                #print(s)
                if 'Unknown' in s:
                    #print('a')
                    final_species[s] = 0
                    continue
                else:
                    split_sp = s.split(';')
                    for t in split_sp:
                        #print(t.strip())
                        this_t = t.strip()
                        if ':' in this_t: #not just a space
                            species_name = t.split(':')[0]
                            species_count = int(t.split(':')[1].strip())
                            final_species[species_name] = species_count

            result_df[index] = result
        final_df = pd.DataFrame.from_dict(result_df)
        final_df = final_df.transpose()
        final_df = final_df[["sample_name", "gene_family", "aro_accession", "drug_class", "high_level_drug_class", "resistance_mechanism", "model_type", "num_contigs",
                             "cutoff", "contig_coverage_breadth", "contig_percent_id", "contig_species", "num_reads", "read_gene_id", "read_coverage_breadth", "read_coverage_depth", "read_species"]]
        final_df.sort_index(inplace=True)
        final_df.dropna(subset=['drug_class'], inplace=True)
        final_df.to_csv("primary_AMR_report.tsv", sep='\t', index_label="gene_name")


        CODE
    >>>
    output {
        File merge_a = "merge_a.tsv"
        File merge_b = "merge_b.tsv"
        File final_summary = "comprehensive_AMR_metrics.tsv"
        File bigtable = "bigtable_report.tsv"
        File synthesized_report = "primary_AMR_report.tsv"
    }
    runtime {
        docker: docker_image_id
    }
}

task RunRgiKmerMain {
    input {
        File main_output_json
        File card_json
        File kmer_db
        File amr_kmer_db
        File wildcard_data 
        File wildcard_index
        String docker_image_id
    }
    command <<< 
        set -exuo pipefail
        time rgi load \
            -i "~{card_json}" \
            --wildcard_annotation "~{wildcard_data}" \
            --wildcard_version 4.0.0 \
            --wildcard_index "~{wildcard_index}" \
            --kmer_database "~{kmer_db}" \
            --amr_kmers "~{amr_kmer_db}" \
            --kmer_size 61
        rgi kmer_query --rgi -k 61 -i "~{main_output_json}" --output contig_species_report 
    >>>
    output { 
        Array[File]+ output_kmer_main = glob("contig_species_report*")
        File species_calling = "contig_species_report_61mer_analysis_rgi_summary.txt"
    }
    runtime {
        docker: docker_image_id
    }

}

task RunRgiKmerBwt { 
    input {
        File output_sorted_length_100
        File card_json
        File kmer_db
        File amr_kmer_db
        File wildcard_data 
        File wildcard_index
        String docker_image_id
    }
    command <<<
        set -exuo pipefail
        
        time rgi load \
            -i "~{card_json}" \
            --wildcard_annotation "~{wildcard_data}" \
            --wildcard_version 4.0.0 \
            --wildcard_index "~{wildcard_index}" \
            --kmer_database "~{kmer_db}" \
            --amr_kmers "~{amr_kmer_db}" \
            --kmer_size 61
        rgi kmer_query --bwt -k 61 -i "~{output_sorted_length_100}" --output sr_species_report
    >>>
    output {
        Array[File]+ output_kmer_bwt = glob("sr_species_report*")
        File sr_species_allele = "sr_species_report_61mer_analysis.allele.txt"
        File kma_species_calling = "sr_species_report_61mer_analysis.gene.txt"
    }
    runtime {
        docker: docker_image_id
    }

}
task RunRgiMain { 
    input { 
        File contigs_fa
        File card_json
        String docker_image_id
    }
    command <<<
        set -exuo pipefail
        if [[ $(head -n 1 "~{contigs_fa}") == ";ASSEMBLY FAILED" ]]; then
            # simulate empty outputs
            echo "{}" > contig_amr_report.json
            cp /tmp/empty-main-header.txt contig_amr_report.txt
        else
            rgi main -i "~{contigs_fa}" -o contig_amr_report -t contig -a BLAST --include_nudge --clean
        fi
    >>>
    output {
        Array[File]+ output_main = glob("contig_amr_report*")
        File output_json = "contig_amr_report.json"
        File main_amr_results = "contig_amr_report.txt"
    }

    runtime {
        docker: docker_image_id
    }
}
task RunRgiBwtKma {
    input {
        Array[File]+ non_host_reads_fa
        File card_json
        String docker_image_id
    }

    command <<<
        set -exuo pipefail
        rgi bwt -1 ~{sep=' -2 ' non_host_reads_fa} -a kma -o sr_amr_report --clean
    >>>

    output {
        Array[File]+ output_kma = glob("sr_amr_report*")
        File kma_amr_results = "sr_amr_report.allele_mapping_data.txt"
        File artifacts_mapping_stats = "sr_amr_report.artifacts_mapping_stats.txt"
        File gene_mapping_data = "sr_amr_report.gene_mapping_data.txt"
        File overall_mapping_stats = "sr_amr_report.overall_mapping_stats.txt"
        File reference_mapping_stats = "sr_amr_report.reference_mapping_stats.txt"
        File output_sorted_length_100 = "sr_amr_report.sorted.length_100.bam"
        File output_sorted_length_100_bai = "sr_amr_report.sorted.length_100.bam.bai"
        File kma_amr_results_json = "sr_amr_report.allele_mapping_data.json"
    }

    runtime {
        docker: docker_image_id
    }
}
task RunSpades { 
    input { 
        Array[File]+ reduplicated_reads
        Int min_contig_length
        String docker_image_id
    }
    command <<< 
        set -euxo pipefail
        function handle_failure() 
        {
            echo ";ASSEMBLY FAILED" > spades/contigs.fasta
            echo ";ASSEMBLY FAILED" > spades/scaffolds.fasta

            python3 <<CODE
    import csv
    import os

    WARNINGS_LOG = "spades/warnings.log"
    SPADES_LOG = "spades/spades.log"

    failure_info = {
        "log_present": False,
        "warnings_present": False,
        "warnings": None,
        "stack_trace": None,
        "errors": None,
    }

    if os.path.isfile(WARNINGS_LOG):
        failure_info["warnings_present"] = True
        with open(WARNINGS_LOG) as warnings_file:
            failure_info["warnings"] = "".join(warnings_file.readlines()).encode("unicode_escape").decode("utf-8")

    if os.path.isfile(SPADES_LOG):
        failure_info["log_present"] = True
        with open(SPADES_LOG) as log_file:
            log = log_file.readlines()
            error_lines = [line.replace("== Error ==", "*") for line in log if line.startswith("== Error ==")]
            if len(error_lines) > 0:
                failure_info["errors"] = "".join(error_lines).encode("unicode_escape").decode("utf-8")
            try:
                stack_trace_line_no = log.index("=== Stack Trace ===\n")
                stack_trace_end = log.index("\n", stack_trace_line_no)
                stack_trace = "".join(log[stack_trace_line_no:stack_trace_end])
                failure_info["stack_trace"] = stack_trace.encode("unicode_escape").decode("utf-8")
            except ValueError:
                pass

    with open("spades_failure.tsv", "w") as tsvfile:
        fieldnames = list(failure_info.keys())
        writer = csv.DictWriter(tsvfile, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(failure_info)
    CODE

            exit 0

        }
        trap handle_failure ERR
        if [[ "~{length(reduplicated_reads)}" -gt 1 ]]; then
            spades.py -1 ~{sep=" -2 " reduplicated_reads} -o "spades/" -m 100 -t 36 --only-assembler 1>&2
        else
            spades.py -s ~{reduplicated_reads[0]} -o "spades/" -m 100 -t 36 --only-assembler 1>&2
        fi

        if [[ $(head -n 1 spades/contigs.fasta) ==  "" ]]; then
            handle_failure
        else
            seqtk seq -L ~{min_contig_length} spades/contigs.fasta > spades/contigs_filtered.fasta
            mv spades/contigs_filtered.fasta spades/contigs.fasta
        fi
    >>>
    output { 
        File contigs = "spades/contigs.fasta"
        File? scaffolds = "spades/scaffolds.fasta"
        File? failure_log = "spades_failure.tsv"
    }
    runtime { 
        docker: docker_image_id
    }
}
task ZipOutputs {
    input {
        File contigs_fa
        Array[File]+ non_host_reads_fa
        Array[File]+ main_reports
        Array[File]+ raw_reports
        Array[File]+ intermediate_files
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        export TMPDIR=${TMPDIR:-/tmp}

        mkdir ${TMPDIR}/outputs
        mkdir ${TMPDIR}/outputs/final_reports
        mkdir ${TMPDIR}/outputs/raw_reports
        mkdir ${TMPDIR}/outputs/intermediate_files

        # copy contigs and interleave non_host_reads
        cp ~{contigs_fa} contigs.fasta

        if [[ "~{length(non_host_reads_fa)}" == 2 ]]; then
            seqfu ilv -1 ~{sep=" -2 " non_host_reads_fa} > non_host_reads_fa_ilv.fasta
            seqkit rename -1 -s "/" non_host_reads_fa_ilv.fasta -o non_host_reads.fasta
        else 
            cat ~{sep=" " non_host_reads_fa} > non_host_reads.fasta
        fi

        cp ~{sep=' ' main_reports} ${TMPDIR}/outputs/final_reports
        cp ~{sep=' ' raw_reports} ${TMPDIR}/outputs/raw_reports
        cp ~{sep=' ' intermediate_files} ${TMPDIR}/outputs/intermediate_files
        export WORK=$(pwd)
        cd ${TMPDIR}/outputs; zip -r ${WORK}/outputs.zip .
    >>>

    output {
        File non_host_reads = "non_host_reads.fasta"
        File contigs = "contigs.fasta"
        File output_zip = "outputs.zip"
    }

    runtime {
        docker: docker_image_id
    }
}

task MakeGeneCoverage {
    input { 
        File main_amr_results 
        File main_output_json
        String docker_image_id
    }
    command <<<
    set -euxo pipefail
    python3 <<CODE
    import pandas as pd
    import numpy as np
    import json
    df = pd.read_csv("~{main_amr_results}", delimiter="\t").loc[:, ["ORF_ID", "ID", "ARO", "Model_ID", "Hit_Start", "Hit_End", "Percentage Length of Reference Sequence"]]

    # create seq length reference map
    with open("~{main_output_json}") as json_file:
        rgi_main_json = json.load(json_file)
    db_seq_length = {}
    for ind, row in df.iterrows():
        db_seq_length[row["Model_ID"]] = len(rgi_main_json[row["ORF_ID"]][row["ID"]]["dna_sequence_from_broadstreet"])

    agg_res = df.groupby(["ARO", "Model_ID"]).agg(lambda x: list(x))

    gene_coverage = []
    # We are measuring gene coverage by keeping track of the coordinate ranges
    # of the reference genome that are covered by the contigs.
    # Here, ind is a tuple: (ARO, Model_ID)
    for ind, row, in agg_res.iterrows():
        covered_coords = set()
        for start, end in zip(row["Hit_Start"], row["Hit_End"]):
            contig_length = end - start
            # Get a positive (1) or negative (-1) step for our range, depending on if
            # we have a forward or reverse contig.
            sign = int(contig_length / abs(contig_length))
            covered_coords = covered_coords.union(range(start, end, sign))
        gene_coverage_bps = len(covered_coords)
        gene_coverage.append({
            "ID": ind[0],
            "gene_coverage_bps": gene_coverage_bps,
            "db_seq_length": db_seq_length[ind[1]],
            "gene_coverage_perc": np.round((gene_coverage_bps/db_seq_length[ind[1]])*100, 4)
        })

    if gene_coverage:
        gene_coverage_df = pd.DataFrame(gene_coverage)
        gene_coverage_df.to_csv("gene_coverage.tsv", index=None, sep="\t")
    else:
        gene_coverage_df = pd.DataFrame(columns=["ID", "gene_coverage_bps", "db_seq_length", "gene_coverage_perc"])
        gene_coverage_df.to_csv("gene_coverage.tsv", index=None, sep="\t")
    CODE
    >>>
    output {
        File output_gene_coverage = "gene_coverage.tsv"
    }

    runtime {
        docker: docker_image_id
    } 
}

task tsvToSam {
    input {
        File contigs_fa
        File final_summary
        String docker_image_id
    }

    command <<<
        set -euxo pipefail
        python3 <<CODE
        import pandas as pd
        import pysam
        import re
        import sys

        # NOTE: Contigs are indexed by ARO accession, not gene_id. This is because contigs and reads
        # come from different sources, so not every contig will have a corresponding gene_id,
        # even if you join the two datasets on the ARO accession. Thus indexing on gene_id is not
        # possible.

        COLUMN_ARO_CONTIG = "ARO_contig_amr"
        COLUMN_CONTIG_NAME = "Contig_contig_amr"
        OUTPUT_BAM = "contig_amr_report.sorted.bam"

        # Create an empty BAM/BAI file if the SPADES assembly failed, then exit
        with open("~{contigs_fa}") as f:
            first_line = f.readline()
            if first_line == ";ASSEMBLY FAILED\n":
                output_bam = pysam.AlignmentFile(OUTPUT_BAM, "wb", reference_names=["NoGenes"], reference_lengths=[100])
                output_bam.close()
                pysam.index(OUTPUT_BAM)
                sys.exit()

        # Create index to enable querying fasta file
        pysam.faidx("~{contigs_fa}")
        contigs_fasta = pysam.Fastafile("~{contigs_fa}")

        # Load columns of interest from CSV and drop rows with at least one NaN
        df = pd.read_csv("~{final_summary}", sep="\t", usecols=[COLUMN_ARO_CONTIG, COLUMN_CONTIG_NAME])
        df = df.dropna(subset=COLUMN_CONTIG_NAME)

        # Format accessions to follow the pattern 'ARO:3000000'; Remove extraneous _* at the end of contig names
        df[COLUMN_ARO_CONTIG] = df[COLUMN_ARO_CONTIG].apply(lambda x: f'ARO:{int(x)}')
        df[COLUMN_CONTIG_NAME] = df[COLUMN_CONTIG_NAME].apply(lambda x: re.sub(r'_\d*$', '', x))

        # Create BAM with mock reference lengths for the header. If contigss are found at all, have a mock accession
        # to make sure we can create the SAM file with no errors (web app will look for that file)
        contig_aros = df[COLUMN_ARO_CONTIG].dropna().unique().tolist() or ["NoGenes"]
        output_bam = pysam.AlignmentFile(OUTPUT_BAM, "wb", reference_names=contig_aros, reference_lengths=[100] * len(contig_aros))

        # Go through each line of the TSV and create a SAM record (https://wckdouglas.github.io/2021/12/pytest-with-pysam)
        for index, row in df.iterrows():
            contig_aro = row[COLUMN_ARO_CONTIG]
            contig_name = row[COLUMN_CONTIG_NAME]
            contig_sequence = contigs_fasta.fetch(contig_name)

            # Create new alignment
            alignment = pysam.AlignedSegment(output_bam.header)
            alignment.reference_name = contig_aro
            alignment.query_name = contig_name
            alignment.query_sequence = contig_sequence
            alignment.reference_start = 1
            output_bam.write(alignment)

        output_bam.close()
        pysam.index(OUTPUT_BAM)
        CODE
    >>>

    output {
        File output_sorted = "contig_amr_report.sorted.bam"
        File output_sorted_bai = "contig_amr_report.sorted.bam.bai"
    }

    runtime {
        docker: docker_image_id
    }
}
