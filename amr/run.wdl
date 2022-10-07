version 1.1

import "../short-read-mngs/host_filter.wdl" as host_filter

workflow amr {
    input {
        Array[File]? non_host_reads
        File? raw_reads_0
        File? raw_reads_1
        File? contigs
        String docker_image_id
        String sample_name
        String host_filtering_docker_image_id = "czid-short-read-mngs" # default local value
        File card_json = "s3://czid-public-references/test/AMRv2/card.json"
        File kmer_db = "s3://czid-public-references/test/AMRv2/61_kmer_db.json"
        File amr_kmer_db = "s3://czid-public-references/test/AMRv2/all_amr_61mers.txt"
        File wildcard_data = "s3://czid-public-references/test/AMRv2/wildcard_database_v3.1.0.fasta"
        File wildcard_index = "s3://czid-public-references/test/AMRv2/index-for-model-sequences.txt"
        Int min_contig_length = 100
        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }
    if (defined(raw_reads_0)) { 
        call host_filter.czid_host_filter as host_filter_stage { 
            input:
            fastqs_0 = select_first([raw_reads_0]),
            fastqs_1 = if defined(raw_reads_1) then select_first([raw_reads_1]) else None,
            docker_image_id = host_filtering_docker_image_id,
            s3_wd_uri = s3_wd_uri
        }
        call RunSpades {
            input:
            non_host_reads = select_all(
                [
                    host_filter_stage.gsnap_filter_out_gsnap_filter_1_fa,
                    host_filter_stage.gsnap_filter_out_gsnap_filter_2_fa
                ]
            ),
            min_contig_length = min_contig_length,
            docker_image_id = host_filtering_docker_image_id,
        }
    }
    call RunRgiBwtKma {
        input:
        non_host_reads = select_first([non_host_reads, 
            select_all(
                [
                    host_filter_stage.gsnap_filter_out_gsnap_filter_1_fa,
                    host_filter_stage.gsnap_filter_out_gsnap_filter_2_fa
                ]
            )]),
        card_json = card_json, 
        docker_image_id = docker_image_id
    }
    call RunRgiMain {
        input:
        contigs = select_first([contigs, RunSpades.contigs]),
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
        docker_image_id = docker_image_id,
        sample_name = sample_name
    }
    call ZipOutputs {
        input:
        nonHostReads = select_first([
                           non_host_reads,
                           select_all(
                               [
                                   host_filter_stage.gsnap_filter_out_gsnap_filter_1_fa,
                                   host_filter_stage.gsnap_filter_out_gsnap_filter_2_fa
                               ]
                           )
                       ]),
        outputFiles = select_all(
            [
                select_first([contigs, RunSpades.contigs]),
            ]
        ),
        mainReports = select_all(
            [
                RunResultsPerSample.final_summary,
                RunResultsPerSample.synthesized_report,
            ]
        ),
        rawReports = select_all(
            [
                RunRgiKmerBwt.kma_species_calling,
                RunRgiKmerBwt.sr_species_allele,
                RunRgiKmerMain.species_calling,
                RunRgiMain.main_amr_results,
                RunRgiBwtKma.kma_amr_results,
                RunRgiBwtKma.gene_mapping_data,
            ]
        ),
        intermediateFiles = select_all(
            [
                RunRgiBwtKma.artifacts_mapping_stats,
                RunRgiBwtKma.overall_mapping_stats,
                RunRgiBwtKma.reference_mapping_stats,
                RunRgiBwtKma.output_sorted_length_100,
                RunRgiBwtKma.output_sorted_length_100_bai,
            ]
        ),
        docker_image_id = docker_image_id
    }

}

task RunResultsPerSample {
    input {
        File main_output
        File main_species_output
        File kma_output
        File kma_species_output
        String docker_image_id
        String sample_name
    }
    command <<< 
        set -euxo pipefail

        python3 <<CODE
        import pandas as pd
        def clean_aro(df, column_name):
            """modifies dataframe inplace to clean the ARO string"""
            df[column_name] = df[column_name].map(lambda x: x.lower().strip())


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
        
        append_suffix_to_colname(kma_output, "_kma_amr")  # ARO Term
        append_suffix_to_colname(kma_species_output, "_kma_sp")  # ARO term

        merge_b = kma_output.merge(kma_species_output, left_on = 'ARO Term_kma_amr', right_on = 'ARO term_kma_sp',
                                   how = 'outer', suffixes = [None, None])

        
        # remove kma results from variant models (because these results are inaccurate)
        merge_b = merge_b[merge_b['Reference Model Type_kma_amr'] != "protein variant model"] # remove protein variant model
        merge_b = merge_b[merge_b['Reference Model Type_kma_amr'] != "rRNA gene variant model"] # remove rRNA variant model
        merge_b["ARO_kma"] = this_or_that(
            merge_b, "ARO Term_kma_amr", "ARO term_kma_sp"
        )
        merge_b.to_csv("merge_b.tsv", index=None, sep="\t")

        # final merge of MAIN and KMA combined results
        merge_x = merge_a.merge(merge_b, left_on = 'ARO_contig',
                                right_on = 'ARO_kma', how='outer',
                                suffixes = [None, None])

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
        big_table = df[['ARO_overall', 'Gene_Family_overall', 'Drug_Class_overall', 'Resistance_Mechanism_overall', 'model_overall', 'All Mapped Reads_kma_amr', 'Percent Coverage_kma_amr','Depth_kma_amr', 'CARD*kmer Prediction_kma_sp', 'Cut_Off_contig_sp', 'Percentage Length of Reference Sequence_contig_amr', 'Best_Identities_contig_amr', 'CARD*kmer Prediction_contig_sp']]
        big_table.sort_values(by='Gene_Family_overall', inplace=True)

        big_table.to_csv("bigtable_report.tsv", sep='\t', index=None)

        if big_table.empty: # if no outputs, simply return
            open("primary_AMR_report.tsv", "a").close()
            exit(0)

        def remove_na(input_set):
            set_list = list(input_set)
            return([i for i in set_list if i == i])

        this_list = list(set(df['ARO_overall']))
        result_df = {}
        for index in this_list:#[0:1]:
            print(index)
            sub_df = df[df['ARO_overall'] == index]
            result = {}
            gf = remove_na(set(sub_df['AMR Gene Family_contig_amr']).union(set(sub_df['AMR Gene Family_kma_amr'])))
            result['sample_name'] = "~{sample_name}"
            result['gene_family'] = ';'.join(gf) if len(gf) > 0 else None
        
            dc = remove_na(set(sub_df['Drug Class_contig_amr']).union(set(sub_df['Drug Class_kma_amr'])))
            result['drug_class'] = ';'.join(dc) if len(dc) > 0 else None
            rm = remove_na(set(sub_df['Resistance Mechanism_contig_amr']).union(set(sub_df['Resistance Mechanism_kma_amr'])))
            result['resistance_mechanism'] = ';'.join(rm) if len(rm) > 0 else None
            
            co = remove_na(set(sub_df['Cut_Off_contig_amr']))
            result['cutoff'] = ';'.join(co) if len(co) > 0 else None
            
            mt = remove_na(set(sub_df['Model_type_contig_amr']).union(set(sub_df['Reference Model Type_kma_amr'])))
            result['model_type'] = ';'.join([' '.join(i.split(' ')[0:-1]) for i in mt]) if len(mt) > 0 else None
            
            cb = remove_na(set(sub_df['Percentage Length of Reference Sequence_contig_amr']).union(
            set(sub_df['Percent Coverage_kma_amr'])))
            result['coverage_breadth'] = max(cb) if len(cb) > 0 else None
            
            cd = remove_na(set(sub_df['Depth_kma_amr']))
            result['coverage_depth'] = max(cd) if len(cd) > 0 else None
            
            result['num_contigs'] = len(remove_na(set(sub_df['Contig_contig_amr'])))
            
            nr = remove_na(set(sub_df['All Mapped Reads_kma_amr']))
            result['num_reads'] = max(nr) if len(nr) > 0 else None
            
            pid = remove_na(set(sub_df['Best_Identities_contig_amr']))
            result['percent_id'] = max(pid) if len(pid) > 0 else None
            
            sp_contig = ' '.join(remove_na(set(sub_df['CARD*kmer Prediction_contig_sp'])))
            result['sp_contig'] = sp_contig
            
            sp_kma = ' '.join(remove_na(set(sub_df['CARD*kmer Prediction_kma_sp'])))
            result['sp_kma'] = sp_kma
            
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

            if(len(final_species)) > 0:
                result['species'] = max(final_species, key = final_species.get)
            else:
                result['species'] = None

            result_df[index] = result
        final_df = pd.DataFrame.from_dict(result_df)
        final_df = final_df.transpose()
        final_df.sort_index(inplace=True)
        final_df.sort_values(by='gene_family', inplace=True)
        #print(final_df.columns)
        final_df.dropna(subset=['drug_class'], inplace=True)
        final_df[['sample_name', 'gene_family', 'drug_class', 'resistance_mechanism', 'model_type', 'num_reads', 'num_contigs', 'coverage_breadth', 'coverage_depth', 'percent_id', 'cutoff', 'species', 'sp_contig', 'sp_kma']]
        final_df.to_csv("primary_AMR_report.tsv", sep='\t', index=None)


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
        File kmer_db
        File amr_kmer_db
        File wildcard_data 
        File wildcard_index
        String docker_image_id
    }
    command <<< 
        set -exuo pipefail
        time rgi load \
            --wildcard_annotation "~{wildcard_data}" \
            --wildcard_version 3.1.0 \
            --wildcard_index "~{wildcard_index}" \
            --kmer_database "~{kmer_db}" \
            --amr_kmers "~{amr_kmer_db}" \
            --kmer_size 61
        rgi kmer_query --rgi -k 61 -i "~{main_output_json}" --output contig_species_report 
    >>>
    output { 
        Array[File] output_kmer_main = glob("contig_species_report*")
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
            --wildcard_annotation "~{wildcard_data}" \
            --wildcard_version 3.1.0 \
            --wildcard_index "~{wildcard_index}" \
            --kmer_database "~{kmer_db}" \
            --amr_kmers "~{amr_kmer_db}" \
            --kmer_size 61
        rgi kmer_query --bwt -k 61 -i "~{output_sorted_length_100}" --output sr_species_report
    >>>
    output {
        Array[File] output_kmer_bwt = glob("sr_species_report*")
        File sr_species_allele = "sr_species_report_61mer_analysis.allele.txt"
        File kma_species_calling = "sr_species_report_61mer_analysis.gene.txt"
    }
    runtime {
        docker: docker_image_id
    }

}
task RunRgiMain { 
    input { 
        File contigs
        File card_json
        String docker_image_id
    }
    command <<<
        set -exuo pipefail
        if [[ $(head -n 1 "~{contigs}") == ";ASSEMBLY FAILED" ]]; then
            # simulate empty outputs
            echo "{}" > contig_amr_report.json
            cp /tmp/empty-main-header.txt contig_amr_report.txt
        else
            rgi main -i "~{contigs}" -o contig_amr_report -t contig -a BLAST --clean --include_nudge
        fi
    >>>
    output {
        Array[File] output_main = glob("contig_amr_report*")
        File output_json = "contig_amr_report.json"
        File main_amr_results = "contig_amr_report.txt"
    }

    runtime {
        docker: docker_image_id
    }
}
task RunRgiBwtKma {
    input {
        Array[File] non_host_reads
        File card_json
        String docker_image_id
    }

    command <<<
        set -exuo pipefail
        rgi bwt -1 ~{sep=' -2 ' non_host_reads} -a kma -o sr_amr_report --clean
    >>>

    output {
        Array[File] output_kma = glob("sr_amr_report*")
        File kma_amr_results = "sr_amr_report.allele_mapping_data.txt"
        File artifacts_mapping_stats = "sr_amr_report.artifacts_mapping_stats.txt"
        File gene_mapping_data = "sr_amr_report.gene_mapping_data.txt"
        File overall_mapping_stats = "sr_amr_report.overall_mapping_stats.txt"
        File reference_mapping_stats = "sr_amr_report.reference_mapping_stats.txt"
        File output_sorted_length_100 = "sr_amr_report.sorted.length_100.bam"
        File output_sorted_length_100_bai = "sr_amr_report.sorted.length_100.bam.bai"
    }

    runtime {
        docker: docker_image_id
    }
}
task RunSpades { 
    input { 
        Array[File] non_host_reads
        Int min_contig_length
        String docker_image_id
    }
    command <<< 
        set -euxo pipefail
        function handle_failure() 
        {
            echo ";ASSEMBLY FAILED" > spades/contigs.fasta
            echo ";ASSEMBLY FAILED" > spades/scaffolds.fasta
            exit 0

        }
        trap handle_failure ERR
        if [[ "~{length(non_host_reads)}" -gt 1 ]]; then 
            spades.py -1 ~{sep=" -2 " non_host_reads} -o "spades/" -m 100 -t 36 --only-assembler 1>&2
        else
            spades.py -s ~{non_host_reads[0]} -o "spades/" -m 100 -t 36 --only-assembler 1>&2
        fi
        seqtk seq -L ~{min_contig_length} spades/contigs.fasta > spades/contigs_filtered.fasta
        mv spades/contigs_filtered.fasta spades/contigs.fasta
    >>>
    output { 
        File contigs = "spades/contigs.fasta"
        File scaffolds = "spades/scaffolds.fasta"
    }
    runtime { 
        docker: docker_image_id
    }
}
task ZipOutputs {
    input {
        Array[File] nonHostReads
        Array[File] outputFiles
        Array[File] mainReports
        Array[File] rawReports
        Array[File] intermediateFiles
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        export TMPDIR=${TMPDIR:-/tmp}

        mkdir ${TMPDIR}/outputs
        mkdir ${TMPDIR}/outputs/final_reports
        mkdir ${TMPDIR}/outputs/raw_reports
        mkdir ${TMPDIR}/outputs/intermediate_files

        counter=1
        for fastx in ~{sep= ' ' nonHostReads}; do 
            cp $fastx ${TMPDIR}/outputs/non_host_reads_R$counter."${fastx#*.}"
            ((counter++))
        done
        cp ~{sep=' ' outputFiles} ${TMPDIR}/outputs/
        cp ~{sep=' ' mainReports} ${TMPDIR}/outputs/final_reports
        cp ~{sep=' ' rawReports} ${TMPDIR}/outputs/raw_reports
        cp ~{sep=' ' intermediateFiles} ${TMPDIR}/outputs/intermediate_files
        export WORK=$(pwd)
        cd ${TMPDIR}/outputs; zip -r ${WORK}/outputs.zip .
    >>>

    output {
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
    df = pd.read_csv("~{main_amr_results}", delimiter="\t").loc[:, ["ORF_ID", "ID", "Model_ID", "Hit_Start", "Hit_End", "Percentage Length of Reference Sequence"]]

    # create seq length reference map
    with open("~{main_output_json}") as json_file:
        rgi_main_json = json.load(json_file)
    db_seq_length = {}
    for ind, row in df.iterrows():
        db_seq_length[row["Model_ID"]] = len(rgi_main_json[row["ORF_ID"]][row["ID"]]["dna_sequence_from_broadstreet"])

    agg_res = df.groupby(["ID", "Model_ID"]).agg(lambda x: list(x))

    gene_coverage = []
    for ind, row, in agg_res.iterrows():
      max_end = -1
      gene_coverage_bps = 0
      for start, end, in sorted(zip(row["Hit_Start"], row["Hit_End"]), key=lambda x: x[0]):
        gene_coverage_bps += max(0, # if end-max(max_end, start) is negative
                      # segment is contained in a previous segment, so don't add
                     end - max(max_end, start) # if max_end > start, don't double count the already added portion of the gene
                    )
        max_end = max(max_end, end)
      gene_coverage.append({
        "ID": ind[0],
        "gene_coverage_bps": gene_coverage_bps,
        "db_seq_length": db_seq_length[ind[1]],
        "gene_coveraage_perc": np.round((gene_coverage_bps/db_seq_length[ind[1]])*100, 4)
      })
    gene_coverage_df = pd.DataFrame(gene_coverage)
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
