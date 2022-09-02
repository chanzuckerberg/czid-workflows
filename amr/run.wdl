version 1.1

import "../short-read-mngs/host_filter.wdl" as host_filter

workflow AMR {
    input {
        Array[File]? non_host_reads
        File? raw_reads_0
        File? raw_reads_1
        File contigs
        String docker_image_id
        File card_json = "s3://czid-public-references/test/AMRv2/card.json"
        File kmer_db = "s3://czid-public-references/test/AMRv2/61_kmer_db.json"
        File amr_kmer_db = "s3://czid-public-references/test/AMRv2/all_amr_61mers.txt"
        File wildcard_data = "s3://czid-public-references/test/AMRv2/wildcard_data.tar.bz2"
        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }
    if (defined(raw_reads_0)) { 
        call host_filter.czid_host_filter as host_filter_stage { 
            input:
            fastqs_0 = select_first([raw_reads_0]),
            fastqs_1 = if defined(raw_reads_1) then select_first([raw_reads_1]) else None,
            docker_image_id = "czid-short-read-mngs",
            s3_wd_uri = s3_wd_uri
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
        contigs = contigs,
        card_json = card_json, 
        docker_image_id = docker_image_id
    }
    call RunRgiKmerBwt { 
        input:
        output_sorted_length_100 = RunRgiBwtKma.output_sorted_length_100,
        card_json = card_json,
        kmer_db = kmer_db,
        amr_kmer_db = amr_kmer_db,
        wildcard_data = wildcard_data,
        docker_image_id = docker_image_id
    }
    call RunRgiKmerMain { 
        input:
        main_output_json = RunRgiMain.output_json,
        card_json = card_json, 
        kmer_db = kmer_db,
        amr_kmer_db = amr_kmer_db,
        wildcard_data = wildcard_data,
        docker_image_id = docker_image_id
    }
    call RunResultsPerSample { 
        input: 
        main_output = RunRgiMain.main_amr_results,
        main_species_output = RunRgiKmerMain.species_calling,
        kma_output = RunRgiBwtKma.kma_amr_results,
        kma_species_output = RunRgiKmerBwt.kma_species_calling,
        docker_image_id = docker_image_id
    }
    call ZipOutputs {
        input:
        outputFiles = select_all(
            flatten([
                RunRgiBwtKma.output_kma,
                RunRgiMain.output_main,
                RunRgiKmerBwt.output_kmer_bwt,
                RunRgiKmerMain.output_kmer_main,
                [
                    RunResultsPerSample.merge_a,
                    RunResultsPerSample.merge_b,
                    RunResultsPerSample.final_summary,
                    RunResultsPerSample.bigtable,
                ]
            ])
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
            return df.apply(
                lambda x: x[this_colname] if not pd.isna(x[this_colname]) else x[that_colname],
                axis=1,
            )


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
        merge_b['ARO_kma'] = [merge_b.iloc[i]['ARO Term_kma_amr'] if str(merge_b.iloc[i]['ARO Term_kma_amr']) != 'nan' else merge_b.iloc[i]['ARO term_kma_sp'] for i in range(len(merge_b.index))]
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
        merge_x.to_csv("final_summary.tsv", index=None, sep='\t')

        df = pd.read_csv("final_summary.tsv", sep='\t')
        big_table = df[['ARO_overall', 'Gene_Family_overall', 'Drug_Class_overall', 'Resistance_Mechanism_overall', 'model_overall', 'All Mapped Reads_kma_amr', 'Percent Coverage_kma_amr','Depth_kma_amr', 'CARD*kmer Prediction_kma_sp', 'Cut_Off_contig_sp', 'Percentage Length of Reference Sequence_contig_amr', 'Best_Identities_contig_amr', 'CARD*kmer Prediction_contig_sp']]
        big_table.sort_values(by='Gene_Family_overall', inplace=True)

        big_table.to_csv("bigtable_report.tsv", sep='\t', index=None)

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
        final_df[['gene_family', 'drug_class', 'resistance_mechanism', 'model_type', 'num_reads', 'num_contigs', 'coverage_breadth', 'coverage_depth', 'percent_id', 'cutoff', 'species', 'sp_contig', 'sp_kma']]
        final_df.to_csv("synthesized_report.tsv", sep='\t', index=None)


        CODE
    >>>
    output {
        File merge_a = "merge_a.tsv"
        File merge_b = "merge_b.tsv"
        File final_summary = "final_summary.tsv"
        File bigtable = "bigtable_report.tsv"
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
        String docker_image_id
    }
    command <<< 
        set -exuo pipefail
        rgi kmer_query --rgi -k 61 -i "~{main_output_json}" --output output.rgi.main.kmerspecies 
    >>>
    output { 
        Array[File] output_kmer_main = glob("output.rgi.main.kmerspecies*")
        File species_calling = "output.rgi.main.kmerspecies_61mer_analysis_rgi_summary.txt"
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
        String docker_image_id
    }
    command <<<
        set -exuo pipefail
        rgi kmer_query --bwt -k 61 -i "~{output_sorted_length_100}" --output output.rgi.kma.kmerspecies
    >>>
    output {
        Array[File] output_kmer_bwt = glob("output.rgi.kma.kmerspecies*")
        File kma_species_calling = "output.rgi.kma.kmerspecies_61mer_analysis.gene.txt"
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
        rgi main -i "~{contigs}" -o output.rgi.main -t contig -a BLAST --clean --include_nudge

    >>>
    output {
        Array[File] output_main = glob("output.rgi.main*")
        File output_json = "output.rgi.main.json"
        File main_amr_results = "output.rgi.main.txt"
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
        rgi bwt -1 ~{sep=' -2 ' non_host_reads} -a kma -o output_kma.rgi.kma --clean
    >>>

    output {
        Array[File] output_kma = glob("output_kma.rgi.kma*")
        File output_sorted_length_100 = "output_kma.rgi.kma.sorted.length_100.bam"
        File kma_amr_results = "output_kma.rgi.kma.allele_mapping_data.txt"
    }

    runtime {
        docker: docker_image_id
    }
}
task ZipOutputs {
    input {
        Array[File] outputFiles
        String docker_image_id
    }

    command <<<
        set -euxo pipefail

        export TMPDIR=${TMPDIR:-/tmp}

        mkdir ${TMPDIR}/outputs
        cp ~{sep=' ' outputFiles} ${TMPDIR}/outputs/
        zip -r -j outputs.zip ${TMPDIR}/outputs/
    >>>

    output {
        File output_zip = "outputs.zip"
    }

    runtime {
        docker: docker_image_id
    }
}
