version 1.1

workflow AMR {
    input {
        Array[File] non_host_reads
        File contigs
        String docker_image_id
        File card_json = "s3://czid-public-references/test/AMRv2/card.json"
        File kmer_db = "s3://czid-public-references/test/AMRv2/61_kmer_db.json"
        File amr_kmer_db = "s3://czid-public-references/test/AMRv2/all_amr_61mers.txt"
        File wildcard_data = "s3://czid-public-references/test/AMRv2/wildcard_data.tar.bz2"
        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }

    call RunRgiBwtKma {
        input:
        non_host_reads = non_host_reads,
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

        main_output = pd.read_csv("~{main_output}", sep="\t")
        main_output['Best_Hit_ARO'] = [i.lower() for i in main_output['Best_Hit_ARO']] 
        main_output.sort_values(by='Best_Hit_ARO', inplace=True)
        main_output.reset_index()

        main_species_output = pd.read_csv("~{main_species_output}", sep="\t")
        main_species_output['Best_Hit_ARO'] = [i.lower() for i in main_species_output['Best_Hit_ARO']] 
        main_species_output.sort_values(by = 'Best_Hit_ARO', inplace=True)
        main_species_output.reset_index()

        kma_output = pd.read_csv("~{kma_output}", sep="\t")
        kma_output['ARO Term'] = [i.lower() for i in kma_output['ARO Term']] 
        kma_output.sort_values(by = 'ARO Term', inplace=True)
        kma_output.reset_index()

        kma_species_output = pd.read_csv("~{kma_species_output}", sep="\t")
        kma_species_output['ARO term'] = [i.lower() for i in kma_species_output['ARO term']] 
        kma_species_output.sort_values(by='ARO term', inplace=True)
        kma_species_output.reset_index()

        # rename main-amr and main-species columns to account for duplciate column names
        main_output.columns = [i+'_contig_amr' for i in main_output.columns]
        main_species_output.columns = [i+'_contig_sp' for i in main_species_output.columns]

        # vv get rid of weird additional space in some of the contig amr fields to make matching work downstream
        main_output['Contig_contig_amr'] = [i.strip() for i in main_output['Contig_contig_amr']] 

        # merge the data where Best_Hit_ARO and Contig name must match
        merge_a = main_output.merge(main_species_output, left_on = ['Best_Hit_ARO_contig_amr', 'Contig_contig_amr'],
                                                                    right_on = ['Best_Hit_ARO_contig_sp', 'Contig_contig_sp'], 
                                                                    how='outer',
                                                                    suffixes = [None, None])
        merge_a['ARO_contig'] = [merge_a.iloc[i]['Best_Hit_ARO_contig_amr'] if str(merge_a.iloc[i]['Best_Hit_ARO_contig_amr']) != 'nan' else merge_a.iloc[i]['Best_Hit_ARO_contig_sp'] for i in range(len(merge_a.index))]
        merge_a.to_csv("merge_a.tsv", index=None, sep="\t")

        kma_output.columns = [i+'_kma_amr' for i in kma_output.columns]  # ARO Term
        kma_species_output.columns = [i+'_kma_sp' for i in kma_species_output.columns]  #ARO term
        merge_b = kma_output.merge(kma_species_output, left_on = 'ARO Term_kma_amr', right_on = 'ARO term_kma_sp',
                                   how = 'outer', suffixes = [None, None])

        
        # remove kma results from variant models (because these results are inaccurate)
        merge_b = merge_b[merge_b['Reference Model Type_kma_amr'] != "protein variant model"] # remove protein variant model
        merge_b = merge_b[merge_b['Reference Model Type_kma_amr'] != "rRNA gene variant model"] # remove rRNA variant model
        merge_b['ARO_kma'] = [merge_b.iloc[i]['ARO Term_kma_amr'] if str(merge_b.iloc[i]['ARO Term_kma_amr']) != 'nan' else merge_b.iloc[i]['ARO term_kma_sp'] for i in range(len(merge_b.index))]

        merge_b.to_csv("merge_b.tsv", index=None, sep="\t")

        # final merge of MAIN and KMA combined results
        merge_x = merge_a.merge(merge_b, left_on = 'ARO_contig',
                                right_on = 'ARO_kma', how='outer',
                                suffixes = [None, None])

        merge_x['ARO_overall'] = [merge_x.iloc[i]['ARO_contig'] if str(merge_x.iloc[i]['ARO_contig']) != 'nan' else merge_x.iloc[i]['ARO_kma'] for i in range(len(merge_x.index))]
        merge_x['Gene_Family_overall'] = [merge_x.iloc[i]['AMR Gene Family_contig_amr'] if not pd.isna(merge_x.iloc[i]['AMR Gene Family_contig_amr']) else merge_x.iloc[i]['AMR Gene Family_kma_amr'] for i in range(len(merge_x.index))]
        merge_x['Drug_Class_overall'] = [merge_x.iloc[i]['Drug Class_contig_amr'] if not pd.isna(merge_x.iloc[i]['Drug Class_contig_amr']) else merge_x.iloc[i]['Drug Class_kma_amr'] for i in range(len(merge_x.index))]
        merge_x['model_overall'] = [merge_x.iloc[i]['Model_type_contig_amr'] if not pd.isna(merge_x.iloc[i]['Model_type_contig_amr']) else merge_x.iloc[i]['Reference Model Type_kma_amr'] for i in range(len(merge_x.index))]
        merge_x['Resistance_Mechanism_overall'] = [merge_x.iloc[i]['Resistance Mechanism_contig_amr'] if not pd.isna(merge_x.iloc[i]['Resistance Mechanism_contig_amr']) else merge_x.iloc[i]['Resistance Mechanism_kma_amr'] for i in range(len(merge_x.index))]

        merge_x.to_csv("final_summary.tsv", index=None, sep='\t')

        CODE
    >>>
    output {
        File merge_a = "merge_a.tsv"
        File merge_b = "merge_b.tsv"
        File final_summary = "final_summary.tsv"

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
        source /usr/local/miniconda/etc/profile.d/conda.sh
        conda activate rgi
        mkdir -p wildcard
        tar -xjf "~{wildcard_data}" -C wildcard
        gunzip wildcard/*.gz
        rgi card_annotation -i "~{card_json}" > card_annotation.log
        rgi load -i "~{card_json}" --card_annotation card_database_*.fasta
        rgi wildcard_annotation -i wildcard/ --card_json "~{card_json}" -v 3.1.0 
        rgi load --wildcard_annotation wildcard_database_v3.1.0.fasta --wildcard_index wildcard/index-for-model-sequences.txt --card_annotation card_database_v3.2.3.fasta --local
        rgi load --kmer_database "~{kmer_db}" --amr_kmers "~{amr_kmer_db}" --kmer_size 61 --debug --local > kmer_load.61.log 2>&1
        rgi kmer_query --rgi -k 61 -i "~{main_output_json}" --output output.rgi.main.kmerspecies --local
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
        source /usr/local/miniconda/etc/profile.d/conda.sh
        conda activate rgi
        mkdir -p wildcard
        tar -xjf "~{wildcard_data}" -C wildcard
        gunzip wildcard/*.gz
        rgi card_annotation -i "~{card_json}" > card_annotation.log
        rgi load -i "~{card_json}" --card_annotation card_database_*.fasta
        rgi wildcard_annotation -i wildcard/ --card_json "~{card_json}" -v 3.1.0 
        rgi load --wildcard_annotation wildcard_database_v3.1.0.fasta --wildcard_index wildcard/index-for-model-sequences.txt --card_annotation card_database_v3.2.3.fasta --local
        rgi load --kmer_database "~{kmer_db}" --amr_kmers "~{amr_kmer_db}" --kmer_size 61 --debug --local 
        rgi kmer_query --bwt -k 61 -i "~{output_sorted_length_100}" --output output.rgi.kma.kmerspecies --local
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
        source /usr/local/miniconda/etc/profile.d/conda.sh
        conda activate rgi
        rgi card_annotation -i "~{card_json}" > card_annotation.log
        rgi load -i "~{card_json}" --card_annotation card_database_*.fasta
        rgi main -i "~{contigs}" -o output.rgi.main -t contig -a BLAST --clean

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
	source /usr/local/miniconda/etc/profile.d/conda.sh
        conda activate rgi
        rgi card_annotation -i "~{card_json}" > card_annotation.log 
        rgi load -i "~{card_json}" --card_annotation card_database_*.fasta  
    	rgi bwt -1 "~{non_host_reads[0]}" -2 "~{non_host_reads[0]}" -a kma -o output_kma.rgi.kma --clean
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
