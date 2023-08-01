amr run.wdl
flowchart TD
    if2[raw_reads_0] --> |true| id1{host_filter_stage}
    id1{host_filter_stage} --> |non_host_reads|id2{RunSpades}
    if2[raw_reads_0] --> |true| id2{RunSpades}
    id1{host_filter_stage} --> |non_host_reads|id3{RunRgiBwtKma}
    id2{RunSpades} --> |contigs|id4{RunRgiMain}
    id4{RunRgiMain} --> |main_amr_results,main_output_json|id5{MakeGeneCoverage}
    id3{RunRgiBwtKma} --> |output_sorted_length_100|id6{RunRgiKmerBwt}
    id4{RunRgiMain} --> |main_output_json|id7{RunRgiKmerMain}
    id4{RunRgiMain} --> |main_output|id8{RunResultsPerSample}
    id7{RunRgiKmerMain} --> |main_species_output|id8{RunResultsPerSample}
    id3{RunRgiBwtKma} --> |kma_output|id8{RunResultsPerSample}
    id6{RunRgiKmerBwt} --> |kma_species_output|id8{RunResultsPerSample}
    id5{MakeGeneCoverage} --> |gene_coverage|id8{RunResultsPerSample}
    id2{RunSpades} --> |contigs|id9{tsvToSam}
    id8{RunResultsPerSample} --> |final_summary|id9{tsvToSam}
    id2{RunSpades} --> |contigs_in|id10{ZipOutputs}
    id1{host_filter_stage} --> |nonHostReads|id10{ZipOutputs}
    id8{RunResultsPerSample} --> |mainReports|id10{ZipOutputs}
    id6{RunRgiKmerBwt} --> |rawReports|id10{ZipOutputs}
    id7{RunRgiKmerMain} --> |rawReports|id10{ZipOutputs}
    id4{RunRgiMain} --> |rawReports,intermediateFiles|id10{ZipOutputs}
    id3{RunRgiBwtKma} --> |rawReports,intermediateFiles|id10{ZipOutputs}

long-read-mngs run.wdl
flowchart TD
    id1{RunValidateInput} --> |input_fastq|id2{RunQualityFilter}
    id2{RunQualityFilter} --> |input_fastq|id3{RunHostFilter}
    id3{RunHostFilter} --> |input_fastq|id4{RunHumanFilter}
    id4{RunHumanFilter} --> |input_fastq|id5{ReadLengthMetrics}
    id4{RunHumanFilter} --> |input_fastq|id6{RunSubsampling}
    id6{RunSubsampling} --> |input_fastq|id7{PreAssemblyFasta}
    id6{RunSubsampling} --> |input_fastq|id8{RunAssembly}
    id6{RunSubsampling} --> |input_fastq|id9{RunReadsToContigs}
    id8{RunAssembly} --> |assembled_reads|id9{RunReadsToContigs}
    id8{RunAssembly} --> |assembled_reads|id10{RemoveUnmappedContigs}
    id9{RunReadsToContigs} --> |reads_to_contig_tsv|id10{RemoveUnmappedContigs}
    id9{RunReadsToContigs} --> |reads_to_contig_tsv,reads_to_contigs_sam|id11{GenerateContigStats}
    id9{RunReadsToContigs} --> |non_contig_reads_fa|id12{PrepareNTAlignmentInputs}
    id10{RemoveUnmappedContigs} --> |assembled_reads_fa|id12{PrepareNTAlignmentInputs}
    id12{PrepareNTAlignmentInputs} --> |all_sequences_to_align|id13{RunNTAlignment}
    id10{RemoveUnmappedContigs} --> |assembled_reads_fa|id14{RunNRAlignment}
    id13{RunNTAlignment} --> |nt_m8|id15{FindTopHitsNT}
    id14{RunNRAlignment} --> |nr_m8|id16{FindTopHitsNR}
    id15{FindTopHitsNT} --> |nt_top_m8|id17{SummarizeHitsNT}
    id9{RunReadsToContigs} --> |read_to_contig_tsv|id17{SummarizeHitsNT}
    id16{FindTopHitsNR} --> |nr_top_m8|id18{SummarizeHitsNR}
    id9{RunReadsToContigs} --> |read_to_contig_tsv|id18{SummarizeHitsNR}
    id10{RemoveUnmappedContigs} --> |contigs_fasta|id19{GenerateCoverageStats}
    id9{RunReadsToContigs} --> |read_contig_sam|id19{GenerateCoverageStats}
    id6{RunSubsampling} --> |input_file|id20{UnmappedReads}
    id17{SummarizeHitsNT} --> |nt_hit_summary|id20{UnmappedReads}
    id18{SummarizeHitsNR} --> |nr_hit_summary|id20{UnmappedReads}
    id9{RunReadsToContigs} --> |reads_to_contigs_tsv|id20{UnmappedReads}
    id9{RunReadsToContigs} --> |read_to_contig_tsv|id21{SummarizeContigsNT}
    id17{SummarizeHitsNT} --> |nt_hit_summary,nt_m8_reassigned|id21{SummarizeContigsNT}
    id9{RunReadsToContigs} --> |read_to_contig_tsv|id22{SummarizeContigsNR}
    id18{SummarizeHitsNR} --> |nr_hit_summary,nr_m8_reassigned|id22{SummarizeContigsNR}
    id17{SummarizeHitsNT} --> |nt_hit_summary,nt_m8_reassigned|id23{ComputeMergedTaxonCounts}
    id21{SummarizeContigsNT} --> |nt_contig_summary_json|id23{ComputeMergedTaxonCounts}
    id18{SummarizeHitsNR} --> |nr_hit_summary,nr_m8_reassigned|id23{ComputeMergedTaxonCounts}
    id22{SummarizeContigsNR} --> |nr_contig_summary_json|id23{ComputeMergedTaxonCounts}
    id21{SummarizeContigsNT} --> |counts_json_files|id24{CombineTaxonCounts}
    id22{SummarizeContigsNR} --> |counts_json_files|id24{CombineTaxonCounts}
    id23{ComputeMergedTaxonCounts} --> |counts_json_files|id24{CombineTaxonCounts}
    id21{SummarizeContigsNT} --> |json_files|id25{CombineJson}
    id22{SummarizeContigsNR} --> |json_files|id25{CombineJson}
    id23{ComputeMergedTaxonCounts} --> |json_files|id25{CombineJson}
    id7{PreAssemblyFasta} --> |pre_alignment_fasta|id26{GenerateAnnotatedFasta}
    id17{SummarizeHitsNT} --> |nt_m8_reassigned|id26{GenerateAnnotatedFasta}
    id18{SummarizeHitsNR} --> |nr_m8_reassigned|id26{GenerateAnnotatedFasta}
    id26{GenerateAnnotatedFasta} --> |annotated_merged_fa|id27{GenerateTaxidFasta}
    id17{SummarizeHitsNT} --> |nt_hit_summary|id27{GenerateTaxidFasta}
    id18{SummarizeHitsNR} --> |nr_hit_summary|id27{GenerateTaxidFasta}
    id27{GenerateTaxidFasta} --> |assembly_refined_taxid_annot_fasta|id28{GenerateTaxidLocator}
    id17{SummarizeHitsNT} --> |nt_hit_summary,nt_m8_reassigned|id29{GenerateCoverageViz}
    id15{FindTopHitsNT} --> |nt_top_m8|id29{GenerateCoverageViz}
    id19{GenerateCoverageStats} --> |contig_in_contig_coverage_json|id29{GenerateCoverageViz}
    id11{GenerateContigStats} --> |contig_in_contig_stats_json|id29{GenerateCoverageViz}
    id10{RemoveUnmappedContigs} --> |contig_in_contigs_fasta|id29{GenerateCoverageViz}

consensus-genome run.wdl
flowchart TD
    if1[ref_accession_id != None] --> |true| id2{FetchSequenceByAccessionId}
    id1{ValidateInput} --> |fastqs|id3{ApplyLengthFilter}
    if2[technology == _ONT_ && apply_length_filter] --> |true| id3{ApplyLengthFilter}
    id3{ApplyLengthFilter} --> |fastqs|id4{RemoveHost}
    id1{ValidateInput} --> |fastqs|id4{RemoveHost}
    id4{RemoveHost} --> |fastqs|id5{QuantifyERCCs}
    if12[technology == _Illumina_] --> |true| id5{QuantifyERCCs}
    id4{RemoveHost} --> |fastqs|id6{FilterReads}
    id2{FetchSequenceByAccessionId} --> |ref_fasta|id6{FilterReads}
    if4[filter_reads] --> |true| id6{FilterReads}
    if12[technology == _Illumina_] --> |true| id6{FilterReads}
    id6{FilterReads} --> |fastqs|id7{TrimReads}
    id4{RemoveHost} --> |fastqs|id7{TrimReads}
    if6[trim_adapters] --> |true| id7{TrimReads}
    if12[technology == _Illumina_] --> |true| id7{TrimReads}
    id7{TrimReads} --> |fastqs|id8{AlignReads}
    id6{FilterReads} --> |fastqs|id8{AlignReads}
    id4{RemoveHost} --> |fastqs|id8{AlignReads}
    id2{FetchSequenceByAccessionId} --> |ref_fasta|id8{AlignReads}
    if12[technology == _Illumina_] --> |true| id8{AlignReads}
    id8{AlignReads} --> |alignments|id9{TrimPrimers}
    if12[technology == _Illumina_] --> |true| id9{TrimPrimers}
    id9{TrimPrimers} --> |bam|id10{MakeConsensus}
    if12[technology == _Illumina_] --> |true| id10{MakeConsensus}
    id9{TrimPrimers} --> |call_variants_bam|id11{CallVariants}
    id2{FetchSequenceByAccessionId} --> |ref_fasta|id11{CallVariants}
    if12[technology == _Illumina_] --> |true| id11{CallVariants}
    id2{FetchSequenceByAccessionId} --> |ref_fasta|id12{RealignConsensus}
    id10{MakeConsensus} --> |consensus|id12{RealignConsensus}
    if12[technology == _Illumina_] --> |true| id12{RealignConsensus}
    id4{RemoveHost} --> |fastqs|id13{RunMinion}
    if13[technology == _ONT_] --> |true| id13{RunMinion}
    id10{MakeConsensus} --> |assembly|id14{Quast}
    id13{RunMinion} --> |assembly,bam|id14{Quast}
    id9{TrimPrimers} --> |bam|id14{Quast}
    id7{TrimReads} --> |fastqs|id14{Quast}
    id6{FilterReads} --> |fastqs|id14{Quast}
    id4{RemoveHost} --> |fastqs|id14{Quast}
    id2{FetchSequenceByAccessionId} --> |ref_fasta|id14{Quast}
    id9{TrimPrimers} --> |cleaned_bam|id15{ComputeStats}
    id13{RunMinion} --> |vcf,assembly,cleaned_bam|id15{ComputeStats}
    id10{MakeConsensus} --> |assembly|id15{ComputeStats}
    id5{QuantifyERCCs} --> |ercc_stats|id15{ComputeStats}
    id11{CallVariants} --> |vcf|id15{ComputeStats}
    id4{RemoveHost} --> |fastqs|id15{ComputeStats}
    id4{RemoveHost} --> |outputFiles|id16{ZipOutputs}
    id10{MakeConsensus} --> |outputFiles|id16{ZipOutputs}
    id13{RunMinion} --> |outputFiles|id16{ZipOutputs}
    id15{ComputeStats} --> |outputFiles|id16{ZipOutputs}
    id9{TrimPrimers} --> |outputFiles|id16{ZipOutputs}
    id14{Quast} --> |outputFiles|id16{ZipOutputs}
    id8{AlignReads} --> |outputFiles|id16{ZipOutputs}
    id5{QuantifyERCCs} --> |outputFiles|id16{ZipOutputs}
    id11{CallVariants} --> |outputFiles|id16{ZipOutputs}
    id12{RealignConsensus} --> |outputFiles|id16{ZipOutputs}

phylotree-ng run.wdl
flowchart TD
    id1{GetSampleContigFastas} --> |sample_and_reference_fastas|id3{RunSKA}
    id2{GetReferenceAccessionFastas} --> |sample_and_reference_fastas|id3{RunSKA}
    id3{RunSKA} --> |ska_distances|id4{ComputeClusters}
    id4{ComputeClusters} --> |clusters_directory|id5{GenerateClusterPhylos}
    id3{RunSKA} --> |ska_hashes|id5{GenerateClusterPhylos}
    id3{RunSKA} --> |distances|id6{AddSampleNamesToDistances}
    id5{GenerateClusterPhylos} --> |variants|id7{AddSampleNamesToVariants}
    id5{GenerateClusterPhylos} --> |newick|id8{AddSampleNamesToNewick}
    if1[GenerateClusterPhylos.phylotree_newick != None] --> |true| id8{AddSampleNamesToNewick}

legacy-host-filter legacy-host-filter.wdl
flowchart TD
    id1{RunValidateInput} --> |validate_input_summary_json,valid_input_fastq|id2{RunStar}
    id2{RunStar} --> |unmapped_fastq|id3{RunTrimmomatic}
    id3{RunTrimmomatic} --> |trimmomatic_fastq|id4{RunPriceSeq}
    id4{RunPriceSeq} --> |priceseq_fa|id5{RunCZIDDedup}
    id5{RunCZIDDedup} --> |duplicate_clusters_csv,duplicate_cluster_sizes_tsv,dedup_fa|id6{RunLZW}
    id6{RunLZW} --> |lzw_fa|id7{RunBowtie2_bowtie2_out}
    id5{RunCZIDDedup} --> |duplicate_clusters_csv,duplicate_cluster_sizes_tsv,dedup_fa|id7{RunBowtie2_bowtie2_out}
    id7{RunBowtie2_bowtie2_out} --> |bowtie2_fa|decl3{RunSubsample}
    id5{RunCZIDDedup} --> |duplicate_clusters_csv,duplicate_cluster_sizes_tsv,dedup_fa|decl3{RunSubsample}
    decl3{RunSubsample} --> |subsampled_fa|id9{RunStarDownstream}
    id1{RunValidateInput} --> |validate_input_summary_json,valid_input_fastq|id9{RunStarDownstream}
    id5{RunCZIDDedup} --> |duplicate_clusters_csv,duplicate_cluster_sizes_tsv,dedup_fa|id9{RunStarDownstream}
    if2[host_genome != _human_] --> |true| id9{RunStarDownstream}
    id9{RunStarDownstream} --> |unmapped_human_fa|decl6{RunBowtie2_bowtie2_human_out}
    id5{RunCZIDDedup} --> |duplicate_clusters_csv,duplicate_cluster_sizes_tsv,dedup_fa|decl6{RunBowtie2_bowtie2_human_out}
    if2[host_genome != _human_] --> |true| decl6{RunBowtie2_bowtie2_human_out}
    decl3{RunSubsample} --> |subsampled_2_fa,subsampled_1_fa,subsampled_merged_fa| decl7((if_then_else))
    decl6{RunBowtie2_bowtie2_human_out} --> |bowtie2_human_merged_fa,bowtie2_human_2_fa,bowtie2_human_1_fa| decl7((if_then_else))
    decl7((if_then_else)) --> |subsampled_fa| id11{RunGsnapFilter}
    id5{RunCZIDDedup} --> |duplicate_clusters_csv,duplicate_cluster_sizes_tsv,dedup_fa|id11{RunGsnapFilter}

index-generation index_generation.wdl
flowchart TD
    id1{DownloadNR} --> |nr|id5{GenerateIndexAccessions}
    id2{DownloadNT} --> |nt|id5{GenerateIndexAccessions}
    id3{DownloadAccession2Taxid} --> |accession2taxid|id5{GenerateIndexAccessions}
    id2{DownloadNT} --> |nt|id6{GenerateNTDB}
    id1{DownloadNR} --> |nr|id7{GenerateNRDB}
    id1{DownloadNR} --> |nr|id8{GenerateIndexDiamond}
    id4{DownloadTaxdump} --> |taxdump|id9{GenerateIndexLineages}
    id2{DownloadNT} --> |nt|id10{GenerateIndexMinimap2}
    id9{GenerateIndexLineages} --> |versioned_taxid_lineages_csv|id11{LoadTaxonLineages}
    if1[write_to_db && defined_environ_ && defined_s3_dir_] --> |true| id11{LoadTaxonLineages}

short-read-mngs non_host_alignment.wdl
flowchart TD
    id1{RunAlignment_minimap2_out} --> |m8_file|id2{RunCallHitsMinimap2}
    id3{RunAlignment_diamond_out} --> |m8_file|id4{RunCallHitsDiamond}
    id2{RunCallHitsMinimap2} --> |counts_json_files|id5{CombineTaxonCounts}
    id4{RunCallHitsDiamond} --> |counts_json_files|id5{CombineTaxonCounts}
    id1{RunAlignment_minimap2_out} --> |gsnap_m8|id6{GenerateAnnotatedFasta}
    id2{RunCallHitsMinimap2} --> |gsnap_hitsummary_tab,gsnap_counts_with_dcr_json,gsnap_deduped_m8|id6{GenerateAnnotatedFasta}
    id3{RunAlignment_diamond_out} --> |rapsearch2_m8|id6{GenerateAnnotatedFasta}
    id4{RunCallHitsDiamond} --> |rapsearch2_deduped_m8,rapsearch2_hitsummary_tab,rapsearch2_counts_with_dcr_json|id6{GenerateAnnotatedFasta}

short-read-mngs host_filter.wdl
flowchart TD
    id1{RunValidateInput} --> |valid_input2_fastq,valid_input1_fastq|id2{ercc_bowtie2_filter}
    id2{ercc_bowtie2_filter} --> |bowtie2_ercc_filtered2_fastq,bowtie2_ercc_filtered1_fastq|id3{fastp_qc}
    id3{fastp_qc} --> |fastp2_fastq,fastp1_fastq|id4{kallisto}
    id3{fastp_qc} --> |fastp2_fastq,fastp1_fastq|id5{bowtie2_filter}
    id5{bowtie2_filter} --> |bowtie2_host_filtered1_fastq,bowtie2_host_filtered2_fastq|decl4{hisat2_filter}
    id5{bowtie2_filter} --> |bam|id7{collect_insert_size_metrics}
    if1[fastqs_1] --> |true| id7{collect_insert_size_metrics}
    decl4{hisat2_filter} --> |hisat2_host_filtered1_fastq,hisat2_host_filtered2_fastq|id8{bowtie2_human_filter}
    if3[host_genome != _human_] --> |true| id8{bowtie2_human_filter}
    id8{bowtie2_human_filter} --> |bowtie2_human_filtered1_fastq,bowtie2_human_filtered2_fastq|decl3{hisat2_human_filter}
    if3[host_genome != _human_] --> |true| decl3{hisat2_human_filter}
    decl3{hisat2_human_filter} --> |hisat2_human_filtered1_fastq| decl5((select_first))
    decl4{hisat2_filter} --> |hisat2_host_filtered1_fastq| decl5((select_first))
    decl3{hisat2_human_filter} --> |hisat2_human_filtered2_fastq| decl6((if_then_else))
    decl4{hisat2_filter} --> |hisat2_host_filtered2_fastq| decl6((if_then_else))
    decl5((select_first)) --> |hisat2_filtered1_fastq| id10{RunCZIDDedup}
    decl6((if_then_else)) --> |hisat2_filtered2_fastq| id10{RunCZIDDedup}
    id10{RunCZIDDedup} --> |dedup2_fastq,dedup1_fastq,duplicate_cluster_sizes_tsv|id11{RunSubsample}

short-read-mngs postprocess.wdl
flowchart TD
    id1{RunAssembly} --> |assembly_read_contig_sam,assembly_contig_stats_json,assembly_scaffolds_fasta,assembly_contigs_fasta|id2{GenerateCoverageStats}
    id1{RunAssembly} --> |assembly_read_contig_sam,assembly_contig_stats_json,assembly_scaffolds_fasta,assembly_contigs_fasta|id5{BlastContigs_refined_gsnap_out}
    id3{DownloadAccessions_gsnap_accessions_out} --> |assembly_nt_refseq_fasta|id5{BlastContigs_refined_gsnap_out}
    id1{RunAssembly} --> |assembly_read_contig_sam,assembly_contig_stats_json,assembly_scaffolds_fasta,assembly_contigs_fasta|id6{BlastContigs_refined_rapsearch2_out}
    id4{DownloadAccessions_rapsearch2_accessions_out} --> |assembly_nr_refseq_fasta|id6{BlastContigs_refined_rapsearch2_out}
    id5{BlastContigs_refined_gsnap_out} --> |nt_m8,nt_hitsummary2_tab,nt_contig_summary_json|id7{ComputeMergedTaxonCounts}
    id6{BlastContigs_refined_rapsearch2_out} --> |nr_m8,nr_hitsummary2_tab,nr_contig_summary_json|id7{ComputeMergedTaxonCounts}
    id5{BlastContigs_refined_gsnap_out} --> |counts_json_files|id8{CombineTaxonCounts}
    id6{BlastContigs_refined_rapsearch2_out} --> |counts_json_files|id8{CombineTaxonCounts}
    id7{ComputeMergedTaxonCounts} --> |counts_json_files|id8{CombineTaxonCounts}
    id5{BlastContigs_refined_gsnap_out} --> |json_files|id9{CombineJson}
    id6{BlastContigs_refined_rapsearch2_out} --> |json_files|id9{CombineJson}
    id7{ComputeMergedTaxonCounts} --> |json_files|id9{CombineJson}
    id5{BlastContigs_refined_gsnap_out} --> |assembly_gsnap_blast_top_m8,assembly_refined_gsnap_counts_with_dcr_json,assembly_gsnap_contig_summary_json,assembly_gsnap_reassigned_m8,assembly_gsnap_blast_m8,assembly_gsnap_hitsummary2_tab|id10{GenerateAnnotatedFasta}
    id6{BlastContigs_refined_rapsearch2_out} --> |assembly_rapsearch2_contig_summary_json,assembly_rapsearch2_reassigned_m8,assembly_rapsearch2_hitsummary2_tab,assembly_rapsearch2_blast_m8,assembly_refined_rapsearch2_counts_with_dcr_json,assembly_rapsearch2_blast_top_m8|id10{GenerateAnnotatedFasta}
    id10{GenerateAnnotatedFasta} --> |assembly_refined_unidentified_fa,assembly_refined_annotated_merged_fa|id11{GenerateTaxidFasta}
    id5{BlastContigs_refined_gsnap_out} --> |assembly_gsnap_blast_top_m8,assembly_refined_gsnap_counts_with_dcr_json,assembly_gsnap_contig_summary_json,assembly_gsnap_reassigned_m8,assembly_gsnap_blast_m8,assembly_gsnap_hitsummary2_tab|id11{GenerateTaxidFasta}
    id6{BlastContigs_refined_rapsearch2_out} --> |assembly_rapsearch2_contig_summary_json,assembly_rapsearch2_reassigned_m8,assembly_rapsearch2_hitsummary2_tab,assembly_rapsearch2_blast_m8,assembly_refined_rapsearch2_counts_with_dcr_json,assembly_rapsearch2_blast_top_m8|id11{GenerateTaxidFasta}
    id11{GenerateTaxidFasta} --> |assembly_refined_taxid_annot_fasta|id12{GenerateTaxidLocator}

short-read-mngs local_driver.wdl
flowchart TD
    id1{host_filter} --> |host_filter_out_gsnap_filter_1_fa,host_filter_out_gsnap_filter_2_fa,host_filter_out_gsnap_filter_merged_fa,czid_dedup_out_duplicate_clusters_csv,duplicate_cluster_sizes_tsv|id2{non_host_alignment}
    id1{host_filter} --> |host_filter_out_gsnap_filter_1_fa,host_filter_out_gsnap_filter_2_fa,host_filter_out_gsnap_filter_merged_fa,czid_dedup_out_duplicate_clusters_csv,duplicate_cluster_sizes_tsv|id3{postprocess}
    id2{non_host_alignment} --> |rapsearch2_out_rapsearch2_deduped_m8,rapsearch2_out_rapsearch2_m8,gsnap_out_gsnap_deduped_m8,rapsearch2_out_rapsearch2_hitsummary_tab,gsnap_out_gsnap_m8,gsnap_out_gsnap_hitsummary_tab,gsnap_out_gsnap_counts_with_dcr_json,rapsearch2_out_rapsearch2_counts_with_dcr_json|id3{postprocess}
    id2{non_host_alignment} --> |gsnap_m8_gsnap_deduped_m8,taxid_fasta_in_rapsearch2_hitsummary_tab,taxid_fasta_in_annotated_merged_fa,taxid_fasta_in_gsnap_hitsummary_tab|id4{experimental}
    id3{postprocess} --> |nonhost_fasta_refined_taxid_annot_fasta,refined_gsnap_in_gsnap_hitsummary2_tab,refined_gsnap_in_gsnap_blast_top_m8,contig_in_contig_coverage_json,refined_gsnap_in_gsnap_reassigned_m8,contig_in_contigs_fasta,contig_in_contig_stats_json|id4{experimental}
    id1{host_filter} --> |duplicate_clusters_csv|id4{experimental}

short-read-mngs experimental.wdl
flowchart TD
    id1{GenerateTaxidFasta} --> |taxid_annot_fasta|id2{GenerateTaxidLocator}
    id2{GenerateTaxidLocator} --> |taxid_annot_sorted_nr_fasta,taxid_locations_genus_nt_json,taxid_locations_combined_json,taxid_locations_nr_json,taxid_locations_nt_json,taxid_locations_family_nr_json,taxid_locations_family_nt_json,taxid_annot_sorted_family_nr_fasta,taxid_annot_sorted_genus_nr_fasta,taxid_locations_genus_nr_json,taxid_annot_sorted_family_nt_fasta,taxid_annot_sorted_genus_nt_fasta,taxid_annot_sorted_nt_fasta|id3{GenerateAlignmentViz}

