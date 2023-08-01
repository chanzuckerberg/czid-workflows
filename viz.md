```mermaid
amr run.wdl
flowchart TD;
    if2[raw_reads_0] --> |true| id1{host_filter_stage};
    id1{host_filter_stage} --> |non_host_reads|id2{RunSpades};
    if2[raw_reads_0] --> |true| id2{RunSpades};
    id1{host_filter_stage} --> |non_host_reads|id3{RunRgiBwtKma};
    id2{RunSpades} --> |contigs|id4{RunRgiMain};
    id4{RunRgiMain} --> |main_amr_results,main_output_json|id5{MakeGeneCoverage};
    id3{RunRgiBwtKma} --> |output_sorted_length_100|id6{RunRgiKmerBwt};
    id4{RunRgiMain} --> |main_output_json|id7{RunRgiKmerMain};
    id4{RunRgiMain} --> |main_output|id8{RunResultsPerSample};
    id7{RunRgiKmerMain} --> |main_species_output|id8{RunResultsPerSample};
    id3{RunRgiBwtKma} --> |kma_output|id8{RunResultsPerSample};
    id6{RunRgiKmerBwt} --> |kma_species_output|id8{RunResultsPerSample};
    id5{MakeGeneCoverage} --> |gene_coverage|id8{RunResultsPerSample};
    id2{RunSpades} --> |contigs|id9{tsvToSam};
    id8{RunResultsPerSample} --> |final_summary|id9{tsvToSam};
    id2{RunSpades} --> |contigs_in|id10{ZipOutputs};
    id1{host_filter_stage} --> |nonHostReads|id10{ZipOutputs};
    id8{RunResultsPerSample} --> |mainReports|id10{ZipOutputs};
    id6{RunRgiKmerBwt} --> |rawReports|id10{ZipOutputs};
    id7{RunRgiKmerMain} --> |rawReports|id10{ZipOutputs};
    id4{RunRgiMain} --> |rawReports,intermediateFiles|id10{ZipOutputs};
    id3{RunRgiBwtKma} --> |rawReports,intermediateFiles|id10{ZipOutputs};
```