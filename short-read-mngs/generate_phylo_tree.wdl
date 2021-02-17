version development

struct TaxonByteRange {
  Int first_byte
  Int last_byte
  String refined_taxid_annotated_sorted_fasta_s3_path
}

struct SampleInfo {
  String sample_name
  String run_id
  File hitsummary2_nr
  File hitsummary2_nt
  TaxonByteRange taxon_byterange_nr
  TaxonByteRange taxon_byterange_nt
}

task PrepareTaxonFasta {
  input {
    String docker_image_id
    String s3_wd_uri
    Array[SampleInfo] samples
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name generate_phylo_tree \
    --step-module idseq_dag.steps.prepare_taxon_fasta \
    --step-class PipelineStepPrepareTaxonFasta \
    --step-name prepare_taxon_fasta \
    --input-files '[["~{write_json(samples)}"]]' \
    --output-files '[]' \
    --output-dir-s3 '~{s3_wd_uri}' \
  >>>
  output {
    Array[File] taxon_fastas = glob("taxon_fastas/*.fasta")
  }
  runtime {
    docker: docker_image_id
  }
}

task GeneratePhyloTree {
  input {
    String docker_image_id
    String s3_wd_uri
    String superkingdom_name
    Array[SampleInfo] samples
    Array[File] taxon_fastas
    Int taxid
    Array[Int] reference_taxids
    File nt_loc_db
    String nt_db
  }
  command<<<
  set -euxo pipefail
  idseq-dag-run-step --workflow-name host_filter \
    --step-module idseq_dag.steps.generate_phylo_tree \
    --step-class PipelineStepGeneratePhyloTree \
    --step-name generate_phylo_tree \
    --input-files '[["~{sep('", "', taxon_fastas)}"], ["~{write_json(samples)}"]]' \
    --output-files '["phylo_tree.newick", "ncbi_metadata.json"]' \
    --output-dir-s3 '~{s3_wd_uri}' \
    --additional-attributes '{"superkingdom_name": "~{superkingdom_name}", "taxid": ~{taxid}, "reference_taxids": [~{sep(", ", prefix("", reference_taxids))}], "nt_loc_db": "~{nt_loc_db}", "nt_db": "~{nt_db}"}'
  >>>
  output {
    File phylo_tree_newick = "phylo_tree.newick"
    File ncbi_metadata_json = "ncbi_metadata.json"
    Directory ksnp3_outputs = "ksnp3_outputs/"
  }
  runtime {
    docker: docker_image_id
  }
}

workflow idseq_generate_phylo_tree {
  input {
    String docker_image_id
    String s3_wd_uri
    String superkingdom_name
    Array[SampleInfo] samples
    Int taxid
    Array[Int] reference_taxids
    File nt_loc_db
    String nt_db
  }

  call PrepareTaxonFasta {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      samples = samples,
  }

  call GeneratePhyloTree {
    input:
      docker_image_id = docker_image_id,
      s3_wd_uri = s3_wd_uri,
      superkingdom_name = superkingdom_name,
      samples = samples,
      taxid = taxid,
      reference_taxids = reference_taxids,
      nt_loc_db = nt_loc_db,
      nt_db = nt_db,
      taxon_fastas = PrepareTaxonFasta.taxon_fastas,
  }

  output {
    Array[File] taxon_fastas = PrepareTaxonFasta.taxon_fastas
    File phylo_tree_newick = GeneratePhyloTree.phylo_tree_newick
    File ncbi_metadata_json = GeneratePhyloTree.ncbi_metadata_json
    Directory ksnp3_outputs = GeneratePhyloTree.ksnp3_outputs
  }
}
