version 1.1

struct SampleInfo {
    String sample_name
    Int workflow_run_id
    File contig_fasta # assembly_out_assembly_contigs_fasta from mngs
    File combined_contig_summary # contig_summary_out_assembly_combined_contig_summary_json from mngs
}

workflow phylotree {
    input {
        Array[SampleInfo] samples
        Int reference_taxon_id

        String superkingdom_name # viruses, bacteria, or eukaryota

        # allow the user to pass specific refaccessions to include with the tree
        Array[String] additional_reference_accession_ids = []

        String cut_height = .16
        String ska_align_p = .9
        String docker_image_id

        # Dummy values - required by SFN interface
        String s3_wd_uri = ""
    }

    call GetSampleContigFastas {
        input:
        reference_taxon_id = reference_taxon_id,
        samples = samples,
        docker_image_id = docker_image_id
    }

    call GetReferenceAccessionFastas {
        input:
        accession_ids = additional_reference_accession_ids,
        docker_image_id = docker_image_id
    }

    call RunSKA {
        input:
        sample_and_reference_fastas = flatten([GetSampleContigFastas.sample_contig_fastas, GetReferenceAccessionFastas.reference_fastas]), 
        superkingdom_name = superkingdom_name,
        docker_image_id = docker_image_id
    }

    call ComputeClusters {
        input:
        ska_distances = RunSKA.distances,
        cut_height = cut_height,
        samples = samples,
        docker_image_id = docker_image_id
    }

    call GenerateClusterPhylos {
        input:
        clusters_directory = ComputeClusters.clusters_directory,
        ska_hashes = RunSKA.ska_hashes,
        ska_align_p = ska_align_p,
        docker_image_id = docker_image_id
    }

    call AddSampleNamesToDistances {
        input:
        distances = RunSKA.distances,
        samples = samples,
        docker_image_id = docker_image_id
    }

    call AddSampleNamesToVariants {
        input:
        variants = GenerateClusterPhylos.variants,
        samples = samples,
        docker_image_id = docker_image_id
    }

    if (GenerateClusterPhylos.phylotree_newick != None) {
      call AddSampleNamesToNewick {
        input:
        newick = select_first([GenerateClusterPhylos.phylotree_newick]),
        samples = samples,
        docker_image_id = docker_image_id
      }
    }

    call FetchNCBIMetadata {
        input:
        reference_accession_ids = additional_reference_accession_ids,
        docker_image_id = docker_image_id
    }

    output {
        File ska_distances = AddSampleNamesToDistances.sample_name_distances
        File clustermap_png = ComputeClusters.clustermap_png
        File clustermap_svg = ComputeClusters.clustermap_svg
        File? phylotree_newick = AddSampleNamesToNewick.phylotree_newick
        File variants = AddSampleNamesToVariants.sample_name_variants
        File ncbi_metadata_json = FetchNCBIMetadata.ncbi_metadata_json
    }
}

task GetReferenceAccessionFastas {
    input {
        Array[String] accession_ids
        String docker_image_id
    }

    command <<<
    set -euxo pipefail
        for accession_id in ~{sep=' ' accession_ids}; do
            taxoniq get-from-s3 --accession-id $accession_id > $accession_id.fasta
            if [[ $? == 4 ]]; then
                export error=AccessionIdNotFound cause="Accession ID $accession_id not found in the index"
                jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
                exit 4
            fi
        done
    >>>

    output {
        Array[File] reference_fastas = glob("*.fasta")
    }

    runtime {
        docker: docker_image_id
    }
}

task GetSampleContigFastas {
    # Given a list of samples, their workflow run IDs, and a reference taxon ID, retrieve contigs assigned to that taxon
    # ID or its descendants within each sample. Emit one fasta file per sample, to be graphed in a phylotree along with
    # references.

    # For each input sample, scan its contig summary to identify the contigs that need to be pulled from contigs.fasta.
    # TODO: For now, we do exact matching of the given reference taxon ID to the taxon IDs in the contig summary.
    #       This implies that the reference is selected at the species level, and that the contig summary is rolled up
    #       to species level as well. In the future we want to remove this assumption and select all contigs that map
    #       at all ranks under the given reference taxon.
    input {
        Int reference_taxon_id
        Array[SampleInfo] samples
        String docker_image_id
    }

    command <<<
    python3 /bin/get_sample_contig_fastas.py --reference-taxid ~{reference_taxon_id} --samples "~{write_json(samples)}"
    >>>

    output {
        Array[File] sample_contig_fastas = glob("*.fasta")
    }

    runtime {
        docker: docker_image_id
    }
}

task RunSKA {
    input {
        Array[File] sample_and_reference_fastas
        String superkingdom_name
        String docker_image_id
    }

    command <<<
    set -euxo pipefail
    k=18
    if [[ "~{superkingdom_name}" == viruses ]]; then
        k=12
    fi

    for i in ~{sep=' ' sample_and_reference_fastas}; do
        ska fasta -o $(basename $i | sed 's/\.fasta//g') -k $k $i
    done

    mkdir ska_hashes
    mv *.skf ska_hashes

    ska distance -o ska ska_hashes/*.skf

    tar -czf ska_hashes.tar.gz ska_hashes
    >>>

    output {
        File distances = "ska.distances.tsv"
        File ska_hashes = "ska_hashes.tar.gz"
    }

    runtime {
        docker: docker_image_id
    }
}

task ComputeClusters {
    input {
        File ska_distances
        String cut_height
        Array[SampleInfo] samples
        String docker_image_id
    }

    command <<<
    set -euxo pipefail
    mkdir cluster_files
    python3 /bin/compute_clusters.py \
        --ska-distances ~{ska_distances} \
        --cut-height ~{cut_height} \
        --samples "~{write_json(samples)}" \
        --output-clusters-dir cluster_files
    tar -czf clusters.tar.gz cluster_files
    >>>

    output {
        File clusters_directory = "clusters.tar.gz"
        File clustermap_png = "clustermap.png"
        File clustermap_svg = "clustermap.svg"
    }

    runtime {
        docker: docker_image_id
    }
}

task GenerateClusterPhylos {
    input {
        File clusters_directory
        File ska_hashes
        String ska_align_p
        String docker_image_id
    }

    command <<<
    # -e omitted because iqtree can exit with a non-zero exit code despite
    #   producing valid output
    set -uxo pipefail
    tar -xzvf "~{clusters_directory}"
    tar -xzvf "~{ska_hashes}"

    ska distance -o ska ska_hashes/*.skf
    ska merge -o ska.merged ska_hashes/*.skf
    ska align -p "~{ska_align_p}" -o ska -v ska.merged.skf
    mv ska_variants.aln ska.variants.aln

    if [[ $(ls cluster_files | wc -l) -gt 1 ]]; then
      # If we have more than one cluster then the samples are too
      #   divergent and we should not generate a tree
      exit 0
    fi


    num_zero_variant_samples=`grep -c "^$" ska.variants.aln`
    if [[ $num_zero_variant_samples -gt 0 ]]; then
      # If we have more than one 0 variant samples then the samples
      #   are too divergent and we should not generate a tree
      exit 0
    fi

    iqtree -s ska.variants.aln
    mv ska.variants.aln.treefile phylotree.nwk
    >>>

    output {
        File? phylotree_newick = "phylotree.nwk"
        File variants = "ska.variants.aln"
    }

    runtime {
        docker: docker_image_id
    }
}

task AddSampleNamesToDistances {
    input {
        File distances
        Array[SampleInfo] samples
        String docker_image_id
    }

    command <<<
    python3 /bin/add_sample_names_to_distances.py \
        --distances ~{distances} \
        --samples "~{write_json(samples)}" \
        --output-distances ska.distances.tsv
    >>>

    output {
        File sample_name_distances = "ska.distances.tsv"
    }

    runtime {
        docker: docker_image_id
    }
}

task AddSampleNamesToVariants {
    input {
        File variants
        Array[SampleInfo] samples
        String docker_image_id
    }

    command <<<
    python3 /bin/add_sample_names_to_variants.py \
        --variants ~{variants} \
        --samples "~{write_json(samples)}" \
        --output-variants ska.variants.aln
    >>>

    output {
        File sample_name_variants = "ska.variants.aln"
    }

    runtime {
        docker: docker_image_id
    }
}

task AddSampleNamesToNewick {
    input {
        File newick
        Array[SampleInfo] samples
        String docker_image_id
    }

    command <<<
    python3 /bin/add_sample_names_to_newick.py \
        --newick ~{newick} \
        --samples "~{write_json(samples)}" \
        --output-newick phylotree.nwk
    >>>

    output {
        File phylotree_newick = "phylotree.nwk"
    }

    runtime {
        docker: docker_image_id
    }
}

task FetchNCBIMetadata {
    input {
        Array[String] reference_accession_ids
        String docker_image_id
    }

    command <<<
    python3 /bin/fetch_ncbi_metadata.py \
        --reference-accession-ids "~{write_json(reference_accession_ids)}" \
        --output-ncbi-metadata ncbi_metadata.json
    >>>

    output {
        File ncbi_metadata_json = "ncbi_metadata.json"
    }

    runtime {
        docker: docker_image_id
    }
}
