# IDseq PhyloTree-NG workflow

version 1.1

struct SampleInfo {
    String sample_name
    Int workflow_run_id
}

struct ReferenceInfo {
    Int? taxon_id
    String? accession_id
}

workflow phylotree {
    input {
        Array[SampleInfo] samples
        ReferenceInfo reference

        # TODO: pass this to the relevant tasks and adjust SKA parameters as appropriate
        # (kSNP3 used a variable kmer length for each superkingdom: Viruses: 13, Bacteria: 19, Eukaryota: 19)
        String superkingdom_name

        # allow the user to pass specific reference taxids/accessions to include with the tree
        Array[ReferenceInfo] additional_references

        String cut_height = .16
        String ska_align_p = .9
        String docker_image_id
        String outgroup = ""
    }

    call GetSampleContigFastas {
        input:
        samples = samples,
        docker_image_id = docker_image_id
    }

    call GetReferenceFastas {
        input:
        references = additional_references,
        docker_image_id = docker_image_id
    }

    call RunSKA {
        input:
        sample_and_reference_fastas = flatten([GetSampleContigFastas.sample_contig_fastas, GetReferenceFastas.reference_fastas]),
        docker_image_id = docker_image_id
    }

    call ComputeClusters {
        input:
        ska_distances = RunSKA.distances,
        cut_height = cut_height,
        docker_image_id = docker_image_id
    }

    call GenerateClusterPhylos {
        input:
        clusters_directory = ComputeClusters.clusters_directory,
        ska_hashes = RunSKA.ska_hashes,
        ska_align_p = ska_align_p,
        docker_image_id = docker_image_id
    }

    call PlotClusterPhylos {
        input:
        ska_results = GenerateClusterPhylos.ska_results,
        outgroup = outgroup,
        docker_image_id = docker_image_id
    }

    output {
        File ska_hashes = RunSKA.ska_hashes
        File ska_distances = RunSKA.distances
        File stats_json = ComputeClusters.stats_json
        File dendrogram_png = ComputeClusters.dendrogram_png
        File clustermap_png = ComputeClusters.clustermap_png
        File clusters_directory = ComputeClusters.clusters_directory
        File dummy_subcluster_output = GenerateClusterPhylos.dummy_subcluster_output
        File ska_results = GenerateClusterPhylos.ska_results
        File dummy_plotphylos = PlotClusterPhylos.dummy_plotphylos
        File phylo_plot_outputs = PlotClusterPhylos.phylo_plot_outputs

        # TODO: These are the output names from the old phylotree. Check which of these we still need to emit
        # Array[File] taxon_fastas = PrepareTaxonFasta.taxon_fastas
        # File phylo_tree_newick = GeneratePhyloTree.phylo_tree_newick
        # File ncbi_metadata_json = GeneratePhyloTree.ncbi_metadata_json
    }
}

task GetSampleContigFastas {
    # Given a list of samples, their workflow run IDs, and a taxon ID, retrieve contigs assigned to that taxon ID or its
    # descendants within each sample. Emit one fasta file per sample, to be graphed in a phylotree along with references
    input {
        Array[SampleInfo] samples
        String docker_image_id
    }

    command <<<

    >>>

    output {
        Array[File] sample_contig_fastas = glob("*.fasta")
    }

    runtime {
        docker: docker_image_id
    }
}

task GetReferenceFastas {
    # Given a list of reference taxon IDs or sequence accession IDs, retrieve the GenBank nt reference sequences
    # associated with them. Emit one fasta file per reference, to be graphed in a phylotree along with samples
    input {
        Array[ReferenceInfo] references
        String docker_image_id
    }

    command <<<
    for accession_id in $(jq -r .[].accession_id "~{write_json(references)}"); do
        taxoniq get-from-s3 --accession-id $accession_id > $accession_id.fasta
    done
    >>>

    output {
        Array[File] reference_fastas = glob("*.fasta")
    }

    runtime {
        docker: docker_image_id
    }
}

task RunSKA {
    input {
        Array[File] sample_and_reference_fastas
        String docker_image_id
    }

    command <<<
    for i in ~{sep=' ' sample_and_reference_fastas}; do
        ska fasta $i
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
        String docker_image_id
    }

    command <<<
    mkdir cluster_files

    python3 <<CODE
    import pandas as pd
    import json
    import numpy as np
    from scipy import cluster
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    import math

    def processed_lines():
        with open("~{ska_distances}", 'r') as distances:
            for i, line in enumerate(distances):
                split_line = line.split()
                if i == 1:
                    # on the first line, add entry to both columns to ensure square distance matrix
                    yield [split_line[0], split_line[0], 0, 0, 0, 0, 0, 0]
                elif i > 1: # 0 is purposefully omitted
                    yield split_line
           # on the last line, add entry to both columns
            yield [split_line[1], split_line[1], 0, 0, 0, 0, 0, 0]

    # create dataframe from the distance data
    df = pd.DataFrame.from_dict(dict(enumerate(processed_lines())), orient='index')
    df.columns=['Sample_1','Sample_2','Matches','Mismatches','Jaccard_Index','Mash-like_distance','SNPs','SNP_distance']
    df["Mash-like_distance"] = pd.to_numeric(df["Mash-like_distance"], downcast="float")

    # long dataframe to wide
    df2 = df.pivot_table(index=['Sample_1'], columns='Sample_2', values='Mash-like_distance')

    # fill lower triangle of the matrix (currently NA)
    npdf = df2.to_numpy()
    i_lower = np.tril_indices(len(df2.index), -1)
    npdf[i_lower] = npdf.T[i_lower]
    df3 = pd.DataFrame(npdf)
    df3.columns = df2.columns
    df3.index = df2.index
    df3.fillna(0, inplace=True)

    # cluster the data, create a dendrogram
    Z = cluster.hierarchy.linkage(1-df3, method='complete')
    dn = cluster.hierarchy.dendrogram(Z, leaf_rotation=90, labels=df3.index)
    plt.savefig('dendrogram.png', bbox_inches='tight')

    trim_height = "~{cut_height}"
    cutree = cluster.hierarchy.cut_tree(Z, height=float(trim_height))
    ordered_clusterids = [i[0] for i in cutree]
    cluster_assignments = dict(zip(df3.index, ordered_clusterids))
    n_clusters = len(set(ordered_clusterids))
    cluster_sets = dict(zip(list(set(ordered_clusterids)), [[] for i in range(n_clusters)]))

    for i in cluster_assignments.keys():
        cluster_sets[cluster_assignments[i]].append(i)

    stats = {"sample_name": "insert_sample_name"}

    # write cluster contents to files for future processing
    for c in cluster_sets.keys():
        stats[c] = ' '.join(cluster_sets[c]) # record cluster IDs in stats file for all clusters
        if(len(cluster_sets[c])) > 2: # only output files where there are > 2 samples
            filenames = '\n'.join(cluster_sets[c])
            with open("./cluster_files/cluster_" + str(c), "w") as text_file:
                text_file.write(filenames)

    color_list = sns.color_palette("Dark2", 8)
    long_color_list = color_list*math.ceil(len(set(ordered_clusterids))/len(color_list))
    col_colors = [long_color_list[i] for i in ordered_clusterids]
    h = sns.clustermap(df3, cmap='coolwarm_r', vmin = 0, vmax = 0.15, col_linkage = Z, col_colors = col_colors, figsize=(15,15))
    plt.savefig('clustermap.png', bbox_inches='tight')

    stats["dataframe_shape_0"] = df.shape[0]
    stats["dataframe_shape_1"] = df.shape[1]

    with open("stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    CODE

    tar -czf clusters.tar.gz cluster_files
    >>>

    output {
        File stats_json = "stats.json"
        File clusters_directory = "clusters.tar.gz"
        File dendrogram_png = "dendrogram.png"
        File clustermap_png = "clustermap.png"
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
    echo "hello world" > output.txt
    tar -xzvf "~{clusters_directory}"
    tar -xzvf "~{ska_hashes}"

    ls . >> output.txt
    echo "unzipped files" >> output.txt
    ls cluster_files/* >> output.txt
    ls ska_hashes/* >> output.txt

    CLUSTER_COUNTER=1
    mkdir ska_outputs

    for i in `ls cluster_files/*`
    do
        rm -r temp_cluster # remove if already existed
        mkdir temp_cluster

        for j in `cat $i`
        do
            cp ska_hashes/$j.skf temp_cluster
        done

        ska distance -o ska temp_cluster/*.skf
        ska merge -o ska.merged temp_cluster/*.skf
        ska align -p "~{ska_align_p}" -o ska -v ska.merged.skf
        mv ska_variants.aln ska.variants.aln
        iqtree -s ska.variants.aln

        mkdir ska_outputs/cluster_$CLUSTER_COUNTER
        mv ska.* ska_outputs/cluster_$CLUSTER_COUNTER

        (( CLUSTER_COUNTER++ ))
    done

    tar -czf ska_outputs.tar.gz ska_outputs
    >>>

    output {
        File dummy_subcluster_output = "output.txt"
        File ska_results = "ska_outputs.tar.gz"
    }

    runtime {
        docker: docker_image_id
    }
}

task PlotClusterPhylos {
    input {
        File ska_results
        String outgroup
        String docker_image_id
    }

    command <<<
    tar -xzvf "~{ska_results}"

    ls > output.txt
    ls ska_outputs >> output.txt

    mkdir phylo_plot_outputs

    python3 <<CODE
        import toytree
        import glob
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import toyplot
        import toyplot.pdf

        tree_file_paths = glob.glob('ska_outputs/*/*.treefile')
        print(tree_file_paths)

        for tree_file in tree_file_paths:
            file = open(tree_file, 'r')
            lines = file.readlines()
            newick = lines[0].strip()

            # set the outgroup to first sample unless otherwise specified
            # TODO: create a more robust method for picking outgroup...
            # ... the current method of specifying it will break / error in cases where there
            # ... is > 1 cluster
            wildcard_value = newick.split('.')[0].split('(')[1]
            if("~{outgroup}" == ""):
                wildcard_value = newick.split('.')[0].split('(')[1]
            else:
                wildcard_value = "~{outgroup}"

            tre0 = toytree.tree(newick, tree_format=0)
            rtre = tre0.root(wildcard=wildcard_value)

            style = {
                "tip_labels_align": False,
                "tip_labels_style": {
                    "font-size": "9px"
                },
            }
            canvas, axes, makr = rtre.draw(tip_labels_colors='indigo', **style);
            axes.show = True
            axes.x.ticks.show = True
            axes.y.ticks.show = False
            print("about to create plot")
            print('phylo_plot_outputs/' + tree_file.split('/')[1] + '.treeplot.pdf')
            toyplot.pdf.render(canvas, 'phylo_plot_outputs/' + tree_file.split('/')[1] + '.treeplot.pdf')
            #plt.savefig('phylo_plot_outputs/' + tree_file.split('/')[1] + '.treeplot.png', bbox_inches='tight')
        CODE

    tar -czf phylo_plot_outputs.tar.gz phylo_plot_outputs
    >>>

    output {
        File dummy_plotphylos = "output.txt"
        File phylo_plot_outputs = "phylo_plot_outputs.tar.gz"
    }

    runtime {
        docker: docker_image_id
    }
}

task FetchSequenceByAccessionId {
    input {
        String accession_id
        String docker_image_id
    }

    command <<<
        taxoniq get-from-s3 --accession-id "~{accession_id}" > sequence.fa
        if [[ $? == 4 ]]; then
            export error=AccessionIdNotFound cause="Accession ID ~{accession_id} not found in the index"
            jq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr
            exit 4
        fi
        exit $?
    >>>

    output {
        File sequence_fa = "sequence.fa"
    }

    runtime {
        docker: docker_image_id
    }
}
