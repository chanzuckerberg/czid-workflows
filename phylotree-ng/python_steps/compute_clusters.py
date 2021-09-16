import argparse
import json
from typing import Iterable, TypedDict

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import rcParams
from scipy import cluster

matplotlib.use('Agg')


class Sample(TypedDict):
    sample_name: str
    workflow_run_id: int
    contig_fasta: str
    combined_contig_summary: str


def main(ska_distances: str, trim_height: float, samples: Iterable[Sample], output_clusters_dir: str):
    # we have observed some strange parsing behavior of this file, ensure it works with end to end testing
    # we may need to use a regex separator
    df = pd.read_csv(ska_distances, sep='\t')
    df.reset_index(drop=True, inplace=True)
    df.fillna(1, inplace=True)  # fill NA values for SNP and Mash-like Distance with 1 (maximum dist)

    # long dataframe to wide
    df2 = df.pivot_table(index=['Sample 1'], columns='Sample 2', values='Mash-like distance')
    df2.columns = [str(i) for i in df2.columns]  # ensure that columns and index values are strings
    df2.index = [str(i) for i in df2.index]
    index = df2.index.union(df2.columns)
    df2 = df2.reindex(index=index, columns=index, fill_value=0)
    df2 = df2[df2.index]  # re-order columns to match index order

    # fill lower triangle of the matrix (currently NA)
    npdf = df2.to_numpy()
    i_lower = np.tril_indices(len(df2.index), -1)
    npdf[i_lower] = npdf.T[i_lower]
    df3 = pd.DataFrame(npdf)
    df3.columns = df2.columns
    df3.index = df2.index
    df3.fillna(0, inplace=True)

    Z = cluster.hierarchy.linkage(1-df3, method='average')
    cutree = cluster.hierarchy.cut_tree(Z, height=float(trim_height))
    ordered_clusterids = [i[0] for i in cutree]
    cluster_assignments = dict(zip(df3.index, ordered_clusterids))
    n_clusters = len(set(ordered_clusterids))

    cluster_sets = dict(zip(list(set(ordered_clusterids)), [[] for _ in range(n_clusters)]))

    for i in cluster_assignments.keys():
        cluster_sets[cluster_assignments[i]].append(i)

    # write cluster contents to files for future processing
    n_clusters_out = 0
    for c, s in cluster_sets.items():
        if(len(s)) > 2:  # only output files where there are > 2 samples
            n_clusters_out += 1
            filenames = '\n'.join(cluster_sets[c])
            with open(f"{output_clusters_dir}/cluster_{str(c)}", "w") as text_file:
                text_file.write(filenames)

    sample_name_by_workflow_id = {str(s["workflow_run_id"]): s["sample_name"] for s in samples}
    df3.columns = [sample_name_by_workflow_id.get(a, a) for a in df3.columns]
    df3.index = [sample_name_by_workflow_id.get(a, a) for a in df3.index]

    # Use Open Sans for the labels
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Open Sans']
    # 1.0 is baseline font size
    sns.set(font_scale=1.15)

    chart = sns.clustermap(
        df3,
        # cbar means colorbar. It's the key/legend for the colors.
        cbar_kws={'orientation': 'horizontal'},
        # (left, bottom, width, height) numbers as ratios of the figure.
        cbar_pos=(1, 1, 0.2, 0.02),
        cmap='YlOrRd_r',
        col_linkage=Z,
        # Higher number means larger dendrograms.
        dendrogram_ratio=0.1,
        figsize=(15, 15),
        # The lines are for the cell borders.
        linecolor='w',
        linewidth=5,
        row_linkage=Z,
        # Set dendogram color.
        tree_kws={'colors': "#767676"},
        vmin=0,
        vmax=0.15,
    )

    heatmap = chart.ax_heatmap

    # make the color bar last tick >.15
    cbar = heatmap.collections[0].colorbar
    cbar.set_ticks([0, 0.05, 0.1, 0.15])
    cbar.set_ticklabels(["0", "0.05", "0.1", ">0.15"])

    # Make the x tick labels slanted:
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
    # Hide the actual tick marks:
    heatmap.tick_params(bottom=False, right=False)

    plt.savefig('clustermap.png', bbox_inches='tight', dpi=300)
    plt.savefig('clustermap.svg', bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ska-distances")
    parser.add_argument("--cut-height", type=float)
    parser.add_argument("--samples")
    parser.add_argument("--output-clusters-dir")
    args = parser.parse_args()

    with open(args.samples) as f:
        samples: Iterable[Sample] = json.load(f)

    main(args.ska_distances, args.cut_height, samples, args.output_clusters_dir)
