import json
import os
from csv import DictReader

from Bio import SeqIO
from Bio.Phylo import NewickIO

from test_util import WDLTestCase


class TestPhylotree(WDLTestCase):
    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    samples_dir = os.path.join(os.path.dirname(__file__), "full_zika_test_data")

    samples = {}
    for i, name in enumerate(os.listdir(samples_dir)):
        first_dot = name.find(".")
        sample_name = name[:first_dot]
        sample = samples.get(sample_name, {
            "workflow_run_id": i,
            "sample_name": sample_name,
        })
        if "contig_summary" in name:
            sample["combined_contig_summary"] = os.path.join(samples_dir, name)
        else:
            sample["contig_fasta"] = os.path.join(samples_dir, name)
        samples[sample_name] = sample

    accession_ids = ["NC_012532.1", "NC_035889.1"]

    common_inputs = {
        "samples": list(samples.values()),
        "reference_taxon_id": 64320,
        "additional_reference_accession_ids": accession_ids,
        "superkingdom_name": "viruses"
    }

    def test_phylotree(self):
        sample_names = [s["sample_name"] for s in self.common_inputs["samples"]]

        res = self.run_miniwdl()
        outputs = res["outputs"]

        self.assertCountEqual(outputs.keys(), [
            "phylotree.clustermap_png",
            "phylotree.clustermap_svg",
            "phylotree.ncbi_metadata_json",
            "phylotree.phylotree_newick",
            "phylotree.ska_distances",
            "phylotree.variants",
        ])

        with open(outputs["phylotree.phylotree_newick"]) as f:
            tree = next(NewickIO.parse(f))
            nodes = [n.name for n in tree.get_terminals() + tree.get_nonterminals() if n.name]
            self.assertCountEqual(nodes, sample_names + self.accession_ids)

        identifiers = sorted(sample_names + self.accession_ids)
        with open(outputs["phylotree.ska_distances"]) as f:
            pairs = [sorted([r["Sample 1"], r["Sample 2"]]) for r in DictReader(f, delimiter="\t")]
            expected = [[a, b] for a in identifiers for b in identifiers if a < b]
            self.assertCountEqual(pairs, expected)

        with open(outputs["phylotree.variants"]) as f:
            self.assertCountEqual(identifiers, [r.id for r in SeqIO.parse(f, "fasta")])

        with open(outputs["phylotree.ncbi_metadata_json"]) as f:
            self.assertEqual(json.load(f), {
                "NC_012532.1": {
                    "name": "Zika virus, complete genome",
                    "country": "Uganda",
                },
                "NC_035889.1": {
                    "name": "Zika virus isolate ZIKV/H. sapiens/Brazil/Natal/2015, complete genome",
                    "country": "Brazil: Rio Grande do Norte, Natal",
                    "collection_date": "2015",
                },
            })

        with open(outputs["phylotree.clustermap_svg"]) as f:
            full_text = "\n".join(f.readlines())
            for name in sample_names + self.accession_ids:
                self.assertEqual(full_text.count(name), 2, name)
