"""
Accession2Taxid
After alignment, IDseq uses the NCBI accession2taxid database to map accessions to taxonomic IDs.
The full database contains billions of entries, but only ~15% of those are found in either NR or NT databases.
Therefore, the full NCBI accession2taxid database is subsetted to include only the relevant entries and then used
to map accessions to taxonomic IDs.
Finally, the taxonomic lineage for each read is computed using the
[ncbitax2lin](https://github.com/chanzuckerberg/ncbitax2lin) script.
For each taxonomic ID this results in the following: taxid → [superkingdom, phylum, class, order ..., species].
__Final Read Assignment__
For any given read, the tax ID is assigned by the following steps:
1. If the read was assembled into a contig and the contig maps to an accession via blast,
    then the read is assigned the taxID of its respective contig.
2. If the read was not assembled into a contig, then it is assigned to the taxID identified
    by short read alignment (NT with GSNAP and NR with Rapsearch2).
The progression of read assignment (from initial short read TaxID to refined assembly-based TaxID) can be identified
in the hitsummary2 files (gsnap.hitsummary2.tab and rapsearch2.hitsummary2.tab)
```
NB501961:211:HGWKCBGX9:1:11205:5935:16320/1 1   155900  GQ881617.1  155900  -200    -300
NODE_1497_length_451_cov_1335.782828    GQ881617.1  155900.0    -200.0  -300.0
```
The .tab file contains 12 columns:
1. Read ID
2. Taxonomy level (all -1)
3. Final TaxID assignment
4. Initial (single short-read) GenBank alignment
5. TaxID (species), from single alignment - obtained from the accession2taxid mapping after GSNAP or Rapsearch2.
    *note: *the pipeline outputs species-level counts as well as genus-level counts. For rows corresponding to
    genus-level counts (tax_level = 2), the species taxID is listed as “-100”.
6. TaxID (genus), from single alignment - obtained by walking the phylogenetic tree backwards. If there is no
    genus-level classification, this value will be “-200”.
7. TaxID (family), from single alignment - obtained by walking the phylogenetic tree backwards. If there is no
    family-level classification, this value will be “-300”.
8. Assembled contig that this read maps to
9. The GenBank ID that the contig mapped to
10. TaxID (species), from assembled contig - obtained from the accession2taxid mapping after BLAST
11. TaxID (genus), from assembled contig
12. TaxID (family), from assembled contig
13. “from_assembly” - this flag indicates whether this read was mapped ONLY through assembly. If this flag is present,
    then fields (4) - (7) are duplicates of the assembly-based TaxID call, and are not single short-read alignments.
__note:__ columns 3 and 10 of the hitsummary2.xxx.tab files should always be identical.
1. Download NT/NR
2. Extract accessions
3. Subsetting accessions to the ones appearing in NT/NR
"""
import argparse
import sys

import marisa_trie


def accession_mapping(accessions: marisa_trie.Trie, *accession_mapping_files: str):
    for accession_mapping_file in accession_mapping_files:
        with open(accession_mapping_file) as mapf:
            for i, _line in enumerate(mapf):
                line = _line

                accession_line = line.split("\t")
                accession = accession_line[0]
                # If using the prot.accession2taxid.FULL file, should add a column with no version
                if len(accession_line) < 3 and accession.split(".")[0] in accessions:
                    accession_no_version = accession.split(".")[0]
                    line = f"{accession_no_version}\t{line}"

                fields = line.rstrip().split("\t")
                yield (fields[0], (int(fields[2]),))

                if i % 1000000 == 0:
                    print(f"{accession_mapping_file} line {i/1000000}M", file=sys.stderr)


def get_accessions(*source_files: str):
    for source_file in source_files:
        with open(source_file) as source:
            for line in source:
                if line[0] == '>':
                    accession = line[1:].split(' ')[0].split(".")[0]
                yield accession


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate accession2taxid.')
    parser.add_argument('accession_mapping_files', nargs=4)
    parser.add_argument('--nt_file')
    parser.add_argument('--nr_file')
    parser.add_argument('--output_gz')
    parser.add_argument('--accession2taxid_db')
    args = parser.parse_args()

    accession_mapping_files = args.accession_mapping_files
    num_partitions = args.parallelism
    nt_file = args.nt_file
    nr_file = args.nr_file

    # Build accession trie
    accessions = marisa_trie.Trie(get_accessions(nt_file, nr_file))
    print("Building trie", file=sys.stderr)
    marisa_trie.RecordTrie(
        "L",
        accession_mapping(accessions, *args.accession_mapping_files),
    ).save(args.accession2taxid_db)
