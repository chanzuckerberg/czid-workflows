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
import dbm
import logging
import os
import shelve
import sys
from multiprocessing.pool import ThreadPool


def output_dicts_to_db(mapping_files, accession2taxid_db):
    # generate the accession2taxid db and file
    accession_dict = shelve.Shelf(dbm.ndbm.open(accession2taxid_db.replace(".db", ""), 'c'))  # type: ignore
    for partition_list in mapping_files:
        for partition in partition_list:
            with open(partition, 'r', encoding="utf-8") as pf:
                for line in pf:
                    if len(line) <= 1:
                        break
                    fields = line.rstrip().split("\t")
                    accession_dict[fields[0]] = fields[2]

    accession_dict.close()


def grab_accession_names(source_file, dest_file):
    with open(source_file) as source:
        with open(dest_file, 'w') as dest:
            for line in source:
                if line[0] == '>':
                    dest.write(line.split(' ')[0] + "\n")


def grab_accession_mapping_list(source, num_partitions, partition_id,
                                accessions, output_file):
    num_lines = 0
    with open(output_file, 'w') as out, open(source, 'r') as mapf:
        for line in mapf:
            if num_lines % num_partitions == partition_id:
                accession_line = line.split("\t")
                accession = accession_line[0]
                # If using the prot.accession2taxid.FULL file, should add a column with no version
                if len(accession_line) < 3 and accession.split(".")[0] in accessions:
                    accession_no_version = accession.split(".")[0]
                    out.write(f"{accession_no_version}\t{line}")
                elif accession in accessions:
                    out.write(line)
            num_lines += 1
            if num_lines % 1000000 == 0:
                print(f"{source} partition {partition_id} line {num_lines/1000000}M", file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate accession2taxid.')
    parser.add_argument('accession_mapping_files', nargs=4)
    parser.add_argument('--parallelism', type=int)
    parser.add_argument('--nt_file')
    parser.add_argument('--nr_file')
    parser.add_argument('--output_gz')
    parser.add_argument('--accession2taxid_db')
    args = parser.parse_args()

    accession_mapping_files = args.accession_mapping_files
    num_partitions = args.parallelism
    nt_file = args.nt_file
    nr_file = args.nr_file
    accession2taxid_db = args.accession2taxid_db

    # Get accession_list
    source_files = [nt_file, nr_file]
    accessions_files = [f"{source_file}.accessions" for source_file in source_files]
    pool = ThreadPool()
    pool.starmap(grab_accession_names, zip(source_files, accessions_files))

    accessions = set()
    for accession_file in accessions_files:
        with open(accession_file, 'r') as acf:
            for line in acf:
                accession = line[1:].split(".")[0]
                accessions.add(accession)

    grab_accession_mapping_list_args = []
    mapping_files = []
    for accession_mapping_file in accession_mapping_files:
        partition_list = []
        for p in range(num_partitions):
            part_file = f"{os.path.basename(accession_mapping_file)}-{p}"
            partition_list.append(part_file)
            grab_accession_mapping_list_args.append([
                accession_mapping_file,
                num_partitions,
                p,
                accessions,
                part_file
            ])
        mapping_files.append(partition_list)

    pool.starmap(grab_accession_mapping_list, grab_accession_mapping_list_args)
    accessions = set()  # reset accessions to release memory

    logging.info("starting writing output dictionaries to db")

    output_dicts_to_db(mapping_files, accession2taxid_db)
