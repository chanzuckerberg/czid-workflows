''' Accession2Taxid'''
import dbm
import gzip
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.log as log
import shelve
import threading

from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.command import run_in_subprocess
from idseq_dag.util.dict import IdSeqDictForUpdate, IdSeqDictValue

BATCH_INSERT_SIZE = 300
NUM_PARTITIONS = 8


class PipelineStepGenerateAccession2Taxid(PipelineStep):
    """ After alignment, IDseq uses the NCBI accession2taxid database to map accessions to taxonomic IDs.
    The full database contains billions of entries, but only ~15% of those are found in either NR or NT databases.
    Therefore, the full NCBI accession2taxid database is subsetted to include only the relevant entries and then used
    to map accessions to taxonomic IDs.
    Finally, the taxonomic lineage for each read is computed using the [ncbitax2lin](https://github.com/chanzuckerberg/ncbitax2lin) script.
    For each taxonomic ID this results in the following: taxid → [superkingdom, phylum, class, order ..., species].

    __Final Read Assignment__
    For any given read, the tax ID is assigned by the following steps:

    1. If the read was assembled into a contig and the contig maps to an accession via blast, then the read is assigned the taxID of its respective contig.
    2. If the read was not assembled into a contig, then it is assigned to the taxID identified by short read alignment (NT with GSNAP and NR with Rapsearch2).

    The progression of read assignment (from initial short read TaxID to refined assembly-based TaxID) can be identified in the hitsummary2 files (gsnap.hitsummary2.tab and rapsearch2.hitsummary2.tab)

    ```
    NB501961:211:HGWKCBGX9:1:11205:5935:16320/1 1   155900  GQ881617.1  155900  -200    -300    NODE_1497_length_451_cov_1335.782828    GQ881617.1  155900.0    -200.0  -300.0
    ```

    The .tab file contains 12 columns:

    1. Read ID
    2. Taxonomy level (all -1)
    3. Final TaxID assignment
    4. Initial (single short-read) GenBank alignment
    5. TaxID (species), from single alignment - obtained from the accession2taxid mapping after GSNAP or Rapsearch2. *note: *the pipeline outputs species-level counts as well as genus-level counts. For rows corresponding to genus-level counts (tax_level = 2), the species taxID is listed as “-100”.
    6. TaxID (genus), from single alignment - obtained by walking the phylogenetic tree backwards. If there is no genus-level classification, this value will be “-200”.
    7. TaxID (family), from single alignment - obtained by walking the phylogenetic tree backwards. If there is no family-level classification, this value will be “-300”.
    8. Assembled contig that this read maps to
    9. The GenBank ID that the contig mapped to
    10. TaxID (species), from assembled contig - obtained from the accession2taxid mapping after BLAST
    11. TaxID (genus), from assembled contig
    12. TaxID (family), from assembled contig
    13. “from_assembly” - this flag indicates whether this read was mapped ONLY through assembly. If this flag is present, then fields (4) - (7) are duplicates of the assembly-based TaxID call, and are not single short-read alignments.

    __note:__ columns 3 and 10 of the hitsummary2.xxx.tab files should always be identical.
    """

    def run(self):
        """
        1. Download NT/NR
        2. Extract accessions
        3. Subsetting accessions to the ones appearing in NT/NR
        Example command:
        aegea batch submit --command="pip install --upgrade git+git://github.com/chanzuckerberg/s3mi.git; cd /mnt; git clone https://github.com/chanzuckerberg/idseq-dag.git; cd idseq-dag;  pip3 install -e .; sed 's/<NCBI_DATE>/2018-12-01/g' examples/accession2taxid_dag.json > examples/accession2taxid_dag.json.now  idseq_dag --no-versioned-output examples/accession2taxid_dag.json.now"  --storage /mnt=1000 --volume-type gp2 --ecr-image idseq_dag --memory 240000 --queue idseq-prod-himem --vcpus 32 --job-role idseq-pipeline

        """
        accession_mapping_files = self.input_files_local[0]
        nt_file = self.input_files_local[1][0]
        nr_file = self.input_files_local[2][0]
        (
            output_gz,
            output_wgs_gz,
            accession2taxid_db_sqlite,
            taxid2wgs_accession_db_sqlite,
            accession2taxid_db,
            taxid2wgs_accession_db,
         ) = self.output_files_local()

        # Get accession_list
        accessions_files = []
        threads = []
        for source_file in [nt_file, nr_file]:
            accession_file = f"{source_file}.accessions"
            thread = threading.Thread(target=self.grab_accession_names,
                                      args=[source_file, accession_file])
            accessions_files.append(accession_file)
            threads.append(thread)
            thread.start()
        wgs_accessions = f"{nt_file}.wgs_acc"
        wgs_thread = threading.Thread(target=self.grab_wgs_accessions,
                                      args=[nt_file, wgs_accessions])
        wgs_thread.start()

        for t in threads:
            t.join()

        accessions = set()
        for accession_file in accessions_files:
            with open(accession_file, 'r') as acf:
                for line in acf:
                    accession = line[1:].split(".")[0]
                    accessions.add(accession)

        threads = []
        mapping_files = []
        for accession_mapping_file in accession_mapping_files:
            partition_list = []
            for p in range(NUM_PARTITIONS):
                part_file = f"{accession_mapping_file}-{p}"
                partition_list.append(part_file)
                thread = threading.Thread(target=self.grab_accession_mapping_list,
                                          args=[accession_mapping_file, NUM_PARTITIONS, p,
                                                accessions, part_file])
                accessions_files.append(accession_file)
                threads.append(thread)
                thread.start()
            mapping_files.append(partition_list)

        for t in threads:
            t.join()
        wgs_thread.join()
        accessions = []  # reset accessions to release memory

        self.output_dicts_to_db_for_sqlite(mapping_files, wgs_accessions,
                                           accession2taxid_db_sqlite, taxid2wgs_accession_db_sqlite,
                                           output_gz, output_wgs_gz)
        self.output_dicts_to_db_for_shelf(mapping_files, wgs_accessions,
                                          accession2taxid_db, taxid2wgs_accession_db,
                                          output_gz, output_wgs_gz)

    def output_dicts_to_db_for_sqlite(self, mapping_files, wgs_accessions,
                                      accession2taxid_db, taxid2wgs_accession_db,
                                      output_gz, output_wgs_gz):
        # generate the accession2taxid db and file
        with IdSeqDictForUpdate(accession2taxid_db, IdSeqDictValue.VALUE_TYPE_SCALAR) as accession_dict:
            batch_list = {}
            with gzip.open(output_gz, "wt") as gzf:
                for partition_list in mapping_files:
                    for partition in partition_list:
                        with open(partition, 'r', encoding="utf-8") as pf:
                            for line in pf:
                                if len(line) <= 1:
                                    break
                                fields = line.rstrip().split("\t")
                                gzf.write(line)
                                batch_list[fields[0]] = fields[2]
                                if len(batch_list) >= BATCH_INSERT_SIZE:
                                    accession_dict.batch_inserts(batch_list.items())
                                    batch_list = {}
                accession_dict.batch_inserts(batch_list.items())

            # generate taxid2 accession
            taxid2accession_dict = {}
            with gzip.open(output_wgs_gz, "wt") as gzf:
                with open(wgs_accessions, 'r', encoding="utf-8") as wgsf:
                    for line in wgsf:
                        accession = line[1:].split(".")[0]
                        taxid = accession_dict.get(accession)
                        if taxid:
                            current_match = taxid2accession_dict.get(taxid, "")
                            taxid2accession_dict[taxid] = f"{current_match},{accession}"
                            gzf.write(line)

        with IdSeqDictForUpdate(taxid2wgs_accession_db, IdSeqDictValue.VALUE_TYPE_SCALAR) as taxid2wgs_accession_dict:
            taxid2wgs_accession_dict.batch_inserts(taxid2accession_dict.items())

    def output_dicts_to_db_for_shelf(self, mapping_files, wgs_accessions,
                                     accession2taxid_db, taxid2wgs_accession_db,
                                     output_gz, output_wgs_gz):
        # generate the accession2taxid db and file
        accession_dict = shelve.Shelf(dbm.ndbm.open(accession2taxid_db.replace(".db", ""), 'c'))
        with gzip.open(output_gz, "wt") as gzf:
            for partition_list in mapping_files:
                for partition in partition_list:
                    with open(partition, 'r', encoding="utf-8") as pf:
                        for line in pf:
                            if len(line) <= 1:
                                break
                            fields = line.rstrip().split("\t")
                            accession_dict[fields[0]] = fields[2]
                            gzf.write(line)

        # generate taxid2 accession
        with shelve.Shelf(dbm.ndbm.open(taxid2wgs_accession_db.replace(".db", ""), 'c')) as taxid2accession_dict:
            with gzip.open(output_wgs_gz, "wt") as gzf:
                with open(wgs_accessions, 'r', encoding="utf-8") as wgsf:
                    for line in wgsf:
                        accession = line[1:].split(".")[0]
                        taxid = accession_dict.get(accession)
                        if taxid:
                            current_match = taxid2accession_dict.get(taxid, "")
                            taxid2accession_dict[taxid] = f"{current_match},{accession}"
                            gzf.write(line)

        accession_dict.close()

    def grab_accession_names(self, source_file, dest_file):
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''grep '^>' "${source_file}" | cut -f 1 -d' ' > "${dest_file}";''',
                named_args={
                    'source_file': source_file,
                    'dest_file': dest_file
                }
            )
        )

    def grab_wgs_accessions(self, source_file, dest_file):
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''grep '^>' "${source_file}" | grep 'complete genome' | cut -f 1 -d' ' > "${dest_file}";''',
                named_args={
                    'source_file': source_file,
                    'dest_file': dest_file
                }
            )
        )

    @run_in_subprocess
    def grab_accession_mapping_list(self, source_gz, num_partitions, partition_id,
                                    accessions, output_file):
        num_lines = 0
        with open(output_file, 'w') as out:
            with gzip.open(source_gz, 'r') as mapf:
                for line in mapf:
                    if num_lines % num_partitions == partition_id:
                        line = line.decode('utf-8')
                        accession = line.split("\t")[0]
                        if accession in accessions:
                            out.write(line)
                    num_lines += 1
                    if num_lines % 1000000 == 0:
                        log.write(f"{source_gz} partition {partition_id} line {num_lines/1000000}M")

    def count_reads(self):
        ''' Count reads '''
        pass
