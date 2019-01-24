''' Accession2Taxid'''
import gzip
import shelve
import dbm
import threading
from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.command import run_in_subprocess
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count

NUM_PARTITIONS = 8

class PipelineStepGenerateAccession2Taxid(PipelineStep):
    ''' Curate Accession2Taxid based on NT/NR Entries '''
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
        (output_gz, output_wgs_gz, accession2taxid_db, taxid2wgs_accession_db) = self.output_files_local()

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
        # generate the accession2taxid db and file
        accessions = [] # reset accessions to release memory
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
        wgs_thread.join()
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
        command.execute(f"grep '^>' {source_file} |cut -f 1 -d' ' > {dest_file}")

    def grab_wgs_accessions(self, source_file, dest_file):
        command.execute(f"grep '^>' {source_file} | grep 'complete genome' | cut -f 1 -d' ' > {dest_file}")

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

