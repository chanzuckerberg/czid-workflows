''' Generate Taxon Lineage DB in sqlite3 '''
import gzip
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.log as log
import idseq_dag.util.count as count
from idseq_dag.util.dict import IdSeqDict, IdSeqDictValue

BATCH_INSERT_SIZE = 300

class PipelineStepGenerateLineageDB(PipelineStep):
    ''' Generate TaxonLineage Table in sqlite3 format. '''
    def run(self):
        """
        Download taxon-lineage.csv.gz and convert it to sqlite3 dict
        Example command:
        aegea batch submit --command="pip install --upgrade git+git://github.com/chanzuckerberg/s3mi.git; cd /mnt; git clone https://github.com/chanzuckerberg/idseq-dag.git; cd idseq-dag;  pip3 install -e .; sed 's/<NCBI_DATE>/2018-12-01/g' templates/lineage_db.json > templates/lineage_db.json.now  idseq_dag --no-versioned-output templates/lineage_db.json.now"  --storage /mnt=1000 --volume-type gp2 --ecr-image idseq_dag --memory 240000 --queue idseq-prod-himem --vcpus 32 --job-role idseq-pipeline

        """
        lineage_csv_gz = self.input_files_local[0][0]
        output_db = self.output_files_local()[0]
        log.write(f"input: {lineage_csv_gz} output: {output_db}")

        lineage_dict = IdSeqDict(output_db, IdSeqDictValue.VALUE_TYPE_ARRAY)
        batch_list = {}
        with gzip.open(lineage_csv_gz, "rt") as gzf:
            for line in gzf:
                fields = line.rstrip().split(",")
                taxid = fields[0]
                species, genus, family = fields[-1:-4:-1]
                batch_list[taxid] = [species, genus, family]
                if len(batch_list) >= BATCH_INSERT_SIZE:
                    lineage_dict.batch_inserts(batch_list.items())
                    batch_list = {}
            lineage_dict.batch_inserts(batch_list.items())

    def count_reads(self):
        ''' Count reads '''
        pass
