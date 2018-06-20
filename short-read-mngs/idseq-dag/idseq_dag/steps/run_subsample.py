import random

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count

class PipelineStepRunSubsample(PipelineStep):
    '''
    Subsample the input data to max_fragments
    '''
    def run(self):
        ''' Invoking subsampling '''
        input_fas = self.input_files_local[0]
        output_fas = self.output_files_local()
        max_fragments = self.additional_attributes["max_fragments"]
        PipelineStepRunSubsample.subsample_fastas(input_fas, output_fas, max_fragments)

    def count_reads(self):
        files_to_count = self.output_files_local()[0:2]
        read_count = count.reads_in_group(files_to_count)
        self.counts_dict[self.name] = read_count
        # If the read count is exactly equal to the maximum allowed number,
        # infer that subsampling occurred:
        max_read_count = len(files_to_count) * self.additional_attributes["max_fragments"]
        if read_count == max_read_count:
            self.counts_dict["subsampled"] = 1

    @staticmethod
    def subsample_fastas(input_fas, output_fas, max_fragments):
        ''' In memory subsampling '''
        paired = len(input_fas) >= 2
        # count lines
        cmd = "wc -l %s | cut -f1 -d ' '" % input_fas[0]
        total_records = int(command.execute_with_output(cmd).decode('utf-8')) // 2
        log.write("total reads: %d" % total_records)
        log.write("target reads: %d" % max_fragments)
        if total_records <= max_fragments:
            for infile, outfile in zip(input_fas, output_fas):
                command.execute("cp %s %s" % (infile, outfile))
            return

        # total_records > max_fragments, sample
        randgen = random.Random(x=hash(input_fas[0]))
        records_to_keep = randgen.sample(range(total_records), max_fragments)
        PipelineStepRunSubsample.subset(input_fas[0], output_fas[0], records_to_keep)
        if paired:
            PipelineStepRunSubsample.subset(input_fas[1], output_fas[1], records_to_keep)
            if len(input_fas) == 3 and len(output_fas) == 3:
                # subset the merged fasta
                records_to_keep_merged = []
                for r in records_to_keep:
                    records_to_keep_merged += [2*r, 2*r+1]
                PipelineStepRunSubsample.subset(input_fas[2], output_fas[2],
                                                records_to_keep_merged)



    @staticmethod
    def subset(input_fa, output_fa, records_to_keep):
        ''' Subset input_fa based on record indice specified in records_to_keep '''
        records_to_keep = set(records_to_keep)
        record_number = 0
        kept_records = 0
        with open(input_fa, 'r', encoding='utf-8') as input_f, open(output_fa, 'w') as output_f:
            # Iterate through the FASTA file records
            sequence_name = input_f.readline()
            sequence_data = input_f.readline()
            while len(sequence_name) > 0 and len(sequence_data) > 0:
                if record_number in records_to_keep:
                    output_f.write(sequence_name)
                    output_f.write(sequence_data)
                    kept_records += 1
                sequence_name = input_f.readline()
                sequence_data = input_f.readline()
                record_number += 1
        if len(records_to_keep) != kept_records:
            raise RuntimeError("subsample subset error: records to keep: %d, records kept: %d" % (
                len(records_to_keep), kept_records))
