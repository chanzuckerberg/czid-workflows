from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count

class PipelineStepRunLZW(PipelineStep):

    def run(self):
        input_fas = self.input_files_local[0]
        output_fas = self.output_files_local()
        cutoff_fractions = self.additional_attributes["thresholds"]
        PipelineStepRunLZW.generate_lzw_filtered(input_fas, output_fas, cutoff_fractions)

    def count_reads(self):
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])

    @staticmethod
    def lzw_fraction(sequence):
        sequence = str(sequence)
        if sequence == "":
            return 0.0
        sequence = sequence.upper()

        dictionary = {}
        dict_size = 0
        for c in sequence:
            if c not in dictionary:
                dict_size += 1
                dictionary[c] = dict_size

        word = ""
        results = []
        for c in sequence:
            wc = word + c
            if dictionary.get(wc):
                word = wc
            else:
                results.append(dictionary[word])
                dict_size += 1
                dictionary[wc] = dict_size
                word = c
        if word != "":
            results.append(dictionary[word])
        return float(len(results)) / len(sequence)

    @staticmethod
    def generate_lzw_filtered(fasta_files, output_files, cutoff_fractions):
        assert(len(fasta_files) == len(output_files))
        cutoff_fractions.sort(reverse=True) # Make sure cutoff is from high to low

        read_streams = []
        readcount_list = [] # one item per cutoff
        outstream_list = [] # one item per cutoff
        outfiles_list = [] # one item per cutoff

        for fasta_file in fasta_files:
            read_streams.append(open(fasta_file, 'rb'))
        paired = len(fasta_files) == 2

        for cutoff in cutoff_fractions:
            readcount_list.append(0)
            outstream = []
            outfiles = []
            for f in output_files:
                outfile_name = "%s-%f" % (f, cutoff)
                outfiles.append(outfile_name)
                outstream.append(open(outfile_name, 'wb'))

            outstream_list.append(outstream)
            outfiles_list.append(outfiles)

        total_reads = 0
        while True:
            line_r1_header = read_streams[0].readline()
            line_r1_sequence = read_streams[0].readline()
            line_r2_header = None
            line_r2_sequence = None

            if line_r1_header and line_r1_sequence:
                total_reads += 1
                if paired:
                    line_r2_header = read_streams[1].readline()
                    line_r2_sequence = read_streams[1].readline()
                fraction_1 = PipelineStepRunLZW.lzw_fraction(line_r1_sequence.rstrip())
                fraction_2 = PipelineStepRunLZW.lzw_fraction(line_r2_sequence.rstrip()) if paired else 1.0
                fraction = min(fraction_1, fraction_2)
                for i in range(len(cutoff_fractions)):
                    cutoff = cutoff_fractions[i]
                    if fraction >= cutoff:
                        readcount_list[i] += 1
                        outstream = outstream_list[i]
                        outstream[0].write(line_r1_header)
                        outstream[0].write(line_r1_sequence)
                        if paired:
                            outstream[1].write(line_r2_header)
                            outstream[1].write(line_r2_sequence)
                        break
            else:
                break

        # closing all the streams
        for os in outstream_list:
            for o in os:
                o.close()
        for r in read_streams:
            r.close()
        # get the right output file and metrics
        cutoff_frac = cutoff_fractions[-1]
        kept_count = 0
        filtered = total_reads
        for i, readcount in enumerate(readcount_list):
            if readcount > 0:
                # found the right bin
                cutoff_frac = cutoff_fractions[i]
                kept_count = readcount
                filtered = total_reads - kept_count
                outfiles = outfiles_list[i]
                # move the output files over
                for outfile, output_file in zip(outfiles, output_files):
                    command.execute("mv %s %s" % (outfile, output_file))
                break

        if kept_count == 0:
            raise RuntimeError("All the reads are filtered by LZW with lowest cutoff: %f" % cutoff_frac)

        kept_ratio = float(kept_count)/float(total_reads)
        msg = "LZW filter: cutoff_frac: %f, total reads: %d, filtered reads: %d, " \
              "kept ratio: %f" % (cutoff_frac, total_reads, filtered, kept_ratio)
        log.write(msg)

