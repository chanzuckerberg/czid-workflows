from multiprocessing import cpu_count
from typing import Iterator
import os
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
from idseq_dag.util.command import run_in_subprocess
import idseq_dag.util.log as log
import idseq_dag.util.count as count
import idseq_dag.util.fasta as fasta
from idseq_dag.util.thread_with_result import mt_map


class PipelineStepRunLZW(PipelineStep):

    MAX_SUBPROCS = 16

    # Core count caveats:
    #
    #   * Due to hyperthreading, the core count is exagerated 2x.
    #
    #   * When running inside a docker container, cpu_count reports the number of
    #     virtual CPU cores on the instance that is hosting the container.  There
    #     could be limits preventing the container from using all those cores.
    REAL_CORES = (cpu_count() + 1) // 2

    NUM_SLICES = min(MAX_SUBPROCS, REAL_CORES)

    def run(self):
        input_fas = self.input_files_local[0]
        output_fas = self.output_files_local()
        cutoff_scores = self.additional_attributes["thresholds"]
        threshold_readlength = self.additional_attributes.get("threshold_readlength", 150)
        PipelineStepRunLZW.generate_lzw_filtered(input_fas, output_fas, cutoff_scores, threshold_readlength)

    def count_reads(self):
        self.should_count_reads = True
        self.counts_dict[self.name] = count.reads_in_group(self.output_files_local()[0:2])

    @staticmethod
    def lzw_score(sequence, threshold_readlength):
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

        seq_length = len(sequence)
        lzw_fraction = float(len(results)) / seq_length

        if seq_length > threshold_readlength:
            # Make sure longer reads don't get excessively penalized
            adjustment_heuristic = (1 + (seq_length - threshold_readlength) / 1000) # TODO: revisit 
            score = lzw_fraction * adjustment_heuristic
        else:
            score = lzw_fraction
        return score

    @staticmethod
    def lzw_compute(input_files, threshold_readlength, slice_step=NUM_SLICES):
        """Spawn subprocesses on NUM_SLICES of the input files, then coalesce the
        scores into a temp file, and return that file's name."""

        temp_file_names = [f"lzwslice_{slice_step}_{slice_start}.txt" for slice_start in range(slice_step + 1)]
        for tfn in temp_file_names:
            assert not os.path.exists(tfn)

        @run_in_subprocess
        def lzw_compute_slice(slice_start):
            """For each read, or read pair, in input_files, such that read_index % slice_step == slice_start,
            output the lzw score for the read, or the min lzw score for the pair."""
            lzw_score = PipelineStepRunLZW.lzw_score
            with open(temp_file_names[slice_start], "a") as slice_output:
                for i, reads in enumerate(fasta.synchronized_iterator(input_files)):
                    if i % slice_step == slice_start:
                        lzw_min_score = min(lzw_score(r.sequence, threshold_readlength) for r in reads)
                        slice_output.write(str(lzw_min_score) + "\n")

        # slices run in parallel
        mt_map(lzw_compute_slice, range(slice_step))

        slice_outputs = temp_file_names[:-1]
        coalesced_score_file = temp_file_names[-1]
        # Paste can insert newlines at the end;  we grep those out.
        command.execute("paste -d '\n' " + " ".join(slice_outputs) + " | grep -v ^$ > " + coalesced_score_file)
        for tfn in slice_outputs:
            os.remove(tfn)
        return coalesced_score_file

    @staticmethod
    def generate_lzw_filtered(fasta_files, output_files, cutoff_scores, threshold_readlength):
        assert len(fasta_files) == len(output_files)

        # This is the bulk of the computation.  Everything else below is just binning by cutoff score.
        coalesced_score_file = PipelineStepRunLZW.lzw_compute(fasta_files, threshold_readlength)

        cutoff_scores.sort(reverse=True) # Make sure cutoff is from high to low

        readcount_list = [] # one item per cutoff
        outstreams_list = [] # one item per cutoff
        outfiles_list = [] # one item per cutoff

        for cutoff in cutoff_scores:
            readcount_list.append(0)
            outstreams = []
            outfiles = []
            for f in output_files:
                outfile_name = "%s-%f" % (f, cutoff)
                outfiles.append(outfile_name)
                outstreams.append(open(outfile_name, 'w'))

            outstreams_list.append(outstreams)
            outfiles_list.append(outfiles)

        outstreams_for_cutoff = list(zip(outstreams_list, cutoff_scores))

        def score_iterator(score_file: str) -> Iterator[float]:
            with open(score_file, "r") as sf:
                for line in sf:
                    yield float(line)

        total_reads = 0
        for reads, score in zip(fasta.synchronized_iterator(fasta_files), score_iterator(coalesced_score_file)):
            total_reads += 1
            for i, (outstreams, cutoff) in enumerate(outstreams_for_cutoff):
                if score > cutoff:
                    readcount_list[i] += 1
                    for ostr, r in zip(outstreams, reads):
                        ostr.write(r.header + "\n")
                        ostr.write(r.sequence + "\n")
                    break
        os.remove(coalesced_score_file)

        # closing all the streams
        for outstreams in outstreams_list:
            for ostr in outstreams:
                ostr.close()

        # get the right output file and metrics
        kept_count = 0
        filtered = total_reads
        cutoff_frac = None
        for cutoff_frac, readcount, outfiles in zip(cutoff_scores, readcount_list, outfiles_list):
            if readcount > 0:
                # found the right bin
                kept_count = readcount
                filtered = total_reads - kept_count
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
