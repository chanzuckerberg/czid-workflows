''' Generate Coverage Statistics '''
import json
import os
import sys
import traceback
from multiprocessing import cpu_count

import pysam

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
from idseq_dag.util.command import run_in_subprocess, LongRunningCodeSection
import idseq_dag.util.command_patterns as command_patterns
from idseq_dag.util.thread_with_result import mt_map


MIN_CONTIG_FILE_SIZE = 50


COVERAGE_STATS_SCHEMA = {
    "contig_name": str,
    "avg": float,
    "min": int,
    "max": int,
    "p25": int,
    "p50": int,
    "p75": int,
    "avg2xcnt": int,
    "cnt0": int,
    "cnt1": int,
    "cnt2": int
}


class PipelineStepGenerateCoverageStats(PipelineStep):
    ''' Generate Coverage Statistics from Assembly Output '''
    def run(self):
        """
          1. extract contigs.fasta and read-contig.sam
          2. run pile up
        """
        contigs, _scaffolds, read_contig_sam, _stats = self.input_files_local[0]
        coverage_json, coverage_summary_csv = self.output_files_local()

        if os.path.getsize(contigs) < MIN_CONTIG_FILE_SIZE:
            command.write_text_to_file('{}', coverage_json)
            command.write_text_to_file('No Contigs', coverage_summary_csv)
            return

        # generate bam files
        bam_file = read_contig_sam.replace(".sam", ".bam")
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''samtools view -S -b "${read_contig_sam}" | samtools sort - -o "${bam_file}";''',
                named_args={
                    'read_contig_sam': read_contig_sam,
                    'bam_file': bam_file
                }
            )
        )
        command.execute(
            command_patterns.SingleCommand(
                cmd="samtools",
                args=[
                    "index",
                    bam_file
                ]
            )
        )
        # run coverage info
        output_csv, output_json = self.calc_contig2coverage(bam_file)
        os.rename(output_csv, coverage_summary_csv)
        os.rename(output_json, coverage_json)

    @staticmethod
    def _process_contig(input_bam, output_csv, output_json, contig_name):
        pileup = input_bam.pileup(contig=contig_name)
        coverage = [pileup_column.get_num_aligned() for pileup_column in pileup]
        if not coverage:
            return False
        contig_len = len(coverage)
        avg = sum(coverage) / contig_len
        avg2xcnt = 0
        cnt0 = 0
        cnt1 = 0
        cnt2 = 0
        for t in coverage:
            if t > 2 * avg:
                avg2xcnt += 1
            if t == 0:
                cnt0 += 1
            elif t == 1:
                cnt1 += 1
            elif t == 2:
                cnt2 += 1
        sorted_coverage = sorted(coverage)
        stats = {
            "coverage": coverage,
            "avg": avg,
            "p0": sorted_coverage[0],
            "p100": sorted_coverage[-1],
            "p25": sorted_coverage[int(0.25 * contig_len)],
            "p50": sorted_coverage[int(0.5 * contig_len)],
            "p75": sorted_coverage[int(0.75 * contig_len)],
            "avg2xcnt": avg2xcnt / contig_len,
            "cnt0": cnt0 / contig_len,
            "cnt1": cnt1 / contig_len,
            "cnt2": cnt2 / contig_len
        }
        output_json.write(json.dumps(contig_name) + ": ")
        output_json.write(json.dumps(stats))
        output_json.write(", ")
        stats["contig_name"] = contig_name
        stats["min"] = stats.pop("p0")
        stats["max"] = stats.pop("p100")
        output_csv.write(",".join(str(stats[col]) for col in COVERAGE_STATS_SCHEMA))
        output_csv.write("\n")
        return True

    @staticmethod
    def calc_contig2coverage(bam_filename):
        # PySAM pileup is CPU-intenstive.  Each CPU core is assigned a slice of the input bam file on which to perform pileup.  The slice contigs are selected by slice_idx modulo num_slices.  Each slice gets its own pair of temporary output files, one in CSV format and one in JSON.  In the end, these slice outputs are concatenated.  This is a similar pattern to run_lzw.
        num_physical_cpu = (cpu_count() + 1) // 2
        num_slices = num_physical_cpu
        output_csv_filenames = [f"tmp_slice_{num_slices}_{slice}.csv" for slice in range(num_slices + 1)]
        output_json_filenames = [f"tmp_slice_{num_slices}_{slice}.json" for slice in range(num_slices + 1)]
        for fn in output_csv_filenames + output_json_filenames:
            if os.path.exists(fn):
                os.remove(fn)

        @run_in_subprocess
        def compute_slice(slice_idx):
            try:
                with open(output_csv_filenames[slice_idx], "w") as output_csv, \
                     open(output_json_filenames[slice_idx], "w") as output_json, \
                     pysam.AlignmentFile(bam_filename, "rb") as input_bam:  # noqa: E126
                    for contig_idx, contig_name in enumerate(input_bam.references):
                        if contig_idx % num_slices == slice_idx:
                            PipelineStepGenerateCoverageStats._process_contig(input_bam, output_csv, output_json, contig_name)
            except Exception as e:
                print(f"Exception in subprocess: {e.__class__} '{e}'")
                print("-" * 60)
                traceback.print_exc(file=sys.stdout)
                print("-" * 60)
                raise e

        # Compute pileup for each slice
        with LongRunningCodeSection("PipelineStepGenerateCoverageStats.calc_contig2coverage.mt_map"):
            mt_map(compute_slice, range(num_slices))
        # Output CSV headers
        with open(output_csv_filenames[-1], "w") as ocsv:
            ocsv.write(",".join(COVERAGE_STATS_SCHEMA))
            ocsv.write("\n")
        # Output JSON dict open paren
        with open(output_json_filenames[-1], "w") as ojson:
            ojson.write("{")
        # Collate CSV slices
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''cat "${individual_slice_outputs[@]}" >> "${collated_csv}";''',  # note >> for appending
                named_args={
                    'collated_csv': output_csv_filenames[-1],
                    'individual_slice_outputs': output_csv_filenames[:-1]
                }
            )
        )
        for tfn in output_csv_filenames[:-1]:
            os.remove(tfn)
        # Collate JSON slices, replacing final ", " with "}"
        command.execute(
            command_patterns.ShellScriptCommand(
                script=r'''cat "${individual_slice_outputs[@]}" | sed 's=, $=}=' >> "${collated_json}";''',  # note >> for appending
                named_args={
                    'collated_json': output_json_filenames[-1],
                    'individual_slice_outputs': output_json_filenames[:-1]
                }
            )
        )
        for tfn in output_json_filenames[:-1]:
            os.remove(tfn)
        return (output_csv_filenames[-1], output_json_filenames[-1])

    def count_reads(self):
        ''' Count reads '''
        pass
