''' Generate Coverage Statistics '''
import json
import os
import pysam
import functools
from multiprocessing import Pool

from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.command_patterns as command_patterns
import idseq_dag.util.log as log
from idseq_dag.util.sequence import chunks
MIN_CONTIG_FILE_SIZE = 50

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
        contig2coverage = self.calc_contig2coverage(bam_file)

        with open(coverage_json, 'w') as csf:
            json.dump(contig2coverage, csf)

        with open(coverage_summary_csv, 'w') as csc:
            csc.write("contig_name,avg,min,max,p25,p50,p75,avg2xcnt,cnt0,cnt1,cnt2\n")
            for contig, stats in contig2coverage.items():
                output_row = [
                    contig, stats['avg'], stats['p0'], stats['p100'],
                    stats['p25'], stats['p50'], stats['p75'],
                    stats['avg2xcnt'], stats['cnt0'], stats['cnt1'], stats['cnt2']
                ]
                output_str = ','.join(str(e) for e in output_row)
                csc.write(output_str + "\n")

    @staticmethod
    def _process_contig_names(bam_file, f, contig_names):
        contig2coverage = {}
        if len(contig_names) > 0:
            with log.log_context(f"calc_contig2coverage_chunk", {"bam_file": bam_file, "from": contig_names[0], "to": contig_names[-1]}):
                for c in contig_names:
                    coverage = []
                    for pileup_column in f.pileup(contig=c):
                        coverage.append(pileup_column.get_num_aligned())
                    sorted_coverage = sorted(coverage)
                    contig_len = len(coverage)
                    if contig_len <= 0:
                        continue

                    avg = sum(coverage)/contig_len
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

                    contig2coverage[c] = {
                        "coverage": coverage,
                        "avg": avg,
                        "p0": sorted_coverage[0],
                        "p100": sorted_coverage[-1],
                        "p25": sorted_coverage[int(0.25*contig_len)],
                        "p50": sorted_coverage[int(0.5*contig_len)],
                        "p75": sorted_coverage[int(0.75*contig_len)],
                        "avg2xcnt": avg2xcnt / contig_len,
                        "cnt0": cnt0 / contig_len,
                        "cnt1": cnt1 / contig_len,
                        "cnt2": cnt2 / contig_len
                    }
        return contig2coverage

    @staticmethod
    def calc_contig2coverage(bam_file):
        CHUNK_SIZE = 100
        with pysam.AlignmentFile(bam_file, "rb") as f:
            all_references = f.references
            with log.log_context("calc_contig2coverage", {"bam_file": bam_file, "all_references_count": len(all_references)}):
                chunked_references = chunks(f.references, CHUNK_SIZE)
                fn = functools.partial(PipelineStepGenerateCoverageStats._process_contig_names, bam_file, f)
                with Pool() as pool:
                    results = pool.map(fn, chunked_references)
        contig2coverage = {}
        for result in results:
            contig2coverage.update(result)
        return contig2coverage


    def count_reads(self):
        ''' Count reads '''
        pass
