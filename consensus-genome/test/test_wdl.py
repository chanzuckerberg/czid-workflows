import os
import json
import zipfile
import gzip
import tempfile
from subprocess import CalledProcessError
import hashlib
import yaml
from test_util import WDLTestCase

from Bio import SeqIO


nt_s3_path = "s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/nt"
nt_loc_db = "s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/nt_loc.marisa"
nr_s3_path = "s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/nr"
nr_loc_db = "s3://czid-public-references/ncbi-indexes-prod/2021-01-22/index-generation-2/nr_loc.marisa"


class TestConsensusGenomes(WDLTestCase):
    wdl = os.path.join(os.path.dirname(__file__), "..", "run.wdl")
    with open(os.path.join(os.path.dirname(__file__), "local_test.yml")) as fh:
        common_inputs = yaml.safe_load(fh)
    sc2_ref_fasta = "s3://czid-public-references/consensus-genome/MN908947.3.fa"

    def test_vadr_error_caught(self):
        # use the long filename error as a proxy for testing VADR error handling
        consensus = os.path.join(os.path.dirname(__file__), "vadr_input", "really-long-name-consensus.fa")
        vadr_opts_string = ("-s -r --nomisc --mkey NC_045512 --lowsim5term 2 --lowsim3term 2 --fstlowthr 0.0 "
                            "--alt_fail lowscore,fsthicnf,fstlocnf")
        args = ["vadr_model=s3://czid-public-references/consensus-genome/vadr-models-sarscov2-1.2-2.tar.gz",
                f"vadr_options={vadr_opts_string}",
                f"assembly={consensus}"]
        res = self.run_miniwdl(args=args, task="Vadr", task_input={"prefix": ""})
        self.assertIsNotNone(res["outputs"]["Vadr.vadr_errors"])
        self.assertEqual(res["outputs"]["Vadr.vadr_alerts"], None)
        self.assertEqual(res["outputs"]["Vadr.vadr_quality"], None)

    def test_vadr_flag_works(self):
        consensus = os.path.join(os.path.dirname(__file__), "vadr_input", "really-long-name-consensus.fa")
        vadr_opts_string = ("-s -r --nomisc --mkey sarscov2 --lowsim5term 2 --lowsim3term 2 "
                            "--fstlowthr 0.0 --alt_fail lowscore,fsthicnf,fstlocnf --noseqnamemax")
        args = ["vadr_model=s3://czid-public-references/consensus-genome/vadr-models-sarscov2-1.2-2.tar.gz",
                f"vadr_options={vadr_opts_string}",
                f"assembly={consensus}"]
        res = self.run_miniwdl(args=args, task="Vadr", task_input={"prefix": ""})
        self.assertIsNotNone(res["outputs"]["Vadr.vadr_alerts"])
        self.assertIsNone(res["outputs"]["Vadr.vadr_errors"])

    # test the depths associated with SNAP ivar trim -x 5
    def test_sars_cov2_illumina_cg_snap(self):
        aligned_reads = os.path.join(os.path.dirname(__file__), "trim_primers_input", "snap_aligned_reads.bam")
        args = [f"alignments={aligned_reads}",
                "primer_bed=s3://czid-public-references/consensus-genome/snap_primers.bed"]
        res = self.run_miniwdl(args, task="TrimPrimers", task_input={"prefix": ""})
        with open(res["outputs"]["TrimPrimers.trimmed_bam_bai"], 'rb') as f:
            hash = hashlib.md5(f.read()).hexdigest()
        self.assertEqual(hash, "3081c5e09cc31194821a84fc24b5685f")
        with open(res["outputs"]["TrimPrimers.trimmed_bam_ch"], 'rb') as f:
            hash = hashlib.md5(f.read()).hexdigest()
        self.assertEqual(hash, "a4d8b1d6e5d6a0bbc4b336919baddffd")

    # test the depths associated with tailedseq protocol, ivar trim -x 2
    def test_sars_cov2_illumina_cg_tailedseq(self):
        aligned_reads = os.path.join(os.path.dirname(__file__), "trim_primers_input", "tailedseq_aligned_reads.bam")
        args = [f"alignments={aligned_reads}",
                "primer_bed=s3://czid-public-references/consensus-genome/artic_v3_short_275_primers.bed"]
        res = self.run_miniwdl(args, task="TrimPrimers", task_input={"prefix": ""})
        with open(res["outputs"]["TrimPrimers.trimmed_bam_bai"], 'rb') as f:
            hash = hashlib.md5(f.read()).hexdigest()
        self.assertEqual(hash, "010dfb9df3c7a45c9bf6556925859ce7")
        with open(res["outputs"]["TrimPrimers.trimmed_bam_ch"], 'rb') as f:
            hash = hashlib.md5(f.read()).hexdigest()
        self.assertEqual(hash, "7028bf8450548391f264f2948b9e19f0")

    def test_sars_cov2_medaka_model(self):
        """
        Test that the pipeline will run a variety of different medaka models
        """
        models = ["r10_min_high_g340", "r941_prom_fast_g303"]
        fastq = os.path.join(os.path.dirname(__file__), "no_host_1.fq.gz")
        for model in models:
            args = ["prefix=''", "sample=test_sample", f"fastqs={fastq}",
                    "normalise=1000", f"medaka_model={model}",
                    "primer_schemes=s3://czid-public-references/consensus-genome/artic-primer-schemes_v2.tar.gz",
                    "primer_set=nCoV-2019/V1200"]
            res = self.run_miniwdl(args, task="RunMinion")
            for filename in res["outputs"].values():
                self.assertGreater(os.path.getsize(filename), 0)
            with open(res["outputs"]["RunMinion.log"]) as f:
                log_output = f.read()
            self.assertIn("primer_schemes/nCoV-2019/V1200", log_output)
            for filename in res["outputs"].values():
                self.assertGreater(os.path.getsize(filename), 0)

    def test_sars_cov2_medaka_fail(self):
        """
        Test that the pipeline will fail if the medaka model is incompatible
        """
        model = "r941_prom_snp_g360"
        fastq = os.path.join(os.path.dirname(__file__), "no_host_1.fq.gz")
        args = ["prefix=''", "sample=test_sample", f"fastqs={fastq}",
                "normalise=1000", f"medaka_model={model}",
                "primer_schemes=s3://czid-public-references/consensus-genome/artic-primer-schemes.tar.gz",
                "primer_set=nCoV-2019/V3"]
        with self.assertRaises(CalledProcessError) as ecm:
            self.run_miniwdl(args, task="RunMinion")
        miniwdl_error = json.loads(ecm.exception.output)
        self.assertEqual(miniwdl_error["error"], "RunFailed")
        self.assertEqual(miniwdl_error["cause"]["error"], "CommandFailed")

    def test_sars_cov2_midnight_primers_quast(self):
        """
        Test that Quast will run with midnight primers
        """
        assembly = os.path.join(os.path.dirname(__file__), "quast_input", "test_sample.consensus.fasta")
        bam = os.path.join(os.path.dirname(__file__), "quast_input", "test_sample.primertrimmed.rg.sorted.bam")
        fastq = os.path.join(os.path.dirname(__file__), "no_host_1.fq.gz")
        args = [
            "prefix=''",
            f"assembly={assembly}",
            f"bam={bam}",
            f"fastqs={fastq}",
            "no_reads_quast=false",
            "technology=ONT",
            "primer_schemes=s3://czid-public-references/consensus-genome/artic-primer-schemes_v2.tar.gz",
            "primer_set=nCoV-2019/V1200"
        ]
        res = self.run_miniwdl(args, task="Quast")
        for filename in res["outputs"]["Quast.quast_dir"]:
            self.assertGreater(os.path.getsize(filename), 0)

    def test_sars_cov2_ont_cg_input_file_format(self):
        """
        Tests that .fq.gz inputs are correctly converted to .fastq.gz formats to ensure they are read by guppyplex in
        ApplyLengthFilter
        """
        fastqs = os.path.join(os.path.dirname(__file__), "Ct20K.fq.gz")
        res = self.run_miniwdl(
            task="ValidateInput",
            args=[
                "prefix=test",
                f"fastqs={fastqs}",
                "technology=ONT",
                "max_reads=75000000",
            ],
        )
        for output_name, output in res["outputs"].items():
            for filename in output:
                if output_name == 'ValidateInput.input_stats' or output_name == 'ValidateInput.seqkit_version':
                    continue
                self.assertTrue(filename.endswith(".fastq.gz"))
                self.assertGreater(os.path.getsize(filename), 0)

    def test_zip_outputs(self):
        res = self.run_miniwdl(task="ZipOutputs", args=[
            "prefix=test",
            f"outputFiles={self.wdl}",
        ])
        with zipfile.ZipFile(res["outputs"]["ZipOutputs.output_zip"]) as fh:
            self.assertEqual(fh.namelist(), ["run.wdl"])

    def test_fetch_sequence_by_accession_id(self):
        res = self.run_miniwdl(task="FetchSequenceByAccessionId", args=[
            "accession_id=HIN71308.1",
            f"nt_s3_path={nt_s3_path}",
            f"nt_loc_db={nt_loc_db}",
            f"nr_s3_path={nr_s3_path}",
            f"nr_loc_db={nr_loc_db}",
        ])
        with open(res["outputs"]["FetchSequenceByAccessionId.sequence_fa"]) as fh:
            self.assertEqual(fh.readline().strip(), ">HIN71308.1 hypothetical protein [Dehalococcoidia bacterium]")
            self.assertTrue(fh.readline().startswith("MSMIIQSGKLTQTRVEAQHMDFEPRYTEEQQTFRTEVRDWLKDNVPDDLANPADSAD"))

        with self.assertRaises(CalledProcessError) as ecm:
            self.run_miniwdl(task="FetchSequenceByAccessionId", args=[
                "accession_id=NO_ACCESSION_ID",
                f"nt_s3_path={nt_s3_path}",
                f"nt_loc_db={nt_loc_db}",
                f"nr_s3_path={nr_s3_path}",
                f"nr_loc_db={nr_loc_db}",
            ])
        self.assertRunFailed(ecm, task="FetchSequenceByAccessionId",
                             error="AccessionIdNotFound",
                             cause=("The Accession ID was not found in the CZ ID database, "
                                    "so a generalized consensus genome could not be run"))

    def test_no_coverage_error(self):
        ref_host_blank = os.path.join(os.path.dirname(__file__), "blank.fastq.gz")
        assembly_blank = os.path.join(os.path.dirname(__file__), "blank.fastq.gz")
        vcf_blank = os.path.join(os.path.dirname(__file__), "blank.fastq.gz")
        fastqs_blank = os.path.join(os.path.dirname(__file__), "blank.fastq.gz")
        cleaned_bam = os.path.join(os.path.dirname(__file__), "empty.bam")
        with self.assertRaises(CalledProcessError) as ecm:
            self.run_miniwdl(task="ComputeStats", args=[
                f"ref_host={ref_host_blank}",
                f"assembly={assembly_blank}",
                f"vcf={vcf_blank}",
                f"fastqs={fastqs_blank}",
                "technology=Illumina",
                f"cleaned_bam={cleaned_bam}"],
                task_input={
                    "sample": "sample",
                    "prefix": "",
                })
        self.assertRunFailed(ecm, task="ComputeStats",
                             error="InsufficientReadsError",
                             cause="There was insufficient coverage so a consensus genome could not be created.")

    def test_max_reads_illumina(self):
        fastq_0 = os.path.join(os.path.dirname(__file__), "SRR11741455_65054_nh_R1.fastq.gz")
        fastq_1 = os.path.join(os.path.dirname(__file__), "SRR11741455_65054_nh_R2.fastq.gz")

        res = self.run_miniwdl(
            task="ValidateInput",
            task_input={
                "max_reads": 100,
                "technology": "Illumina",
                "fastqs": [fastq_0, fastq_1],
                "prefix": "",
            },
        )
        for output in res["outputs"]["ValidateInput.validated_fastqs"]:
            with gzip.open(output, 'rt') as f:
                self.assertEqual(sum(1 for _ in SeqIO.parse(f, "fastq")), 100)

    def test_max_reads_ont(self):
        fastq = os.path.join(os.path.dirname(__file__), "SRR11741455_65054_nh_R1.fastq.gz")

        res = self.run_miniwdl(
            task="ValidateInput",
            task_input={
                "max_reads": 100,
                "technology": "ONT",
                "fastqs": [fastq],
                "prefix": "",
            },
        )
        for output in res["outputs"]["ValidateInput.validated_fastqs"]:
            with gzip.open(output, 'rt') as f:
                self.assertEqual(sum(1 for _ in SeqIO.parse(f, "fastq")), 100)

    def test_max_reads_uncompressed_input(self):
        fastq = os.path.join(os.path.dirname(__file__), "SRR11741455_65054_nh_R1.fastq.gz")
        with tempfile.NamedTemporaryFile('wb') as f, gzip.open(fastq) as gzipped_f:
            f.write(gzipped_f.read())

            res = self.run_miniwdl(
                task="ValidateInput",
                task_input={
                    "max_reads": 100,
                    "technology": "ONT",
                    "fastqs": [f.name],
                    "prefix": "",
                }
            )
        for output in res["outputs"]["ValidateInput.validated_fastqs"]:
            with gzip.open(output, 'rt') as f:
                self.assertEqual(sum(1 for _ in SeqIO.parse(f, "fastq")), 100)
