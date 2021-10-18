# IDSeq Workflows Changelog

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
and is based on recommendations from [Keep a Changelog](https://keepachangelog.com/en/1.0.0/). 

To add clarity for users it is [recommended](https://keepachangelog.com/en/1.0.0/) to differentiate between the types of changes that have been made: 

`Added` for new features.
`Changed` for changes in existing functionality.
`Deprecated` for soon-to-be removed features.
`Removed` for now removed features.
`Fixed` for any bug fixes.
`Security` in case of vulnerabilities.

In addition, if a modification is made that may affect the results of a pipeline run, we suggest using noting the `[Pipeline Change]` as well as including the specific change that was made, the predicted result to the output. 

### Unreleased

### consensus-genome-v3.4.3
- fix a bug to print out errors properly 

### consensus-genome-v3.4.2
- reduce size of tests
- retry accession id with updated version
- refactor validation input to use seqkit

### short-read-mngs-v6.8.7
- don't fail if picard fails, emit warning instead

### short-read-mngs-v6.8.6
- add step-name in STAR comment


### consensus-genome-v3.4.1
- [Pipeline Change] Modify parameters to reduce stringency in consensus-genome pipeline
	- Change: reduce `ivarFreqThreshold` parameter from 0.9 to 0.75 and reduce `minDepth` parameter from 10 to 5 when using generalized consensus genome
	- Predicted result: fewer ambiguous bases
- Fixed broken consensus-genome tests
- Fixed bitrot error in benchmarking notebook

### short-read-mngs-v6.8.6
- add a comment into STAR filter step to allow miniwdl to parse step name. Fixes a bug where the "Reads Lost" chart included STAR under Trimmomatic counts. 

### phylotree-ng-v1.2.3
- changed the phylo-heatmap scale bar from '0.15' to '>0.15' to remove ambiguity

### short-read-mngs-v6.8.5
- mitigation for intermittent GenerateCoverageStats error (#161)

### consensus-genome-v3.4.0
- truncate consensus genome inputs (#134)
- [Pipeline Change] Fix input min/max length when using midnight primers (#145)
	- Change: Relax length filter from 350-700 to 250-1500 when  using midnight primers
	- Predicted result: outputs for consensus genomes run with midnight primers will probably make more sense 

### phylotree-ng-v1.2.2
- ignore divergent iqtree error (#162)
- output ska.variants.aln if no tree is produced (#160)

### short-read-mngs-v6.8.4
- fix output gene file name (#159)

### short-read-mngs-v6.8.3
- fix bug with reading back idseq clusters (#158)

### short-read-mngs-v6.8.2 
- force path to genome to standard name (#157)

### phylotree-ng-v1.2.1
- fix issue with clusters in phylotree (#156)
- remove another instance of TooDivergentError from phylo tree (#153)

### short-read-mngs-v6.8.1
- enhance csv security (#154)
- Add coverage breadth to coverage_viz_summary output (#155)
- add error logging to generate coverage stats (#151)
- add short-read-mngs/auto_benchmark/harvest2.py
- Replace RunStar in host filter (#139)
- fix sample swap in short-read-mngs benchmark library

### phylotree-ng-v1.2.0
- fix select_first for newick (#150)
- Tweak phylo heatmap style (#149)
- Handle windows line endings (#144)
- [CH-143966] Adjust phylo heatmap styling (#147)
- Add names to newick (#148)
- updates to README prior to v0 launch (#146)
- Update README.md

### phylotree-ng-v1.1.0
- succeed without phylo tree (#143)
(refs/heads/handle-no-lines
- CI: push docker images to ghcr.io (#141)

### consensus-genome-v3.3.0
- add new RealignConsensus step (#138)
- Add option for primer set  (#133)

### phylotree-ng-v1.0.1
- add s3_wd_uri input to phylotree-ng (#136)

### phylotree-ng-v1.0.0
- fix bug with NA values in ska.distances.tsv and add known user error for divergent samples (#131)
- Phylo Tree (#123)
- fix reading ska tsv (#107)
- Expand phylotree-ng description
- PhyloTree-NG workflow skeleton (#104)

### short-read-mngs-v6.8.0
- re-enable threading in bowtie2 (#132)

### consensus-genome-v3.2.1
- add ability to select medaka model(#130)

### consensus-genome-v3.2.0
- Upgrade vadr model to vadr-models-sarscov2-1.2-2.tar.gz (#129)

### short-read-mngs-v6.7.0
- Update samtools parameters to improve indel calling (#122)
- short-read-mngs RunAssembly: contig length filter (#127)
- don't emit valid input files (#128)

### short-read-mngs-v6.6.0
- remove nr-refseq output (#124)

### short-read-mngs-v6.5.1
- handle null output_log_file (#121)
- switch to initial space safety (#120)
- quote clusters csv (#119)
- optimize github actions jobs to avoid redundant runs (#118)
- remove nt.refseq output (#117)

### short-read-mngs-v6.5.0
- emit star logfile (#115)

### short-read-mngs-v6.4.4
- fix batch waiting bug (#114)

### short-read-mngs-v6.4.3
- remove redundant manual timeout (#113)

### short-read-mngs-v6.4.2
- raise proper insufficient reads error (#112)
- Upgrade VADR to v1.2 (#108)
- Add 2bp wiggle room to ivar primer trimming for tailedseq protocol (#105)

### consensus-genome-v3.1.1
- Bump taxoniq to v0.6.0 (#106)
- Update README.md (#100)
- Call VADR only on SC2 samples (#103)
- Avoid querying GitHub API in Dockerfile
- Emit legible error on missing accession id (#101)
- Use flake8 in lint (#102)
- Ensure recent miniwdl in dev requirements

### consensus-genome-v3.1.0
- CG: Make length filter optional (#99)
- Remove references to idseq-database (#98)
- adding provenance documentation for cg files (#92)

### consensus-genome-v3.0.0
- Migrate readme comments from #74 (#95)
- Add missing refactored default input
- Add support for consensus genome generation on any reference sequence (#87)

### consensus-genome-v2.2.0
- update to catch vadr errors and associated unit tests (#94)

### short-read-mngs-v6.4.1
- fix fasta writing bug (#93)
- fix short-read-mngs viral benchmarks
- relocate tests/short-read-mngs/auto_benchmark to short-read-mngs/auto_benchmark
- automating short-read-mngs full-size benchmarks (#91)
- retire idseq-dag step status JSON (#85)
- simplify pytest fixture scheme (#88)
- Update ivar trim to address issues with trimming in SWIFT/SNAP protocol (#90)

### consensus-genome-v2.1.0
- Updates to address gaps identified in initial ont validation (#89)

### consensus-genome-v2.0.0
- Update Consensus Genomes workflow to handle Oxford Nanopore data (#86)
- Use automake and libhts-dev
- Install ivar 1.3.1
- Add tabix
- Add libhts3
- Add trailing newline
- Update CG Dockerfile


### short-read-mngs-v6.4.0 short-read-mngs-v6.3.0
- short-read-mngs auto_benchmark: measure precision/recall wrt taxa truth sets (#83)
- Update README.md
- short-read-mngs/auto_benchmark/ref_libs: add three CAMI results
- short-read-mngs/auto_benchmark/run_dev.py: make workflow_version more consistent
- Fix typo in RemoveHost no reads left check (#81)
- Prepare to migrate JSON files out of idseq-dag (#80)
- Add truth file for idseq-bench-5
- Add benchmark truth files to manifest (#78)
- Update Entrez contact address (#72)

### short-read-mngs-v6.2.0
- filter taxids out of hit summaries (#73)
- increase sj collapsed limit (#75)

### consensus-genome-v1.5.0
- enable the ability to create generic consensus genomes (#71)
- Hitsummary + BlastOutput6 Parsing Refactor (#68)

### consensus-genome-v1.4.2
- [CH-27791] Add check for "No reads after MakeConsensus" (#67)

### short-read-mngs-v6.1.6
- fail if kSNP3 produces empty output (#66)

### short-read-mngs-v6.1.5
- Low hanging fruit optimization of cluster file parsing (#65)
- idseq_dag/util/fasta.py:sort_fastx_by_entry_id yet another fix for sorting with various quality score encodings
- add CAMI challenge benchmark inputs
- Added logs to deduplication step (#64)

### short-read-mngs-v6.1.4
- idseq_dag/util/fasta.py: don't mangle FASTQ read names with pipe characters (#63)
- Revert "idseq_dag/util/fasta.py: don't mangle FASTQ read names with pipe characters (#62)"
- idseq_dag/util/fasta.py: don't mangle FASTQ read names with pipe characters (#62)
- auto_benchmark: surface build warning if deviations detected
- flesh out auto_benchmark README

### short-read-mngs-v6.1.3
- add detail to tests/README.md
- Fix count on short check (#61)

### short-read-mngs-v6.1.2
- Stop filtering short reads in validate input and make count more reliable (#60)
- auto_benchmark: fix reference to deleted branch
- short-read-mngs auto benchmarking framework (#59)
- Add instructions to run CG locally (#58)

### consensus-genome-v1.4.1
- Fix glob on trim reads that would not catch all output files (#57)

### short-read-mngs-v6.1.1
- Fix logical condition on length fix (#56)
- local alignment: extract db to temp space instead of output dir (#54)
- Fix file ref (#55)
- github CI: upgrade deprecated set-env (#53)

### short-read-mngs-v6.1.0
- Compute and export alignment results from combining NT and NR hits (#51)

### short-read-mngs-v6.0.1
- idseq dedup cluster header bugfix (#52)

### short-read-mngs-v6.0.0
- Replace cdhit-dup with idseq-dedup (#50)

### consensus-genome-v1.4.0
- Add samtools_depth.txt to the output zip (#49)

### consensus-genome-v1.3.0
- Add support for single-end fastqs (#48)

### consensus-genome-v1.2.5
- Remove dag_branch from consensus-genome (#47)

### short-read-mngs-v5.2.1
- short-read-mngs unit tests v2 (#46)
- Remove dag_branch parameter (#40)
- add short-read-mngs end-to-end test to CI (using viral databases) (#44)
- Catch shell pipeline failures in run_srst2 (#43)

### consensus-genome-v1.2.4
- Emit depths file (#45)

### short-read-mngs-v5.2.0
- RunSRST2 unconditionally again (#42)

### short-read-mngs-v5.1.1
- Catch shell pipeline failures in run_srst2
- switch to a warning for multiple emitted reads from the same cluster (#41)

### short-read-mngs-v5.1.0
- local driver wdl for short-read-mngs (#39)
- Rename main to short-read-mngs in Readme

### consensus-genome-v1.2.3
- Fix issues with consensus genome pipeline (#38)

### short-read-mngs-v4.12.2
- Remove unused logo and badges
- Fix copy command
- Install bundled idseq-dag instead of repo archive
- Merge idseq-dag repository contents
- Rename main workflow to short-read-mngs

### vconsensus-genome-1.2.1
- Fix issues with error reporting (#32)
- Use unshallow checkout; pass tag to deployment
- upgrade idseq-dag to v4.11.8 (#31)
- fix version 4.11.8 (#324)
- support local alignment (#321)

### vconsensus-genome-1.1.0
- support local non host alignment (#25)
- Add "no reads left" errors (#28)
- Create CONTRIBUTING.md (#24)
- s/idseq-database/idseq-public-references (#323)
- Switch to public references (#19)
- Warn only about the first failure to upload step status JSON (#322)
- Simplify CG docker image and update a few packages (#27)
- add host_filter RunStar test case using smaller ERCC database (#22)
- add miniwdl download cache cleaner script to main docker (#23)
- miniwdl v0.8.1 handles SwarmContainer.global_init (#26)

### v1-consensus-genome
- Add consensus genome workflow (#21)
- Readd STARlong (#20)
- Re-add LICENSE file
- Run small/medium S3 downloads through miniwdl instead of idseq-dag (#17)
- task test skeletons for postprocess & experimental (#18)
- Bump idseq-dag to 4.11.6
- Raise InsufficientReadsError on empty input file; add test (#320)
- add auto-untar to 'local additional files' workaround (#318)
- Alignment unify (#319)
- Update head branch in CI
- Keep curl in idseq-wdl image (#16)
- Bump idseq-dag to v4.11.4
- count_reads: correctly count reads in empty file (#317)
- Add rapsearch2 to docker image (#11)

### v1.0.2
- Fix S3 path parameterization for nt_loc.db reference file
- add e-value threshold to require internal alignments have e-value < 1 (#309)
- Simplify WDL input (#6)
- run `miniwdl check` through pre-commit (#2)
- basic host_filter test (#4)

### v1
- Fix passing of boolean JSON literals
- Raise error on pip install failure
- Port workflow conditional logic from idseq-web (#3)
- new version for index update (#307)
- Bump version (#306)
- IDSEQ-2372: Fix contig read count (#304)
- Fix error in counting of unidentified.fa (#305)
- Begin IDseq WDL workflows
- Add step entry point (#303)
- Add nonunique reads to unmapped fasta for download (#300)
- Add nonunique reads to nonhost fastqs for download (#299)
- Determinism (#302)
- Disable reusing previous chunks when re-running alignment (#301)
- Smarter contig assignment to accessions (#297)
- Increase Alignment Chunk Size and Decrese Chunks in Flight (#298)
- Remove index generation + sqlite (#296)
- Remove use_taxon_whitelist requirement for betacoronavirus fastqs (#288)
- Added spot attempt with on demand fallback to alignment (#294)
- Check for filename collisions (#287)
- Refactored m8 summary to improve memory usage (#278)
- Update status by merging directly to S3 file (#291)
- Remove [Datadog] from logs (#292)
- Prod hotfix (#289)
- IDSEQ-2600: fastqs for coronavirus (#282)
- Logging + s3 download race condition fix (#285)
- Catch invalid gzip files in run_validate_input step (#284)
- Set superkingdom_name to None if the dag passes an empty string in Phylo Tree creation (#286)
- Alignment batch Unrevert (#281)
- Revert "Alignment batch (#276)" (#279)
- Alignment batch (#276)
- Apply deuterostome_db, taxon_whitelist, taxon_blacklist, human filter to contigs (#272)
- Add note about git tags for versioning (#273)
- remove dead dag hacks (#268)
- Move taxon whitelisting to taxon count method (#270)
- Upgrade to 4.0.0, switch on cdhitdup change (#271)
- fix v4 counts for nonhuman hosts human filter steps (#267)
- tolerate assert failures in SPADES assembler (#266)
- Revert "new read counting mode for v4, cd-hit-dup fix for unpaired reads (#252)" (#265)
- new read counting mode for v4, cd-hit-dup fix for unpaired reads (#252)
- Switch gsnap to s3parcp (#255)
- Add support for a custom taxon whitelist mode + fix blacklist as well (#263)
- IDSEQ-2336: Lz4 compress more files (#262)
- Removed final instance of optional_files_to_upload (#261)
- STAR title switch (#260)
- Insert Metrics Copy Update + Failure Resiliancy (#257)
- Handle null nucleotide type (#256)
- Put cdhitdup clstr file in additional_output_files_visible (#254)
- Upload clstr files with cdhitdup output (#251)
- Compute insert size metrics for humans only (#248)
- Append [Datadog] to AmrAlleleMismatchError log (#253)
- Tag idseq 3.18.1
- steps.run_validate_input: do not clobber input files (#250)
- Update version marker for index update to 2020-02-10 (#249)
- Update NCBI index databases to those downloaded on 2020-02-03 (#247)
- Alignment isolation and cleanup (#243)
- Fix rapsearch2_index generation (#245)
- Fail gracefully when no reads are left after LZW filter (#244)
- Add checksum uploads for rapsearch (#241)
- Fix boto import + concurrency issue (#242)
- complete the implementation of LRU cache policy for reference downloads (#240)
- Add check for organism specific gtf (#236)
- Raise default Rapsearch timeout from 1 hr to 3 hrs (#239)
- okay_if_missing should not mean okay_if_corrupt (#238)
- remove excess concurrent logging (#237)
- ensure subprocesses write to separate files and remove excess logging to prevent hangs (#235)
- All Hosts + RNA insert size metrics (#234)
- v3.14.5
- Allow keyword arguments to be passed through in step constructor (#233)
- idseq_dag.util.log: Use named logger (#232)
- Update logo on README.md (#231)
- reduce memory use of blast reranking algorithm using generators (#230)
- reduce BATCH_SIZE for blastn to prevent OOM when the input contains many short low complexity contigs (#229)
- fix coverage extent bug seen in coverage viz (#228)
- Make the contents of idseq-dag installable via pip/setuptools (#227)
- cache AWS credentials for aws s3 cp subcommands (#224)
- Update pull_request_template.md (#226)
- Create pull_request_template.md and explain versioning (#225)
- Generate info DBs in shelve format (#221)
- Average insert size (#215)
- amendment to the fix in PR 222 to cover an additional clobbering scenario (#223)
- fix clobbering of references that occurs when the etag mismatches an earlier download from the same job (for example if the same reference is downloaded twice, from different source paths) (#222)
* remove code paths that execute frequent concurrent s3 get-range operations to fetch reference sequences a-la-carte (#219)
- Make phylo tree and alignment viz steps more robust to index errors. (#220)
- Initialize boto3 S3 client from scratch every time
- IDSEQ-1759 Check reference cache against S3 object etag (#218)
- rerank blastn (NT) results to improve query coverage for longer contigs (#199)
- Fix loose read counting bug with coverage viz. (#217)
- Add GitHub Actions, flake8, fix lint errors (#214)
- Bump version (#213)
- Add missing dependency (#212)
- Initial commit
- Add explicit output file to lz4 (#208)
- Flexible Host Genome Index Generation (#203)
- Modify LZW score adjustment heuristic for long reads to reduce the number of long, low-complexity sequences (#207)
- Log error on amr allele mismatch. (#206)
- Generate LZ4 files (#204)
- Revert "Increment pipeline version for new NCBI indexes"
- Increment pipeline version for new NCBI indexes
- Add back berkeley db file generation after sqlite (#202)
- Increase gsnapl threads from 36 to 48 for better utilization of r5d.metal machines (#201)
- Remove deprecated NCBI files from input (#200)
- Handle no 'name' in a dag file (#198)
- IDSEQ-1326 - [Stabilization] Performance optimization for coverage_out (#197)
- IDSEQ-1152 - Avoid corrupted upload files (idseq-dag) (#196)
- Add a 36-nt filter for nucleotide alignments (#194)
- Use sorted mode of bedtools coverage to reduce RAM usage. (#195)
- Fix typo in phylo tree pipeline step. (#193)
- Revert file merged by mistake (#192)
- Fix an issue with ncbi index regarding oddly annotated descriptions (#191)
- Fix experimental step failure for unpaired sample inputs (#189)
- Fix an issue with bedtools command invocation (#188)
- IDSEQ-1132 - Command patterns refactoring (#187)
- Switch from using out target name to service (#185)
- Boris/blacklist project 549 (#184)
- disable project 549 (#183)
- Fix use of wrong method to add to Python list (#182)
- Fix typo in run_srst2 (#181)
- Add total reads, reads per million, and depth per million to SRST2 output (#175)
- Record start time not at instantiation, but at running start time (#179)
- Update Docstrings (#178)
- Reapply "Upload step-level updates during pipeline run" (#174)
- Revert "Upload step-level updates during pipeline run (#172)" (#173)
- Upload step-level updates during pipeline run (#172)
- Fail gracefully if read pairing is broken (#171)
- Add support to a new header pattern (#170)
- improve rapsearch and gsnap chunk retry logic
- mutex around sqlite open implicated in hangs
- change description from unbiased to hypothesis-free (#166)
- use lz4, limit s3mi concurrency, switch from dash to bash (#167)
- minimum updates to SQLite wrapper to ensure stability (#157)
- catch invalid fastq during read counting (#165)
- catch gunzip errors (#164)
- catch gunzip errors (#163)
- RunAlignmentRemotely retry pattern (#161)
- Add input file validation for pipeline steps. (#153)
- Rename `timestamp` field name. (#160)
- Log output in JSON format. Change TraceLock log level to DEBUG. (#159)
- Upgrade to python 3.7.3 (#155)
- Add progress logs to GenerateCoverageStats step (#156)
- Remove db_hack. Standardize db_open/db_assert_table/db_close log entries. (#152)
- Fix division by zero error in generate_coverage_viz. (#151)
- modify Trimmomatic command; allow small insert-size (rev. complement) reads through QC (#149)
- Another attempt to fix the sqlite3 multiprocessing lock issue (#150)
- Fix an issue with the log event when trying to log non json serializable fields (#148)
- Open sqlite3 dbs in readonly mode whenever possible (#147)
- Add logs to m8.generate_taxon_count_json_from_m8 and remove redundant locks. (#146)
- Address bug in rounding error code. (#145)
- Logging RLocks (#144)
- [Version 3.6] Add GenerateCoverageViz pipeline step to dag. (#133)
- Log improvements (#143)
- fix inconsistent behavior in populate_reference_sequences (#141)
- Optimize downloads and uploads in alignment viz step (#139)
- Sanitize user-facing error messages (#138)
- Generate nt_info and nr_info dbs that map accession id to name and length. (#137)
- get reference accessions from hitsummary2 instead of align_viz files (#134)
- handle empty m8 (#136)
- fix unbound variable (#132)
- Indexing steps (#131)
- 32 to 16 threads to avoid blastx flakiness (#129)
- implement upstream and downstream star classes (#127)
- optimize align viz opt (#128)
- Shelve 2 sqlite (#126)
- Generate host genome (#125)
- fixing newline issue w/ phylo trees (#123)
- bug fix (#122)
- Check read count between two files (#121)
- output
- Adding "validate input" step to host filtering and touching up pipeline to accept FASTA and PacBio files (#118)
- use nmdb specifically for shelve db (#120)
- Merge branch 'master' of github.com:chanzuckerberg/idseq-dag
- move dag
- Generate nonhost fastq files. (#119)
- Merge branch 'master' of github.com:chanzuckerberg/idseq-dag
- add suffix
- deploy 2018-10-26 gsnap to prod and bump version to 3.3.0 (#114)
- gsnap 2018-10-26 pre-release testing on staging only (#113)
- revert 3.2.3 (#112)
- pre-release testing of gsnap 2018-10-26 on staging (#111)
- Example running accession2taxid (#110)
- Custom db index building  (#107)
- Accession2taxid (#106)
- revert 3.2.1 (new gsnap for staging) (#105)
- handle null case (#103)
- update gsnap to 2018-10-20 prerelease by twu (#102)
- Assembly improvement (#101)
- Change username for gsnap running on Amazon Linux 2 (#98)
- skip fasta for trimming (#100)
- set defaults for chunk scheduling properly (#95)
- scale-in / do not dispatch to draining servers, add job-tag to selected server (#91)
- readme
- Generic assembly (#86)
- fix example (#94)
- remove quality trimming but keep adapter trimming (#93)
- add threshold_readlength argument (#89)
- remove genbank (#92)
- Fix issue with PriceSeqFilter outputing multiline FASTA files (#90)
- Merge branch 'master' of github.com:chanzuckerberg/idseq-dag
- update docker file
- require minimum LZW vocabulary length (JIRA-625) (#88)
- Update README.md
- Update README links
- replace pipeline_run_ids by sample names (#83)
- quality and adapter trimming (#82)
- parallelize adapter trimming (#85)
- fix bug
- generate and upload vcf (#81)
- Better wiki coverage (#78)
- drop reads only if EVERY accession is on the blacklist (#74)
- Update README.md
- Update README.md
- Create SECURITY.md
- Download taxid descriptions (#73)
- unpack align viz dict (#72)
- different k for bacteria and viruses (#71)
- get NCBI metadata (#70)
- generate SNP annotations (#68)
- Reclassify skeleton (#67)
- use both NT and NR hits for phylo trees (#64)
- restore NCBI reference nodes (#66)
- Update generate_phylo_tree.py (#65)
- Add blacklist (#63)
- Indicate where indexing scripts live (#62)
- Suppress s3api calls output (#59)
- acquire lock before fork (#61)
- Removing amr test script
- Removed amr test dags
- Meera add amr (#44)
- add comment for the genus-level taxid case (#60)
- Prototype for outbreak phylogenetic trees (#54)
- Fail at the right stage (#58)
- i3metal instances have 36 physical cores (#56)
- parallelize lzw (#52)
- Make fasta work for pipeline (#51)
- Tweak (#50)
- Setting up integration tests (#47)
- Version alignment index (#49)
- fix readme typo
- Upload with retries (#48)
- subprocess (#46)
- Use gsnap and rapsearch's private instead of public IPs  (#45)
- Update folder (#43)
- Idseq dag documentation (#42)
- try this (#41)
- Better error handling (#40)
- Revert align viz to using Unix commands (#39)
- Copy over some release notes (#38)
- WIP-ish - Align viz bug fix (#37)
- up version fix lzw
- change file type
- master
- exit when a step failed
- Merge branch 'master' of github.com:chanzuckerberg/idseq-dag
- truncation
- Tweak NT downloading (#36)
- docker file update. upgrade to python 3.6
- docker file move from idseq-pipeline
- add version
- Index changes for the m8 (#35)
- Clean changes (#34)
- bug fix (#33)
- Non host alignment (#27)
- Style tweaks and delete empty file remnant (#31)
- Move test JSONS to new files (#32)
- counts files (#28)
- Alignment viz (#26)
- count reads after (#30)
- Copy postprocess steps (#29)
- counts file (#25)
- Create LICENSE (#24)
- Generate taxid locator (#23)
- Allow subdir structure (#22)
- Generate taxid fasta (#20)
- Subsampling (#21)
- Bowtie2 & GsnapFilter (#19)
- Integrate run STAR step (#18)
- Cdhit dup lzw step (#17)
- Add PipelineStepRunPriceSeq (#16)
- Additional files (#15)
- Node to target (#14)
- Fix download bug (#13)
- fix logger (#12)
- quick diff with retry (#11)
- Minor log changes (#8)
- Move downloading uploading over (#7)
- Add logging and command.execute functions (#6)
- Merge pull request #3 from chanzuckerberg/tests_bowtie2
- some quick fix
- change head_nodes to hash format
- bug fixes
- idseq_step_setup
- step setup
- Merge branch 'master' into tests_bowtie2
- fix dag json to the right order
- add a base class for all the unittests
- a test that runs
- Merge pull request #2 from chanzuckerberg/add_name_to_step
- Merge pull request #1 from chanzuckerberg/boris/2018_06_05_a
- add name
- better conform with python3 style
- test bowtie2
- test_bowtie2
- change version
- add validation step. fix input file bug
- make font smaller
- quick readme
- fix typo
- dummy steps
- linting
- lint
- fix some syntax errors
- implementing pipeline_step.py
- pipeline_step
- some further improvement
- move conf validation to pipeline_flow
- change from step to class
- fix some bugs
- hook into main
- main pipeline_flow function
- make sure all the nodes are covered by the steps
- idseq_dag
- steps
- add a realistic example
- build up
- more detailed example
- more realistic sample
- add a new dag json
- wip pipeline flow/step
- add unittests
- python
- initial checkin


Previous releases to the short-read-mngs pipeline can be found [here](https://github.com/chanzuckerberg/idseq-dag#release-notes)
