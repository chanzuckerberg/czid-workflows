#!/bin/bash
#!/bin/bash
idseq_dag examples/amr_test_dags/fastq.gz_paired_rds/amr_human_sample_paired_rds_dag.json --no-lazy-run 
echo "====Done with fastq.gz_paired_rds/amr_human_sample_paired_rds_dag.json===="
idseq_dag examples/amr_test_dags/fastq.gz_paired_rds/amr_water_sample_paired_rds_dag.json --no-lazy-run
echo "====Done with fastq.gz_paired_rds/amr_water_sample_paired_rds_dag.json===="
idseq_dag examples/amr_test_dags/fastq.gz_single_rds/amr_human_sample_single_rd_dag.json --no-lazy-run
echo "====Done with fastq.gz_single_rds/amr_human_sample_single_rd_dag.json===="
idseq_dag examples/amr_test_dags/fastq.gz_single_rds/amr_water_sample_single_rd_dag.json --no-lazy-run
echo "====Done with fastq.gz_single_rds/amr_human_water_single_rd_dag.json===="
idseq_dag examples/amr_test_dags/fastq_paired_rds/amr_human_sample_paired_rds_dag.json --no-lazy-run
echo "====Done with fastq_paired_rds/amr_human_sample_paired_rds_dag.json===="
idseq_dag examples/amr_test_dags/fastq_paired_rds/amr_water_sample_paired_rds_dag.json --no-lazy-run
echo "====Done with fastq_paired_rds/amr_water_sample_paired_rds_dag.json===="
idseq_dag examples/amr_test_dags/fastq_single_rds/amr_human_sample_single_rd_dag.json --no-lazy-run
echo "=====Done with fastq_single_rds/amr_human_sample_single_rd_dag.json======"
idseq_dag examples/amr_test_dags/fastq_single_rds/amr_water_sample_single_rd_dag.json --no-lazy-run
echo "=====Done with fastq_single_rds/amr_water_sample_single_rd_dag.json======"
idseq_dag examples/amr_test_dags/fasta_rds/amr_fasta_sample_dag.json --no-lazy-run
echo "=====Done with fasta_rds/amr_fasta_sample_dag.json======"
idseq_dag examples/amr_test_dags/fasta.gz_rds/amr_fasta.gz_sample_dag.json --no-lazy-run
echo "=====Done with fasta.gz_rds/amr_fasta.gz_sample_dag.json======"	
 
 
 
 

