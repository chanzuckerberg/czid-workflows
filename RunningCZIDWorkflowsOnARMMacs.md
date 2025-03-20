## Running czid-workflows on Silicon (M1/M2) Macs

Creating the docker containers for the workflows is challenging on machines that don't use the `x86` architecture. Many of the dependencies that we use right now don't have options for building for other architectures. 

The following is how to set up a virtual machine using `colima` where you can build and run docker containers locally in machines that use the ARM CPU architecture


### Install Colima

`brew install qemu colima`

Go to the czid-workflows repo

`cd ~/czid-workflows`

### Boot and mount the colima VM
``` 
# uses 4GB of memory & 4 CPUs, change this to suit your needs 
colima start --mount $PWD:/work:w -m 4 -c 4 -a x86_64

```

### Set up the VM

SSH into the VM

`colima ssh`

In the VM run:

`sudo apt-get update &&  sudo apt install python3-virtualenv gcc python3-dev musl-dev make qemu-user-static`

Go to where the `czid-workflows` repo is mounted

`cd /work`

Then you should be able to build and run the docker containers as normal

`make build WORKFLOW=minimap2`

or 

`make pull WORKFLOW=consensus-genome`

Then try running an environment with: 

```
virtualenv -p python3 .venv
source .venv/bin/activate

pip install -r requirements-dev.txt
sudo ./.venv/bin/miniwdl  run workflows/consensus-genome/run.wdl docker_image_id=czid-consensus-genome fastqs_0=workflows/consensus-genome/test/sample_sars-cov-2_paired_r1.fastq.gz fastqs_1=workflows/consensus-genome/test/sample_sars-cov-2_paired_r2.fastq.gz technology=Illumina sample="test" ref_fasta=s3://czid-public-references/consensus-genome/MN908947.3.fa -i workflows/consensus-genome/test/local_test.yml -v

```