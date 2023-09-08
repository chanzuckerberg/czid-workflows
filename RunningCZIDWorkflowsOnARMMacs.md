## Running czid-workflows on Silicon (M1/M2) Macs

Creating the docker containers for the workflows is challenging on machines that don't use the `x86` architecture. Many of the dependencies that we use right now don't have options for building for other architectures. 

The following is how to set up a virtual machine using `colima` where you can build and run docker containers locally in machines that use the ARM CPU architecture


### Install Colima

`brew install colima`

Go to the czid-workflows repo

`cd ~/czid-workflows`

### Boot and mount the colima VM
``` 
# uses 4GB of memory & 4 CPUs, change this to suit your needs 
colima start --mount $PWD:/work:w -m 4 -c 4

```

### Set up the VM

SSH into the VM

`colima ssh`

In the VM run:

`sudo apk add py3-virtualenv gcc python3-dev musl-dev linux-headers make`

Go to where the `czid-workflows` repo is mounted

`cd /work`

Then you should be able to build and run the docker containers as normal

`make build WORKFLOW=minimap2`
