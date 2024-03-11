# Diamond Scatter

`diamond_scatter.py` works with a [modified version of diamond](https://github.com/morsecodist/diamond) to use diamond's parallelization to work in separate jobs that can be run at any time instead of on an HPC cluster with a shared file system. This allows us to run alignment on batch with spot instances while dynamically scaling to meet uneven demand more easily.


The approach was to change the diamond code base as little as possible. That is why the `diamond_scatter.py` script exists, to control some aspects of the diamond parallelization management without adding it to diamond. If we want to fold this into diamond we would want to add what is in this script in addition to the fork but I am not sure this sort of parallelization is a goal of the diamond project so they may not be receptive.

# How Diamond Parallelization Works

Refer to [the docs](https://github.com/bbuchfink/diamond/wiki/6.-Distributed-computing) for a detailed explaination.

Diamond expects to run in different processes with a file system accessible by all processes. The directory contains a record of what work needs to be done that the different processes can update to communicate with one another. The diretory also allows the sharing of files between the processes. A parallel run starts with running an initialization command that creates the list of work needed for a particular alignment. Then, you can run as many processes as you want. Each processes will pick up a chunk of alignment work based on the work in the list, and will update the list when it begins. Each alignment chunk produces some intermediate files that they store in the directory. The final chunk, is able to detect it is the final chunk based on the work list and it has access to all of the intermediate files so it joins them.

# How Diamond was Modified

To make diamond run in parallel without a shared directory you need two main elements:

1. We need to split the joining from the alignment because the alignment chunk won't have all of the intermediate files. We want to be able to run the join later after we've gathered up all the intermediate files.
2. We need diamond to be able to run a single, specific, chunk and output intermediate files without deleting them.

To achieve 1 diamond was modified slightly to split the alignment command in two by adding flags: `--single-chunk` and `--join-chunks`. The `single-chunk` skips the checking for joining and cleaning up intermediate files and `join-chunks` skips alignment and goes right to that check.

Now all you need is to make sure the directory is in the correct state at each step of the processes. That is where the `diamond_scatter.py` script comes in. The script fills the role of initializing the parallel directory.

To do a single chunk of alignment work, instead of specifying all of the work in the work list the script just specifies that one chunk needs to be completed, but a specific chunk with a particular number. This is necessary because the intermediate file naming scheme uses chunk numbers so they will be needed for the final join. The modified diamond won't clean these up so you can copy them into remote storage and join at any time.

To join the chunks you just need to copy all of the intermediate files onto the node you are using to join. Thie script initializes the directory with the correct format and runs with the `join-chunks` flag.
