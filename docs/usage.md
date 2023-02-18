# Usage

This page contains information on using PySEP in non-standard environments

## Running PySEP on UAF Chinook

Chinook is University of Alaska Fairbanks' (UAF) high performance computer. 
We can run PySEP on Chinook using Docker containers through 
Singularity/Apptainer. 

PySEP has been containerized directly, and any changes pushed to the repo will
trigger the container to rebuild, keeping everything up-to-date.
The following code snippet downloads the correct 
[SCOPED container](https://github.com/SeisSCOPED/pysep/pkgs/container/pysep). 
and runs the PySEP help message on Chinook.

```bash
module load singularity
singularity pull ghcr.io/seisscoped/pysep:centos7
singularity exec -c pysep_centos7.sif pysep -h
```

To run a data download we will need to mount the local filesystem into the
container using the ``--bind`` command. Using the Anchorage example event:

```bash
singularity exec -c --bind $(pwd):/home1 pysep_centos7.sif \
    bash -c "cd /home1/; pysep -p mtuq_workshop_2022 -e 2009-04-07T201255_ANCHORAGE.yaml"
```

In the above example, the `-c/--contain` flag preserves the internal container
filesystem, the `-B/--bind` flag binds the current working directory on Chinook
(i.e., pwd)to a directory called */home1* within the container, and then the 
`bash -c` command changes to this */home1* directory and runs PySEP. Files are 
subsequently saved to the local filesystem in our current working directory.

`RecSec`, the record section plotting tool, can be run from the command line 
using a similar format. With the Anchorange example files we just generated:

```bash
cd 2009-04-07T201255_SOUTHERN_ALASKA/
singularity exec -c --bind $(pwd):/home1 ../pysep_centos7.sif \
    bash -c "cd /home1; recsec --pysep_path SAC/ --min_period 10 --save record_section_tmin10.png"
```

---------------------

## Accessing LLNL Waveform Database

`PySEP` interfaces with the databases of:

* W. Walter et al. (2006)
  An assembled western United States dataset for regional seismic analysis
  ISSO 9660 CD, LLNL release UCRL-MI-222502

  Download link: https://ds.iris.edu/mda/18-001

To use the LLNL database, set the input parameter `client` as:

```yaml
client: 'LLNL'
```

