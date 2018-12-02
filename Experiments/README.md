# Reproducing Liu, et al. (_Nature Communications_, 2017)

Source code for reproducing signature discovery from Liu, et al. [1]. Briefly, our goal is to reproduce Figure 2a-2c.
Additional details of the experiments, data, and conclusions can be found in the [`report/`](report/).

## Setup

The source code is written in Python 3. We use `snakemake` to manage the workflow. We suggest using Conda to install dependencies, which you can do directly from the provided [`environment.yml`](environment.yml) file.

    conda env create -f environment.yml
    source activate reproducing-liu2017

We have included git modules that you will need to download before you can run `snakemake all`

    git submodule init
    git submodule update

## Usage

To download the input data files and reproduce the experiments and..., simply run:

    snakemake all

Since the cohort considered is quite small (30 samples), the runtime required to produce all figures is less than 5 minutes.

### Configuration

Additional configuration options are detailed at the beginning of the [`Snakefile`](Snakefile).
There are also various command-line options for the provided scripts in [`src/`](src).

## References
1. Liu, et al. (2017) "Mutational patterns in chemotherapy resistant muscle-invasive bladder cancer." _Nature Communications_ **8**, Article number: 2193. [doi: 10.1038/s41467-017-02320-7](https://doi.org/10.1038/s41467-017-02320-7)
