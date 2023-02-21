# Introduction

This repository contains code and supplementary materials for paper: *Methyl ...*.

Authors: Jose Maria Zamora-Fuentes, Jesus Espinal-Enriquez, Enrique Hernandez-Lemus

## Pre-requisites

### R

Considerations:

- R (3.6.3)

Pre-requisites to run scripts in theses paheses are obtained with:

`$ Rscript pkgs-requiremnts.R`

### CONDA

```
bash Miniconda3-py38_4.12.0-Linux-x86_64.sh -b -p ~/bioconda38
export PATH=~/bioconda38/bin:$PATH
conda install -y numpy scipy matplotlib pandas ipython
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
# mamba init
conda activate snakemake
```

## Running and building whole project

`$ snakemake -c all`

## Data to start 

Download from link provided file: `data-tostart.tgz` and untar into the project.

`$ tar xzvf data-tostart.tgz`
`$ mv data-tostart  data`