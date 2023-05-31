# Introduction

This repository contains code and supplementary materials for paper _Methylation-driven genes involved in renal carcinoma progression_

Authors: Jose Maria Zamora-Fuentes, Jesus Espinal-Enriquez, Enrique Hernandez-Lemus

## Abstract

Renal carcinomas are a group of malignant tumors often originating in the cells lining the small tubes in the kidney responsible for waste filtering from the blood and urine production.  Kidney tumors arise from the uncontrolled growth of cells in the kidneys and are responsible for a large share of global cancer-related morbidity and mortality. Understanding the molecular mechanisms driving renal carcinoma progression results crucial for the development of targeted therapies leading to an improvement of patient outcomes. Epigenetic mechanisms such as DNA methylation are known factors underlying the development of several cancer types. There is solid experimental evidence of relevant biological functions modulated by methylation-driven genes, associated with the progression of different carcinomas. Those mechanisms can be often associated to different epigenetic marks, such as DNA methylation sites or chromatin conformation patterns. Currently, there is no definitive method to establish clear relations between genetic and epigenetic factors that influence the progression of Cancer. Here we have developed a data-driven method to find methylation-driven genes, so we could find relevant bonds between gene coexpression and methylation-wide-genome regulation patterns able to drive biological processes during the progression of Clear Cell Renal Carcinoma (ccRC). With this approach we found out genes such as ITK and TSFRN9 that appear hypomethylated during all four stages of ccRC progression, and are strongly involved in immune response functions. Also we found out relevant tumor suppressor genes such as RAB25 hypermethylated, thus potentially avoiding repressed functions in the AKT signaling pathway during the evolution of ccRC. Our results have relevant implications to further understand some epigenetic-genetic affected roles underlying the progression of renal cancer.

## Pre-requisites

### R

Considerations:

- R (4.2.2)

Don't forget setting R_LIBS_USER variable to install personal packages. 

To install pre-requisites to run scripts:

`$ Rscript R/install-pkgs.R`

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

`$ snakemake -c all rule_name`

## Dataset to start 

Download file from link provided: `https://zenodo.org/record/7988316#.ZHbSrOzMKrx` and untar into the project.

```
tar xzvf data-tostart.tgz
mv data-tostart  data
```
