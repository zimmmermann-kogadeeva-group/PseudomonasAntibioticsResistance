
[![DOI](https://zenodo.org/badge/908679580.svg)](https://doi.org/10.5281/zenodo.14559333)

# Pseudomonas antibiotic resistance project

All the scripts to reproduce all the analysis and figures in the paper.

Data is stored using git-lfs in this repo.

Alternative download [link](https://oc.embl.de/index.php/s/1IwrFS6Khohg3p1).

## Setup

The following tools are needed:
- conda
- snakemake >= 8.0
- apptainer >= 1.3
- git
- git-lfs

## Running the pipeline

First download the transcriptomics and WGS data from [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJEB76120).

After placing the data in the `Data` directory, run the pipeline with the following command:
```
snakemake --workflow-profile Profiles/Slurm
```

Then create conda environment needed to run the notebook `main.ipynb` with:
```
conda env create -n omics -f conda_env.yml
```
Activate it with
```
conda activate omics
```
Next, open jupyter-lab with 
```
jupyter-lab
```
command and open the `main.ipynb` notebook.


