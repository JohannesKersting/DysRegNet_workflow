# DysRegNet analysis workflow  

This repository contains the analysis workflow for the DysRegNet manuscript: [DysRegNet: Patient-specific and confounder-aware dysregulated network inference](https://doi.org/10.1101/2022.04.29.490015). 
The workflow is implemented using the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system.
The workflow will perform the following steps:
- Data preprocessing
- Reference network inference
- Patient-specific network inference with [DysRegNet](https://github.com/biomedbigdata/DysRegNet_package) and [SSN](https://academic.oup.com/nar/article/44/22/e164/2691334?login=false)
- Evaluation and plot generation
- Build the database for the [DysRegNet web application](https://exbio.wzw.tum.de/dysregnet/)

The main workflow is orchestrated in the `workflow/Snakefile`.
The `profiles` directory contains configuration files for different execution environments (e.g., slurm, local).
The `workflow\scripts` directory contains all analysis scripts used by the workflow. 

## Usage  
  
After cloning the repository, install the required packages using [conda](https://docs.anaconda.com/free/anacondaorg/user-guide/index.html) and the environment.yaml file: 
```bash  
conda env create -f environment.yaml
```  

Activate the environment:
```bash  
conda activate drn
```  

If necessary, adjust data (expected to contain all input files) and result directories specified by the `DATA_DIR` and `RESULT_DIR` variables in `workflow/Snakefile`.

Run the workflow via slurm:
``` bash
snakemake --profile profiles/slurm
```

The slurm  profile contains cluster-specific configuration. For development and pipeline testing, a local executor and/or the dry-run feature can be used:
``` bash
snakemake --profile profiles/local -n
```
