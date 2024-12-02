
# Network-Based Method Benchmarking for Cancer Research

This repository contains the code, data, and instructions for replicating analyses from our study, *[Title of Your Paper]*. In this paper, we benchmark network-based methods for identifying functionally important genes and molecular pathways in cancer by leveraging protein-protein interaction (PPI) networks.

## Table of Contents
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Demo](#demo)
- [TCGA Case Study](#tcga-case-study)
- [Citing this Work](#citing-this-work)
- [Contact](#contact)

## Overview
Network-based methods have become essential tools in computational biology, particularly for cancer research. This repository provides all resources necessary to reproduce our benchmarking results, compare network-based methods, and analyze how they perform in the presence of topological biases at both gene and subnetwork levels. All the code and tutorials here are designed for guide uses to run benchmark study for these network-based methods on high-performance computing systems.

## Repository Structure

- `Data_Generation/`: Jupyter notebooks for generation of synthetic and pre-processing of real data

- `Running_Scripts/`: Scripts for all tested methods. Each method needs to be set up at high-performance computing system. The provided scripts are for Slurm workload manager.
- `TCGA/`: Data and scripts for nine cancer datasets used in the benchmark study. `TCGA/data` directory has the score files in different formats, such as p-value, q-value, z-values and `TCGA/scripts` directory has the code for loading results and analysing results.

## Installation

Each method need to be installed following their own guide. We recommend use Easybuild to set up the environment for each method on high-performance computing system. See details in table `network_methods_comparison.xlsx`.


## Usage

See Jupyter notebooks in `Data_Generation/` for step-by-step tutorials on how to generate synthetic datasets and pre-process protein-protein networks.

See Jupyter notebooks in `Running_Scripts/` for scripts to run each method. The provided scripts are for Slurm workload manager and did not include how to set up running environment for each method.

## Data

Due to size constraints, raw data are not included directly in this repository. Please follow the Jupyter notebook survey_data_generation.ipynb to generate synthetic datasets.

Once generated, place the generated scores and processed networks in the `data/` and `network/` folders respectively in each method's folder in `Running_Scripts` in order to benchmark them.

## Demo

Here is a Demo for running FDRnet algorithm under our benchmark framework:
 - Step 1: Follow the Jupyter notebook survey_data_generation.ipynb to generate synthetic datasets. Make `data/`, `network/`, `result/` folders in FDRnet's folder in `Running_Scripts`. Put generated lfdr files into `data/` and processed network files to `network/`.
 - Step 2: Upload all the data and code to high-performance computing system.
 - Step 3: Install the environment of FDRnet following the guideline of `https://github.com/yangle293/FDRnet`.
 - Step 4: If using Slurm, run `sbatch sbatch_fdrnet` to run batch jobs.

## TCGA Case Study

We provide the data and code for nine TCGA datasets (BLCA, LUAD, LUSC, COADREAD, PRAD, HNSC,
                 UCEC, KIPAN, BRCA) as case studies: 
- The raw datasets were downloaded from the TCGA firehose website `https://gdac.broadinstitute.org/` and preprocessed to assign a p-value for each gene.
- These studies can be submitted to high-performance computing platform using the same scripts in `Running_Scripts`.
- The resulting subnetwork can be analyzed by the scripts in the `TCGA/scripts` directory.
- You can use the `analysis_script_real.ipynb to reproduce the figures in the paper.

## Citing this Work

If you use this repository in your research, please cite our paper:

```
@article{your_paper,
  author    = {Le Yang, Runpu Chen, Steve Goodison and Yijun Sun},
  title     = {A Comprehensive Benchmark Study of Methods for Identifying Significantly Perturbed Subnetworks in Cancer},
  journal   = {TBD},
  year      = {2024},
  volume    = {TBD},
  number    = {TBD},
  pages     = {TBD},
  doi       = {DOI link}
}
```

## Contact

For questions or issues, please open an issue on GitHub or contact lyang25@buffalo.edu.

