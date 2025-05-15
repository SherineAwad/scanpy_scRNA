# Zebrafish scRNA-seq Analysis using Snakemake and Scanpy

This repository contains a Snakemake-based pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data from **zebrafish retina** using [Scanpy](https://scanpy.readthedocs.io/en/stable/).  
The data used in this analysis is derived from the study:

> **Title**: Molecular architecture of the zebrafish retina  
> **Authors**: Thomas J. Pandolfi et al.  
> **Journal**: *Nature Communications* (2024)  
> **Link**: [https://www.nature.com/articles/s41467-023-44142-w](https://www.nature.com/articles/s41467-023-44142-w)

## 🔧 Tools Used

- [Snakemake](https://snakemake.readthedocs.io/en/stable/): Workflow management
- [Scanpy](https://scanpy.readthedocs.io/en/stable/): Single-cell RNA-seq analysis in Python
- [Anndata](https://anndata.readthedocs.io/): Annotated data object
- Python (≥3.8)


## 🚀 Getting Started

1. **Set up the environment**:
    ```bash
    conda env create -f envs/scRNAseq.yaml
    conda activate scRNAseq
    ```

2. **Run the pipeline**:
    ```bash
    snakemake --cores 4
    ```

## 📊 Analysis Overview

The pipeline performs:
- Quality control
- Normalization
- Feature selection
- Dimensionality reduction (PCA, UMAP)
- Clustering
- Marker gene visualization

## 📌 Notes

- The markers used in the analysis were derived from literature and curated to match zebrafish-specific genes.
- The pipeline is modular and easy to extend for other datasets or species.

---



