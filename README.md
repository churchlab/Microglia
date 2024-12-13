# Microglia

# Gene-TF Linear Regression Analysis

This repository contains an R script for analyzing gene expression data to identify relationships between target genes and transcription factors (TFs) using linear regression models. The script evaluates TF interactions at various levels and outputs significant findings for further analysis.

## Features

- **Multilevel Linear Regression**: Fits models with up to 4th-order interactions among TFs.
- **Significance Testing**: Filters results based on p-value thresholds (< 0.01).
- **Automated Workflow**: Processes multiple genes in a single run.
- **Customizable Output**: Saves results in organized tab-separated files.

---

## Prerequisites

- **R**: [Download R](https://www.r-project.org/) if not already installed.
- **Required Files**:
  - `input_files/tf_allgenes_normalized_matrix_2020-01-10.rds`: Normalized matrix of gene and TF expression data.
  - `SCT_genes_index.rds`: List of gene indices to analyze.

  *(Optional)* Additional input files include:
  - `PGP1_QCed_raw_counts_sparse.rds`: Sparse gene expression matrix.
  - `PGP1_QCed_raw_counts_df.rds`: Dense gene expression matrix.
  - `pearson_cor_pgp1_2019-12-26.rds`: Precomputed Pearson correlation matrix.

---

Input Files

Normalized Gene-TF Expression Matrix:
File: tf_allgenes_normalized_matrix_2020-01-10.rds
Content: A normalized matrix of gene and TF expression levels.
Gene Index:
File: SCT_genes_index.rds
Content: A list of gene indices to analyze.
(Optional) Additional Input Files:
Sparse gene matrix: PGP1_QCed_raw_counts_sparse.rds
Dense gene matrix: PGP1_QCed_raw_counts_df.rds
Precomputed correlation matrix: pearson_cor_pgp1_2019-12-26.rds
Output

For each gene analyzed, the script generates:

Regression Results with Significant TFs:
Format: <gene_name>_lm_res_fit[1-4].txt
Content: Coefficients and p-values for TFs with p-value < 0.01.
Single TF Regression Results:
Format: <gene_name>_lm_res_singletf_fit[1-4].txt
Content: Results of single TF contributions for each model.
Outputs are saved in the output_files_3000/ directory.

Script Details

Gene Selection:
Iterates through gene indices specified in SCT_genes_index.rds.
Data Preparation:
Extracts normalized gene-TF matrix for each gene.
Sets the gene of interest as gene1.
Regression Models:
Model 1: Linear regression without interaction terms.
Model 2: Includes 2nd-order TF interactions.
Model 3: Includes 3rd-order TF interactions.
Model 4: Includes 4th-order TF interactions.
Significance Filtering:
Extracts TFs with p-values < 0.01 for each model.
File Outputs:
Saves significant results and single TF results as tab-separated files.
Usage

Run the script in R:

R
Copy code
source("<script_name>.R")
Ensure that all input files are placed in the input_files/ directory and the output directory exists.

Notes

Adjust the index variable in the script to specify a custom range of gene indices.
Uncomment optional input files in the script if additional data processing is required.
Ensure sufficient computational resources for large datasets.



## Directory Structure

Ensure the following directory structure is in place:

```plaintext
project_root/
├── input_files/
│   ├── tf_allgenes_normalized_matrix_2020-01-10.rds
│   ├── SCT_genes_index.rds
│   └── (optional) other input files
├── output_files_3000/   # Directory for output files
├── analysis_script.R    # The script provided
