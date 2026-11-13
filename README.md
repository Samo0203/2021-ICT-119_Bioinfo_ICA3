# Leukemia Gene Expression Analysis: ALL vs AML

## 1. Project Overview

This project analyzes gene expression data from patients diagnosed with Acute Lymphoblastic Leukemia (ALL) and Acute Myeloid Leukemia (AML).

The main objective is to identify gene-expression probes that show statistically significant differences between these two leukemia types using bioinformatics techniques in R.

## 2. Research Question

Which gene-expression probes show significant differences between ALL and AML patients?

## 3. Objectives

- Load and explore a publicly available leukemia gene expression dataset.
- Examine the distribution of gene expression across patient samples.
- Visualize similarities and differences between patient groups.
- Perform differential expression analysis using the limma package.
- Identify significantly differentially expressed probes.
- Visualize the analysis results using volcano plots and heatmaps.

## 4. Dataset Description

This project uses the Golub leukemia dataset available through the Bioconductor `golubEsets` package.

| Description | Value |
|---|---:|
| Total expression probes | 7,129 |
| Total patient samples | 72 |
| ALL samples | 47 |
| AML samples | 25 |

The dataset contains gene expression measurements from patients diagnosed with ALL and AML.

## 5. Technologies and Packages

**Programming language:** R

**Development environment:** RStudio

**R packages:**

- `golubEsets` – Accessing the Golub dataset.
- `Biobase` – Working with Bioconductor expression data.
- `limma` – Differential expression analysis.
- `ggplot2` – Data visualization.
- `dplyr` – Data manipulation.
- `tidyr` – Data organization.
- `pheatmap` – Heatmap visualization.
- `ggrepel` – Plot labeling.
- `knitr` – Report generation.
- `rmarkdown` – Reproducible reporting.

## 6. Methodology

### 6.1 Data Loading and Preprocessing

- Loaded the Golub leukemia dataset.
- Extracted the gene expression matrix and patient information.
- Identified ALL and AML patient groups.
- Checked missing and infinite expression values.
- Saved the processed data for further analysis.

### 6.2 Exploratory Data Analysis

Performed exploratory analysis using:

- Patient distribution bar chart.
- Gene expression distribution boxplot.
- Principal Component Analysis (PCA).
- Gene expression heatmap.

### 6.3 Differential Expression Analysis

The `limma` package was used to compare gene expression between ALL and AML groups.

The comparison was defined as:

`AML - ALL`

The Benjamini-Hochberg method was used to adjust p-values for multiple testing.

**Significance criteria:**

- Adjusted p-value (adj.P.Val) less than 0.05.
- Absolute log fold change greater than 1.

## 7. Results

The analysis identified **1,072 significantly differentially expressed probes**.

| Result | Number of probes |
|---|---:|
| Total significant probes | 1,072 |
| Higher expression in AML | 531 |
| Higher expression in ALL | 541 |

Positive log fold change values indicate higher expression in AML, while negative values indicate higher expression in ALL.

These results represent expression probes, not necessarily unique genes.

## 8. Visualizations

The project generates the following visualizations:

- Patient distribution.
- Expression distribution.
- PCA plot.
- Exploratory gene expression heatmap.
- Differential expression volcano plot.
- Differential expression heatmap of the top 50 significant probes.

## 9. Project Structure

```text
Leukemia_Gene_Expression_Analysis/
│
├── code/
│   ├── 02_data_preprocessing.R
│   ├── 03_exploratory_analysis.R
│   └── 04_differential_expression.R
│
├── data/
│   ├── expression_data.rds
│   ├── sample_information.csv
│   └── patient_groups.csv
│
├── results/
│   ├── figures/
│   └── tables/
│
├── report/
│
├── README.md
└── Leukemia_Gene_Expression_Analysis.Rproj
```

## 10. How to Run the Project

1. Install R and RStudio.
2. Clone or download this repository.
3. Open `Leukemia_Gene_Expression_Analysis.Rproj` in RStudio.
4. Install the required R packages if they are not already installed.
5. Run the scripts in the `code/` directory in numerical order.
6. View the generated tables in `results/tables/` and visualizations in `results/figures/`.

## 11. Conclusion

This project demonstrates how bioinformatics techniques can be used to analyze leukemia gene expression data and identify differences between ALL and AML patient groups.

Using the limma statistical framework, 1,072 probes met the selected differential expression criteria. The findings provide a basis for further investigation of expression differences between the two leukemia subtypes.

## 12. Data Source

Golub leukemia dataset, accessed through the Bioconductor `golubEsets` package.
