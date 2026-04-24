# Analysis - BeatAML Drug Combo Manuscript
Methylation analysis & Biomarker identification

## Setup - Methylation Analysis

R version: 4.5.3, renv 1.0.11

Refer to `init.R` for specific versions of packages used. To run the pipeline open the `methylation_pipeline.Rproj` file and initiate as follows

```r
# Restoring R environments
renv::activate()
renv::restore()

# Run the pipeline
source('methylation_pipeline.R')
```

### Data


### Functions - Methylation Analysis

| Scripts | Remarks | Files |
| --- | --- | --- |
| `scripts/01_preprocess_counts.R` | Perform standard RNA-seq preprocessing by converting raw sequencing counts into a normalized, log-transformed expression matrix | Inputs: <br/> `input/beataml_waves1to4_counts_dbgap.txt` <br/><br/> Outputs: <br/> `output/Expression_BAML_707_log2RPM.txt` |
| `scripts/02_CIBERSORT_pipeline.R` | A cellular deconvolution pipeline to estimate the proportions of specific AML (Acute Myeloid Leukemia) cell states within bulk RNA-seq samples. <br/> Modified version that return QC | Inputs: <br/> `output/BAML_Normalized.txt` <br/> `input/AML_signature_matrix_M.txt` <br/><br/> Outputs: <br/> `output/CIBERSORT_mixture_RPM.txt` <br/> `output/CIBERSORT_AML_results.csv` <br/> `output/CIBERSORT_AML_malignant_normalized.csv` <br/> `output/AML_dominant_cell_state.csv` <br/><br/> Plots: <br/> `plots/CIBERSORT_AML_heatmap.pdf` <br/> `plots/CIBERSORT_AML_PCA.pdf` |
| `scripts/03_Processing_RNA_Methylation_datasets.R` | Annotates methylation data and filters by promoters (TSS2000/TSS1500). Collapses to gene level and filters based on Figure 6F programs. | Inputs: <br/> `output/CIBERSORT_mixture_RPM.txt` <br/> `input/Methylation Beta Values_clean.csv` <br/><br/> Outputs: <br/> `output/RNA_Expression.txt` <br/> `output/Meth_Expression.txt` |
| `scripts/04_ModuleScore.R` | Calculates AML differentiation module scores and generates ordered trajectory samples. Uses Seurat for AddModuleScore and PCA. | Inputs: <br/> `output/CIBERSORT_mixture_RPM.txt` <br/> `output/AML_dominant_cell_state.csv` <br/> `input/AMLCellType_Genesets.gmt` <br/><br/> Outputs: <br/> `output/AML_ModuleScores.csv` <br/> `output/AML_ModuleScore_PCA.csv` <br/> `output/AML_Trajectory_Ordering.csv` |
| `scripts/05_RNA_heatmap.R` | Generates RNA differentiation trajectory heatmap | Inputs: <br/> `output/RNA_Expression.txt` <br/> `output/AML_Trajectory_Ordering.csv` <br/><br/> Plots: <br/> `plots/RNA_Trajectory_Heatmap.pdf` <br/> Figure 6E |
| `scripts/06_Methylation_heatmap.R` | Performs promoter methylation processing to generate a DNA methylation trajectory heatmap. | Inputs: <br/> `output/Meth_Expression.txt` <br/> `output/AML_Trajectory_Ordering.csv` <br/><br/> Outputs: <br/> `Methylation_expression.csv` <br/><br/> Plots: <br/> `plots/Methylation_Trajectory_Heatmap.pdf` <br/> Figure 6E |
| `scripts/07_Correlation_plots.R` | Generates RNA vs methylation correlation Ridge and Mirror plots based on specific gene programs. | Inputs: <br/> `output/RNA_Expression.txt` <br/> `output/Meth_Expression.txt` <br/><br/> Plots: <br/> `plots/Correlation_RidgePlot.pdf` <br/> `plots/Correlation_MirrorPlot.pdf` <br/> SuppFig 7 f-g |

<br/>
<br/>


## Setup - Biomarker Identification


### Data

| Script | Remarks | Files |
| --- | --- | --- |
| `BiomarkerIdentification.R` | Processes AML expression data using Seurat (PCA, clustering, UMAP). Integrates clinical response metadata, calculates an "Inflammation" module score, and generates a targeted Z-scored heatmap comparing Complete Responders (CR) and Non-Responders (NR). | Inputs: <br/> `Merged_AML_data.csv` <br/> `response.csv` <br/><br/> Plots: <br/> `plots/Inflammatory_Genes.pdf` <br/> Figure 7 g |
