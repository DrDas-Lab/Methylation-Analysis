# Trajectory Analysis Pipeline - BeatAML Drug Combo Manuscript
Trajectory Analysis for the Paper

# Setup

R version: 4.5.3, renv 1.0.11

To reproduce fully, refer to `src/init.R` for specific versions of packages used

```r
## restoring renv
renv::activate()
renv::restore()
```

Note: Seurat versions (4.3.0 vs 5.1.0) used to process single cell is different from the loaded renv and not included.

# Data

Abbreviations for datasets involved

- Beat AML drug combo - bamlcombo; bamlcombi; oshu
- Beat AML 2.0 ([Tyner et al 2018](https://doi.org/10.1038/s41586-018-0623-z), [Bottomly et al 2022](https://doi.org/10.1016/j.ccell.2022.07.002)) - baml
- Beat AML Ven Combinations ([Eide et al 2023](https://doi.org/10.1158/2643-3230.bcd-23-0014)) - bamlvencombo
- Helsinki ([Malani et al 2022](https://doi.org/10.1158/2159-8290.cd-21-0410)) - fpmtb

# Workflow

## Functions

| Functions | Remarks |
| --- | --- |
| `scripts/01_preprocess_counts_no.R` | Perform standard RNA-seq preprocessing by converting raw sequencing counts into a normalized, log-transformed expression matrix |
| `scripts/02_CIBERSORT_pipeline_no.R` | A cellular deconvolution pipeline to estimate the proportions of specific AML (Acute Myeloid Leukemia) cell states within bulk RNA-seq samples. <br/> Modified version that return QC |
| `scripts/03_ModuleScore_trajectory_no.R` | Calculates differentiation module scores to quantify stemness or mature-like for each AML sample and then orders them into a linear differentiation trajectory. |
| `scripts/04_RNA_heatmap.R` | |
| `scripts/Methylation_processing.R` |  |

## Plotting

| Code for Plotting | Remarks | Files | 
| --- | --- | --- |
| `scripts/01_preprocess_counts_no.R`  |  |  |
| `scripts/02_CIBERSORT_pipeline_no.R`  |  |  | |
| `src/analysis/3.catboost.Rmd`  |  |  |  |
| `src/analysis/4.RNA_clustering.Rmd`  |  |  |  |
| `src/analysis/5.WGCNA.Rmd`  |  |  | Figure 5e <br/> SuppFig e-f |
| `src/analysis/6.scrna.Rmd`  |  | | Figure 6 <br/> SuppFig 7 |
| `src/analysis/7.clusterassoc.Rmd`  |  |  | Figure 7a-c |
| `src/analysis/8.integrating_rna_prot.Rmd` |  |  | Figure 7d |
| `src/analysis/9.biomarker_elastic.Rmd` |  | Supp Tables <br/> `output/oshu_biomarkers.csv` <br/> `output/biomarkers_elasticnet_summary.csv` <br/> Other Tables <br/> `output/supp/biomarkers_elasticnet_long.csv` <br/> `output/supp/elasticnet_perf_summary.csv` <br/> Other rds <br/> `output/rds/biomarkers_elasticnet.rds` | Figure 7e-g <br/> SuppFig 8 |
