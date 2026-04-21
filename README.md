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
| `scripts/04_RNA_heatmap.R` | runs ida prediction, but also calculates the actual DSS, HR, IDAComboscore |
| `scripts/Methylation_processing.R` | Standard parameters used for plots |

## Analysis

R objects, python objects

| Analysis/Code | Remarks | Files | 
| --- | --- | --- |
| `scripts/01_preprocess_counts_no.R`  | Perform standard RNA-seq preprocessing by converting raw sequencing counts into a normalized, log-transformed expression matrix | Supp Tables <br/> `output/Expression_LAML_707_log2RPM.txt` <br/> | Fig 1 <br/> SuppFig 1 <br/> Fig 2 <br/> SuppFig 2 |
| `scripts/02_CIBERSORT_pipeline_no.R`  | A cellular deconvolution pipeline to estimate the proportions of specific AML (Acute Myeloid Leukemia) cell states within bulk RNA-seq samples. | Supp Tables <br/> `output/baml_IDA_prediction.csv` <br/> `output/combine_baml_fp_IDA_prediction.csv` <br/> Other tables <br/> `output/supp/baml_woOSHU_IDA_prediction.csv` <br/> `output/supp/fp_full_IDA_prediction.csv` <br/> `output/supp/fp_commondrug_IDA_prediction.csv` <br/> Other rds <br/> (folder) `output/supp/synergy/` | Figure 3 <br/> SuppFig 3 <br/> SuppFig 4 |
| `src/analysis/3.catboost.Rmd`  | Refer to `src/analysis/catboost_py.ipynb` for codes modeling | Supp Tables <br/> `output/all_feature_summary_df.csv` <br/> `output/model_perf_summary_perfold.csv` <br/> `output/shap_summary_df.csv` <br/> `output/top_feature_summary_df.csv` <br/> `output/wilcox_pwc_combi.csv` <br/> `output/spear_numeric_combi.csv` <br/> Other tables <br/> (folder) `output/supp/ML/` <br/> <br/> Other rds <br/> `output/rds/cv_results_all170.pkl` | Figure 4 <br/> SuppFig 5 |
| `src/analysis/4.RNA_clustering.Rmd`  | Refer to `src/preprocess/4.rnaseq_preprocess.Rmd` for VST normalization and cellular hierarchy gene signature calculation using PLAGE-like method | Supp Tables <br/> `output/oshu_GSEA_hall_cluster.csv` <br/> `output/oshu_GSEA_hall_gobp_cluster.csv` <br/> `output/oshu_GSEA_gobp_cluster.csv` <br/> Other rds <br/> `output/rds/oshu_cluster_info.rds` <br/> `output/rds/oshu_seurat_obj.rds` | Figure 5a-d <br/> SuppFig 6a-d |
| `src/analysis/5.WGCNA.Rmd`  |  | Supp Tables  <br/> `output/oshu_wgcna_module_genes.csv` <br/> `output/oshu_wgcna_ora_hall_q25_go_q05.csv` <br/> `output/oshu_corr_ME_trait_results.csv` <br/> Other rds <br/> `output/rds/oshu_wgcna.rds` | Figure 5e <br/> SuppFig e-f |
| `src/analysis/6.scrna.Rmd`  | Refer to `src/preprocess/6.scrna_preprocessing.Rmd` for preprocessing <br/> Refer to `src/preprocess/6.scrna_mast.R` for DEG | Supp Tables <br/> `output/scrna_deg_all.csv` <br/> Other Tables <br/> (folder) `output/supp/scrna/` <br/> Other rds <br/> `output/rds/scrna_selected_obj.rds` <br/> `output/rds/scrna_selected_pseudotime_fitGAM.rds` <br/> `output/rds/scrna_selected_pseudotime_slingshot.rds` | Figure 6 <br/> SuppFig 7 |
| `src/analysis/7.clusterassoc.Rmd`  |  | Supp Tables <br/> `output/oshu_feature_cluster_fisher.csv` <br/> `output/oshu_feature_cluster_pwc.csv` <br/> `output/oshu_drug_cluster_kw.csv` <br/> `output/oshu_drug_cluster_wilcox_inflam.csv` <br/> `output/oshu_drug_cluster_posthoc_pwc_ref.csv` <br/> `output/oshu_genesig_drugresp.csv` | Figure 7a-c |
| `src/analysis/8.integrating_rna_prot.Rmd` |  | Supp Tables <br/> `output/oshu_corr_global_prot_rna_per_gene.csv` <br/> `output/oshu_rna_prot_overlap.csv` <br/> Other Tables <br/> `output/supp/corr_prot_combi_all.csv` <br/> `output/supp/corr_prot_rna_combi.csv` <br/> `output/supp/corr_rna_combi_all.csv` <br/> `output/supp/inflam_DEG_protein.csv` <br/> `output/supp/inflam_DEG_transcript.csv` | Figure 7d |
| `src/analysis/9.biomarker_elastic.Rmd` |  | Supp Tables <br/> `output/oshu_biomarkers.csv` <br/> `output/biomarkers_elasticnet_summary.csv` <br/> Other Tables <br/> `output/supp/biomarkers_elasticnet_long.csv` <br/> `output/supp/elasticnet_perf_summary.csv` <br/> Other rds <br/> `output/rds/biomarkers_elasticnet.rds` | Figure 7e-g <br/> SuppFig 8 |
