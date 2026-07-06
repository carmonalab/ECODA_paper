# ============================================================
# CONSTANTS AND LOOKUP TABLES
# ============================================================

# Method label mappings (centralized to avoid repeated recode calls)
method_label_map_main <- c(
  "ECODA_authors_HR" = "ECODA_authors",
  "ECODA_authors_HR_NULL" = "Shuffled_baseline",
  "ECODA_HiTME_HR_layer2" = "ECODA_automated",
  "ECODA_scATOMIC_HR" = "ECODA_scATOMIC",
  "ECODA_seuratres_2" = "ECODA_unsupervised",
  "Pseudobulk_hvg2000" = "Pseudobulk",
  "MOFA_hvg2000_factors15" = "MOFA",
  "scITD_hvg2000_factors5" = "scITD",
  "MrVI_hvg2000" = "MrVI",
  "scPoli_hvg2000_dims15_highres" = "scPoli",
  "PILOT_hvg2000" = "PILOT",
  "GloScope_hvg2000_pcadims30_sqrtmat" = "GloScope",
  "GloProp" = "GloProp",
  "Freq_highres" = "CellType_Frequency%"
)

# Method label mappings for annotation method comparison figures
method_label_map_annotation <- c(
  "ECODA_seuratres_0.1" = "ECODA_Leiden_res_0.1",
  "ECODA_seuratres_0.4" = "ECODA_Leiden_res_0.4",
  "ECODA_seuratres_2" = "ECODA_Leiden_res_2",
  "ECODA_seuratres_5" = "ECODA_Leiden_res_5",
  "ECODA_seuratres_20" = "ECODA_Leiden_res_20",
  "ECODA_HiTME_HR_layer2" = "ECODA_HiTME",
  "ECODA_scATOMIC_HR" = "ECODA_scATOMIC"
)

# Score label mappings
score_label_map <- c(
  "anosim_score" = "ANOSIM",
  "mod_knn3_score" = "Modularity",
  "cluster_score" = "ARI"
)

# Score facet labels
score_facet_labels <- c(
  "anosim_score" = "ANOSIM score",
  "cluster_score" = "Adjusted rand index",
  "mod_knn3_score" = "Modularity score"
)

# GongSharma dataset facet labels
gs_dataset_facet_labels <- c(
  "GongSharma_cmv_females"         = "Gong & Sharma\nSubset: females\nSeparation by: CMV infection",
  "GongSharma_age_females_cmvneg"  = "Gong & Sharma\nSubset: CMV negative females\nSeparation by: age"
)

# Dataset display label mappings
dataset_label_map <- c(
  "Adams" = "Adams (pulmonary fibrosis)",
  "Bassez" = "Bassez (ICB treatment response)",
  "GongSharma" = "GongSharma (healthy, CMV)",
  "GongSharma_all_age" = "GongSharma (healthy, age)",
  "GongSharma_all_cmv" = "GongSharma (healthy, CMV)",
  "GongSharma_all_sex" = "GongSharma (healthy, sex)",
  "Kfoury" = "Kfoury (prostate metastasis)",
  "Kim" = "Kim (metastatic lung adenocarcinoma)",
  "Lee" = "Lee (ICB treatment vs naive)",
  "Pelka" = "Pelka (colon cancer vs control)",
  "Smillie" = "Smillie (ulcerative colitis)",
  "Stephenson" = "Stephenson (COVID-19)",
  "Wu" = "Wu (breast cancer subtype)",
  "Zhang" = "Zhang (tissue)"
)