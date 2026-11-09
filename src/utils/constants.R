# ============================================================
# CONSTANTS AND LOOKUP TABLES
# ============================================================

# Method label mappings (centralized to avoid repeated recode calls)
benchmark_main_method_keys <- c(
  "ECODA_authors_HR", "ECODA_authors_HR_NULL", "Pseudobulk_hvg2000",
  "MOFA_hvg2000_factors15", "scITD_hvg2000_factors5", "MrVI_hvg2000",
  "scPoli_hvg2000_dims15_highres", "PILOT_hvg2000", "GloScope_hvg2000_pcadims30",
  "GloProp", "QOT_hvg2000", "PILOT-GM-VAE_hvg2000"
)
benchmark_extended_method_keys <- c(
  "ECODA_authors_HR", "ECODA_authors_HR_NULL", "ECODA_HiTME_HR_layer2", "ECODA_scATOMIC_HR",
  "ECODA_seuratres_2", "Freq_highres", "Pseudobulk_hvg2000", "Pseudobulk_hvg2000_bl",
  "Pseudobulk_CT_HR_hvg2000", "Pseudobulk_CT_HR_hvg500", "Pseudobulk_CT_LR_hvg500",
  "Pseudobulk_CT_LR_hvg2000", "MOFA_hvg2000_factors15", "scITD_hvg2000_factors5",
  "MrVI_hvg2000", "scPoli_hvg2000_dims15_highres", "PILOT_hvg2000",
  "GloScope_hvg2000_pcadims30", "GloProp", "QOT_hvg2000", "PILOT-GM-VAE_hvg2000"
)
method_label_map_main <- c(
  "ECODA_authors_HR" = "ECODA_authors",
  "ECODA_authors_HR_NULL" = "Shuffled_baseline",
  "ECODA_HiTME_HR_layer2" = "ECODA_automated",
  "ECODA_scATOMIC_HR" = "ECODA_scATOMIC",
  "ECODA_seuratres_2" = "ECODA_unsupervised",
  "ecoda_unsupervised_harmony" = "ECODA_Harmony_Leiden_CLR",
  "ecoda_unsupervised_harmony_limma" = "ECODA_Harmony_Leiden_CLR_limma",
  "Pseudobulk_hvg2000" = "Pseudobulk",
  "Pseudobulk_hvg2000_bl" = "Pseudobulk_BL",
  "MOFA_hvg2000_factors15" = "MOFA",
  "scITD_hvg2000_factors5" = "scITD",
  "MrVI_hvg2000" = "MrVI",
  "scPoli_hvg2000_dims15_highres" = "scPoli",
  "PILOT_hvg2000" = "PILOT",
  "PILOT-GM-VAE_hvg2000" = "PILOT-GM-VAE",
  "QOT_hvg2000" = "QOT",
  "GloScope_hvg2000_pcadims30" = "GloScope",
  "GloProp" = "GloProp",
  "Freq_highres" = "CellType_Frequency%"
)

# Method label mappings for annotation method comparison figures
method_label_map_annotation <- c(
  "ECODA_seuratres_0.1" = "ECODA_Leiden_res_0.1",
  "ECODA_seuratres_0.4" = "ECODA_Leiden_res_0.4",
  "ECODA_seuratres_2" = "ECODA_Leiden_res_2",
  "ECODA_seuratres_2_harmony" = "ECODA_Leiden_res_2_harmony",
  "ECODA_seuratres_5" = "ECODA_Leiden_res_5",
  "ECODA_seuratres_20" = "ECODA_Leiden_res_20",
  "ECODA_seuratres_50" = "ECODA_Leiden_res_50",
  "ECODA_HiTME_HR_layer2" = "ECODA_HiTME",
  "ECODA_scATOMIC_HR" = "ECODA_scATOMIC"
)

# Score label mappings
score_label_map <- c(
  "anosim_score" = "ANOSIM",
  "mod_knn3_score" = "Modularity",
  "cluster_score" = "ARI",
  "sil_score" = "Silhouette",
  "lisi_score" = "LISI"
)

# Score facet labels
score_facet_labels <- c(
  "anosim_score" = "ANOSIM score",
  "cluster_score" = "Adjusted rand index",
  "mod_knn3_score" = "Modularity score",
  "sil_score" = "Silhouette score",
  "lisi_score" = "LISI separation (higher is better)"
)

method_color_map <- c(
  "ecoda_unsupervised_harmony" = "#0072B2",
  "ecoda_unsupervised_harmony_limma" = "#D55E00"
)
