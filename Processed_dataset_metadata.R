
# Define datasets

datasets <- list(
  "Adams" = list(
    ds_name = "AdamsT_2020_32832599",
    label_col = "Disease_Identity",
    low_res_ct_col = "OriginalAnnotationLevel1",
    hi_res_ct_col = "OriginalAnnotationLevel2"
  ),
  
  "Bassez" = list(
    ds_name = "BassezA_2021_33958794whole",
    label_col = "expansion",
    low_res_ct_col = "cellType",
    hi_res_ct_col = "cellSubType"
  ),
  
  "GongSharma_age_females_cmvneg" = list(
    ds_name = "GongSharma_age_females_cmvneg",
    label_col = "age",
    low_res_ct_col = "AIFI_L1",
    hi_res_ct_col = "AIFI_L3"
  ),
  
  "GongSharma_cmv_females" = list(
    ds_name = "GongSharma_cmv_females",
    label_col = "subject.cmv",
    low_res_ct_col = "AIFI_L1",
    hi_res_ct_col = "AIFI_L3"
  ),
  
  "Kfoury" = list(
    ds_name = "Kfoury_2021_34719426",
    label_col = "Status",
    low_res_ct_col = "cells_lowres",
    hi_res_ct_col = "cells"
  ),
  
  "Lee" = list(
    ds_name = "LeeA_2021_34836966p_tumor_seurat",
    label_col = "condition",
    low_res_ct_col = "layer1",
    hi_res_ct_col = "layer2"
  ),
  
  "Pelka" = list(
    ds_name = "PelkaK_2022_34450029whole",
    label_col = "Tissue",
    low_res_ct_col = "OriginalAnnotationLevel1",
    hi_res_ct_col = "OriginalAnnotationLevel3"
  ),
  
  "Smillie" = list(
    ds_name = "SmillieC_2019_31348891",
    label_col = "Condition",
    low_res_ct_col = "OriginalAnnotationLevel2",
    hi_res_ct_col = "OriginalAnnotationLevel3"
  ),
  
  "Stephenson" = list(
    ds_name = "MitchelJ_2023_PrePrintTBD_Stephenson_preprocessed",
    label_col = "Status",
    low_res_ct_col = "initial_clustering",
    hi_res_ct_col = "full_clustering"
  ),
  
  "Wu" = list(
    ds_name = "WuS_2021_34493872",
    label_col = "subtype",
    low_res_ct_col = "celltype_major",
    hi_res_ct_col = "celltype_subset"
  ),
  
  "Zhang" = list(
    ds_name = "ZhangY_2022_34653365",
    label_col = "Tissue",
    low_res_ct_col = "layer1",
    hi_res_ct_col = "layer2"
  )
)