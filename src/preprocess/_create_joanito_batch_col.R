library(Seurat)

input <- "data/JoaI_2022_35773407_Nofilt_whole.rds"
seurat <- readRDS(input)

seurat$seqtec <- ifelse(
  seurat$dataset %in% c("CRC-SG1", "KUL5"),
  "5' seq",
  "3' seq"
)

saveRDS(seurat, input)
