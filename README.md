[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2026.03.27.714811-b31b1b.svg)](https://doi.org/10.64898/2026.03.27.714811)

# ECODA: Exploratory Compositional Data Analysis for scRNA-seq Cohorts

This repository contains the code to reproduce the analyses, benchmark, and figures from the paper:  
**"Cell type composition drives patient stratification in single-cell RNA-seq cohorts"**.

---

## Overview

Single-cell RNA sequencing (scRNA-seq) enables high-resolution characterization of cellular heterogeneity, but summarizing single-cell data for cohort-level patient stratification remains a challenge. We benchmarked nine state-of-the-art sample representation methods—**MOFA+, scITD, GloScope, GloProp, MrVI, PILOT, PILOT-GM-VAE, QOT, and scPoli**—alongside baseline approaches (**Pseudobulk** and cell-type composition via **ECODA**) across standardized clinical cohorts (see [`datasets.json`](datasets.json)).

### Key Findings
- **Performance:** Centered log-ratio (CLR)-transformed cell-type proportions (**ECODA**) consistently match or outperform complex gene-expression-based embeddings in recovering known biological groupings in an unsupervised setting.
- **Efficiency:** ECODA produces patient embeddings in seconds with minimal computational overhead compared to deep-learning or optimal-transport methods.
- **Robustness:** Highly robust to technical batch effects and across various cell-type annotation granularities.
- **Interpretability:** Stratification is often driven by a small subset of **Highly Variable Cell Types (HVCs)**, providing direct mechanistic insights.

The standalone **scECODA** R package for cohort-level compositional analysis is available at [github.com/carmonalab/scECODA](https://github.com/carmonalab/scECODA).

---

## Installation & Environment Setup

This project uses [Pixi](https://pixi.sh) for reproducible multi-language (R & Python) environment management.

### 1. Install Pixi
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### 2. Set Up the Environment

- **Local macOS / Workstation (Lightweight analysis & figure generation):**
  ```bash
  pixi install && pixi run setup
  ```
  *Note:* Only the interactive analysis notebooks (`notebooks/benchmark_analysis.rmd` and `notebooks/batch_effect_analysis.rmd`) are designed to run locally on precomputed results.

- **HPC Cluster (SLURM worker build for full pipeline execution):**
  ```bash
  sbatch src/utils/bash/setup_env_sbatch.sh
  ```

---

## Workflow Overview

The reproducible pipeline is structured into four main stages. Detailed execution instructions, SLURM array commands, and folder architectures are documented in [**`docs/ARCHITECTURE.md`**](docs/ARCHITECTURE.md).

```
Stage 1: QC Filtering           -> notebooks/QC_filtering/
Stage 2: Preprocessing & Annot  -> src/1_stage_data/ -> src/2_dataset_specific_preprocessing/
                                   -> src/3_scrnaseq_preprocessing/ -> src/4_cell_type_annotation/
Stage 3: Benchmark Methods      -> src/5_run_benchmark_methods/ (Python & R arrays)
                                   -> notebooks/benchmark_analysis.rmd (Figures & Evaluation)
Stage 4: Batch Effect Analysis  -> notebooks/batch_effect_analysis.rmd
```

---

## Repository Structure

```
├── datasets.json               # Authoritative dataset metadata, columns, and views
├── docs/                       # Architecture documentation and hotfix notes
│   └── ARCHITECTURE.md         # Full pipeline architecture, HPC layout, and call flow
├── notebooks/                  # R Markdown analysis notebooks & dataset onboarding
│   ├── benchmark_analysis.rmd  # Main benchmark figures, separation metrics, and tables
│   ├── batch_effect_analysis.rmd # Batch effect quantification & mitigation
│   ├── QC_filtering/           # Per-dataset raw QC filtering notebooks
│   └── dataset_onboarding/     # Onboarding scripts for new cohorts
├── src/                        # Pipeline source code
│   ├── 1_stage_data/           # Raw data staging from storage to compute scratch
│   ├── 2_dataset_specific_preprocessing/ # Dataset-specific cohort adaptations
│   ├── 3_scrnaseq_preprocessing/         # Standardized Scanpy preprocessing pipeline
│   ├── 4_cell_type_annotation/           # Parallel scATOMIC + HiTME annotation
│   ├── 5_run_benchmark_methods/          # Python & R benchmark method arrays
│   ├── slurm_config.sh         # Centralized HPC environment and path configuration
│   └── utils/                  # Core math, scoring, plotting, and I/O utilities
└── aux/                        # Reference maps, gene blocklists, and gene annotation tables
```

---

## Citation

If you use ECODA or this benchmark suite in your research, please cite:

> **Cell type composition drives patient stratification in single-cell RNA-seq cohorts.**  
> Halter, C., Andreatta, M., & Carmona, S. J. (2026). *bioRxiv*.  
> doi: [10.64898/2026.03.27.714811v1](https://doi.org/10.64898/2026.03.27.714811v1)

---

## Reference Methods

Sample representation and embedding benchmark methods:

- **MOFA+**  
  Argelaguet, R., Arnol, D., Bredikhin, D., Deloro, Y., Velten, B., Bonn, M. Y., & Stegle, O. (2020). MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biology*, 21(1), 111. doi: [10.1186/s13059-020-02015-1](https://doi.org/10.1186/s13059-020-02015-1)

- **MrVI**  
  Boyeau, P., Hong, J., Gayoso, A., Jordan, M. I., Azizi, E., Ergen, C., & Yosef, N. (2024). Deep generative modeling of sample-level heterogeneity in single-cell genomics. *bioRxiv*, 2022.10.04.510898. doi: [10.1101/2022.10.04.510898](https://doi.org/10.1101/2022.10.04.510898)

- **scPoli**  
  De Donno, C., Hediyeh-Zadeh, S., Moinfar, A. A., Xu, P., & Theis, F. J. (2023). Population-level integration of single-cell datasets enables multi-scale analysis across samples. *Nature Methods*, 20(11), 1683–1692. doi: [10.1038/s41592-023-02035-2](https://doi.org/10.1038/s41592-023-02035-2)

- **PILOT**  
  Joodaki, M., Shaigan, M., Parra, V., Bülow, R. D., Kuppe, C., Hölscher, D. L., Cheng, M., Nagai, J. S., Goedertier, M., Bouteldja, N., Tesar, V., Barratt, J., Roberts, I. S., Coppo, R., Kramann, R., Boor, P., & Costa, I. G. (2024). Detection of PatIent-Level distances from single cell genomics and pathomics data with Optimal Transport (PILOT). *Molecular Systems Biology*, 20(2), 57–74. doi: [10.1038/s44320-023-00003-8](https://doi.org/10.1038/s44320-023-00003-8)

- **PILOT-GM-VAE**  
  Joodaki, M., Shaigan, M., Samiei, S., Nagai, J., Maié, T., Kuppe, C., & Costa, I. G. (2025). PILOT-GM-VAE: patient-level analysis of single-cell disease atlas with optimal transport of Gaussian mixture variational autoencoders. *Briefings in Bioinformatics*, 26(5), bbaf547. doi: [10.1093/bib/bbaf547](https://doi.org/10.1093/bib/bbaf547)

- **Pseudobulk (DESeq2)**  
  Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. doi: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

- **scITD**  
  Mitchel, J., Gordon, M. G., Perez, R. K., Biederstedt, E., Bueno, R., Ye, C. J., & Kharchenko, P. V. (2024). Coordinated, multicellular patterns of transcriptional variation that stratify patient cohorts are revealed by tensor decomposition. *Nature Biotechnology*, 43(7), 1192–1201. doi: [10.1038/s41587-024-02411-z](https://doi.org/10.1038/s41587-024-02411-z)

- **GloScope & GloProp**  
  Wang, H., Torous, W., Gong, B., & Purdom, E. (2024). Visualizing scRNA-Seq data at population scale with GloScope. *Genome Biology*, 25(1), 259. doi: [10.1186/s13059-024-03398-1](https://doi.org/10.1186/s13059-024-03398-1)

- **QOT**  
  Wang, Z., Zhan, Q., Yang, S., Mu, S., Chen, J., Garai, S., Orzechowski, P., Wagenaar, J., & Shen, L. (2025). QOT: Quantized Optimal Transport for sample-level distance matrix in single-cell omics. *Briefings in Bioinformatics*, 26(1), bbae713. doi: [10.1093/bib/bbae713](https://doi.org/10.1093/bib/bbae713)

---

## Reference Datasets

Single-cell RNA-seq cohorts used across the benchmark and batch effect analyses:

- **Adams (pulmonary fibrosis)**  
  Adams, T. S., Schupp, J. C., Poli, S., Ayaub, E. A., Neumark, N., Ahangari, F., Chu, S. G., Raby, B. A., DeIuliis, G., Januszyk, M., Duan, Q., Arnett, H. A., Siddiqui, A., Washko, G. R., Homer, R., Yan, X., Rosas, I. O., & Kaminski, N. (2020). Single-cell RNA-seq reveals ectopic and aberrant lung-resident cell populations in idiopathic pulmonary fibrosis. *Science Advances*, 6(28), eaba1983. doi: [10.1126/sciadv.aba1983](https://doi.org/10.1126/sciadv.aba1983)

- **Bassez (breast cancer ICB response)**  
  Bassez, A., Vos, H., Van Dyck, L., Floris, G., Arijs, I., Desmedt, C., Boeckx, B., Vanden Bempt, M., Neven, P., Smeets, A., Wildiers, H., Vandekerkhove, P., Nevelsteen, I., Punie, K., Pochet, R., Keupers, M., Schepers, M., De Schepper, M., Salihi, R., ... Lambrechts, D. (2021). A single-cell map of intratumoral changes during anti-PD1 treatment of patients with breast cancer. *Nature Medicine*, 27(5), 820–832. doi: [10.1038/s41591-021-01323-8](https://doi.org/10.1038/s41591-021-01323-8)

- **Alzheimer (SEA-AD)**  
  Gabitto, M. I., Travaglini, K. J., Rachleff, V. M., Kaplan, E. S., Long, B., Ariza, J., Ding, Y., Mahoney, J. T., Dee, N., Goldy, J., ... Lein, E. S. (2024). Integrated multimodal cell atlas of Alzheimer’s disease. *Nature Neuroscience* / *bioRxiv*, 2023.05.08.539485. doi: [10.1101/2023.05.08.539485](https://doi.org/10.1101/2023.05.08.539485)

- **GongSharma (healthy young males, CMV)**  
  Gong, Q., Sharma, M., Kuan, E. L., ... Davis, M. M. (2024). Longitudinal multi-omic immune profiling reveals age-related immune cell dynamics in healthy adults. *bioRxiv*, 2024.09.10.612119. doi: [10.1101/2024.09.10.612119](https://doi.org/10.1101/2024.09.10.612119)

- **Diabetes (mouse pancreas)**  
  Hrovatin, K., Bastidas-Ponce, A., Bakhti, M., Zunder, E., Lickert, H., & Theis, F. J. (2023). Delineating mouse β-cell identity during lifetime and in diabetes with a single cell atlas. *Nature Metabolism*, 5(9), 1615–1637. doi: [10.1038/s42255-023-00876-x](https://doi.org/10.1038/s42255-023-00876-x)

- **Joanito (colorectal cancer)**  
  Joanito, I., Wirapati, P., Zhao, N., Nawaz, Z., Yeo, G., Teng, F., DasGupta, R., Skanderup, A. J., Tan, P. B. O., Rozen, S. G., & Tejpar, S. (2022). Single-cell and bulk transcriptome sequencing identifies two epithelial tumor cell states and refines the consensus molecular classification of colorectal cancer. *Nature Genetics*, 54(7), 963–975. doi: [10.1038/s41588-022-01100-4](https://doi.org/10.1038/s41588-022-01100-4)

- **Parkinson (dopamine neurons)**  
  Kamath, T., Abdulraouf, A., Burris, S. J., Langlieb, J., Gazestani, V., Nadaf, N. M., Balderrama, K., Vanderburg, C., & Macosko, E. Z. (2022). Single-cell genomic profiling of human dopamine neurons identifies a population that selectively degenerates in Parkinson's disease. *Nature Neuroscience*, 25(5), 588–595. doi: [10.1038/s41593-022-01061-1](https://doi.org/10.1038/s41593-022-01061-1)

- **Kfoury (prostate metastasis)**  
  Kfoury, Y., Baryawno, N., Severe, N., Mei, S., Gustafsson, K., Hirz, T., Brouse, T., Scadden, E. W., Iglesias-Bartolome, R., Wu, C., Haber, D. A., Sykes, D. B., & Scadden, D. T. (2021). Human prostate cancer bone metastases have an actionable immunosuppressive microenvironment. *Cancer Cell*, 39(11), 1464–1478.e8. doi: [10.1016/j.ccell.2021.09.005](https://doi.org/10.1016/j.ccell.2021.09.005)

- **Kim (metastatic lung adenocarcinoma)**  
  Kim, N., Kim, H. K., Lee, K., Hong, Y., Cho, J. H., Choi, Y. L., Lee, J., Suh, Y. L., Ku, B. M., Eum, H. H., Choi, S., Choi, J. W., Ahn, M. J., Park, K., & Lee, H. O. (2020). Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma. *Nature Communications*, 11(1), 2285. doi: [10.1038/s41467-020-16164-1](https://doi.org/10.1038/s41467-020-16164-1)

- **Breast cancer (Kumar atlas)**  
  Kumar, T., Nee, K., Wei, R., He, S., Nguyen, Q. H., Bai, S., Blake, K., Pe’er, D., Rozenblatt-Rosen, O., Regev, A., Werb, Z., Kessenbrock, K., & Navin, N. E. (2023). A spatially resolved single-cell genomic atlas of the adult human breast. *Nature*, 620(7972), 181–191. doi: [10.1038/s41586-023-06252-9](https://doi.org/10.1038/s41586-023-06252-9)

- **Myocardial infarction (Kuppe)**  
  Kuppe, C., Ramirez Flores, R. O., Li, Z., Hayat, S., Levinson, R. T., Liao, M. C., ... Kramann, R. (2022). Spatial multi-omic map of human myocardial infarction. *Nature*, 608(7924), 766–777. doi: [10.1038/s41586-022-05060-x](https://doi.org/10.1038/s41586-022-05060-x)

- **Kidney (KPMP / Lake)**  
  Lake, B. B., Menon, R., Winfree, S., Hu, Q., Ferreira, R. M., Kalhor, K., ... Jain, S. (2023). An atlas of healthy and injured cell states and niches in the human kidney. *Nature*, 619(7970), 585–594. doi: [10.1038/s41586-023-05769-3](https://doi.org/10.1038/s41586-023-05769-3)

- **Lee (glioblastoma ICB)**  
  Lee, A. H., Sun, L., Mochizuki, A. Y., Reynoso, J. G., Orpilla, J., Chow, F., Kienzler, J. C., Everson, R. G., Nathanson, D. A., Bensinger, S. J., Liau, L. M., & Prins, R. M. (2021). Neoadjuvant PD-1 blockade induces T cell and cDC1 activation but fails to overcome the immunosuppressive tumor associated macrophages in recurrent glioblastoma. *Nature Communications*, 12(1), 6938. doi: [10.1038/s41467-021-26940-2](https://doi.org/10.1038/s41467-021-26940-2)

- **Pelka (colon cancer)**  
  Pelka, K., Hofree, M., Chen, J. H., Sarkizova, S., Pirl, J. D., Jorgji, V., Bejnood, A., Dionne, D., Ge, W. H., Xu, K. H., Chao, K. X., Trotman, G., Pelham, S. J., Simmons, D. P., Zalmas, L. P., Aguiar, B. R., Smith, S. M., Wang, Y., ... Hacohen, N. (2021). Spatially organized multicellular immune hubs in human colorectal cancer. *Cell*, 184(18), 4734–4752.e20. doi: [10.1016/j.cell.2021.08.003](https://doi.org/10.1016/j.cell.2021.08.003)

- **Lupus PBMC (Perez)**  
  Perez, R. K., Gordon, M. G., Subramaniam, M., Kim, M. C., Hartoularos, G. C., Targ, S., Sun, Y., Chesnut, J. D., Fiorentino, D., Bao, D. X., Lu, A., ... Ye, C. J. (2022). Single-cell RNA-seq reveals cell type-specific molecular and genetic associations to lupus. *Science*, 376(6589), eabf1970. doi: [10.1126/science.abf1970](https://doi.org/10.1126/science.abf1970)

- **Covid-19 PBMC (Ren)**  
  Ren, X., Wen, W., Fan, X., Hou, W., Su, B., Cai, P., ... Zhang, Z. (2021). COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas. *Cell*, 184(7), 1895–1913.e19. doi: [10.1016/j.cell.2021.01.053](https://doi.org/10.1016/j.cell.2021.01.053)

- **Lung atlas (HLCA / Sikkema)**  
  Sikkema, L., Ramírez-Suástegui, C., Strobl, D. C., Gillett, T. E., Sauler, M., ... Theis, F. J. (2023). An integrated cell atlas of the lung in health and disease. *Nature Medicine*, 29(6), 1563–1577. doi: [10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2)

- **Smillie (ulcerative colitis)**  
  Smillie, C. S., Biton, M., Ordovas-Montanes, J., Sullivan, K. M., Burgin, G., Graham, D. B., Herbst, R. H., Rogel, N., Slyper, M., Waldman, J., Sud, M., Andrews, E., Velonias, G., Haber, A. L., Jagadeesh, K., Vickovic, S., Yao, J., ... Regev, A. (2019). Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis. *Cell*, 178(3), 714–730.e22. doi: [10.1016/j.cell.2019.06.029](https://doi.org/10.1016/j.cell.2019.06.029)

- **Stephenson (COVID-19 PBMC)**  
  Stephenson, E., Reynolds, G., Botting, R. A., Calero-Nieto, F. J., Morgan, M. D., Tuong, Z. K., Bach, K., Sungnak, W., Worlock, K. B., Yoshida, M., Kleshchevnikov, V., Chaddock, N., Fortune, M. D., Innes, A., Raghavan, S., Gershlick, D. C., Goh, I., ... Haniffa, M. (2021). Single-cell multi-omics analysis of the immune response in COVID-19. *Nature Medicine*, 27(5), 904–916. doi: [10.1038/s41591-021-01329-2](https://doi.org/10.1038/s41591-021-01329-2)

- **Wu (breast cancer subtypes)**  
  Wu, S. Z., Al-Eryani, G., Roden, D. L., Junankar, S., Harvey, K., Andersson, A., Baker, L. A., Claridge, L. D., Lu, J., Liu, D. L., Galli, M., ... Swarbrick, A. (2021). A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics*, 53(9), 1334–1347. doi: [10.1038/s41588-021-00911-1](https://doi.org/10.1038/s41588-021-00911-1)

- **Zhang (triple-negative breast cancer)**  
  Zhang, Y., Chen, H., Mo, H., Hu, X., Gao, R., Zhao, Y., Liu, B., Niu, L., Sun, X., Yu, X., Wang, Y., ... Zhang, Z. (2021). Single-cell analyses reveal key immune cell subsets associated with response to PD-L1 blockade in triple-negative breast cancer. *Cancer Cell*, 39(12), 1578–1593.e8. doi: [10.1016/j.ccell.2021.09.010](https://doi.org/10.1016/j.ccell.2021.09.010)

- **Zhu (PBMC aging clocks)**  
  Zhu, H., Chen, J., Liu, K., ... Zhang, Z. (2023). Human PBMC scRNA-seq–based aging clocks reveal ribosome to inflammation balance as a single-cell aging hallmark and super longevity. *Science Advances*, 9(26), eabq7599. doi: [10.1126/sciadv.abq7599](https://doi.org/10.1126/sciadv.abq7599)

---

## Reference Data

- **Ensembl 105 (GRCh38.p13):** Human gene reference for name standardization (`aux/EnsemblGenes105_Hsa_GRCh38.p13.txt.gz`).
- **Carmona Lab Reference Atlases:** ScRNA-seq reference atlases for cell-type mapping (Garnica et al., 2024; Figshare DOI: [10.6084/m9.figshare.26310994](https://doi.org/10.6084/m9.figshare.26310994); repo: [`carmonalab/Reference_maps`](https://github.com/carmonalab/Reference_maps)).