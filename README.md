# Analysis Scripts for "Defined human tri-lineage brain microtissues" Manuscript

This repository contains the single-cell and single-nuclei RNA-seq analysis scripts that produced the results shown in the manuscript "Defined human tri-lineage brain microtissues" by Uenaka T. et al. (currently under review). A link to the manuscript will be included in this repository once it is published.

The scripts process the scRNA-seq and snRNA-seq data of the tri-culture, integrate them with three public reference datasets (adult human brain and brain organoids), perform differential expression and GO enrichment analyses against these references, infer cell–cell interactions, and generate the corresponding panels of Figures 4, 5 and S4 as well as Supplementary Tables 4.1–4.6. The computational environments in which the results were produced are documented as renv lockfiles (see [Computational environments](#computational-environments)).

## Repository structure

```
.
├── Process_tri-culture_scRNAseq.R     # Processing and annotation of the tri-culture scRNA-seq data
├── Process_tri-culture_snRNAseq.R     # Processing and annotation of the tri-culture snRNA-seq data
├── IntegrationWithOtherStudies/      # Processing of public datasets and integration with the tri-culture
│   ├── Process_GSE179590.R
│   ├── Process_Linnarsson.R
│   └── Process_scOrganoids.R
├── Analyses/                         # Differential expression and GO enrichment vs. public datasets
│   ├── DiffExp_w_GSE179590.R
│   ├── DiffExp_w_Linnarsson.R
│   ├── DiffExp_w_scOrganoids.R
│   ├── GOEnrichment_w_GSE179590.R
│   ├── GOEnrichment_w_Linnarsson.R
│   └── GOEnrichment_w_scOrganoids.R
├── Figures/                          # One script per figure panel (Fig4A–J, Fig5A–E, FigS4G–H)
├── aux/                              # Color definitions and scripts creating figure input data
├── InputData/                        # Gene and ontology lists used by the scripts
├── lockfiles/                        # renv lockfiles and Python environment files
│   ├── env2023.lock
│   ├── env2025.lock
│   ├── env2026.lock
│   ├── env2023_densvis_python.txt
│   ├── env2026_python_linnarsson_explicit.txt
│   ├── env2026_python_linnarsson_pip.txt
│   ├── env2026_python_linnarsson.yml
├── LICENSE
└── README.md
```

## Workflow

All scripts use paths relative to the repository root and must be run with the working directory set to the repository root, in a UTF-8 locale (see [Running the analysis](#running-the-analysis)). Each step requires the output of the previous steps. To reproduce the results, run each step in the environment in which it was computed.

| Step | Scripts | Environment | Requires | Produces |
|---|---|---|---|---|
| 1. Processing of the tri-culture data | `Process_tri-culture_scRNAseq.R`, `Process_tri-culture_snRNAseq.R` | 2023 | 10x count matrices | Annotated Seurat objects |
| 2. Integration with public datasets | `IntegrationWithOtherStudies/Process_GSE179590.R`, `Process_scOrganoids.R` | 2023, then 2025 | Step 1, public datasets | Integrated Seurat objects |
| | `IntegrationWithOtherStudies/Process_Linnarsson.R` | 2026 (Python preprocessing, then R) | Step 1, public datasets | Integrated Seurat object |
| 3. Differential expression | `Analyses/DiffExp_w_*.R` | 2026 | Step 2 | Supplementary Tables 4.1, 4.3, 4.5 |
| 4. GO enrichment | `Analyses/GOEnrichment_w_*.R` | 2026 | Step 3 | Supplementary Tables 4.2, 4.4, 4.6 |
| 5. Figures | `Figures/*.R` | 2026 | Steps 1–4 | Figure panels in `PaperPlots/` |

Figure scripts check whether their input data exist and otherwise source the corresponding script in `aux/` to create them.

The objects produced by Steps 1 and 2 are available as processed data (see [Processed data](#processed-data)). Steps 3–5 can be run from them directly.

## 1. Processing of the tri-culture data

| Script | Input | Output | Manuscript |
|---|---|---|---|
| `Process_tri-culture_scRNAseq.R` | `scRNA-seq/CellRanger-Count-10X/` | `scRNA-seq/Y1_SeuratObj.rds` (processed), `scRNA-seq/Y1_final.rds` (annotated) | Figures S4A–C |
| `Process_tri-culture_snRNAseq.R` | `sNuc-seq/CellRanger-Count-10X/` | `sNuc-seq/sNuc_SeuratObj.rds` (processed, then overwritten with the annotated object) | Figures S4D–F |

Each script has two parts, separated by the first `saveRDS()`. The first part processes the data:

1. Quality control: Cells are removed if their fraction of mitochondrial reads or their log10 number of detected genes deviates by more than 3 MADs from the median, if their log10 UMI count exceeds the median by more than 3 MADs, or if they fall below the linear fit of log10 detected genes on log10 UMI counts (intercept offset −0.09).
2. Normalization: Cell-cycle scores are computed with the genes in `InputData/CellCycleGenes_Human.csv`, whose Ensembl IDs are mapped to gene symbols through Ensembl BioMart at run time. The data are normalized with SCTransform (v2), regressing out the difference between S and G2/M scores.
3. Doublet removal: DoubletFinder (pN = 0.25, pK maximizing the BC metric), with the expected doublet rate interpolated from the 10x Genomics estimates for 3' v3.1 chemistry. Cells with pANN ≥ 0.5 are removed. The number of principal components is selected automatically (`findNumPCs`).
4. Embedding: UMAP, t-SNE and densMAP.
5. Automated annotation with SCINA using brain markers from CellMarker 2.0 (metadata column `SCINA_annot`).

The second part reads the processed object back and annotates it: Leiden clustering on the UMAP (monocle3) and assignment of cell types to clusters through hard-coded maps keyed by cluster number (`ct_map`, `subct_map`; metadata columns `CTs` and, for snRNA-seq, `SubCTs`). In the scRNA-seq data, cells in a fixed region of the UMAP are additionally labelled as neurons. These assignments are valid for the original UMAP and clusters only (see [Reproducibility of the results](#reproducibility-of-the-results)).

Of 5,648 scRNA-seq and 11,781 snRNA-seq barcodes, 4,551 cells and 11,201 nuclei are retained. Annotated cell types are Astrocyte, Radial Glia, OPC, OLIG1+GFAP+ Cells, Neuron, Fibroblast-like and Microglia. The scRNA-seq data additionally contain Activated Microglia and Proliferating Microglia. QC and embedding plots are written to `scRNA-seq/` and `sNuc-seq/`. Label sizes in Figures S4A–F were adjusted in BioRender.

## 2. Integration with public datasets (`IntegrationWithOtherStudies/`)

| Script | Reference dataset | Integrated with | Output |
|---|---|---|---|
| `Process_GSE179590.R` | snRNA-seq of normal-appearing prefrontal cortex from patients with relapsing-remitting (RRMS) or secondary progressive (SPMS) multiple sclerosis, 8 samples (GEO: GSE179590) | tri-culture snRNA-seq | `GSE179590/Integrated.rds` (reference only), `GSE179590/Combined_SeuratObj.RData` |
| `Process_Linnarsson.R` | Human Brain Cell Atlas v1.0, snRNA-seq of the adult human brain (Siletti et al., 2023) | tri-culture snRNA-seq | `Linnarsson/Combined_SeuratObj.rds` |
| `Process_scOrganoids.R` | scRNA-seq of 6-month dorsal forebrain organoids from the 11a, GM08330 and PGP1 lines (Velasco et al., 2019) | tri-culture scRNA-seq | `scOrganoids/Combined_SeuratObj.RData` |

Each sample (GSE179590), donor (Linnarsson) or organoid (scOrganoids) and the tri-culture dataset are normalized separately with SCTransform (v2) and integrated using Seurat's reciprocal PCA anchors (3,000 integration features, 30 dimensions), followed by PCA, UMAP and t-SNE on the integrated assay. Downstream scripts use the metadata columns `Origin` (`This study` vs. `Brain` or `Organoid`) and `CTOrigin` (cell type and origin).

Dataset-specific steps:

- GSE179590. Each sample is processed with the pipeline of Step 1 (QC threshold 10 MADs, 30 PCs). The samples are integrated and clustered, and clusters are labelled (`FinalLabel`) using the cell-level metadata published with the dataset, through a map keyed by cluster number. The integration with the tri-culture uses 50 PCs.
- Linnarsson. The Python code included as comments at the top of the script reads all `.h5ad` files in `~/Downloads/`, downsamples them to at most 1,000 cells per cell type (random seed 42) and writes one `.h5ad` file per cell type to `IntegrationWithOtherStudies/Linnarsson/`. tri-culture gene symbols are converted to Ensembl IDs (org.Hs.eg.db; unmapped and duplicated IDs are removed), and atlas cell-type labels are mapped to broader classes (`ctMap`). The integration is performed per donor and uses 50 PCs.
- scOrganoids. Each organoid is processed separately (QC threshold 5 MADs, at most 5,000 cells per organoid, genes detected in more than 3 cells) and keeps its published cell-type label. The integration uses 30 PCs.

`Process_GSE179590.R` and `Process_scOrganoids.R` were run in two environments: everything before the line `# Integration with tri-culture` in the 2023 environment, and the section starting at that line in the 2025 environment (see [Running the analysis](#running-the-analysis)).

## 3. Differential expression and GO enrichment (`Analyses/`)

| Reference | Compared cell types (tri-culture vs. reference) | Differential expression | GO enrichment |
|---|---|---|---|
| GSE179590 | Microglia; Astrocyte; Neuron vs. Excitatory Neuron; Neuron vs. Inhibitory Neuron; OPC | `DiffExp_w_GSE179590.R` (Suppl. Table 4.1) | `GOEnrichment_w_GSE179590.R` (Suppl. Table 4.2) |
| Linnarsson | Microglia; Astrocyte; Neuron vs. Glutamatergic Neuron; Neuron vs. GABAergic Neuron; OPC | `DiffExp_w_Linnarsson.R` (Suppl. Table 4.3) | `GOEnrichment_w_Linnarsson.R` (Suppl. Table 4.4) |
| scOrganoids | Astrocyte vs. Astroglia; Neuron vs. CPNs; Radial Glia vs. RG | `DiffExp_w_scOrganoids.R` (Suppl. Table 4.5) | `GOEnrichment_w_scOrganoids.R` (Suppl. Table 4.6) |

Differential expression is computed with Libra as a pseudobulk test (edgeR, quasi-likelihood F-test) between tri-culture and reference cells of the same cell type. Tri-culture cells are randomly assigned to pseudo-replicates (`set.seed(1234)`), and the number of pseudo-replicates is defined per comparison in each script. In the organoid comparison, cells of the 11a and PGP1 organoids are split in the same way. The Linnarsson astrocyte and OPC comparisons are restricted to donors H18.30.002, H19.30.001 and H19.30.002.

Genes with |avg_logFC| > 0.5 and adjusted P ≤ 0.01 are reported. GO enrichment is computed with WebGestaltR (over-representation analysis, GO Biological Process, Benjamini–Hochberg FDR ≤ 0.01), separately for genes higher in the tri-culture and genes higher in the reference. WebGestaltR retrieves the gene sets from the WebGestalt server at run time.

Output files (in `IntegrationWithOtherStudies/<dataset>/`):

| File | Content |
|---|---|
| `DifferentialExpressionResults.RData` (GSE179590, Linnarsson), `DiffExp.RData` (scOrganoids) | Complete differential expression results |
| `DiffExp_<CellType>.txt` | Differentially expressed genes |
| `GOResults.RData` (GSE179590, Linnarsson), `GOEnrichment.RData` (scOrganoids) | Complete GO enrichment results |
| `GO_<CellType><tri-culture/Brain/Organoid>.txt` | Enriched GO terms among genes higher in the tri-culture or in the reference |

## 4. Figures (`Figures/`)

Each script generates one panel and writes it to `PaperPlots/<Panel>.png`.

| Script | Panel | Data |
|---|---|---|
| `Fig4A.R` | UMAP of annotated cell types | snRNA-seq |
| `Fig4B.R` | UMAP of annotated cell types | scRNA-seq |
| `Fig4C.R` | Clustered dot plot of marker genes | snRNA-seq |
| `Fig4D.R` | Clustered dot plot of marker genes | scRNA-seq |
| `Fig4E.R` | Heatmap of the top 10 cell type-specific genes per cell type (z-scored average expression) | scRNA-seq |
| `Fig4F.R` | UMAP of the integration with GSE179590 | snRNA-seq |
| `Fig4G.R` | UMAP of the integration with organoids | scRNA-seq |
| `Fig4H.R` | SynGO terms enriched among genes higher in tri-culture neurons than in organoid CPNs | scRNA-seq |
| `Fig4I.R` | Correlation of pseudobulk profiles between tri-culture and GSE179590 cell types | snRNA-seq |
| `Fig4J.R` | Correlation of pseudobulk profiles between tri-culture and organoid cell types | scRNA-seq |
| `Fig5A.R` | Number of inferred interactions between astrocytes, neurons and microglia | snRNA-seq |
| `Fig5B.R` | Number of inferred interactions between astrocytes, neurons and microglia (microglia states merged) | scRNA-seq |
| `Fig5C.R` | Receptors ranked by number of interactions, per target cell type | snRNA-seq |
| `Fig5D.R` | Receptors ranked by number of interactions, per target cell type | scRNA-seq |
| `Fig5E.R` | Network of interactions inferred in both datasets | scRNA-seq, snRNA-seq |
| `FigS4G.R` | UMAP of the integration with the Linnarsson data | snRNA-seq |
| `FigS4H.R` | Correlation of pseudobulk profiles between tri-culture and Linnarsson cell types | snRNA-seq |

The label "Inhibitory Neuron (Brain)" of the lower-left cluster in Figure 4F was added manually and is not reproduced by `Fig4F.R`. Label positions placed by ggrepel (Figures 4F, 4G and 5E) can vary between runs whereas the plotted data do not.

## 5. Auxiliary scripts (`aux/`)

| Script | Content | Output | Used by |
|---|---|---|---|
| `colorDefinitions.R` | Cell-type colors and color palettes | – | All figure scripts |
| `CreateInputData_Fig4E.R` | Selects the 10 genes per cell type with the largest difference in scaled average expression to all other cell types (average expression > 1) | `InputData/InputData_Fig4E.rds` | `Fig4E.R` |
| `CreateInputData_Fig4H.R` | Filters the GO enrichment results of tri-culture neurons to SynGO terms | `InputData/InputData_Fig4H.rds` | `Fig4H.R` |
| `CreateInputData_Fig4I.R` | Pearson correlation between pseudobulk profiles of tri-culture and GSE179590 cell types | `InputData/InputData_Fig4I.rds` | `Fig4I.R` |
| `CreateInputData_Fig4J.R` | Pearson correlation between pseudobulk profiles of tri-culture and organoid cell types | `InputData/InputData_Fig4J.rds` | `Fig4J.R` |
| `CreateInputData_FigS4H.R` | Pearson correlation between pseudobulk profiles of tri-culture and Linnarsson cell types. Donor H18.30.001 is excluded due to its low cell number | `InputData/InputData_FigS4H.rds` | `FigS4H.R` |
| `CreateInputData_Fig5ACE.R` | Infers cell–cell interactions in the snRNA-seq data with InterCom and keeps interactions with score ≥ 0.9 | `sNuc-seq/Intercom_forPaper_RNA/FinalInteractions.txt` | `Fig5A.R`, `Fig5C.R`, `Fig5E.R` |
| `CreateInputData_Fig5BDE.R` | Infers cell–cell interactions in the scRNA-seq data with InterCom and keeps interactions with score ≥ 0.9 | `scRNA-seq/Y1_forPaperFinal/FinalInteractions.txt` | `Fig5B.R`, `Fig5D.R`, `Fig5E.R` |

Pseudobulk profiles are computed by summing SCT counts per cell type over the genes shared by both datasets, normalizing to counts per million and z-scoring each gene across cell types.

## Input data (`InputData/`)

| File | Content |
|---|---|
| `CellCycleGenes_Human.csv` | Human S-phase (43) and G2/M-phase (54) genes as Ensembl IDs, identical to `cell_cycle/Homo_sapiens.csv` of [hbc/tinyatlas](https://github.com/hbc/tinyatlas) |
| `syngo_ontologies_BP_GO.txt` | Identifiers of the [SynGO](https://syngoportal.org) biological process terms (GO IDs and SynGO-specific IDs) |
| `InputData/Cell_marker_Human.txt` | List of human marker genes used for automated annotation of the data from CellMarker DB 2.0 |

The `InputData_*.rds` files generated by the scripts in `aux/` are written to the same folder.

## External data (not included)

The following data are required but not part of this repository. The paths are relative to the repository root.

| Location | Content | Source |
|---|---|---|
| `scRNA-seq/CellRanger-Count-10X/` | `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` of the tri-culture scRNA-seq data | GEO sample GSM9089763 (series GSE301787). Remove the prefix `GSM9089763_scRNA_` from the file names. |
| `sNuc-seq/CellRanger-Count-10X/` | `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` of the tri-culture snRNA-seq data | GEO sample GSM9089764 (series GSE301787). Remove the prefix `GSM9089764_snRNA_` from the file names. |
| `IntegrationWithOtherStudies/GSE179590/<sample>/` | 10x count matrix per sample, one folder per sample | GEO: GSE179590 |
| `IntegrationWithOtherStudies/GSE179590/GSE179590_MS_brains_8samples_cell_lvl_metadata.csv` | Cell-level metadata | GEO: GSE179590 |
| `~/Downloads/*.h5ad` | The 32 `.h5ad` files of the Human Brain Cell Atlas v1.0 used for the analysis: `14a02885-c168-4b3d-8775-b1750c828b9f.h5ad`, `2dd6eef4-2bbb-4633-b820-c03bfc0bbf94.h5ad`, `304a31c3-8e07-4325-99d0-f5c4a1948f95.h5ad`, `30cac3e5-8618-4132-9a77-bceeb8233c04.h5ad`, `3d605018-ad96-456a-848b-2841bf0dc338.h5ad`, `49371d50-9cb9-4699-868b-d8f1c9dbf127.h5ad`, `5059cc27-fc8a-448c-8029-16ce5f982567.h5ad`, `55fb42f4-a579-4d92-9eea-76d4bdb7c24d.h5ad`, `59248c6e-7b3a-4d77-9540-e9d3c9f05473.h5ad`, `5d75259a-9a43-45b2-853b-0b2134b99866.h5ad`, `79061388-1453-4f75-a9e0-2dbfe4dce4f6.h5ad`, `7b05d261-5e36-422a-90e6-2b3903861907.h5ad`, `867a02bc-e6e2-48f7-9b26-8c24c734efff.h5ad`, `8d65d563-c55e-4c81-b20d-352384357df5.h5ad`, `98f94301-6b7c-4dfa-a7a1-c6e27a62dfaf.h5ad`, `9a3a5e91-a599-407a-9cad-a4293c0bdf37.h5ad`, `9a4ba995-e83e-409f-afb4-3803ce189f9b.h5ad`, `9d67ef2e-09d9-4145-9f93-ceef377d13c4.h5ad`, `a0733204-201e-47d8-bb92-b2ab95f98938.h5ad`, `ad4355ed-4d8c-47ad-8c26-9d7b72eba69d.h5ad`, `bb9249e4-5515-48b1-b5fa-426876552252.h5ad`, `c1306d34-6880-4657-9432-62648cbed037.h5ad`, `c2118ef1-aafe-44a8-8902-44d4fbbe486d.h5ad`, `c81ddf79-09ce-4d75-a991-f9275484878f.h5ad`, `ce84a389-7b2a-43e5-b599-26e53ea99e6f.h5ad`, `d0fe7a36-36d5-4072-b040-0138d3e024f5.h5ad`, `d9c855b5-f7ed-4b04-8f7e-421aa6c08347.h5ad`, `e7157553-e797-4900-aad5-027ad2b3b047.h5ad`, `e9f3c636-da96-47b4-be7c-d13cf62cc6f7.h5ad`, `ea543fda-b0f0-4b92-91e5-a90293699bf2.h5ad`, `eb33715b-4b80-471f-b2d8-d5189e6d33a3.h5ad`, `f70dbbcf-8d03-40e5-acf3-4562e32e3ad2.h5ad` | CZ CELLxGENE |
| `IntegrationWithOtherStudies/scOrganoids/expression_11a.6mon.txt.gz`, `expression_GM.6mon.txt.gz`, `expression_PGP1.6mon.txt.gz`, `meta.txt` | Expression matrices and cell metadata of the 6-month organoids | Single Cell Portal, study "reproducible-brain-organoids" (GEO: GSE129519) |

Notes:

- `Read10X()` in Seurat 4.3.0 reads only files named `barcodes.tsv.gz`, `features.tsv.gz` and `matrix.mtx.gz` directly inside the given folder.
- `Process_GSE179590.R` processes every subfolder of `IntegrationWithOtherStudies/GSE179590/` and extracts the sample ID with the pattern `S[0-9]+`. Each folder name must therefore contain its sample ID (`S1`, `S2`, …), and each folder must contain only the three 10x files. Per-sample objects are written into these folders. A sample whose object already exists is skipped.
- The Linnarsson preprocessing reads every `.h5ad` file in `~/Downloads/`. The folder must contain exactly the 32 files listed above, otherwise the downsampled data differ.

Seurat objects, result tables and figure panels generated by the scripts are excluded from version control (see `.gitignore`).

## Processed data

The objects that Steps 3-5 start from are available on Zenodo (10.5281/zenodo.24542750). They are the reference for the cell-type annotations used in the manuscript.

| File | Content | Place at |
|---|---|---|
| `Y1_final.rds` | Annotated tri-culture scRNA-seq data | `scRNA-seq/Y1_final.rds` |
| `sNuc_SeuratObj.rds` | Annotated tri-culture snRNA-seq data | `sNuc-seq/sNuc_SeuratObj.rds` |
| `Combined_SeuratObj.RData` (GSE179590) | tri-culture snRNA-seq integrated with GSE179590 | `IntegrationWithOtherStudies/GSE179590/Combined_SeuratObj.RData` |
| `Combined_SeuratObj.RData` (scOrganoids) | tri-culture scRNA-seq integrated with the organoid data | `IntegrationWithOtherStudies/scOrganoids/Combined_SeuratObj.RData` |
| `Combined_SeuratObj.rds` (Linnarsson) | tri-culture snRNA-seq integrated with the Human Brain Cell Atlas | `IntegrationWithOtherStudies/Linnarsson/Combined_SeuratObj.rds` |

## Computational environments

The analysis was carried out in three R environments between 2023 and 2026. Each script, or script section, was run in exactly one of them, and each environment is documented by one renv lockfile in `lockfiles/`.

| Environment | Lockfile | R | Bioconductor | CRAN | Scripts |
|---|---|---|---|---|---|
| 2023 | `env2023.lock` (305 packages) | 4.2.1 | 3.16 | [Posit Package Manager snapshot 2023-03-12](https://packagemanager.posit.co/cran/2023-03-12) | `Process_tri-culture_scRNAseq.R`, `Process_tri-culture_snRNAseq.R`; `Process_GSE179590.R` and `Process_scOrganoids.R` up to the line `# Integration with tri-culture` |
| 2025 | `env2025.lock` (326 packages) | 4.2.1 | 3.16 | [Posit Package Manager snapshot 2024-03-12](https://packagemanager.posit.co/cran/2024-03-12) | `Process_GSE179590.R` and `Process_scOrganoids.R` from the line `# Integration with tri-culture` to the end |
| 2026 | `env2026.lock` (310 packages) | 4.5.3 | 3.22 | recorded in the lockfile | `Process_Linnarsson.R`, `Analyses/`, `aux/`, `Figures/` |

Key packages:

| Package | 2023 | 2025 | 2026 |
|---|---|---|---|
| Seurat / SeuratObject | 4.3.0 / 4.1.3 | 5.0.2 / 5.0.1 | 5.5.0 / 5.4.0 |
| sctransform | 0.3.5 | 0.4.1 | 0.4.3 |
| glmGamPoi | 1.10.2 | 1.10.2 | see lockfile |
| Matrix | 1.5-3 | 1.6-5 | 1.7-5 |
| matrixStats | 0.63.0 | 1.1.0 | see lockfile |
| GitHub packages (commit) | DoubletFinder (`67fb8b5`), monocle3 (`b545460`), SeuratWrappers (`d28512f`), singleseqgset (`3e3d9f0`) | DoubletFinder (`853e2de`), monocle3 (`e658515`), SeuratWrappers (`d9594f6`), singleseqgset (`3e3d9f0`) | InterCom (`f1ac21a`), Libra (`c4cb38c`), ggbump (`fe6d5c7`), ggheatmap (`5c373eb`) |

The lockfiles contain every package the scripts load and every namespace loaded during the runs.

Python environments:

| File | Environment | Used by |
|---|---|---|
| `env2023_densvis_python.txt` | Python environment of densvis 1.8.3 (umap-learn 0.5.0, numba 0.52.0), package lists for x86_64 and arm64 with the commands to create them | densMAP embedding in the 2023 processing |
| `env2026_python_linnarsson_explicit.txt`, `env2026_python_linnarsson_pip.txt` | Python 3.13.9 with scanpy and numpy (566 conda packages, 14 pip packages) | Linnarsson preprocessing |
| `env2026_python_linnarsson.yml` | Readable record of the same environment (`conda env export`); cannot be solved as is | – |

### Restoring the environments (macOS)

The lockfiles were restored and tested on macOS (Apple silicon) with R 4.2.1 (arm64 and x86_64) and R 4.5.3 (arm64), installed with [rig](https://github.com/r-lib/rig) (`rig add 4.2.1`, `rig add -a x86_64 4.2.1`). The x86_64 build runs under Rosetta 2 and is recommended for the 2023 and 2025 environments (see [Reproducibility of the results](#reproducibility-of-the-results)). Restore each lockfile into its own library, with renv 1.3.1 or later installed in a separate library, for example:

```r
renv::restore(lockfile = "lockfiles/env2023.lock", library = "<library for 2023>",
              prompt = FALSE, transactional = FALSE)
```

**2023 and 2025 (R 4.2.1).**

1. System requirements:
   - Xcode command line tools.
   - The gfortran build and CRAN recipe libraries (GDAL, PROJ, GEOS, OpenSSL, FreeType, etc.) listed for R 4.2 at [mac.r-project.org/tools](https://mac.r-project.org/tools/): gfortran 12.0.1 in `/opt/R/arm64` (arm64), gfortran 8.2 in `/opt/R/x86_64` (x86_64). On x86_64, zstd, pcre2 and lz4 are built from source.
   - XQuartz (monocle3 needs Cairo).
2. Build in a clean environment: no conda variables, no `~/.R/Makevars`, only the CRAN toolchain on the `PATH`, `LANG=en_US.UTF-8`; on arm64 also `MACOSX_DEPLOYMENT_TARGET=11.0`.
3. Install CRAN binaries of exactly the lockfile versions from CRAN's frozen R 4.2 binary repositories (`https://cran.r-project.org/bin/macosx/big-sur-arm64/contrib/4.2` or `https://cran.r-project.org/bin/macosx/contrib/4.2`). Build irlba from source, because the CRAN binaries require a newer Matrix package version. leidenbase (arm64, 2025) and lme4 (x86_64, 2025) are also built from source.
4. Build the packages that need specific settings first, one at a time:
   - arm64: Rhdf5lib (`CONFIG_SITE`), ragg (FreeType 2.13.2), igraph (GLPK 5.0), xml2 and fs (C++ standard headers).
   - x86_64: sf and terra (CRAN's static GDAL with all its dependencies), xml2 and fs.
5. Run `renv::restore(transactional = FALSE)` for the remaining packages (source builds from the Posit Package Manager snapshot, Bioconductor 3.16 and the GitHub commits), and repeat it until nothing fails.
6. Create the densvis Python environment as described in `env2023_densvis_python.txt` (`conda create --file`. On arm64 additionally `pip install --no-deps umap-learn==0.5.0`) and make it available to densvis through basilisk's external directory.

**2026 (R 4.5.3).**

1. System requirements for source builds: gfortran 14.2 in `/opt/gfortran` and gettext 0.22.5 in `/opt/R/arm64` (CRAN's R 4.5 toolchain).
2. Run `renv::restore()`. Packages are installed as binaries where available and otherwise built from source.
3. Create the Linnarsson Python environment:
   ```sh
   conda create -p <env> --file lockfiles/env2026_python_linnarsson_explicit.txt
   <env>/bin/python -m pip install --no-deps -r lockfiles/env2026_python_linnarsson_pip.txt
   ```

## Running the analysis

- Run all scripts from the repository root, with R started in a UTF-8 locale (e.g. `LANG=en_US.UTF-8`). In other locales, `list.files()` and `sort()` order files and cell types differently, which changes the merge order in `Process_Linnarsson.R` and the row and column order of Figures 4I, 4J and S4H.
- Run each script with the library restored from its lockfile, e.g. `.libPaths("<library for 2023>")` before `source("Process_Triculture_scRNAseq.R")`.
- Run the Linnarsson preprocessing first: execute the Python code in the comments at the top of `Process_Linnarsson.R` in the Linnarsson Python environment, then run the script in R.
- `Process_GSE179590.R` and `Process_scOrganoids.R` span two environments. The section from `# Integration with tri-culture` on also needs the `library()` calls, maps and functions at the top of the script (the header), which ends before `paths <- list.dirs(` in `Process_GSE179590.R` and before `meta <- read.table(` in `Process_scOrganoids.R`:

```r
f      <- "IntegrationWithOtherStudies/Process_GSE179590.R"
lines  <- readLines(f)
header <- seq_len(grep("^paths <- list.dirs", lines)[1] - 1)    # Process_scOrganoids.R: "^meta <- read.table"
split  <- grep("^# Integration with tri-culture", lines)

# R 4.2.1 with the 2023 library
source(exprs = parse(text = lines[seq_len(split - 1)]))

# New R 4.2.1 session with the 2025 library
source(exprs = parse(text = lines[c(header, split:length(lines))]))
```

- The cell-cycle mapping (Ensembl BioMart) and the GO enrichment (WebGestalt) require internet access.
- The integration steps require substantial memory. On a laptop, the GSE179590 processing takes about 1.5 hours, the other scripts minutes to an hour.

## Reproducibility of the results

The lockfiles were tested by restoring each of them into a fresh library (2023 and 2025 on arm64 and x86_64, 2026 on arm64), re-running all scripts and comparing every output with the original results. The original results were computed on an Intel Mac (2023 and 2025 environments) and on an Apple silicon Mac (2026 environment).

**Restored on the same processor architecture, the 2023 and 2025 lockfiles give bit-identical results:** independent installations produced identical outputs. The processor architecture and operating system can influence certain aspects of the calculations. For instance, running the scripts used within the 2023 and 2025 environments on an Intel Mac results in:

| Result | x86_64 (Rosetta 2) | arm64 |
|---|---|---|
| Retained cells, SCT counts, doublet calls, SCINA labels | identical | identical for the tri-culture scRNA-seq and the organoid data. Where 1-ulp differences in `exp()` change sctransform's seeded gene subsample, individual SCT counts, doublet calls and SCINA labels differ (tri-culture snRNA-seq: 1 doublet call, 90 SCINA labels; GSE179590: 7 doublet calls). |
| PCA, integrated assay | agree within 1e-5 | agree within 1e-5 where the SCT counts are identical |
| UMAP, t-SNE, densMAP, Leiden clusters on the UMAP | differ, except the UMAP and clusters of 7 of the 9 organoids | differ |

UMAP and the Leiden clustering at very low resolution amplify last-bit differences in floating-point arithmetic and in the system math library. Their results can therefore only be reproduced exactly on the machine that produced them, a known limitation of Seurat workflows across machines ([satijalab/seurat#5358](https://github.com/satijalab/seurat/issues/5358)).

**Rerunning the processing from the raw data may re-number the clusters, and the hard-coded maps then assign wrong labels or none.** Assigned by cluster content instead, 91–97% of the tri-culture cells and 99.5% of the GSE179590 cells recover their original label. We recommend the use of the processed data for the cell-type annotations of the manuscript. Given the original UMAP, contained in the processed objects, the annotation part of the processing scripts reproduces the original clusters and labels exactly.

**From the processed data, the 2026 environment reproduces the downstream results:**

| Result | Restored `env2026.lock` |
|---|---|
| GO enrichment results, figure input data (`InputData_*.rds`), InterCom interactions (`FinalInteractions.txt`), Figures 4A–E, 4H–J, 5A, 5B, 5D and S4H | identical |
| Differential expression results | identical within floating-point precision (log fold changes ≤ 1e-14, P values ≤ 8e-10) |
| Linnarsson integration: cells, counts, SCT data, metadata and cell types | identical |

Online services that a lockfile cannot fix:

- **Ensembl BioMart:** maps the Ensembl IDs of the cell-cycle genes to gene symbols. The original analysis used Ensembl release 109. However, releases 109 and 116 give identical mappings for these genes.
- **WebGestalt:** provides the GO gene sets. With WebGestaltR 0.4.6, the results were identical to the original ones in November 2026. Changes to the WebGestalt database may change them.

## Reference datasets

- GSE179590: single-nucleus RNA-seq of normal-appearing brains from RRMS and SPMS patients. Gene Expression Omnibus, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179590
- Siletti K. et al. Transcriptomic diversity of cell types across the adult human brain. Science (2023). https://doi.org/10.1126/science.add7046
- Velasco S. et al. Individual brain organoids reproducibly form cell diversity of the human cerebral cortex. Nature 570, 523–527 (2019). https://doi.org/10.1038/s41586-019-1289-x

## Citation

Uenaka T. et al. Defined human tri-lineage brain microtissues. Currently under review.

## License

This repository is released under the MIT License (see `LICENSE`).
