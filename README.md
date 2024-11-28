# Spontan_AF_project
AF_vs_Control
This directory contains analyses comparing natural AF cases against control samples, with subdirectories organized by analysis type.

### `analysis/`
- **`01_dge/`**: Differential gene expression analysis.
  - **`output/`**: Results include PCA plots, overlap analysis, and volcano plots.
    - Subdirectories:
      - `Overlap_analysis/`: Overlap analysis results.
      - `PCA_all/`: PCA plots for all samples.
      - `PCA_healthy/`: PCA plots for healthy controls.
      - `Volcano_plots/`: Volcano plots of differential expression.
  - **`script/`**: Contains scripts to run the DGE pipeline.

- **`02_gene-enrichment/`**: Gene enrichment analysis.
  - **`output/`**: Results of gene enrichment analyses.
  - **`script/`**: Containst scripts to run the enrichment pipeline

- **`03_ion-channel/`**: Ion channel-specific analyses.
  - **`output/`**: Results of ion channel-specific expression analyses.
  - **`script/`**: Scripts for ion channel analysis.

- **`04_chord/`**: Chord diagram visualizations.
  - **`output/`**: Chord diagrams and related results.
  - **`script/`**: Scripts to generate chord diagrams.

- **`05_sankey/`**: Sankey diagram visualizations.
  - **`output/`**: Sankey diagrams and related results.
  - **`script/`**: Scripts to generate Sankey diagrams.

- **`06_Venn_Diagrams & MitoXplorer/`**: Venn diagram and mitochondrial pathway analyses.
  - **`output/`**: Results of Venn diagram and MitoXplorer analyses.
  - **`script/`**: Scripts to generate Venn-diagrams and Mitoxplorer pibeline

- **`07_proteomic-integration/`**: Integration with proteomics data from the same horses (previous study)
  - **`output/`**: Results of proteomic data integration.
  - **`script/`**: Scripts for integrating proteomics with other data.

- **`08_comparison_to_induced_model/`**: Comparison between natural and induced AF models (from previous study)
  - **`output/`**: Results of comparative analyses.
  - **`script/`**: Scripts to perform comparative analyses.

### `data/`
- **`count/`**: Count matrices for gene expression analyses.
- **`gene_annotation/`**: Gene annotation files, including GO terms.
- **`induced_AF/`**: Datasets specific to induced AF models.
- **`metadata/`**: Metadata associated with the samples.
- **`proteomics/`**: Proteomics datasets from same horses.

### `preprocessing/`
Contains files related to raw data analysis
