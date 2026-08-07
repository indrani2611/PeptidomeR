# Multivariate Analysis

This module visualises sample-level structure using one selected method: PCA, t-SNE, or UMAP.

## Inputs and preprocessing

- QC-filtered peptide table
- MaxLFQ Intensity, Intensity, or Spectral Count
- Experiment names and replicate counts
- Optional log2 transformation
- A valid-value rule requiring at least 50% valid values in at least one experiment or in every experiment

Quantitative columns are assigned to experiments in their uploaded order. Missing values are handled separately within each experiment:

1. partially missing rows are processed with KNN where possible;
2. rows missing across an entire experiment are processed with QRILC; and
3. any remaining missing cells use a down-shifted normal-distribution fallback.

The imputation seed is fixed at 123. Constant peptide features are excluded before multivariate analysis, and the analysis matrix is centred and scaled where required.

## Methods

### PCA

- Requires at least two samples.
- Displays PC1 and PC2 with their explained variance.
- Ranks peptides by their combined PC1 and PC2 loading magnitude.
- Displays a row-scaled heatmap for the requested number of top-loading peptides, up to the number available.

### t-SNE

- Requires at least five samples.
- Automatically chooses perplexity from the sample count.
- Uses a fixed seed of 123 and 1,000 iterations.

### UMAP

- Requires at least three samples.
- Automatically chooses up to 15 neighbours.
- Uses Euclidean distance, `min_dist = 0.1`, and a fixed seed of 123.

Only the selected method is calculated and plotted.

## Download

Every ZIP includes the uploaded QC-filtered table, selected quantitative values, valid-value-filtered rows, imputed data, imputation summary, sample mapping, and analysis parameters.

Method-specific outputs include:

- PCA sample scores, explained variance, all peptide loadings, the top-peptide table, PCA plot, and heatmap;
- t-SNE sample coordinates and plot; or
- UMAP sample coordinates and plot.

[Back to the main README](../README.md)

