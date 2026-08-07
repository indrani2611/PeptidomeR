# Differential Peptide Analysis

This module compares peptide abundance between two selected experiments using the limma empirical-Bayes workflow.

## Before starting

The module assumes that contaminant removal and general quality-control filtering have already been completed. It does not repeat those steps. Each selected comparison group must contain at least two replicates.

Required metadata include peptide sequence, protein, protein identifier/accession, and protein description. See [Input Data Format](INPUT_DATA_FORMAT.md).

## Analysis steps

1. Select MaxLFQ Intensity, Intensity, or Spectral Count.
2. Enter the experiment names and replicate counts. Columns are assigned in file order.
3. Select the control and test groups.
4. Indicate whether the selected values require log2 transformation.
5. Retain peptides with valid values in at least 50% of replicates in each selected group. The required count is rounded up, so two of four and three of five replicates are required.
6. Impute remaining missing cells separately within each group using K-nearest neighbours. The random seed is fixed at 123, and up to 10 neighbours are used depending on the number of retained rows.
7. Fit a limma linear model and empirical-Bayes test for Test minus Control.
8. Apply Benjamini-Hochberg correction or leave p-values unadjusted, according to the selected setting.

Positive `log2FC` indicates higher abundance in the test group; negative `log2FC` indicates higher abundance in the control group. Significance is determined by the selected p-value threshold and the fold-change direction. No absolute fold-change cutoff is applied.

## On-screen results

- Analysis and imputation summary
- Volcano plot with the five most significant peptides labelled on each side, when available
- A table of the top-ranked differential peptides

## Download

The ZIP file contains:

- uploaded and pairwise-filtered intermediate tables;
- KNN-imputed data and an imputation summary;
- sample-to-experiment mapping and analysis parameters;
- the full differential table;
- separate tables for peptides higher in each group; and
- the volcano plot in PNG and PDF formats.

[Back to the main README](../README.md)

