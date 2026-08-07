# Differential Protein Analysis

This module compares protein abundance between two selected experiments using the limma empirical-Bayes workflow.

## Before starting

The uploaded file should already have passed the intended quality-control and contaminant-filtering steps. Each selected comparison group must contain at least two replicates.

The analysis uses protein name, protein identifier/accession, and description metadata. In the uploaded version, the initial file validator also expects a column with the exact header `Peptide Sequence`; retain this column if it is present in the protein summary exported by the upstream workflow.

## Analysis steps

1. Select MaxLFQ Intensity, Intensity, or Spectral Count.
2. Enter experiment names and replicate counts. Columns are assigned in file order.
3. Select the control and test groups.
4. Indicate whether log2 transformation is needed.
5. Retain proteins with valid values in at least 50% of replicates in each selected group.
6. Impute remaining missing cells separately within each group using K-nearest neighbours with a fixed seed of 123.
7. Fit the Test-minus-Control contrast using limma and empirical Bayes.
8. Apply Benjamini-Hochberg correction or no p-value correction, as selected.

Positive `log2FC` means higher abundance in the test group. Negative `log2FC` means higher abundance in the control group. No absolute fold-change cutoff is used.

## Results and download

The app displays an analysis summary, volcano plot, and top-ranked protein table. The ZIP download contains the uploaded and filtered intermediate data, imputed values, imputation summary, experiment mapping, analysis parameters, full differential results, separate direction-specific result tables, and the volcano plot in PNG and PDF formats.

[Back to the main README](../README.md)

