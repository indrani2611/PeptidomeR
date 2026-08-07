# Summary Plots

The Summary Plots module describes peptide detection and peptide-length patterns across experiments.

## Inputs

- Peptide table in `.csv`, `.tsv`, or delimited `.txt` format
- Quantitative metric: MaxLFQ Intensity, Intensity, or Spectral Count
- Number of experiments, experiment names, and replicate counts

A peptide is considered present in an experiment when at least one assigned replicate has a value greater than zero. Quantitative columns are assigned to experiments in their file order.

## Peptide identity

When the required columns are available, a unique peptide key is constructed from peptide sequence, protein identifier, start position, and end position. Peptide length is calculated from the peptide sequence. If those metadata are unavailable, the first column is used as the fallback identifier.

## Results

The module reports the number and percentage of peptides found:

- only in each individual experiment;
- in each possible combination of experiments; and
- in all experiments.

It displays:

- a Venn diagram of peptide overlap; and
- peptide-length histograms faceted by detection compartment.

## Download

The ZIP file contains a CSV subset for every detection compartment together with 300 dpi PNG files for the Venn diagram and peptide-length distribution.

Venn diagrams become difficult to interpret when many experiments are included. In that situation, use the exported compartment tables as the primary numerical result.

[Back to the main README](../README.md)

