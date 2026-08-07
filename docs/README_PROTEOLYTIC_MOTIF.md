# Proteolytic Motif Analysis

This module extracts eight-residue windows around peptide N- and C-terminal cleavage positions and produces sequence-logo representations for each experiment.

## Required inputs

### Peptide table

The peptide table must contain:

- `Protein`
- `Start`
- `End`
- quantitative columns for the selected metric

`Peptide Sequence` is used when available. Rows matching common contaminant labels are excluded.

### FASTA file

A FASTA file is required for both analysis methods. The first token in each FASTA header is used as the protein identifier and must match the peptide table's `Protein` values.

### Protein summary table

This file is required only for the protein-background method. It must contain:

- `Protein`
- `Combined Total Peptides`

## Motif calculation choices

### Protein background

Choose either the top *n* proteins ranked by `Combined Total Peptides` or all proteins. PeptidomeR calculates amino-acid frequencies from the selected FASTA sequences and uses them as the background probabilities for binomial log-odds scoring.

### Frequency of counts

This method calculates the amino-acid count and frequency at each motif position directly from the extracted peptide cleavage windows. A protein summary table is not used.

## Experiment filtering

Quantitative columns are assigned to experiments in file order. A peptide is included for an experiment when at least one of its assigned replicate values is numeric, finite, and greater than zero. Experiments with no valid peptide rows are skipped.

## Cleavage windows

- N-terminal window: four residues before the peptide start through three residues after it
- C-terminal window: three residues before the peptide end through four residues after it

Only complete eight-residue windows without `X` are retained. Positions are displayed as P4, P3, P2, P1, P1', P2', P3', and P4'.

## Results

The app displays N- and C-terminal motif plots for each experiment. The ZIP download contains:

- experiment setup, column mapping, and valid-value counts;
- selected protein background and its amino-acid frequencies, when used;
- valid peptide subsets and peptide-to-FASTA merges;
- extracted terminal sequences and position indicator tables;
- position-specific counts/frequencies or binomial log-odds matrices; and
- individual and combined N- and C-terminal plots as 500 dpi TIFF files.

[Back to the main README](../README.md)

