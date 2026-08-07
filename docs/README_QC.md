# Quality Control

The Quality Control module prepares peptide-level quantitative data for downstream analysis and records the number of rows retained at each stage.

## Inputs

- Peptide table in `.csv`, `.tsv`, or delimited `.txt` format
- Number of experimental conditions
- A unique name and replicate count for each condition
- Quantitative metric: Intensity, MaxLFQ Intensity, or Spectral Count

The selected quantitative columns are assigned to experiments in file order. The total replicate count must match the number of detected columns.

## Missing-value filters

The user can retain:

- features with at least 50% valid values in at least one group;
- features with at least 50% valid values in all groups;
- features with no more than 50% missing values across all samples; or
- complete rows with no missing values.

Zero quantitative values are treated as missing. If a protein-related column is present, rows whose protein value begins with `contam` or `CON__`, or contains `contaminant`, are removed.

## Transformation and missing-value handling

- **Log2 Transformation:** applies `log2` to the selected quantitative columns.
- **None:** leaves the selected values on their uploaded scale.
- **Fixed Constant Value:** replaces missing cells with a user-specified value.
- **Combined Strategy:** in the uploaded version, fully missing group rows are replaced with 95% of the minimum observed value and partially missing group rows are completed using that row's observed mean.

The two on-screen histograms show the selected metric before and after the optional transformation. Imputation is applied after the second plotted stage.

## Outputs

The downloadable ZIP contains:

- the uploaded table;
- the table after contaminant removal;
- the table after missingness filtering;
- the table after transformation;
- the final table after missing-value handling; and
- before- and after-transformation distribution plots as 300 dpi PNG files.

## Recommended use

Download `peptide_final_imputed_output.csv` for downstream work when imputation was selected. Always verify that the chosen missingness rule and transformation are appropriate for the data-generation and statistical workflow.

[Back to the main README](../README.md)

