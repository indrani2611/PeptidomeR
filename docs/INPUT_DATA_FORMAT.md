# Input Data Format

## Supported files

PeptidomeR accepts:

- comma-separated `.csv` files;
- tab-separated `.tsv` files; and
- comma- or tab-delimited `.txt` files.

The first row must contain column names. Keep quantitative cells numeric; blanks and common missing-value strings such as `NA`, `NaN`, `NULL`, and `null` are treated as missing in most modules.

## Quantitative columns

Choose one quantitative metric in the app:

- `MaxLFQ Intensity`
- `Intensity`
- `Spectral Count`

Column names may use spaces, dots, hyphens, or underscores around the metric words. Plain Intensity detection excludes MaxLFQ and Spectral Count columns. In the protein and multivariate modules, combined, unique, and total spectral-count summary columns are excluded from the experiment-level Spectral Count selection.

## Column order and study design

Quantitative columns are assigned to experiments from left to right. For example, if the user enters two control replicates followed by three treatment replicates, the first two detected metric columns are assigned to Control and the next three to Treatment.

```text
CTL_1 MaxLFQ Intensity
CTL_2 MaxLFQ Intensity
TRT_1 MaxLFQ Intensity
TRT_2 MaxLFQ Intensity
TRT_3 MaxLFQ Intensity
```

The total replicate count entered in the app must match the number of detected columns for the selected metric. Experiment names must be unique.

## Metadata requirements

| Module | Required or expected metadata |
| --- | --- |
| Quality Control and Summary Plots | Exact header `Peptide Sequence`; a protein-related column is used for contaminant removal when present |
| Differential Peptide | `Peptide Sequence`/`Sequence`, `Protein`/`Protein Name`, `Protein ID`/accession, and protein description |
| Differential Protein | `Protein`/`Protein Name`, `Protein ID`/accession, and description. The uploaded version's initial file validator also expects an exact `Peptide Sequence` header |
| Multivariate Analysis | Same peptide metadata as Differential Peptide |
| Proteolytic Motifs | `Protein`, `Start`, and `End`; `Peptide Sequence` is used when available |

Header aliases containing spaces, dots, or similar punctuation are accepted in the differential and multivariate analysis helpers. The initial upload validator in the main tabular modules specifically checks for `Peptide Sequence`.

## Example peptide-table header

```text
Peptide Sequence,Protein,Protein ID,Protein Description,Start,End,CTL_1 MaxLFQ Intensity,CTL_2 MaxLFQ Intensity,TRT_1 MaxLFQ Intensity,TRT_2 MaxLFQ Intensity
```

## Additional motif-analysis files

Proteolytic motif analysis always requires a FASTA file. The first whitespace-delimited token in each FASTA header is used as the protein identifier and must match the values in the peptide table's `Protein` column.

The protein-background method additionally requires a protein summary table containing:

- `Protein`
- `Combined Total Peptides`

The frequency-of-counts method does not use the protein summary table.

## Recommended checks before upload

1. Confirm that sample columns are ordered by experiment.
2. Confirm that the entered replicate counts equal the detected quantitative columns.
3. Remove units and thousands separators that prevent numeric conversion.
4. Ensure metadata identifiers are consistent between the peptide, protein-summary, and FASTA files.
5. Keep a copy of the original input file and download the output ZIP from each completed analysis.

