# PeptidomeR

<p align="center">
  <img src="www/peptidomerlogo.png" alt="PeptidomeR logo" width="320">
</p>

PeptidomeR is a browser-based R Shiny application for common peptidomics data-analysis tasks. It brings preprocessing, summary visualisation, differential analysis, multivariate analysis, and proteolytic motif analysis into one interface. The hosted application can be used from a desktop, tablet, or phone without installing R or writing code.

## Use PeptidomeR online

- [UCD deployment](https://peptide.ucd.ie/peptidomer/)
- [shinyapps.io deployment](https://peptigram.shinyapps.io/peptidomer/)

## Analysis modules

| Module | Main purpose | Documentation |
| --- | --- | --- |
| Quality Control | Remove contaminants, filter missing values, transform data, handle missing values, and export intermediate datasets | [Quality Control](docs/README_QC.md) |
| Summary Plots | Summarise peptide detection across experiments using overlap and peptide-length plots | [Summary Plots](docs/README_SUMMARY_PLOTS.md) |
| Differential Peptide | Compare peptide abundance between two experimental groups using limma | [Differential Peptide](docs/README_DIFFERENTIAL_PEPTIDE.md) |
| Differential Protein | Compare protein abundance between two experimental groups using limma | [Differential Protein](docs/README_DIFFERENTIAL_PROTEIN.md) |
| Multivariate Analysis | Run PCA, t-SNE, or UMAP; PCA also provides a loading-based peptide heatmap | [Multivariate Analysis](docs/README_MULTIVARIATE.md) |
| Proteolytic Motifs | Generate N- and C-terminal cleavage logos using amino-acid frequencies or a protein background | [Proteolytic Motifs](docs/README_PROTEOLYTIC_MOTIF.md) |

## Typical workflow

1. Prepare a comma- or tab-delimited peptide table.
2. Run **Quality Control** and download the final processed table.
3. Explore peptide overlap and length distributions in **Summary Plots**.
4. Use the processed table for differential and multivariate analyses.
5. If protein sequences and peptide coordinates are available, examine N- and C-terminal cleavage patterns in **Proteolytic Motifs**.
6. Download the ZIP file from each module to retain the tables, parameters, intermediate data, and figures used in the analysis.

## Input data

The main analysis modules accept `.csv`, `.tsv`, and comma- or tab-delimited `.txt` files. Quantitative columns may contain **MaxLFQ Intensity**, **Intensity**, or **Spectral Count** in their headers.

PeptidomeR assigns quantitative columns to experiments in their order in the uploaded file. The sum of the replicate counts entered in the app must therefore equal the number of columns detected for the selected metric.

For detailed column requirements and an example header, see [Input Data Format](docs/INPUT_DATA_FORMAT.md).

## Run locally with R

### Requirements

- R 4.5.3 is used by the supplied Docker image.
- The repository includes `renv.lock` for restoring the R package environment.

### Installation

```r
install.packages("renv")
renv::restore()
shiny::runApp()
```

Run these commands from the repository root. The app opens in the default web browser.

## Run with Docker

Build the image from the repository root:

```bash
docker build -t peptidomer:latest .
```

Start the container:

```bash
docker run --rm -p 3838:3838 peptidomer:latest
```

Open [http://localhost:3838](http://localhost:3838). If port 3838 is already in use, map a different host port:

```bash
docker run --rm -p 3839:3838 peptidomer:latest
```

Then open [http://localhost:3839](http://localhost:3839).

## Repository structure

```text
PeptidomeR/
├── app.R
├── R/
│   ├── mod_qc.R
│   ├── mod_eda.R
│   ├── mod_differential_peptide.R
│   ├── mod_differential_protein.R
│   ├── mod_multivariate.R
│   └── mod_proteolytic_motif.R
├── www/
│   └── peptidomerlogo.png
├── docs/
├── Dockerfile
└── renv.lock
```

## Reproducibility and interpretation

Each module records key settings and provides downloadable intermediate and final results. Users remain responsible for checking their study design, column assignments, transformations, and statistical assumptions before interpreting results.

## Citation
to be added

## Contact

indrani.bera@ucd.ie

