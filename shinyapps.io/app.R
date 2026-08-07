suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(ggrepel)
  library(ggVennDiagram)
  library(limma)
  library(impute)
  library(imputeLCMD)
  library(Rtsne)
  library(uwot)
  library(pheatmap)
  library(ggseqlogo)
  library(patchwork)
  library(grid)
})

options(
  shiny.maxRequestSize = 200 * 1024^2,
  shiny.fullstacktrace = FALSE
)

source(file.path("R", "mod_qc.R"), local = TRUE)
source(file.path("R", "mod_eda.R"), local = TRUE)
source(file.path("R", "mod_differential_peptide.R"), local = TRUE)
source(file.path("R", "mod_differential_protein.R"), local = TRUE)
source(file.path("R", "mod_multivariate.R"), local = TRUE)
source(file.path("R", "mod_proteolytic_motif.R"), local = TRUE)

app_theme <- bs_theme(
  version = 5,
  bg = "#FFFFFF",
  fg = "#2C3E50",
  primary = "#8D6B9D",
  success = "#8D6B9D",
  secondary = "#4C3E50"
)

workflow_card <- function(number, title, description) {
  card(
    card_header(paste(number, title)),
    card_body(p(description))
  )
}

ui <- page_navbar(
  title = div(
    style = "display: flex; align-items: center;",
    tags$img(
      src = "peptidomerlogo.png?v=3",
      alt = "PeptidomeR logo",
      style = paste(
        "height: 42px;",
        "width: auto;",
        "display: inline-block;",
        "background-color: white;",
        "padding: 3px;",
        "border-radius: 5px;",
        "margin-right: 10px;"
      )
    ),
    span("PeptidomeR")
  ),
  id = "main_navigation",
  selected = "Home",
  theme = app_theme,
  window_title = "PeptidomeR",

  nav_panel(
    "Home",
    div(
      class = "container py-4",
      tags$img(
        src = "peptidomerlogo.png?v=3",
        alt = "PeptidomeR logo",
        width = "300px",
        style = "display: block; margin: 0 auto 20px auto;"
      ),
      h1("PeptidomeR", style = "text-align: center;"),
      p(
        class = "lead",
        paste(
          "An integrated workflow for preprocessing, exploratory analysis,",
          "differential analysis, multivariate analysis, and proteolytic motif",
          "analysis of peptidomics data."
        )
      ),
      layout_column_wrap(
        width = 1 / 3,
        workflow_card(
          "1",
          "Quality control",
          "Filter missing values, normalize, impute, and export intermediate datasets."
        ),
        workflow_card(
          "2",
          "Exploratory analysis",
          "Summarize experiment coverage, peptide overlap, and length distributions."
        ),
        workflow_card(
          "3",
          "Differential peptides",
          "Compare peptide abundance between two selected experimental groups."
        ),
        workflow_card(
          "4",
          "Differential proteins",
          "Compare peptide abundance between two selected experimental groups."
        ),
        workflow_card(
          "5",
          "Multivariate analysis",
          "Select PCA, t-SNE, or UMAP, and visualize the top PCA-loading peptides in a heatmap."
        ),
        workflow_card(
          "6",
          "Proteolytic motifs",
          "Generate N- and C-terminal motif analyses from peptides based on either binomial probability distribution with a backgroumd of top proteins or frequenct counts."
        )
      )
    )
  ),
  nav_panel("Quality Control", mod_qc_ui("qc")),
  nav_panel("Summary Plots", mod_eda_ui("eda")),
  nav_panel(
    "Differential Peptide",
    mod_differential_peptide_ui("differential_peptide")
  ),
  nav_panel(
    "Differential Protein",
    mod_differential_protein_ui("differential_protein")
  ),
  nav_panel("Multivariate Analysis", mod_multivariate_ui("multivariate")),
  nav_panel(
    "Proteolytic Motifs",
    mod_proteolytic_motif_ui("proteolytic_motif")
  )
)

server <- function(input, output, session) {
  mod_qc_server("qc")
  mod_eda_server("eda")
  mod_differential_peptide_server("differential_peptide")
  mod_differential_protein_server("differential_protein")
  mod_multivariate_server("multivariate")
  mod_proteolytic_motif_server("proteolytic_motif")
}

shinyApp(ui = ui, server = server)
