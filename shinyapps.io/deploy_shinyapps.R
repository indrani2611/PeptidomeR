if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("Install rsconnect first: install.packages('rsconnect')", call. = FALSE)
}

app_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

app_files <- c(
  "app.R",
  "R/mod_qc.R",
  "R/mod_eda.R",
  "R/mod_differential_peptide.R",
  "R/mod_differential_protein.R",
  "R/mod_multivariate.R",
  "R/mod_proteolytic_motif.R"
)

missing_files <- app_files[
  !file.exists(file.path(app_dir, app_files))
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Run this script from the PeptidomeR project directory. Missing: ",
      paste(missing_files, collapse = ", ")
    ),
    call. = FALSE
  )
}

accounts <- rsconnect::accounts(server = "shinyapps.io")
if (nrow(accounts) == 0) {
  stop(
    paste(
      "No shinyapps.io account is configured.",
      "Run rsconnect::setAccountInfo(name, token, secret) once, then retry."
    ),
    call. = FALSE
  )
}

rsconnect::deployApp(
  appDir = app_dir,
  appFiles = app_files,
  appPrimaryDoc = "app.R",
  appMode = "shiny",
  quarto = FALSE,
  appName = "peptidomer",
  appTitle = "PeptidomeR",
  account = "peptigram",
  server = "shinyapps.io",
  forceUpdate = TRUE,
  launch.browser = FALSE
)
