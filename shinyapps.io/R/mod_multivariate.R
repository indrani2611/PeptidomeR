# mod_multivariate: converted from peptidomics_multivariate.R

mod_multivariate_ui <- function(id) {
  ns <- NS(id)

  tagList(
    h2("Peptidomics Multivariate Analysis"),
    layout_sidebar(

  sidebar = sidebar(
    width = 400,
    title = "Configuration Panel",

    card(
      card_header("1. File Input"),
      fileInput(
        ns("uploaded_file"),
        "Choose Peptide Data File:",
        accept = c(".csv", ".tsv", ".txt")
      )
    ),

    card(
      card_header("2. Quantitative Metric"),
      selectInput(
        ns("target_metric"),
        "Columns to analyse:",
        choices = c("MaxLFQ Intensity", "Intensity", "Spectral Count"),
        selected = "MaxLFQ Intensity"
      ),
      helpText(
        "The app finds columns containing the selected metric and keeps their original file order."
      ),
      verbatimTextOutput(ns("metric_column_check"))
    ),

    card(
      card_header("3. Study Design"),
      numericInput(
        ns("num_experiments"),
        "Number of experiments:",
        value = 3,
        min = 2,
        step = 1
      ),
      uiOutput(ns("experiment_inputs"))
    ),

    card(
      card_header("4. Preprocessing"),
      helpText(
        paste(
          "The uploaded file is assumed to have passed quality control. Missing values are
          imputed separately within each experiment using KNN and QRILC"
        )
      ),
      radioButtons(
        ns("apply_log2"),
        "Convert selected quantitative columns to log2?",
        choices = c(
          "Yes" = "yes",
          "No" = "no"
        ),
        selected = "yes"
      ),
      radioButtons(
        ns("valid_rule"),
        "Valid-value rule before imputation:",
        choices = c(
          "At least 50% valid in one experiment" = "one",
          "At least 50% valid in every experiment" = "every"
        ),
        selected = "every"
      )
    ),

    card(
      card_header("5. Multivariate Settings"),
      selectInput(
        ns("analysis_method"),
        "Select one analysis method:",
        choices = c(
          "Principal Component Analysis (PCA)" = "pca",
          "t-distributed Stochastic Neighbor Embedding (t-SNE)" = "tsne",
          "Uniform Manifold Approximation and Projection (UMAP)" = "umap"
        ),
        selected = "pca"
      ),
      helpText(
        paste(
          "Only the selected method is calculated and plotted. t-SNE perplexity",
          "and UMAP neighbours are selected automatically from sample count."
        )
      )
    ),
    conditionalPanel(
      condition = "input.analysis_method == 'pca'",
      numericInput(
        ns("top_n_heatmap"),
        "Top PCA-loading peptides in heatmap:",
        value = 100,
        min = 5,
        max = 500,
        step = 5
      ),
      helpText(
        "If fewer peptides are available, the app uses all available peptides."
      )
    ,
      ns = ns),
    actionButton(
      ns("run_analysis"),
      "Run Selected Analysis",
      class = "btn-success w-100 fw-bold mt-2"
    )
  ),

  card(
    card_header("Analysis Summary"),
    verbatimTextOutput(ns("analysis_summary")),
    card_footer(
      downloadButton(
        ns("download_results_zip"),
        "Download All Tables and Plots (.zip)",
        class = "btn-primary w-100 fw-bold"
      )
    )
  ),

  uiOutput(ns("analysis_display_ui"))
,
      fillable = FALSE,
      fill = FALSE,
      border = FALSE
    )
  )
}

mod_multivariate_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns


  # ------------------------------------------------------------------------
  # STEP 2A. BUILD EXPERIMENT NAME AND REPLICATE INPUTS
  # WHAT THIS DOES: Adds one experiment-name and replicate-count box per group.
  # CHANGE HERE: Change default experiment names or default replicate number.
  # ------------------------------------------------------------------------
  output$experiment_inputs <- renderUI({
    req(input$num_experiments)

    tagList(lapply(seq_len(input$num_experiments), function(i) {
      div(
        class = "p-2 mb-2 border rounded bg-light",
        textInput(
          ns(paste0("experiment_name_", i)),
          paste("Experiment", i, "name:"),
          value = paste0("Experiment", i)
        ),
        numericInput(
          ns(paste0("experiment_reps_", i)),
          paste("Experiment", i, "replicates:"),
          value = 2,
          min = 1,
          step = 1
        )
      )
    }))
  })

  # Collect the dynamic experiment inputs into one data frame.
  experiment_info <- reactive({
    req(input$num_experiments)

    names_entered <- vapply(seq_len(input$num_experiments), function(i) {
      value <- input[[paste0("experiment_name_", i)]]
      if (is.null(value)) "" else trimws(value)
    }, character(1))

    replicates_entered <- vapply(seq_len(input$num_experiments), function(i) {
      value <- input[[paste0("experiment_reps_", i)]]
      if (is.null(value)) NA_integer_ else as.integer(value)
    }, integer(1))

    data.frame(
      Experiment = names_entered,
      Replicates = replicates_entered,
      stringsAsFactors = FALSE
    )
  })

  # ------------------------------------------------------------------------
  # STEP 2C. READ THE UPLOADED FILE
  # WHAT THIS DOES: Reads CSV, TSV, or TXT without changing original headers.
  # CHANGE HERE: Add extra read.delim() arguments if your export needs them.
  # ------------------------------------------------------------------------
  raw_data <- reactive({
    req(input$uploaded_file)

    extension <- tolower(tools::file_ext(input$uploaded_file$name))

    if (!extension %in% c("csv", "tsv", "txt")) {
      stop("input file is not in either .csv/.tsv/.txt", call. = FALSE)
    }

    separator_to_use <- if (extension == "csv") {
      ","
    } else if (extension == "tsv") {
      "\t"
    } else {
      header_lines <- tryCatch(
        readLines(
          input$uploaded_file$datapath,
          n = 20,
          warn = FALSE
        ),
        error = function(e) character(0)
      )
      header_lines <- header_lines[nzchar(trimws(header_lines))]

      if (length(header_lines) == 0) {
        stop("input file doesn't match required columns", call. = FALSE)
      }

      header_line <- header_lines[1]
      tab_count <- lengths(regmatches(
        header_line,
        gregexpr("\t", header_line, fixed = TRUE)
      ))
      comma_count <- lengths(regmatches(
        header_line,
        gregexpr(",", header_line, fixed = TRUE)
      ))

      if (tab_count > comma_count) "\t" else ","
    }

    data <- tryCatch(
      read.delim(
        file = input$uploaded_file$datapath,
        header = TRUE,
        sep = separator_to_use,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        quote = "\"",
        comment.char = "",
        na.strings = c("", "NA", "NaN", "NULL", "null")
      ),
      error = function(e) {
        stop("input file doesn't match required columns", call. = FALSE)
      }
    )

    cleaned_names <- trimws(names(data))
    has_peptide_sequence <- any(
      tolower(cleaned_names) == "peptide sequence"
    )

    quantitative_columns <- which(grepl(
      "MaxLFQ\\s+Intensity|Intensity|Spectral",
      cleaned_names,
      ignore.case = TRUE,
      perl = TRUE
    ))
    has_numeric_quantitative_column <-
      length(quantitative_columns) > 0 &&
      any(vapply(
        data[quantitative_columns],
        is.numeric,
        logical(1)
      ))

    if (!has_peptide_sequence || !has_numeric_quantitative_column) {
      stop("input file doesn't match required columns", call. = FALSE)
    }

    data
  })

  observeEvent(input$uploaded_file, {
    tryCatch(
      raw_data(),
      error = function(e) {
        showNotification(e$message, type = "error", duration = NULL)
      }
    )
  }, ignoreInit = TRUE)

  # ------------------------------------------------------------------------
  # STEP 2D. HELPER FUNCTIONS FOR COLUMN DETECTION AND CLEANING
  # WHAT THIS DOES: Accepts common header variants with spaces, dots, or underscores.
  # CHANGE HERE: Add aliases if your metadata columns use different names.
  # ------------------------------------------------------------------------
  normalise_header <- function(x) {
    tolower(gsub("[^a-z0-9]", "", x))
  }

  find_metadata_column <- function(data, aliases, required = TRUE) {
    headers_normalised <- normalise_header(names(data))
    alias_normalised <- normalise_header(aliases)
    match_index <- match(alias_normalised, headers_normalised, nomatch = 0)
    match_index <- match_index[match_index > 0]

    if (length(match_index) > 0) {
      return(names(data)[match_index[1]])
    }

    if (required) {
      stop(paste0(
        "Missing required metadata column. Accepted names include: ",
        paste(aliases, collapse = ", "), "."
      ))
    }

    NULL
  }

  find_metric_columns <- function(data, metric) {
    column_names <- names(data)
    metric_key <- normalise_header(metric)

    pattern <- switch(
      metric_key,
      "spectralcount" = "Spectral[[:space:]._-]*Count|SpectralCount",
      "intensity" = "(^|[[:space:]._-])Intensity($|[[:space:]._-])|^Intensity$",
      "maxlfqintensity" = "MaxLFQ[[:space:]._-]*Intensity|MaxLFQIntensity",
      "MaxLFQ[[:space:]._-]*Intensity|MaxLFQIntensity"
    )

    matched_columns <- grep(
      pattern,
      column_names,
      ignore.case = TRUE,
      value = TRUE
    )

    # Plain Intensity must not also collect MaxLFQ or Spectral Count columns.
    if (identical(metric_key, "intensity")) {
      matched_columns <- matched_columns[
        !grepl(
          "MaxLFQ|Spectral[[:space:]._-]*Count|SpectralCount",
          matched_columns,
          ignore.case = TRUE
        )
      ]
    }

    # Keep experiment-level Spectral Count columns separate from combined,
    # unique, and total spectral-count summary columns when they are present.
    if (identical(metric_key, "spectralcount")) {
      matched_columns <- matched_columns[
        !grepl(
          "Combined|Unique[[:space:]._-]*Spectral|Total[[:space:]._-]*Spectral",
          matched_columns,
          ignore.case = TRUE
        )
      ]
    }

    matched_columns
  }

  # Show detected columns immediately, before the user runs the analysis.
  output$metric_column_check <- renderPrint({
    if (is.null(input$uploaded_file)) {
      cat("Upload a peptide file to detect quantitative columns.")
      return(invisible(NULL))
    }

    data <- raw_data()
    detected_columns <- find_metric_columns(data, input$target_metric)

    cat("Detected", length(detected_columns), input$target_metric, "columns.\n")
    if (length(detected_columns) > 0) {
      cat(paste0(seq_along(detected_columns), ". ", detected_columns), sep = "\n")
    } else {
      cat("No matching columns found. Check the file delimiter and column headers.")
    }
  })

  to_numeric <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
    suppressWarnings(as.numeric(gsub(",", "", x, fixed = TRUE)))
  }

  safe_file_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_-]+", "_", x)
    gsub("_+", "_", x)
  }

  # ------------------------------------------------------------------------
  # STEP 3. RUN THE COMPLETE ANALYSIS ONLY WHEN THE BUTTON IS CLICKED
  # WHAT THIS DOES: Freezes the current settings and returns all result objects.
  # ------------------------------------------------------------------------
  # ------------------------------------------------------------------------
  # STEP 3. PREPROCESS THE DATA AND RUN THE SELECTED ANALYSIS
  # ------------------------------------------------------------------------
  analysis_results <- eventReactive(input$run_analysis, {
    data <- raw_data()
    info <- experiment_info()

    # STEP 3A. Validate the file, experiment design, and user settings.
    validate(
      need(nrow(data) > 0, "The uploaded file contains no peptide rows."),
      need(all(nzchar(info$Experiment)), "Please give every experiment a name."),
      need(!anyDuplicated(info$Experiment), "Experiment names must be unique."),
      need(all(!is.na(info$Replicates) & info$Replicates >= 1),
           "Every experiment must contain at least one replicate."),
      need(input$apply_log2 %in% c("yes", "no"),
           "Please choose whether log2 transformation is required."),
      need(input$valid_rule %in% c("one", "every", "none"),
           "Please choose a valid-value rule."),
      need(input$analysis_method %in% c("pca", "tsne", "umap"),
           "Please select PCA, t-SNE, or UMAP."),
      need(input$top_n_heatmap >= 1,
           "The requested number of heatmap peptides must be at least one.")
    )

    # STEP 3B. Find peptide metadata and selected quantitative columns.
    sequence_col <- find_metadata_column(
      data,
      c("Peptide Sequence", "Peptide.Sequence", "Sequence")
    )
    protein_col <- find_metadata_column(
      data,
      c("Protein", "Protein Name", "Protein.Name")
    )
    protein_id_col <- find_metadata_column(
      data,
      c("Protein ID", "Protein.ID", "Protein Accession", "Protein.Accession")
    )
    description_col <- find_metadata_column(
      data,
      c("Protein Description", "Protein.Description", "Description")
    )

    metric <- input$target_metric
    metric_cols <- find_metric_columns(data, metric)
    expected_replicates <- sum(info$Replicates)

    validate(
      need(length(metric_cols) > 0,
           paste0("No columns matching '", metric, "' were found.")),
      need(!anyDuplicated(metric_cols),
           "Detected quantitative column names must be unique."),
      need(length(metric_cols) == expected_replicates,
           paste0(
             "The study design contains ", expected_replicates,
             " total replicates, but ", length(metric_cols),
             " '", metric, "' columns were detected. Detected columns: ",
             paste(metric_cols, collapse = ", ")
           ))
    )

    numeric_data <- as.data.frame(
      lapply(data[, metric_cols, drop = FALSE], to_numeric),
      check.names = FALSE
    )
    names(numeric_data) <- metric_cols

    # STEP 3C. Keep the qc.R-filtered rows and optionally apply log2.
    qc_filtered_input <- data
    apply_log2 <- identical(input$apply_log2, "yes")
    if (apply_log2) {
      numeric_data[numeric_data <= 0] <- NA_real_
    }
    analysis_data <- if (apply_log2) log2(numeric_data) else numeric_data
    transformation_label <- if (apply_log2) {
      "Log2 transformation applied by this app"
    } else {
      "No transformation applied; uploaded values treated as already log2"
    }

    metadata <- data.frame(
      Peptide.Sequence = as.character(data[[sequence_col]]),
      Protein = as.character(data[[protein_col]]),
      Protein.ID = as.character(data[[protein_id_col]]),
      Protein.Description = as.character(data[[description_col]]),
      stringsAsFactors = FALSE
    )

    # STEP 3D. Assign quantitative columns to experiments in file order.
    group_columns <- list()
    column_design <- data.frame(
      Original.Column = character(),
      Experiment = character(),
      Replicate = integer(),
      stringsAsFactors = FALSE
    )
    start_index <- 1
    for (i in seq_len(nrow(info))) {
      end_index <- start_index + info$Replicates[i] - 1
      selected_columns <- metric_cols[start_index:end_index]
      group_columns[[info$Experiment[i]]] <- selected_columns
      column_design <- rbind(
        column_design,
        data.frame(
          Original.Column = selected_columns,
          Experiment = info$Experiment[i],
          Replicate = seq_len(info$Replicates[i]),
          stringsAsFactors = FALSE
        )
      )
      start_index <- end_index + 1
    }
    validate(need(
      identical(unname(unlist(group_columns, use.names = FALSE)), metric_cols),
      "Experiment column assignment does not reproduce the detected-column order."
    ))

    selected_data_before_filter <- cbind(metadata, analysis_data)

    # STEP 3E. Apply the selected valid-value rule before imputation.
    required_values <- ceiling(info$Replicates * 0.50)
    valid_count_matrix <- vapply(
      info$Experiment,
      function(group_name) {
        rowSums(!is.na(analysis_data[, group_columns[[group_name]], drop = FALSE]))
      },
      numeric(nrow(analysis_data))
    )
    valid_count_matrix <- as.matrix(valid_count_matrix)
    colnames(valid_count_matrix) <- info$Experiment
    validate(
      need(nrow(valid_count_matrix) == nrow(analysis_data),
           "Valid-count matrix has an incorrect number of peptide rows."),
      need(ncol(valid_count_matrix) == nrow(info),
           "Valid-count matrix has an incorrect number of experiments.")
    )

    passes_matrix <- sweep(valid_count_matrix, 2, required_values, FUN = ">=")
    if (input$valid_rule == "one") {
      keep_rows <- rowSums(passes_matrix) >= 1
      valid_rule_label <- "At least 50% valid in at least one experiment"
    } else if (input$valid_rule == "every") {
      keep_rows <- rowSums(passes_matrix) == ncol(passes_matrix)
      valid_rule_label <- "At least 50% valid in every experiment"
    } else {
      keep_rows <- rep(TRUE, nrow(analysis_data))
      valid_rule_label <- "No additional valid-value filtering"
    }

    filtered_metadata <- metadata[keep_rows, , drop = FALSE]
    filtered_analysis_data <- analysis_data[keep_rows, , drop = FALSE]
    filtered_valid_counts <- valid_count_matrix[keep_rows, , drop = FALSE]
    colnames(filtered_valid_counts) <- paste0(
      "Valid.Values.",
      make.names(colnames(filtered_valid_counts), unique = TRUE)
    )
    rows_after_valid_filter <- cbind(
      filtered_metadata,
      as.data.frame(filtered_valid_counts, check.names = FALSE),
      filtered_analysis_data
    )

    validate(need(
      nrow(filtered_analysis_data) > 0,
      paste0("No peptides passed the selected valid-value rule: ", valid_rule_label)
    ))

    # STEP 3F. Impute separately within each experiment.
    # Partial missing rows use KNN; all-missing experiment rows use QRILC;
    # any values still missing use the current down-shifted-normal fallback.
    imputation_seed <- 123L
    set.seed(imputation_seed)
    imputed_data <- as.data.frame(filtered_analysis_data, check.names = FALSE)
    imputation_log <- character()
    imputation_summary <- data.frame(
      Experiment = character(),
      Replicates = integer(),
      Rows = integer(),
      Partial.Missing.Rows = integer(),
      All.Missing.Rows = integer(),
      Missing.Cells.Before = integer(),
      Missing.Cells.After = integer(),
      KNN.Used = logical(),
      QRILC.Used = logical(),
      Fallback.Used = logical(),
      stringsAsFactors = FALSE
    )

    for (group_name in info$Experiment) {
      group_cols <- group_columns[[group_name]]
      group_matrix <- as.matrix(imputed_data[, group_cols, drop = FALSE])
      storage.mode(group_matrix) <- "numeric"
      group_matrix[is.infinite(group_matrix)] <- NA_real_

      n_replicates <- ncol(group_matrix)
      missing_per_row <- rowSums(is.na(group_matrix))
      partial_rows <- missing_per_row > 0 & missing_per_row < n_replicates
      all_missing_rows <- missing_per_row == n_replicates
      missing_before <- sum(is.na(group_matrix))
      knn_used <- FALSE
      qrilc_used <- FALSE
      fallback_used <- FALSE

      if (any(partial_rows)) {
        knn_rows <- !all_missing_rows
        knn_input <- group_matrix[knn_rows, , drop = FALSE]
        if (sum(is.na(knn_input)) > 0 && nrow(knn_input) >= 2 && ncol(knn_input) >= 2) {
          knn_result <- tryCatch(
            impute::impute.knn(knn_input, rng.seed = imputation_seed),
            error = function(e) NULL
          )
          if (!is.null(knn_result)) {
            group_matrix[knn_rows, ] <- knn_result$data
            knn_used <- TRUE
          }
        }
      }

      if (any(all_missing_rows)) {
        qrilc_result <- tryCatch(
          imputeLCMD::impute.QRILC(group_matrix, tune.sigma = 1),
          error = function(e) NULL
        )
        if (!is.null(qrilc_result)) {
          qrilc_matrix <- if (is.list(qrilc_result)) qrilc_result[[1]] else qrilc_result
          qrilc_matrix <- as.matrix(qrilc_matrix)
          storage.mode(qrilc_matrix) <- "numeric"
          if (identical(dim(qrilc_matrix), dim(group_matrix))) {
            group_matrix[all_missing_rows, ] <-
              qrilc_matrix[all_missing_rows, , drop = FALSE]
            qrilc_used <- TRUE
          }
        }
      }

      still_missing <- is.na(group_matrix)
      if (any(still_missing)) {
        observed <- as.numeric(group_matrix[!is.na(group_matrix)])
        if (length(observed) < 2) {
          observed <- as.numeric(filtered_analysis_data)
          observed <- observed[is.finite(observed)]
        }
        validate(need(
          length(observed) >= 2,
          paste0("Not enough observed values to complete imputation for ", group_name)
        ))
        observed_sd <- stats::sd(observed, na.rm = TRUE)
        if (!is.finite(observed_sd) || observed_sd == 0) observed_sd <- 1
        group_matrix[still_missing] <- stats::rnorm(
          sum(still_missing),
          mean = mean(observed, na.rm = TRUE) - 1.8 * observed_sd,
          sd = 0.3 * observed_sd
        )
        fallback_used <- TRUE
      }

      missing_after <- sum(is.na(group_matrix))
      validate(need(
        missing_after == 0,
        paste0("Missing values remain after imputation in ", group_name)
      ))
      imputed_data[, group_cols] <- group_matrix
      imputation_summary <- rbind(
        imputation_summary,
        data.frame(
          Experiment = group_name,
          Replicates = n_replicates,
          Rows = nrow(group_matrix),
          Partial.Missing.Rows = sum(partial_rows),
          All.Missing.Rows = sum(all_missing_rows),
          Missing.Cells.Before = missing_before,
          Missing.Cells.After = missing_after,
          KNN.Used = knn_used,
          QRILC.Used = qrilc_used,
          Fallback.Used = fallback_used,
          stringsAsFactors = FALSE
        )
      )
      imputation_log <- c(
        imputation_log,
        paste0(
          group_name,
          ": partial rows = ", sum(partial_rows),
          "; all-missing rows = ", sum(all_missing_rows),
          "; KNN = ", knn_used,
          "; QRILC = ", qrilc_used,
          "; fallback = ", fallback_used,
          "; final missing cells = ", missing_after
        )
      )
    }

    imputed_output <- cbind(filtered_metadata, imputed_data)
    validate(need(
      all(is.finite(as.matrix(imputed_data[, metric_cols, drop = FALSE]))),
      "The completed imputed matrix contains missing or non-finite values."
    ))

    # STEP 3G. Build the shared samples-by-peptides matrix.
    sample_matrix_raw <- t(as.matrix(imputed_data[, metric_cols, drop = FALSE]))
    storage.mode(sample_matrix_raw) <- "numeric"
    feature_sd <- apply(sample_matrix_raw, 2, stats::sd, na.rm = TRUE)
    features_used <- is.finite(feature_sd) & feature_sd > 0 &
      colSums(is.na(sample_matrix_raw)) == 0

    validate(need(
      sum(features_used) >= 2,
      "At least two non-constant peptides are required after imputation."
    ))

    method_matrix <- sample_matrix_raw[, features_used, drop = FALSE]
    validate(
      need(nrow(method_matrix) == length(metric_cols),
           "The analysis matrix has an incorrect number of samples."),
      need(ncol(method_matrix) == sum(features_used),
           "The analysis matrix has an incorrect number of peptide features.")
    )

    scaled_method_matrix <- scale(method_matrix, center = TRUE, scale = TRUE)
    validate(need(
      all(is.finite(scaled_method_matrix)),
      "The scaled analysis matrix contains non-finite values."
    ))

    sample_info <- column_design[
      match(rownames(method_matrix), column_design$Original.Column),
      ,
      drop = FALSE
    ]
    validate(
      need(!any(is.na(sample_info$Experiment)),
           "One or more samples could not be mapped to an experiment."),
      need(identical(as.character(sample_info$Original.Column), rownames(method_matrix)),
           "Sample metadata order does not match the analysis matrix.")
    )

    selected_method <- input$analysis_method
    sample_count <- nrow(method_matrix)
    feature_metadata <- filtered_metadata[features_used, , drop = FALSE]
    method_tables <- list()
    method_plot <- NULL
    display_table <- NULL
    heatmap_object <- NULL
    screen_heatmap_object <- NULL
    heatmap_top_n <- 0L
    method_details <- character()
    method_seed <- NA_integer_

    # STEP 3H. Run only the method selected by the user.
    if (selected_method == "pca") {
      validate(need(sample_count >= 2, "PCA requires at least two samples."))

      pca_result <- stats::prcomp(method_matrix, center = TRUE, scale. = TRUE)
      validate(
        need(ncol(pca_result$x) >= 2, "PCA returned fewer than two components."),
        need(nrow(pca_result$x) == sample_count,
             "PCA returned an incorrect number of sample scores."),
        need(nrow(pca_result$rotation) == ncol(method_matrix),
             "PCA returned an incorrect number of peptide loadings.")
      )

      pca_variance <- summary(pca_result)$importance[2, ] * 100
      pca_variance_table <- data.frame(
        Component = names(pca_variance),
        Variance.Percent = as.numeric(pca_variance),
        stringsAsFactors = FALSE
      )
      pca_scores <- data.frame(
        Sample = rownames(pca_result$x),
        Experiment = sample_info$Experiment,
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2],
        stringsAsFactors = FALSE
      )
      pca_loadings <- data.frame(
        Row.Index = which(features_used),
        feature_metadata,
        PC1.Loading = pca_result$rotation[, 1],
        PC2.Loading = pca_result$rotation[, 2],
        stringsAsFactors = FALSE
      )
      pca_loadings$Contribution.PC1.PC2 <- sqrt(
        pca_loadings$PC1.Loading^2 + pca_loadings$PC2.Loading^2
      )
      pca_loadings <- pca_loadings[
        order(-pca_loadings$Contribution.PC1.PC2),
        ,
        drop = FALSE
      ]

      # PCA-specific heatmap of the user-selected number of strongest
      # PC1/PC2 loading peptides.
      heatmap_top_n <- min(as.integer(input$top_n_heatmap), nrow(pca_loadings))
      top_heatmap_peptides <- head(pca_loadings, heatmap_top_n)
      heatmap_matrix <- as.matrix(
        imputed_data[top_heatmap_peptides$Row.Index, metric_cols, drop = FALSE]
      )
      storage.mode(heatmap_matrix) <- "numeric"
      heatmap_labels <- trimws(as.character(top_heatmap_peptides$Peptide.Sequence))
      missing_labels <- is.na(heatmap_labels) | !nzchar(heatmap_labels)
      heatmap_labels[missing_labels] <- paste0(
        "Peptide_",
        top_heatmap_peptides$Row.Index[missing_labels]
      )
      rownames(heatmap_matrix) <- make.unique(heatmap_labels)
      sample_annotation <- data.frame(
        Experiment = sample_info$Experiment,
        row.names = sample_info$Original.Column,
        stringsAsFactors = FALSE
      )
      validate(
        need(nrow(heatmap_matrix) == heatmap_top_n,
             "PCA heatmap has an incorrect number of peptide rows."),
        need(ncol(heatmap_matrix) == sample_count,
             "PCA heatmap has an incorrect number of sample columns."),
        need(identical(rownames(sample_annotation), colnames(heatmap_matrix)),
             "PCA heatmap annotations do not match the sample columns.")
      )
      heatmap_colours <- colorRampPalette(
        c("#B2DF8A", "#FFFFFF", "#FB9A99")
      )(100)
      heatmap_object <- pheatmap::pheatmap(
        heatmap_matrix,
        scale = "row",
        color = heatmap_colours,
        annotation_col = sample_annotation,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize_row = 6,
        fontsize_col = 8,
        main = paste("Top", heatmap_top_n, "Peptides Driving PCA Separation"),
        silent = TRUE
      )

        screen_cell_width <- max(35, min(90, floor(900 / sample_count)))
        screen_heatmap_object <- pheatmap::pheatmap(
        heatmap_matrix,
        scale = "row",
        color = heatmap_colours,
        annotation_col = sample_annotation,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize_row = 8,
        fontsize_col = 10,
        cellheight = 12,
        cellwidth = screen_cell_width,
        border_color = NA,
        main = paste("Top", heatmap_top_n, "Peptides Driving PCA Separation"),
        silent = TRUE
      )

      method_plot <- ggplot(
        pca_scores,
        aes(PC1, PC2, colour = Experiment, label = Sample)
      ) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
        theme_bw() +
        labs(
          title = "PCA of Imputed Peptide Data",
          x = paste0("PC1 (", round(pca_variance[1], 2), "%)"),
          y = paste0("PC2 (", round(pca_variance[2], 2), "%)")
        ) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))

      method_tables <- list(
        PCA_sample_scores = pca_scores,
        PCA_explained_variance = pca_variance_table,
        PCA_all_peptide_loadings = pca_loadings
      )
      method_tables[[paste0("PCA_top_", heatmap_top_n, "_peptides")]] <-
        top_heatmap_peptides
      display_table <- pca_loadings
      method_name <- "Principal Component Analysis (PCA)"
      method_stub <- "PCA"
      plot_header <- "PCA Plot"
      table_header <- "Top Peptides Driving PCA Separation"
      method_details <- c(
        paste0("PC1 variance: ", round(pca_variance[1], 2), "%"),
        paste0("PC2 variance: ", round(pca_variance[2], 2), "%"),
        paste0("PCA heatmap peptides: ", heatmap_top_n)
      )

    } else if (selected_method == "tsne") {
      validate(need(
        sample_count >= 5,
        "t-SNE requires at least five samples."
      ))
      tsne_perplexity <- max(1L, min(30L, floor((sample_count - 2L) / 3L)))
      validate(need(
        3 * tsne_perplexity < sample_count - 1,
        "The sample count cannot produce a valid t-SNE perplexity."
      ))

      method_seed <- 123L
      set.seed(method_seed)
      tsne_result <- Rtsne::Rtsne(
        as.matrix(scaled_method_matrix),
        dims = 2,
        perplexity = tsne_perplexity,
        initial_dims = min(50L, ncol(scaled_method_matrix)),
        check_duplicates = FALSE,
        pca = TRUE,
        max_iter = 1000,
        verbose = FALSE
      )
      validate(
        need(nrow(tsne_result$Y) == sample_count,
             "t-SNE returned an incorrect number of sample coordinates."),
        need(ncol(tsne_result$Y) >= 2,
             "t-SNE returned fewer than two dimensions."),
        need(all(is.finite(tsne_result$Y[, 1:2, drop = FALSE])),
             "t-SNE returned non-finite coordinates.")
      )

      tsne_coordinates <- data.frame(
        Sample = rownames(scaled_method_matrix),
        Experiment = sample_info$Experiment,
        tSNE1 = tsne_result$Y[, 1],
        tSNE2 = tsne_result$Y[, 2],
        stringsAsFactors = FALSE
      )
      method_plot <- ggplot(
        tsne_coordinates,
        aes(tSNE1, tSNE2, colour = Experiment, label = Sample)
      ) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
        theme_bw() +
        labs(
          title = paste0("t-SNE of Imputed Peptide Data (perplexity = ",
                         tsne_perplexity, ")"),
          x = "t-SNE 1",
          y = "t-SNE 2"
        ) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))

      method_tables <- list(tSNE_sample_coordinates = tsne_coordinates)
      display_table <- tsne_coordinates
      method_name <- "t-distributed Stochastic Neighbor Embedding (t-SNE)"
      method_stub <- "tSNE"
      plot_header <- "t-SNE Plot"
      table_header <- "t-SNE Sample Coordinates"
      method_details <- c(
        paste0("Automatic perplexity: ", tsne_perplexity),
        paste0("Fixed seed: ", method_seed)
      )

    } else {
      validate(need(
        sample_count >= 3,
        "UMAP requires at least three samples in this app."
      ))
      umap_neighbours <- min(15L, sample_count - 1L)
      method_seed <- 123L
      set.seed(method_seed)
      umap_coordinates_matrix <- uwot::umap(
        as.matrix(scaled_method_matrix),
        n_neighbors = umap_neighbours,
        n_components = 2,
        metric = "euclidean",
        min_dist = 0.1,
        spread = 1,
        n_threads = 1,
        verbose = FALSE,
        ret_model = FALSE
      )
      validate(
        need(nrow(umap_coordinates_matrix) == sample_count,
             "UMAP returned an incorrect number of sample coordinates."),
        need(ncol(umap_coordinates_matrix) >= 2,
             "UMAP returned fewer than two dimensions."),
        need(all(is.finite(umap_coordinates_matrix[, 1:2, drop = FALSE])),
             "UMAP returned non-finite coordinates.")
      )

      umap_coordinates <- data.frame(
        Sample = rownames(scaled_method_matrix),
        Experiment = sample_info$Experiment,
        UMAP1 = umap_coordinates_matrix[, 1],
        UMAP2 = umap_coordinates_matrix[, 2],
        stringsAsFactors = FALSE
      )
      method_plot <- ggplot(
        umap_coordinates,
        aes(UMAP1, UMAP2, colour = Experiment, label = Sample)
      ) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
        theme_bw() +
        labs(
          title = paste0("UMAP of Imputed Peptide Data (neighbours = ",
                         umap_neighbours, ")"),
          x = "UMAP 1",
          y = "UMAP 2"
        ) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))

      method_tables <- list(UMAP_sample_coordinates = umap_coordinates)
      display_table <- umap_coordinates
      method_name <- "Uniform Manifold Approximation and Projection (UMAP)"
      method_stub <- "UMAP"
      plot_header <- "UMAP Plot"
      table_header <- "UMAP Sample Coordinates"
      method_details <- c(
        paste0("Automatic neighbours: ", umap_neighbours),
        "Metric: Euclidean",
        "Minimum distance: 0.1",
        paste0("Fixed seed: ", method_seed)
      )
    }

    validate(
      need(inherits(method_plot, "ggplot"),
           "The selected analysis did not create a ggplot object."),
      need(is.data.frame(display_table) && nrow(display_table) > 0,
           "The selected analysis did not create a results table.")
    )

    parameter_table <- data.frame(
      Parameter = c(
        "Input file", "Metric", "Experiments", "Total samples",
        "Transformation", "Valid-value rule", "Imputation sequence",
        "Imputation method", "Selected analysis", "Non-constant peptides",
        "Method details", "Method", "Requested heatmap peptides",
        "Application version"
      ),
      Value = c(
        input$uploaded_file$name, metric, paste(info$Experiment, collapse = "; "),
        sample_count, transformation_label, valid_rule_label,
        "KNN partial; QRILC all-missing experiment; down-shifted fallback",
        imputation_seed, method_name, sum(features_used),
        paste(method_details, collapse = "; "), method_seed,
        if (selected_method == "pca") input$top_n_heatmap else "Not applicable",
        "6"
      ),
      stringsAsFactors = FALSE
    )

    list(
      qc_filtered_input = qc_filtered_input,
      selected_data_before_filter = selected_data_before_filter,
      rows_after_valid_filter = rows_after_valid_filter,
      imputed_output = imputed_output,
      imputation_summary = imputation_summary,
      imputation_log = imputation_log,
      column_design = column_design,
      parameters = parameter_table,
      method_tables = method_tables,
      method_plot = method_plot,
      display_table = display_table,
      heatmap_object = heatmap_object,
      screen_heatmap_object = screen_heatmap_object,
      heatmap_available = !is.null(heatmap_object),
      heatmap_top_n = heatmap_top_n,
      method_name = method_name,
      method_stub = method_stub,
      method_details = method_details,
      plot_header = plot_header,
      metric_columns = metric_cols,
      metric = metric,
      transformation_label = transformation_label,
      valid_rule_label = valid_rule_label,
      total_input_rows = nrow(data),
      total_retained_rows = nrow(filtered_analysis_data),
      features_used = sum(features_used),
      sample_count = sample_count
    )
  })

  # STEP 4A. Build the results area for the selected method.
  # PCA shows the score plot and its peptide heatmap side by side. t-SNE and
  # UMAP show only their selected plot. Peptide ranking tables are download-only.
  output$analysis_display_ui <- renderUI({
    if (is.null(input$uploaded_file) || is.null(input$run_analysis) ||
        input$run_analysis < 1) {
      return(
        card(
          card_header("Selected Analysis Plot"),
          
        )
      )
    }

    res <- analysis_results()
    plot_card <- card(
      full_screen = TRUE,
      card_header(res$plot_header),
      plotOutput(ns("analysis_plot"), height = "700px")
    )

    if (!isTRUE(res$heatmap_available)) {
      return(plot_card)
    }

    screen_height <- max(900L, 260L + res$heatmap_top_n * 13L)
    screen_width <- max(900L, 450L + res$sample_count * 70L)
    heatmap_card <- card(
      full_screen = TRUE,
      card_header(
        paste("Top", res$heatmap_top_n, "PCA-Loading Peptide Heatmap")
      ),
      helpText(
        paste(
          "Scroll vertically and horizontally to inspect all peptide and",
          "sample labels. The peptide lists are available in the ZIP download."
        )
      ),
      div(
        style = paste(
          "overflow: auto; height: 700px;",
          "border: 1px solid #DEE2E6; background: white; padding: 8px;"
        ),
        plotOutput(
          ns("heatmap_plot"),
          width = paste0(screen_width, "px"),
          height = paste0(screen_height, "px")
        )
      )
    )

    layout_column_wrap(
      width = 1 / 2,
      plot_card,
      heatmap_card
    )
  })

  output$analysis_summary <- renderPrint({
    if (is.null(input$uploaded_file)) {
      
      return(invisible(NULL))
    }
    if (is.null(input$run_analysis) || input$run_analysis < 1) {
      cat("File uploaded. Configure the analysis and click Run Selected Analysis.")
      return(invisible(NULL))
    }
    res <- analysis_results()
    cat("===============================================================\n")
    cat("              PEPTIDOMICS MULTIVARIATE SUMMARY                 \n")
    cat("===============================================================\n")
    cat("Selected method:", res$method_name, "\n")
    cat("Metric:", res$metric, "\n")
    cat("Transformation:", res$transformation_label, "\n")
    cat("Detected quantitative columns:", length(res$metric_columns), "\n")
    cat("Samples:", res$sample_count, "\n")
    cat("QC-filtered input peptide rows:", res$total_input_rows, "\n")
    cat("Valid-value rule:", res$valid_rule_label, "\n")
    cat("Peptide rows retained before imputation:", res$total_retained_rows, "\n")
    cat("Non-constant peptides used:", res$features_used, "\n")
    cat("Method details:", paste(res$method_details, collapse = "; "), "\n\n")
    cat("IMPUTATION LOG\n")
    cat(paste(res$imputation_log, collapse = "\n"), "\n")
    cat("===============================================================\n")
  })
  output$analysis_plot <- renderPlot({
    req(input$uploaded_file)
    req(input$run_analysis > 0)
    res <- analysis_results()
    req(res$method_plot)
    res$method_plot
  })
  output$heatmap_plot <- renderPlot({
    req(input$uploaded_file)
    req(input$run_analysis > 0)
    res <- analysis_results()
    req(isTRUE(res$heatmap_available))
    req(res$screen_heatmap_object)
    grid::grid.newpage()
    grid::grid.draw(res$screen_heatmap_object$gtable)
  })

  # STEP 4B. Download shared intermediates plus only the selected method.
  output$download_results_zip <- downloadHandler(
    filename = function() {
      res <- analysis_results()
      paste0("peptidomics_", res$method_stub, "_", Sys.Date(), ".zip")
    },
    content = function(file) {
      res <- analysis_results()
      output_directory <- tempfile("peptidomics_multivariate_")
      on.exit(unlink(output_directory, recursive = TRUE, force = TRUE), add = TRUE)
      intermediate_directory <- file.path(output_directory, "intermediate_csv")
      results_directory <- file.path(output_directory, "results_csv")
      plots_directory <- file.path(output_directory, "plots")
      dir.create(intermediate_directory, recursive = TRUE)
      dir.create(results_directory, recursive = TRUE)
      dir.create(plots_directory, recursive = TRUE)

      intermediate_files <- c(
        qc_input = file.path(intermediate_directory, "01_qc_filtered_input_as_uploaded.csv"),
        selected = file.path(intermediate_directory, "02_selected_quantitative_data_before_filter.csv"),
        filtered = file.path(intermediate_directory, "03_rows_after_valid_value_filter.csv"),
        imputed = file.path(intermediate_directory, "04_imputed_analysis_data.csv"),
        imputation = file.path(intermediate_directory, "05_imputation_summary.csv"),
        mapping = file.path(intermediate_directory, "06_column_to_experiment_mapping.csv"),
        parameters = file.path(intermediate_directory, "07_analysis_parameters.csv")
      )
      write.csv(res$qc_filtered_input, intermediate_files[["qc_input"]], row.names = FALSE)
      write.csv(res$selected_data_before_filter, intermediate_files[["selected"]], row.names = FALSE)
      write.csv(res$rows_after_valid_filter, intermediate_files[["filtered"]], row.names = FALSE)
      write.csv(res$imputed_output, intermediate_files[["imputed"]], row.names = FALSE)
      write.csv(res$imputation_summary, intermediate_files[["imputation"]], row.names = FALSE)
      write.csv(res$column_design, intermediate_files[["mapping"]], row.names = FALSE)
      write.csv(res$parameters, intermediate_files[["parameters"]], row.names = FALSE)

      for (table_name in names(res$method_tables)) {
        write.csv(
          res$method_tables[[table_name]],
          file.path(results_directory, paste0(table_name, ".csv")),
          row.names = FALSE
        )
      }

      ggsave(
        file.path(plots_directory, paste0(res$method_stub, "_plot.png")),
        res$method_plot,
        width = 8,
        height = 6,
        dpi = 300
      )

      if (isTRUE(res$heatmap_available)) {
        heatmap_file_stub <- paste0(
          "PCA_top_", res$heatmap_top_n, "_peptide_heatmap"
        )
        png(
          file.path(plots_directory, paste0(heatmap_file_stub, ".png")),
          width = 10,
          height = 12,
          units = "in",
          res = 300
        )
        grid::grid.newpage()
        grid::grid.draw(res$heatmap_object$gtable)
        dev.off()
      }

      zip_contents <- list.files(output_directory, recursive = TRUE, full.names = FALSE)
      zip::zipr(zipfile = file, files = zip_contents, root = output_directory)
    }
  )

  })
}
