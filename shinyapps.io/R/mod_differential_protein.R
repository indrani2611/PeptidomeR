# mod_differential_protein: converted from differential_protein.R

mod_differential_protein_ui <- function(id) {
  ns <- NS(id)

  tagList(
    h2("Proteomics Differential Protein Analysis"),
    layout_sidebar(

  sidebar = sidebar(
    width = 400,
    title = "Configuration Panel",

    card(
      card_header("1. File Input"),
      fileInput(
        ns("uploaded_file"),
        "Choose Protein Data File (.csv/.tsv/.txt):",
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
      card_header("4. Comparison"),
      uiOutput(ns("comparison_inputs"))
    ),

    card(
      card_header("5. Statistics and Imputation"),
      helpText(
        paste(
          "For the two selected groups, proteins must have valid values in at least",
          "50% of replicates in EACH group. Remaining missing cells are imputed by",
          "KNN separately within each group."
        )
      ),
      radioButtons(
        ns("apply_log2"),
        "Convert the selected quantitative columns to log2?",
        choices = c(
          "Yes - values are not yet log2" = "yes",
          "No - values are already log2" = "no"
        ),
        selected = "yes"
      ),
      selectInput(
        ns("correction_method"),
        "P-value correction:",
        choices = c(
          "Benjamini-Hochberg (BH)" = "BH",
          "No correction" = "none"
        ),
        selected = "BH"
      ),
      numericInput(
        ns("p_threshold"),
        "Adjusted p-value threshold:",
        value = 0.05,
        min = 0,
        max = 1,
        step = 0.01
      )
    ),

    actionButton(
      ns("run_analysis"),
      "Run Differential Analysis",
      class = "btn-success w-100 fw-bold mt-2"
    )
  ),

  card(
    card_header("Analysis Summary"),
    verbatimTextOutput(ns("analysis_summary")),
    card_footer(
      downloadButton(
        ns("download_results_zip"),
        "Download All Results (.zip)",
        class = "btn-primary w-100 fw-bold"
      )
    )
  ),

  layout_column_wrap(
    width = 1 / 2,
    card(
      card_header("Volcano Plot"),
      plotOutput(ns("volcano_plot"), height = "520px")
    ),
    card(
      card_header("Top 25 Differential Proteins"),
      div(
        style = "overflow-x: auto;",
        tableOutput(ns("results_table"))
      )
    )
  )
,
      fillable = FALSE,
      fill = FALSE,
      border = FALSE
    )
  )
}

mod_differential_protein_server <- function(id) {
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
  # STEP 2B. BUILD CONTROL AND TEST GROUP MENUS
  # WHAT THIS DOES: Uses the experiment names entered above as comparison choices.
  # IMPORTANT: Fold change is always calculated as Test minus Control on log2 data.
  # ------------------------------------------------------------------------
  output$comparison_inputs <- renderUI({
    info <- experiment_info()
    choices <- info$Experiment[nzchar(info$Experiment)]

    if (length(choices) < 2) {
      return(helpText("Enter at least two experiment names."))
    }

    tagList(
      selectInput(
        ns("control_group"),
        "Group 1 / control group:",
        choices = choices,
        selected = choices[1]
      ),
      selectInput(
        ns("test_group"),
        "Group 2 / test group:",
        choices = choices,
        selected = choices[2]
      )
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

    # Protein reports also contain combined, unique, and total spectral-count
    # summary columns. Keep only the ordinary experiment-level Spectral Count
    # columns, such as "Lentil_ctl_1 Spectral Count".
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
      cat("Upload a protein file to detect quantitative columns.")
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
  analysis_results <- eventReactive(input$run_analysis, {
    data <- raw_data()
    info <- experiment_info()

    # ----------------------------------------------------------------------
    # STEP 3A. VALIDATE USER INPUTS
    # WHAT THIS DOES: Stops early with a useful message if the design is invalid.
    # ----------------------------------------------------------------------
    validate(
      need(nrow(data) > 0, "The uploaded file contains no protein rows."),
      need(all(nzchar(info$Experiment)), "Please give every experiment a name."),
      need(!anyDuplicated(info$Experiment), "Experiment names must be unique."),
      need(all(!is.na(info$Replicates) & info$Replicates >= 1),
           "Every experiment must have at least one replicate."),
      need(!is.null(input$control_group), "Please choose a control group."),
      need(!is.null(input$test_group), "Please choose a test group."),
      need(input$control_group != input$test_group,
           "Control and test groups must be different."),
      need(input$p_threshold > 0 && input$p_threshold <= 1,
           "The p-value threshold must be greater than 0 and at most 1."),
      need(input$apply_log2 %in% c("yes", "no"),
           "Please choose whether the selected columns require log2 transformation.")
    )

    control_group <- input$control_group
    test_group <- input$test_group
    metric <- input$target_metric

    selected_replicates <- info$Replicates[
      match(c(control_group, test_group), info$Experiment)
    ]
    validate(need(
      all(selected_replicates >= 2),
      "The control and test groups must each contain at least two replicates for limma testing."
    ))

    # ----------------------------------------------------------------------
    # STEP 3B. FIND METADATA AND QUANTITATIVE COLUMNS
    # WHAT THIS DOES: Finds three protein identifiers and all selected metric columns.
    # IMPORTANT: Metric columns are assigned to experiments in their FILE ORDER.
    # ----------------------------------------------------------------------
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
      c("Description", "Protein Description", "Protein.Description")
    )

    metric_cols <- find_metric_columns(data, metric)
    expected_replicates <- sum(info$Replicates)

    validate(
      need(length(metric_cols) > 0,
           paste0("No columns matching '", metric, "' were found.")),
      need(length(metric_cols) == expected_replicates,
           paste0(
             "The study design contains ", expected_replicates,
             " total replicates, but ", length(metric_cols),
             " '", metric, "' columns were detected. Detected columns: ",
             paste(metric_cols, collapse = ", ")
           ))
    )

    # Convert only selected quantitative columns to numeric.
    numeric_data <- as.data.frame(
      lapply(data[, metric_cols, drop = FALSE], to_numeric),
      check.names = FALSE
    )
    names(numeric_data) <- metric_cols

    # ----------------------------------------------------------------------
    # STEP 3C. PREPARE THE QC-FILTERED INPUT WITHOUT ADDITIONAL ROW FILTERING
    # WHAT THIS DOES: Keeps every row supplied by the user because missing-value
    # filtering and contaminant filtering have already been completed in qc.R.
    # If the user requests log2 transformation, zero/negative raw values are marked
    # NA because log2 is undefined for them. If the input is already log2, zero and
    # negative values are preserved as valid measurements.
    # ----------------------------------------------------------------------
    qc_filtered_input <- data
    apply_log2 <- identical(input$apply_log2, "yes")
    if (apply_log2) {
      numeric_data[numeric_data <= 0] <- NA_real_
    }
    total_initial <- nrow(data)

    # ----------------------------------------------------------------------
    # STEP 3D. OPTIONALLY LOG2 TRANSFORM AND ASSIGN COLUMNS TO EXPERIMENTS
    # WHAT THIS DOES: Applies log2 only when requested, then maps columns in file order.
    # CHANGE HERE: The user's Yes/No choice is handled in the next statement.
    # ----------------------------------------------------------------------
    analysis_data <- if (apply_log2) log2(numeric_data) else numeric_data
    transformation_label <- if (apply_log2) {
      "Log2 transformation applied by this app"
    } else {
      "No transformation applied; uploaded values treated as already log2"
    }
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

    control_cols <- group_columns[[control_group]]
    test_cols <- group_columns[[test_group]]
    comparison_cols <- c(control_cols, test_cols)

    total_control_detected <- sum(
      rowSums(!is.na(analysis_data[, control_cols, drop = FALSE])) > 0
    )
    total_test_detected <- sum(
      rowSums(!is.na(analysis_data[, test_cols, drop = FALSE])) > 0
    )

    # ----------------------------------------------------------------------
    # STEP 3E. APPLY THE FIXED 50% VALID-VALUE RULE TO BOTH SELECTED GROUPS
    # WHAT THIS DOES: Keeps a protein only when at least half of the replicates
    # are valid in EACH selected group. ceiling() gives the required behaviour:
    # 4 replicates -> 2 valid; 5 replicates -> 3 valid.
    # CHANGE HERE: Replace 0.50 below only if the pairwise rule changes later.
    # ----------------------------------------------------------------------
    valid_fraction <- 0.50
    control_needed <- ceiling(length(control_cols) * valid_fraction)
    test_needed <- ceiling(length(test_cols) * valid_fraction)

    valid_control <- rowSums(!is.na(analysis_data[, control_cols, drop = FALSE]))
    valid_test <- rowSums(!is.na(analysis_data[, test_cols, drop = FALSE]))

    passes_control <- valid_control >=
      control_needed
    passes_test <- valid_test >=
      test_needed
    keep_for_testing <- passes_control & passes_test

    metadata <- data.frame(
      Protein = as.character(data[[protein_col]]),
      Protein.ID = as.character(data[[protein_id_col]]),
      Description = as.character(data[[description_col]]),
      stringsAsFactors = FALSE
    )

    selected_groups_before_filter <- cbind(
      metadata,
      analysis_data[, comparison_cols, drop = FALSE]
    )

    filtered_metadata <- metadata[keep_for_testing, , drop = FALSE]
    filtered_analysis_data <- analysis_data[
      keep_for_testing,
      comparison_cols,
      drop = FALSE
    ]
    pairwise_filter_table <- cbind(
      filtered_metadata,
      Valid.Values.Control = valid_control[keep_for_testing],
      Required.Values.Control = control_needed,
      Valid.Values.Test = valid_test[keep_for_testing],
      Required.Values.Test = test_needed,
      filtered_analysis_data
    )
    total_after_valid_filter <- nrow(filtered_metadata)

    validate(need(
      total_after_valid_filter > 0,
      paste0(
        "No proteins have valid values in at least 50% of replicates in both ",
        control_group, " and ", test_group, "."
      )
    ))

    # ----------------------------------------------------------------------
    # STEP 3F. APPLY KNN IMPUTATION SEPARATELY WITHIN EACH SELECTED GROUP
    # WHAT THIS DOES: Runs KNN once for the control replicate matrix and once for
    # the test replicate matrix. Therefore, donors from one experiment are never
    # used to impute missing cells in the other experiment.
    # EXAMPLE: For four replicates, a retained 2/4 row has two cells imputed;
    # a retained 3/4 row has one cell imputed.
    # CHANGE HERE: Change k_max below to use a different maximum KNN neighbour count.
    # ----------------------------------------------------------------------
    knn_random_seed <- 123L
    set.seed(knn_random_seed)
    imputed_data <- as.data.frame(filtered_analysis_data, check.names = FALSE)
    imputation_log <- character()
    imputation_summary <- data.frame(
      Experiment = character(),
      Replicates = integer(),
      Rows = integer(),
      Missing.Cells.Before.KNN = integer(),
      KNN.Neighbours = integer(),
      Missing.Cells.After.KNN = integer(),
      stringsAsFactors = FALSE
    )
    k_max <- 10L

    for (group_name in c(control_group, test_group)) {
      group_cols <- group_columns[[group_name]]
      group_matrix <- as.matrix(imputed_data[, group_cols, drop = FALSE])
      storage.mode(group_matrix) <- "numeric"
      group_matrix[is.infinite(group_matrix)] <- NA_real_

      missing_before <- sum(is.na(group_matrix))
      k_used <- 0L

      if (missing_before > 0) {
        validate(need(
          nrow(group_matrix) >= 2,
          paste0(
            "KNN imputation for ", group_name,
            " requires at least two retained protein rows."
          )
        ))

        k_used <- min(k_max, nrow(group_matrix) - 1L)
        knn_error <- NULL
        knn_result <- tryCatch(
          impute::impute.knn(
            group_matrix,
            k = k_used,
            rowmax = 0.50,
            colmax = 1,
            rng.seed = knn_random_seed
          ),
          error = function(e) {
            knn_error <<- conditionMessage(e)
            NULL
          }
        )

        validate(need(
          !is.null(knn_result),
          paste0("KNN imputation failed for ", group_name, ": ", knn_error)
        ))
        group_matrix <- knn_result$data
      }

      missing_after <- sum(is.na(group_matrix))
      validate(need(
        missing_after == 0,
        paste0(
          "KNN imputation left ", missing_after, " missing cells in ",
          group_name, ". Check whether an entire replicate column is missing."
        )
      ))

      imputed_data[, group_cols] <- group_matrix
      imputation_summary <- rbind(
        imputation_summary,
        data.frame(
          Experiment = group_name,
          Replicates = ncol(group_matrix),
          Rows = nrow(group_matrix),
          Missing.Cells.Before.KNN = missing_before,
          KNN.Neighbours = k_used,
          Missing.Cells.After.KNN = missing_after,
          stringsAsFactors = FALSE
        )
      )
      imputation_log <- c(
        imputation_log,
        paste0(
          group_name,
          ": missing cells before KNN = ", missing_before,
          "; KNN neighbours = ", k_used,
          "; missing cells after KNN = ", missing_after
        )
      )
    }

    validate(need(
      sum(is.na(imputed_data)) == 0,
      "Missing values remain after imputation."
    ))

    # ----------------------------------------------------------------------
    # STEP 3G. FIT THE LIMMA MODEL
    # WHAT THIS DOES: Tests Test group versus Control group using empirical Bayes.
    # IMPORTANT: Positive log2FC means higher in Test; negative means higher in Control.
    # CHANGE HERE: Replace this block if a different statistical model is required.
    # ----------------------------------------------------------------------
    group_factor <- factor(
      c(
        rep(control_group, length(control_cols)),
        rep(test_group, length(test_cols))
      ),
      levels = c(control_group, test_group)
    )

    design <- stats::model.matrix(~ 0 + group_factor)
    safe_group_names <- make.names(c(control_group, test_group), unique = TRUE)
    colnames(design) <- safe_group_names

    fit <- limma::lmFit(as.matrix(imputed_data), design)
    contrast <- limma::makeContrasts(
      contrasts = paste0(safe_group_names[2], "-", safe_group_names[1]),
      levels = design
    )
    fit <- limma::contrasts.fit(fit, contrast)
    fit <- limma::eBayes(fit)

    limma_table <- limma::topTable(
      fit,
      number = Inf,
      adjust.method = input$correction_method,
      sort.by = "none"
    )

    control_mean <- rowMeans(imputed_data[, control_cols, drop = FALSE])
    test_mean <- rowMeans(imputed_data[, test_cols, drop = FALSE])
    log2_fold_change <- test_mean - control_mean

    differential_table <- data.frame(
      filtered_metadata,
      Control.Group = control_group,
      Test.Group = test_group,
      Mean.Control.Log2 = control_mean,
      Mean.Test.Log2 = test_mean,
      log2FC = log2_fold_change,
      P.Value = limma_table$P.Value,
      Adjusted.P.Value = limma_table$adj.P.Val,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    # ----------------------------------------------------------------------
    # STEP 3H. CLASSIFY SIGNIFICANT PROTEINS
    # WHAT THIS DOES: Uses the adjusted p-value cutoff for significance and the
    # sign of log2FC only to determine which selected group is higher.
    # No absolute fold-change cutoff is applied.
    # CHANGE HERE: The adjusted p-value cutoff is set in STEP 1.
    # ----------------------------------------------------------------------
    not_significant <- "Not significant"
    higher_test <- paste0("Higher in ", test_group)
    higher_control <- paste0("Higher in ", control_group)

    differential_table$Significance <- not_significant
    differential_table$Significance[which(
      !is.na(differential_table$Adjusted.P.Value) &
      differential_table$Adjusted.P.Value < input$p_threshold &
        differential_table$log2FC > 0
    )] <- higher_test
    differential_table$Significance[which(
      !is.na(differential_table$Adjusted.P.Value) &
      differential_table$Adjusted.P.Value < input$p_threshold &
        differential_table$log2FC < 0
    )] <- higher_control

    plot_p_values <- pmax(
      differential_table$Adjusted.P.Value,
      .Machine$double.xmin
    )
    differential_table$Minus.Log10.Adjusted.P <- -log10(plot_p_values)

    higher_in_test <- differential_table[
      differential_table$Significance == higher_test,
      ,
      drop = FALSE
    ]
    higher_in_control <- differential_table[
      differential_table$Significance == higher_control,
      ,
      drop = FALSE
    ]

    # ----------------------------------------------------------------------
    # STEP 3I. CREATE THE VOLCANO PLOT
    # WHAT THIS DOES: Colours significant proteins and labels the top five per side.
    # CHANGE HERE: Edit colours, label count, sizes, or theme in this block.
    # ----------------------------------------------------------------------
    top_labels <- rbind(
      head(higher_in_test[order(higher_in_test$Adjusted.P.Value), , drop = FALSE], 5),
      head(higher_in_control[order(higher_in_control$Adjusted.P.Value), , drop = FALSE], 5)
    )

    significance_levels <- c(not_significant, higher_test, higher_control)
    differential_table$Significance <- factor(
      differential_table$Significance,
      levels = significance_levels
    )

    colour_values <- setNames(
      c("#6C757D", "#2E8B57", "#C0392B"),
      significance_levels
    )

    volcano <- ggplot(
      differential_table,
      aes(x = log2FC, y = Minus.Log10.Adjusted.P, colour = Significance)
    ) +
      geom_point(alpha = 0.78, size = 1.8, na.rm = TRUE) +
      geom_hline(
        yintercept = -log10(input$p_threshold),
        linetype = "dashed",
        colour = "#4C3E50"
      ) +
      ggrepel::geom_text_repel(
        data = top_labels,
        aes(label = Protein.ID),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.4,
        point.padding = 0.3,
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      scale_colour_manual(values = colour_values, drop = FALSE) +
      theme_bw() +
      labs(
        title = paste("Volcano Plot:", test_group, "vs", control_group),
        subtitle = paste0(
          "Positive log2FC = higher in ", test_group,
          "; negative log2FC = higher in ", control_group
        ),
        x = paste0("log2 fold change (", test_group, " - ", control_group, ")"),
        y = "-log10 adjusted p-value",
        colour = "Result"
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom"
      )

    comparison_name <- safe_file_name(paste0(test_group, "_vs_", control_group))
    imputed_output <- cbind(filtered_metadata, imputed_data)

    parameter_table <- data.frame(
      Parameter = c(
        "Input file", "Metric", "Control group", "Test group",
        "Quantitative data transformation",
        "Pairwise valid-value percentage", "Pairwise valid-value rule",
        "Imputation method", "P-value correction",
        "Adjusted p-value threshold", "KNN random seed fixed in code"
      ),
      Value = c(
        input$uploaded_file$name, metric, control_group, test_group,
        transformation_label,
        50, "At least 50% valid in each selected group", "KNN within each group",
        input$correction_method,
        input$p_threshold, knn_random_seed
      ),
      stringsAsFactors = FALSE
    )

    list(
      results = differential_table,
      qc_filtered_input = qc_filtered_input,
      selected_groups_before_filter = selected_groups_before_filter,
      pairwise_filter_table = pairwise_filter_table,
      imputed_output = imputed_output,
      higher_in_test = higher_in_test,
      higher_in_control = higher_in_control,
      volcano = volcano,
      column_design = column_design,
      parameters = parameter_table,
      imputation_log = imputation_log,
      imputation_summary = imputation_summary,
      metric_columns = metric_cols,
      metric = metric,
      control_group = control_group,
      test_group = test_group,
      comparison_name = comparison_name,
      transformation_label = transformation_label,
      correction_method = input$correction_method,
      p_threshold = input$p_threshold,
      total_initial = total_initial,
      total_control_detected = total_control_detected,
      total_test_detected = total_test_detected,
      total_after_valid_filter = total_after_valid_filter,
      control_needed = control_needed,
      test_needed = test_needed
    )
  })

  # Keep result outputs silent during app startup. The analysis object does not
  # exist until a file has been uploaded and the Run button has been clicked.
  displayed_results <- reactive({
    req(input$uploaded_file)
    req(input$run_analysis > 0)
    analysis_results()
  })

  # ------------------------------------------------------------------------
  # STEP 4A. PRINT THE ANALYSIS SUMMARY
  # WHAT THIS DOES: Reports filtering, design, significance, and imputation counts.
  # CHANGE HERE: Add/remove cat() lines to change the on-screen report.
  # ------------------------------------------------------------------------
  output$analysis_summary <- renderPrint({
    res <- displayed_results()

    cat("===============================================================\n")
    cat("             DIFFERENTIAL PROTEIN ANALYSIS SUMMARY             \n")
    cat("===============================================================\n")
    cat("Comparison:", res$test_group, "vs", res$control_group, "\n")
    cat("Metric analysed:", res$metric, "\n")
    cat("Quantitative data:", res$transformation_label, "\n")
    cat("Detected quantitative columns:", length(res$metric_columns), "\n")
    cat("  ", paste(res$metric_columns, collapse = " | "), "\n\n")

    cat("PAIRWISE ROW SELECTION\n")
    cat("Rows received from the QC-filtered input:", res$total_initial, "\n")
    cat("No additional QC or contaminant row filtering was performed.\n")
    cat("Proteins with at least one valid value in", res$control_group, ":",
        res$total_control_detected, "\n")
    cat("Proteins with at least one valid value in", res$test_group, ":",
        res$total_test_detected, "\n")
    cat("Valid-value rule: at least 50% of replicates in EACH selected group.\n")
    cat("Required valid replicates:", res$control_group, "=",
        res$control_needed, ";", res$test_group, "=", res$test_needed, "\n")
    cat("Proteins retained for testing:", res$total_after_valid_filter, "\n\n")

    cat("STATISTICAL RESULTS\n")
    cat("P-value correction:", res$correction_method, "\n")
    cat("Significance cutoff: adjusted p-value <", res$p_threshold, "\n")
    cat("No absolute log2 fold-change cutoff was applied.\n")
    cat("Higher in", res$test_group, ":", nrow(res$higher_in_test), "\n")
    cat("Higher in", res$control_group, ":", nrow(res$higher_in_control), "\n\n")

    cat("IMPUTATION LOG\n")
    cat(paste(res$imputation_log, collapse = "\n"), "\n")
    cat("===============================================================\n")
  })

  # ------------------------------------------------------------------------
  # STEP 4B. DISPLAY THE VOLCANO PLOT
  # ------------------------------------------------------------------------
  output$volcano_plot <- renderPlot({
    res <- displayed_results()
    print(res$volcano)
  })

  # ------------------------------------------------------------------------
  # STEP 4C. DISPLAY THE TOP 20 PROTEINS BY ADJUSTED P-VALUE
  # CHANGE HERE: Replace 20 below to show a different number of rows.
  # ------------------------------------------------------------------------
  output$results_table <- renderTable({
    res <- displayed_results()$results
    res <- res[order(res$Adjusted.P.Value, -abs(res$log2FC)), , drop = FALSE]
    display_columns <- c(
      "Protein", "Protein.ID", "Description", "log2FC",
      "P.Value", "Adjusted.P.Value", "Significance"
    )
    head(res[, display_columns, drop = FALSE], 20)
  }, digits = 4, striped = TRUE, bordered = TRUE, spacing = "s")

  # ------------------------------------------------------------------------
  # STEP 4D. DOWNLOAD ALL INTERMEDIATE TABLES, RESULTS, AND PLOTS
  # WHAT THIS DOES: Builds a ZIP with intermediate_csv, results_csv, and plots
  # subfolders so every stage of the pairwise analysis can be inspected.
  # CHANGE HERE: Add or remove output files in this block.
  # ------------------------------------------------------------------------
  output$download_results_zip <- downloadHandler(
    filename = function() {
      res <- analysis_results()
      paste0("differential_proteins_", res$comparison_name, "_", Sys.Date(), ".zip")
    },
    content = function(file) {
      res <- analysis_results()
      output_directory <- tempfile("differential_proteins_")
      dir.create(output_directory, recursive = TRUE)
      on.exit(unlink(output_directory, recursive = TRUE, force = TRUE), add = TRUE)
      intermediate_directory <- file.path(output_directory, "intermediate_csv")
      results_directory <- file.path(output_directory, "results_csv")
      plots_directory <- file.path(output_directory, "plots")
      dir.create(intermediate_directory)
      dir.create(results_directory)
      dir.create(plots_directory)

      files <- c(
        qc_input = file.path(
          intermediate_directory,
          "01_qc_filtered_input_as_uploaded.csv"
        ),
        selected_groups = file.path(
          intermediate_directory,
          paste0("02_selected_groups_before_50_percent_filter_", res$comparison_name, ".csv")
        ),
        pairwise_filtered = file.path(
          intermediate_directory,
          paste0("03_rows_passing_50_percent_each_group_", res$comparison_name, ".csv")
        ),
        imputed_data = file.path(
          intermediate_directory,
          paste0("04_knn_imputed_analysis_data_", res$comparison_name, ".csv")
        ),
        imputation_summary = file.path(
          intermediate_directory,
          "05_knn_imputation_summary.csv"
        ),
        column_design = file.path(
          intermediate_directory,
          "06_column_to_experiment_mapping.csv"
        ),
        parameters = file.path(
          intermediate_directory,
          "07_analysis_parameters.csv"
        ),
        full_results = file.path(
          results_directory,
          paste0("differential_proteins_", res$comparison_name, ".csv")
        ),
        higher_test = file.path(
          results_directory,
          paste0("proteins_higher_in_", safe_file_name(res$test_group), ".csv")
        ),
        higher_control = file.path(
          results_directory,
          paste0("proteins_higher_in_", safe_file_name(res$control_group), ".csv")
        ),
        volcano_png = file.path(
          plots_directory,
          paste0("volcano_", res$comparison_name, ".png")
        ),
        volcano_pdf = file.path(
          plots_directory,
          paste0("volcano_", res$comparison_name, ".pdf")
        )
      )

      write.csv(res$qc_filtered_input, files[["qc_input"]], row.names = FALSE)
      write.csv(
        res$selected_groups_before_filter,
        files[["selected_groups"]],
        row.names = FALSE
      )
      write.csv(
        res$pairwise_filter_table,
        files[["pairwise_filtered"]],
        row.names = FALSE
      )
      write.csv(res$imputed_output, files[["imputed_data"]], row.names = FALSE)
      write.csv(
        res$imputation_summary,
        files[["imputation_summary"]],
        row.names = FALSE
      )
      write.csv(res$column_design, files[["column_design"]], row.names = FALSE)
      write.csv(res$parameters, files[["parameters"]], row.names = FALSE)
      write.csv(res$results, files[["full_results"]], row.names = FALSE)
      write.csv(res$higher_in_test, files[["higher_test"]], row.names = FALSE)
      write.csv(res$higher_in_control, files[["higher_control"]], row.names = FALSE)

      ggsave(
        filename = files[["volcano_png"]],
        plot = res$volcano,
        width = 8,
        height = 6,
        dpi = 300
      )
      ggsave(
        filename = files[["volcano_pdf"]],
        plot = res$volcano,
        width = 8,
        height = 6
      )

      zip_contents <- list.files(
        output_directory,
        recursive = TRUE,
        full.names = FALSE
      )
      zip::zipr(zipfile = file, files = zip_contents, root = output_directory)
    }
  )

  })
}
