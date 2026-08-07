# mod_qc: converted from qc.R

mod_qc_ui <- function(id) {
  ns <- NS(id)

  tagList(
    h2("Preprocessing and Quality Control of Peptidomics Data"),
    layout_sidebar(
  
  sidebar = sidebar(
    width = 400,
    title = "Configuration Panel",
    
    card(
      card_header("1. File Input"),
      fileInput(ns("uploaded_file"), "Choose File (.csv, .tsv, .txt)", 
                accept = c(".csv", ".tsv", ".txt"))
    ),
    card(
      card_header("2. Study Design"),
      numericInput(ns("num_conditions"), "Number of Experimental Conditions", value = 3, min = 1),
      uiOutput(ns("condition_setup_ui")),
      selectInput(ns("target_metric"), "Columns to use for Analysis",
                  choices = c(
                    "Intensity" = "(?<!MaxLFQ )Intensity", 
                    "MaxLFQ Intensity" = "MaxLFQ Intensity", 
                    "Spectral Count" = "Spectral Count"
                  ))
    ),
    card(
      card_header("3. Missing Value Filtering"),
      selectInput(ns("filter_rule"), "Filtering Method",
                  choices = c(
                    "Keep features with >= 50% valid values in at least ONE group" = "at_least_one",
                    "Keep features with >= 50% valid values in ALL groups" = "all_groups",
                    "Remove features with > 50% missing values globally" = "global_half",
                    "Complete Rows (No missing values allowed)" = "complete_rows"
                  ))
    ),
    
    card(
      card_header("4. Normalization & Imputation"),
      selectInput(ns("norm_method"), "Normalization Method",
                  choices = c("Log2 Transformation" = "log2", "None" = "none")),
      selectInput(ns("impute_method"), "Imputation Method",
                  choices = c("None" = "none", 
                              "Fixed Constant Value" = "fixed", 
                              "Combined Strategy (KNN + QRILC)" = "combined")),
      # Conditional placeholder for the user-defined imputation number
      uiOutput(ns("fixed_val_ui"))
    ),

    actionButton(ns("click_run"), "Run", class = "btn-success w-100 mt-2")
  ),
  
  # ---------------------------------------------------------------------
  # Main Display Viewport (Stacked Top-To-Bottom)
  # ---------------------------------------------------------------------
  card(
    card_header("Summary"),
    verbatimTextOutput(ns("status_summary"))
  ),
  
  card(
    card_header("Phase 1: Before Normalization"),
    plotOutput(ns("plot_raw"), height = "350px")
  ),
  
  card(
    card_header("Phase 2: After Normalization"),
    plotOutput(ns("plot_norm"), height = "350px"),
    card_footer(
      downloadButton(ns("download_pipeline_zip"), "Download All Intermediate Files & Plots (.zip)", class = "btn-primary")
    )
  )
,
      fillable = FALSE,
      fill = FALSE,
      border = FALSE
    )
  )
}

mod_qc_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

  
  # Auto-detect and validate CSV, TSV, and TXT input files
  raw_dataframe <- reactive({
    req(input$uploaded_file)

    ext <- tolower(tools::file_ext(input$uploaded_file$name))

    if (!ext %in% c("csv", "tsv", "txt")) {
      stop("input file is not in either .csv/.tsv/.txt", call. = FALSE)
    }

    actual_sep <- if (ext == "csv") {
      ","
    } else if (ext == "tsv") {
      "\t"
    } else {
      # A .txt file may be either comma- or tab-separated. Detect the
      # delimiter from its first non-empty line (normally the header).
      header_lines <- readLines(
        input$uploaded_file$datapath,
        n = 20,
        warn = FALSE
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

    df <- tryCatch(
      read.delim(
        file = input$uploaded_file$datapath,
        header = TRUE,
        sep = actual_sep,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        na.strings = c("", "NA", "NaN", "NULL", "null"),
        comment.char = ""
      ),
      error = function(e) {
        stop("input file doesn't match required columns", call. = FALSE)
      }
    )

    cleaned_names <- trimws(colnames(df))
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
        df[quantitative_columns],
        is.numeric,
        logical(1)
      ))

    if (!has_peptide_sequence || !has_numeric_quantitative_column) {
      stop("input file doesn't match required columns", call. = FALSE)
    }

    df
  })

  # Validate immediately after upload, before the user runs the analysis.
  observeEvent(input$uploaded_file, {
    tryCatch(
      raw_dataframe(),
      error = function(e) {
        showNotification(e$message, type = "error", duration = NULL)
      }
    )
  }, ignoreInit = TRUE)
  
  # Dynamic Experimental Groups naming
  output$condition_setup_ui <- renderUI({
    req(input$num_conditions)
    lapply(1:input$num_conditions, function(i) {
      div(
        class = "p-2 mb-2 border rounded bg-light",
        textInput(ns(paste0("cond_name_", i)), paste0("Experiment ", i, " Name:"), value = paste0("Group_", i)),
        numericInput(ns(paste0("cond_rep_", i)), "Number of Replicates:", value = 2, min = 1)
      )
    })
  })
  
  # Dynamic Imputation
  output$fixed_val_ui <- renderUI({
    req(input$impute_method)
    if (input$impute_method == "fixed") {
      numericInput(ns("custom_impute_constant"), "Specify Imputation Value:", value = 0)
    }
  })
  
  processed_output <- eventReactive(input$click_run, {
    req(raw_dataframe())
    
    df <- raw_dataframe()
    metric_pattern <- input$target_metric
    
    matched_columns <- grep(metric_pattern, colnames(df), value = TRUE,perl = TRUE)
    total_found <- length(matched_columns)
    
    total_user_replicates <- 0
    group_assignments <- list()
    pointer <- 1
    
    for (i in 1:input$num_conditions) {
      name <- input[[paste0("cond_name_", i)]]
      reps <- input[[paste0("cond_rep_", i)]]
      total_user_replicates <- total_user_replicates + reps
      
      end_pointer <- pointer + reps - 1
      if (end_pointer <= total_found) {
        group_assignments[[name]] <- matched_columns[pointer:end_pointer]
      }
      pointer <- pointer + reps
    }
    
    if (total_user_replicates != total_found) {
      stop(paste0("Shape Error! Your configuration expects ", total_user_replicates, 
                  " columns, but we found ", total_found, " matching columns."))
    }
    
    initial_rows <- nrow(df)
    
    # Contaminant Filter
    protein_col_idx <- grep("Protein", colnames(df), ignore.case = TRUE)
    if (length(protein_col_idx) > 0) {
      protein_vector <- df[[protein_col_idx[1]]]
      non_contam_rows <- !grepl("^contam|^CON__|contaminant", protein_vector, ignore.case = TRUE)
      df_nocontam <- df[non_contam_rows, , drop = FALSE]
      post_contam_rows <- nrow(df_nocontam)
    } else {
      df_nocontam <- df
      post_contam_rows <- initial_rows
    }
    
    working_matrix <- df_nocontam[, matched_columns, drop = FALSE]
    working_matrix[working_matrix == 0] <- NA
    
    # Missing Value Filtering
    rows_to_keep <- rep(TRUE, nrow(working_matrix))
    if (input$filter_rule == "at_least_one") {
      group_passes <- sapply(names(group_assignments), function(grp) {
        cols <- group_assignments[[grp]]
        rowSums(!is.na(working_matrix[, cols, drop = FALSE])) >= (length(cols) * 0.50)
      })
      rows_to_keep <- rowSums(group_passes) >= 1
    } else if (input$filter_rule == "all_groups") {
      group_passes <- sapply(names(group_assignments), function(grp) {
        cols <- group_assignments[[grp]]
        rowSums(!is.na(working_matrix[, cols, drop = FALSE])) >= (length(cols) * 0.50)
      })
      rows_to_keep <- rowSums(group_passes) == length(group_assignments)
    } else if (input$filter_rule == "global_half") {
      rows_to_keep <- rowSums(is.na(working_matrix)) <= (total_found * 0.50)
    } else if (input$filter_rule == "complete_rows") {
      rows_to_keep <- rowSums(is.na(working_matrix)) == 0
    }
    
    df_filtered <- df_nocontam[rows_to_keep, , drop = FALSE]
    post_filter_rows <- nrow(df_filtered)
    
    calc_matrix <- df_filtered[, matched_columns, drop = FALSE]
    calc_matrix[calc_matrix == 0] <- NA
    
    phase_a_matrix <- calc_matrix
    
    # Normalization Step
    if (input$norm_method == "log2") {
      calc_matrix <- log2(calc_matrix)
    }
    
    df_normalized <- df_filtered
    df_normalized[, matched_columns] <- calc_matrix
    phase_b_matrix <- calc_matrix
    
    # Imputation Step
    if (input$impute_method == "fixed") {
      # Use custom value selected by user, default back to 0 if empty
      replacement_val <- if (!is.null(input$custom_impute_constant)) input$custom_impute_constant else 0
      calc_matrix[is.na(calc_matrix)] <- replacement_val
    } else if (input$impute_method == "combined") {
      for (grp in names(group_assignments)) {
        cols <- group_assignments[[grp]]
        sub_mat <- calc_matrix[, cols, drop = FALSE]
        for (r in 1:nrow(sub_mat)) {
          missing_count <- sum(is.na(sub_mat[r, ]))
          if (missing_count == length(cols)) {
            sub_mat[r, ] <- min(calc_matrix, na.rm = TRUE) * 0.95
          } else if (missing_count > 0) {
            row_mean <- mean(as.numeric(sub_mat[r, ]), na.rm = TRUE)
            sub_mat[r, is.na(sub_mat[r, ])] <- row_mean
          }
        }
        calc_matrix[, cols] <- sub_mat
      }
    }
    
    df_imputed <- df_filtered
    df_imputed[, matched_columns] <- calc_matrix
    
    list(
      initial = initial_rows,
      post_contam = post_contam_rows,
      filtered = post_filter_rows,
      group_mappings = group_assignments,
      metric_used = metric_pattern,
      norm_used = input$norm_method,
      
      file_step1_raw      = df,
      file_step2_nocontam = df_nocontam,
      file_step3_filtered = df_filtered,
      file_step4_norm     = df_normalized,
      file_step5_imputed  = df_imputed,
      
      plot_mat_a = phase_a_matrix,
      plot_mat_b = phase_b_matrix
    )
  })
  
  output$status_summary <- renderPrint({
    res <- tryCatch({ processed_output() }, error = function(e) { e$message })
    if (is.character(res)) {
      cat(" ", res)
    } else if (!is.null(res)) {
      cat("Number of rows in Initial File:         ", res$initial, " rows\n")
      cat("Number of rows Remaining Post-Contaminant Filter:", res$post_contam, " rows\n")
      cat("Dropped Contaminants rows Count:       ", (res$initial - res$post_contam), " rows\n")
      cat("Number of rows after filtering for missingness:", res$filtered, " rows\n")
      cat("Total Rows Deleted:    ", (res$initial - res$filtered), " rows\n")
      cat("Overall % Data Retention:   ", round((res$filtered / res$initial) * 100, 2), "%\n")
    } else {
      cat("Awaiting file uploading pipeline parameters...")
    }
  })
  
  get_clean_long_data <- function(matrix_input) {
    res <- processed_output()
    flat <- stack(as.data.frame(matrix_input))
    flat <- flat[!is.na(flat$values), ]
    flat$Group <- "Unknown"
    for (grp_name in names(res$group_mappings)) {
      flat$Group[flat$ind %in% res$group_mappings[[grp_name]]] <- grp_name
    }
    return(flat)
  }
  
  get_group_colors <- function() {
    res <- processed_output()
    unique_groups <- names(res$group_mappings)
    palette_pool <- c("#34495E", "#16A085", "#E67E22", "#2980B9", "#8E44AD", "#D35400")
    group_colors <- palette_pool[1:length(unique_groups)]
    names(group_colors) <- unique_groups
    return(group_colors)
  }
  
  make_raw_plot <- reactive({
    res <- processed_output()
    plot_df <- get_clean_long_data(res$plot_mat_a)
    
    ggplot(plot_df, aes(x = values, fill = Group)) +
      geom_histogram(bins = 30, color = "black", alpha = 0.8) +
      scale_fill_manual(values = get_group_colors()) +
      facet_wrap(~Group, nrow = 1, scales = "free") + 
      theme_bw() +
      labs(title = paste(res$metric_used, "Distribution"), x = res$metric_used, y = "Count") +
      theme(
        legend.position = "none", 
        text = element_text(face = "bold"), 
        strip.text = element_text(size = 11, face = "bold")
      )
  })
  
  make_norm_plot <- reactive({
    res <- processed_output()
    plot_df <- get_clean_long_data(res$plot_mat_b)
    
    x_title  <- if(res$norm_used == "log2") paste("log2 mean", res$metric_used) else res$metric_used
    p_title  <- if(res$norm_used == "log2") paste("log2", res$metric_used, "distribution") else paste(res$metric_used, "distribution")
    
    ggplot(plot_df, aes(x = values, fill = Group)) +
      geom_histogram(bins = 30, color = "black", alpha = 0.8) +
      scale_fill_manual(values = get_group_colors()) +
      facet_wrap(~Group, nrow = 1, scales = "free") + 
      theme_bw() +
      labs(title = paste(p_title, "(Normalized)"), x = x_title, y = "Count") +
      theme(
        legend.position = "none", 
        text = element_text(face = "bold"), 
        strip.text = element_text(size = 11, face = "bold")
      )
  })
  
  output$plot_raw  <- renderPlot({ make_raw_plot() })
  output$plot_norm <- renderPlot({ make_norm_plot() })
  
  output$download_pipeline_zip <- downloadHandler(
    filename = function() { paste0("peptidomics_qc_results_", Sys.Date(), ".zip") },
    content = function(file) {
      res <- processed_output()
      tmp_dir <- tempfile("peptidomer_qc_")
      dir.create(tmp_dir, recursive = TRUE)
      on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
      
      write.csv(res$file_step1_raw, file.path(tmp_dir, "peptide_uploaded_raw.csv"), row.names = FALSE)
      write.csv(res$file_step2_nocontam, file.path(tmp_dir, "peptide_nocontam.csv"), row.names = FALSE)
      write.csv(res$file_step3_filtered, file.path(tmp_dir, "peptide_post_missingness_filter.csv"), row.names = FALSE)
      write.csv(res$file_step4_norm, file.path(tmp_dir, "peptide_post_normalization.csv"), row.names = FALSE)
      write.csv(res$file_step5_imputed, file.path(tmp_dir, "peptide_final_imputed_output.csv"), row.names = FALSE)
      
      ggsave(file.path(tmp_dir, "plot_phase1_raw_distribution.png"), plot = make_raw_plot(), 
             device = "png", width = 11, height = 5, dpi = 300)
      
      ggsave(file.path(tmp_dir, "plot_phase2_normalized_distribution.png"), plot = make_norm_plot(), 
             device = "png", width = 11, height = 5, dpi = 300)
      
      files_to_zip <- c(
        "peptide_uploaded_raw.csv", 
        "peptide_nocontam.csv",
        "peptide_post_missingness_filter.csv", 
        "peptide_post_normalization.csv",
        "peptide_final_imputed_output.csv",
        "plot_phase1_raw_distribution.png",
        "plot_phase2_normalized_distribution.png"
      )
      
      zip::zipr(zipfile = file, files = files_to_zip, root = tmp_dir)
    }
  )

  })
}
