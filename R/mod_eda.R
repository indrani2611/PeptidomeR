# mod_eda: converted from eda.R

mod_eda_ui <- function(id) {
  ns <- NS(id)

  tagList(
    h2("Peptidomics Summary Plots"),
    layout_sidebar(
  
  sidebar = sidebar(
    width = 400,
    title = "Configuration Panel",
    
    card(
      card_header("1. File Input"),
      fileInput(ns("uploaded_file"), "Choose Peptide Data File (.csv/.tsv/.txt):", 
                accept = c(".csv", ".tsv", ".txt"))
    ),
    
    card(
      card_header("2. Metric Selection"),
      selectInput(ns("target_metric"), "Quantitative Columns to Extract:",
                  choices = c("MaxLFQ Intensity", "Intensity", "Spectral Count"))
    ),
    
    card(
      card_header("3. Study design"),
      numericInput(ns("num_conditions"), "Number of Experiments:", value = 3, min = 1),
      uiOutput(ns("condition_setup_ui"))
    ),
    
    actionButton(ns("click_run"), "Run Explorative Data Analysis", class = "btn-success w-100 fw-bold mt-2")
  ),
  
  # ---------------------------------------------------------------------
  # Main Dashboard Viewport
  # ---------------------------------------------------------------------
  card(
    card_header("Statistical Summary"),
    verbatimTextOutput(ns("percentage_summary")),
    card_footer(
      downloadButton(ns("download_pipeline_zip"), "Download All Subsets & Figures (.zip)", class = "btn-primary w-100 fw-bold")
    )
  ),
  
  layout_column_wrap(
    width = 1/2,
    card(
      card_header("Venn Diagram Distribution"),
      plotOutput(ns("plot_venn"), height = "460px")
    ),
    card(
      card_header("Peptide Length Frequency Distribution"),
      plotOutput(ns("plot_length"), height = "460px")
    )
  )
,
      fillable = FALSE,
      fill = FALSE,
      border = FALSE
    )
  )
}

mod_eda_server <- function(id) {
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

  # Dynamically create experiment name and replicate inputs.
  output$condition_setup_ui <- renderUI({
    req(input$num_conditions)

    lapply(seq_len(input$num_conditions), function(i) {
      div(
        class = "p-2 mb-2 border rounded bg-light",
        textInput(
          ns(paste0("cond_name_", i)),
          paste0("Experiment ", i, " Name:"),
          value = ""
        ),
        numericInput(
          ns(paste0("cond_rep_", i)),
          paste0("Experiment ", i, " Replicates:"),
          value = 2,
          min = 1
        )
      )
    })
  })
  
  # Analytical Pipeline
  pipeline_results <- eventReactive(input$click_run, {
    req(raw_dataframe())

    num_conditions <- as.integer(input$num_conditions)

    # Validation to ensure user provided custom experimental labels
    labels <- vapply(seq_len(num_conditions), function(i) {
      value <- input[[paste0("cond_name_", i)]]
      if (is.null(value)) "" else trimws(value)
    }, character(1))

    req_replicates <- vapply(seq_len(num_conditions), function(i) {
      value <- input[[paste0("cond_rep_", i)]]
      if (is.null(value)) NA_real_ else as.numeric(value)
    }, numeric(1))

    if (any(labels == "") || any(is.na(req_replicates))) {
      stop("Configuration Error: Please type text label names for all experimental groups before running.")
    }

    if (anyDuplicated(labels)) {
      stop("Configuration Error: Experiment names must be unique.")
    }
    
    df <- raw_dataframe()
    
    # Clean string-matching for Columns to minimize regex issues
    if (input$target_metric == "MaxLFQ Intensity") {
      matched_cols <- grep("MaxLFQ", colnames(df), value = TRUE)
      matched_cols <- grep("Intensity", matched_cols, value = TRUE)
    } else if (input$target_metric == "Spectral Count") {
      matched_cols <- grep("Spectral", colnames(df), value = TRUE)
    } else {
      matched_cols <- grep("Intensity", colnames(df), value = TRUE)
      matched_cols <- matched_cols[!grepl("MaxLFQ", matched_cols)]
    }
    
    if (sum(req_replicates) != length(matched_cols)) {
      stop(paste0("Shape Error: Total replicates (", sum(req_replicates), 
                  ") does not match columns found for '", input$target_metric, "' (", length(matched_cols), ")."))
    }
    
    # Filter rows containing real signal metrics and strip out contaminant lines
    nonzero_rows <- rowSums(df[, matched_cols, drop = FALSE] > 0, na.rm = TRUE) >= 1
    prot_idx <- grep("Protein", colnames(df), ignore.case = TRUE)
    non_contam <- if (length(prot_idx) > 0) !grepl("^contam", df[[prot_idx[1]]], ignore.case = TRUE) else rep(TRUE, nrow(df))
    
    df_filtered <- df[nonzero_rows & non_contam, , drop = FALSE]
    
    # Locate tracking attributes safely to formulate a definitive identifier key
    find_col <- function(choices, dataset) {
      for (c in choices) {
        if (c %in% colnames(dataset)) return(c)
        dot_c <- gsub(" ", ".", c)
        if (dot_c %in% colnames(dataset)) return(dot_c)
      }
      for (c in choices) {
        idx <- grep(c, colnames(dataset), ignore.case = TRUE)
        if (length(idx) > 0) return(colnames(dataset)[idx[1]])
      }
      return(NULL)
    }
    
    seq_col   <- find_col(c("Peptide Sequence", "Sequence"), df_filtered)
    id_col    <- find_col(c("Protein ID", "Protein"), df_filtered)
    start_col <- find_col(c("Start"), df_filtered)
    end_col   <- find_col(c("End"), df_filtered)
    
    if (!is.null(seq_col) && !is.null(id_col) && !is.null(start_col) && !is.null(end_col)) {
      df_filtered$Peptide.Unique.ID <- paste(df_filtered[[seq_col]], df_filtered[[id_col]], df_filtered[[start_col]], df_filtered[[end_col]], sep = "_")
      df_filtered$Peptide.Length <- nchar(as.character(df_filtered[[seq_col]]))
    } else {
      df_filtered$Peptide.Unique.ID <- as.character(df_filtered[[1]])
      df_filtered$Peptide.Length <- nchar(as.character(df_filtered[[1]]))
    }
    
    # Disassemble rows into explicit experiment lists
    peptide_sets <- list()
    
    idx_tracker <- 1
    for (i in seq_len(num_conditions)) {
      n_rep <- req_replicates[i]
      grp_cols <- matched_cols[idx_tracker:(idx_tracker + n_rep - 1)]
      idx_tracker <- idx_tracker + n_rep
      
      active_rows <- rowSums(df_filtered[, grp_cols, drop = FALSE] > 0, na.rm = TRUE) >= 1
      peptide_sets[[labels[i]]] <- unique(df_filtered$Peptide.Unique.ID[active_rows])
    }

    all_peptides <- unique(unlist(peptide_sets, use.names = FALSE))
    total_unique_peptides <- length(all_peptides)

    # Create ordered labels for every possible experiment combination.
    compartment_levels <- unlist(
      lapply(seq_len(num_conditions), function(group_size) {
        combinations <- combn(labels, group_size, simplify = FALSE)
        vapply(combinations, function(combination) {
          if (length(combination) == num_conditions) {
            "All groups"
          } else if (length(combination) == 1) {
            paste(combination, "only")
          } else {
            paste(combination, collapse = " + ")
          }
        }, character(1))
      }),
      use.names = FALSE
    )

    # Map back definitions for dynamic length groupings
    map_compartment <- function(id) {
      present_in <- labels[vapply(
        peptide_sets,
        function(peptide_set) id %in% peptide_set,
        logical(1)
      )]

      if (length(present_in) == 0) return("None")
      if (length(present_in) == num_conditions) return("All groups")
      if (length(present_in) == 1) return(paste(present_in, "only"))
      paste(present_in, collapse = " + ")
    }
    
    # Build complete subset dataframe mapping (retaining ALL original columns)
    plot_df <- df_filtered[df_filtered$Peptide.Unique.ID %in% all_peptides, ]
    plot_df$Compartment <- factor(
      vapply(
        as.character(plot_df$Peptide.Unique.ID),
        map_compartment,
        character(1)
      ),
      levels = compartment_levels
    )

    unique_compartments <- unique(plot_df[, c("Peptide.Unique.ID", "Compartment")])
    compartment_counts <- table(factor(
      unique_compartments$Compartment,
      levels = compartment_levels
    ))
    compartment_counts <- setNames(
      as.integer(compartment_counts),
      compartment_levels
    )
    
    list(
      peptide_sets = peptide_sets,
      plot_df = plot_df,
      total = total_unique_peptides,
      labels = labels,
      counts = compartment_counts
    )
  })
  
  # Helper plotting configurations for rendering and downloading uniformly
  generate_venn_plot <- function(res) {
    ggVennDiagram(res$peptide_sets, label_alpha = 0.4) +
      scale_fill_gradientn(colors = c("#F7FCF0", "#E0F3DB", "#A8DDB5", "#4EB3D3", "#0868AC")) +
      theme(legend.position = "right") 
  }
  
  generate_length_plot <- function(res) {
    ggplot(res$plot_df, aes(x = Peptide.Length, fill = Compartment)) +
      geom_histogram(bins = 20, color = "black", alpha = 0.8) +
      scale_fill_brewer(palette = "Set2") +
      facet_wrap(~ Compartment, scales = "free_y") +
      theme_bw() +
      labs(title = "Peptide Length Distribution for each group",
           x = "Peptide Length (No. of Amino Acids)", y = "Frequency (Count)") +
      theme(legend.position = "none",
            strip.text = element_text(size = 11, face = "bold"),
            plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  }
  
  # 1. Output Text Summary & Percentage Calculations
  output$percentage_summary <- renderPrint({
    res <- tryCatch({ pipeline_results() }, error = function(e) { e$message })
    if (is.character(res)) {
      cat("Input file", res)
    } else {
      tot <- res$total
      cnt <- res$counts
      pct <- function(x) {
        if (tot == 0) 0 else round((x / tot) * 100, 2)
      }
      
      cat("===============================================================\n")
      cat("       COMPARTMENTAL PERCENTAGE BREAKDOWN CALCULATIONS        \n")
      cat("===============================================================\n")
      cat(sprintf("Total Consolidated Peptides Evaluated: %d\n\n", tot))

      for (i in seq_along(cnt)) {
        cat(sprintf("%d) %s:\n", i, names(cnt)[i]))
        cat(sprintf(
          "   Count: %d | Percentage: %s%%\n\n",
          cnt[i],
          pct(cnt[i])
        ))
      }

      cat("===============================================================\n")
    }
  })
  
  # Comprehensive Pipeline Exporter (.zip containing all CSV subsets + figures)
  output$download_pipeline_zip <- downloadHandler(
    filename = function() {
      paste0("peptidomics_eda_", Sys.Date(), ".zip")
    },
    content = function(file) {
      res <- pipeline_results()
      plot_df <- res$plot_df
      
      # A separate directory avoids filename clashes between Shiny sessions.
      tmp_dir <- tempfile("peptidomer_eda_")
      dir.create(tmp_dir, recursive = TRUE)
      on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
      
      csv_files <- c()
      comp_levels <- levels(plot_df$Compartment)
      
      # Extract individual data compartments containing ALL original parent columns
      for (i in seq_along(comp_levels)) {
        lvl <- comp_levels[i]
        safe_name <- gsub(" \\+ ", "_and_", lvl)
        safe_name <- gsub(" ", "_", safe_name)
        csv_filename <- paste0("subset_compartment_", i, "_", safe_name, ".csv")
        
        subset_isolated <- plot_df[plot_df$Compartment == lvl, ]
        
        # Strip out working unique columns to preserve pure raw/original table column structures
        subset_isolated$Peptide.Unique.ID <- NULL
        subset_isolated$Peptide.Length <- NULL
        subset_isolated$Compartment <- NULL
        
        write.csv(
          subset_isolated,
          file.path(tmp_dir, csv_filename),
          row.names = FALSE
        )
        csv_files <- c(csv_files, csv_filename)
      }
      
      # Export publication-ready images inside the bundle
      png(file.path(tmp_dir, "venn_diagram_distribution.png"),
          width = 8, height = 6, units = "in", res = 300)
      print(generate_venn_plot(res))
      dev.off()
      
      ggsave(file.path(tmp_dir, "peptide_length_distribution.png"),
             plot = generate_length_plot(res), width = 11, height = 7, dpi = 300)
      
      # Compile all items into a unified pipeline zip payload
      all_bundle_items <- c(csv_files, "venn_diagram_distribution.png", "peptide_length_distribution.png")
      zip::zipr(zipfile = file, files = all_bundle_items, root = tmp_dir)
    }
  )
  
  # 2. Output Venn Diagram display viewport
  output$plot_venn <- renderPlot({
    res <- pipeline_results()
    generate_venn_plot(res)
  })
  
  # 3. Output Peptide Length Histograms display viewport
  output$plot_length <- renderPlot({
    res <- pipeline_results()
    req(nrow(res$plot_df) > 0)
    generate_length_plot(res)
  })

  })
}
