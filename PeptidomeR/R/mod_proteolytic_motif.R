# mod_proteolytic_motif: converted from proteolyticmotif.R

mod_proteolytic_motif_ui <- function(id) {
  ns <- NS(id)

  tagList(
    h2("Proteolytic Motif Analysis"),
    layout_sidebar(

  sidebar = sidebar(
    width = 400,
    title = "Configuration Panel",

    card(
      card_header("1. File Input "),
      fileInput(
        inputId = ns("peptide_file"),
        label = "Choose Peptide Data File (.csv/.tsv/.txt):",
        accept = c(".csv", ".tsv", ".txt")
      ),
      uiOutput(ns("motif_source_ui")),
      uiOutput(ns("protein_background_ui")),
      uiOutput(ns("fasta_file_ui"))
    ),

    card(
      card_header("2. Study design"),
      numericInput(
        inputId = ns("n_experiments"),
        label = "Number of Experiments:",
        value = 3,
        min = 1,
        step = 1
      ),
      selectInput(
        inputId = ns("quantity_column_type"),
        label = "Quantitative column type to detect:",
        choices = c(
          "MaxLFQ Intensity" = "MaxLFQ.Intensity",
          "Intensity" = "Intensity",
          "Spectral Count" = "Spectral Count"
        ),
        selected = "MaxLFQ.Intensity"
      ),
      helpText("Columns matching the selected type are assigned to experiments by input order."),
      uiOutput(ns("experiment_setup_ui"))
    ),

    actionButton(
      inputId = ns("run_analysis"),
      label = "Run Proteolytic Motif Analysis",
      class = "btn-success w-100 fw-bold mt-2"
    )
  ),

  card(
    card_header("Input Check"),
    verbatimTextOutput(ns("input_check")),
    card_footer(
      downloadButton(
        outputId = ns("download_results_zip"),
        label = "Download All CSVs & Figures (.zip)",
        class = "btn-primary w-100 fw-bold"
      )
    )
  ),

  layout_column_wrap(
    width = 1/2,
    card(
      card_header(textOutput(ns("experiment_mapping_header"))),
      tableOutput(ns("experiment_mapping"))
    ),
    card(
      card_header(textOutput(ns("valid_value_counts_header"))),
      tableOutput(ns("valid_value_counts"))
    )
  ),

  card(
    card_header("Analysis Status"),
    verbatimTextOutput(ns("analysis_status"))
  ),

  card(
    card_header("Motif Plots"),
    uiOutput(ns("motif_plot_ui"))
  )
,
      fillable = FALSE,
      fill = FALSE,
      border = FALSE
    )
  )
}

mod_proteolytic_motif_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || is.na(x)) y else x
  }

  AA20 <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
            "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  POS8 <- c("P4", "P3", "P2", "P1", "p1", "p2", "p3", "p4")
  DISPLAY_X <- c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")

  read_input_table <- function(path) {
    first_line <- readLines(path, n = 1, warn = FALSE)
    sep <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","

    utils::read.table(
      path,
      header = TRUE,
      sep = sep,
      quote = "\"",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  detect_quantitative_columns <- function(data, column_type = "MaxLFQ.Intensity") {
    column_type <- column_type %||% "MaxLFQ.Intensity"
    column_names <- names(data)

    pattern <- switch(
      column_type,
      "Spectral Count" = "Spectral[[:space:]._-]*Count|SpectralCount",
      "Intensity" = "(^|[[:space:]._-])Intensity($|[[:space:]._-])|^Intensity$",
      "MaxLFQ.Intensity" = "MaxLFQ[[:space:]._-]*Intensity|MaxLFQIntensity",
      "MaxLFQ[[:space:]._-]*Intensity|MaxLFQIntensity"
    )

    matched_columns <- grep(pattern, column_names, ignore.case = TRUE, value = TRUE)

    # Keep plain Intensity separate from MaxLFQ.Intensity and Spectral Count.
    if (identical(column_type, "Intensity")) {
      matched_columns <- matched_columns[
        !grepl("MaxLFQ|Spectral[[:space:]._-]*Count|SpectralCount", matched_columns, ignore.case = TRUE)
      ]
    }

    matched_columns
  }

  safe_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_-]+", "_", trimws(as.character(x)))
    x <- gsub("^_+|_+$", "", x)
    ifelse(nzchar(x), x, "Experiment")
  }

  to_numeric_safe <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "NaN", "NULL", "Inf", "-Inf")] <- NA_character_
    suppressWarnings(as.numeric(gsub(",", "", x, fixed = TRUE)))
  }

  has_valid_numeric_value <- function(x) {
    y <- to_numeric_safe(x)
    is.finite(y) & !is.na(y) & y > 0
  }

  row_has_any_valid_value <- function(data, columns) {
    if (length(columns) == 0 || nrow(data) == 0) {
      return(rep(FALSE, nrow(data)))
    }

    value_matrix <- do.call(cbind, lapply(columns, function(col) {
      has_valid_numeric_value(data[[col]])
    }))

    rowSums(value_matrix, na.rm = TRUE) > 0
  }

  assign_experiment_columns <- function(quantitative_cols, experiment_info) {
    start_col <- 1
    out <- vector("list", nrow(experiment_info))
    names(out) <- experiment_info$ExperimentName

    for (i in seq_len(nrow(experiment_info))) {
      end_col <- start_col + experiment_info$Replicates[i] - 1
      out[[i]] <- quantitative_cols[start_col:end_col]
      out[[i]] <- out[[i]][!is.na(out[[i]])]
      start_col <- end_col + 1
    }

    out
  }

  prepare_peptide_table <- function(peptide_df) {
    protein_col <- names(peptide_df)[tolower(names(peptide_df)) == "protein"][1]
    start_col <- names(peptide_df)[tolower(names(peptide_df)) == "start"][1]
    end_col <- names(peptide_df)[tolower(names(peptide_df)) == "end"][1]
    peptide_col <- names(peptide_df)[grepl("^Peptide[ ._]?Sequence$|Peptide.*Sequence",
                                           names(peptide_df),
                                           ignore.case = TRUE)][1]

    if (is.na(protein_col) || is.na(start_col) || is.na(end_col)) {
      stop("Peptide file must contain Protein, Start, and End columns.")
    }

    peptide_df$Protein <- trimws(as.character(peptide_df[[protein_col]]))
    peptide_df$Start <- as.integer(to_numeric_safe(peptide_df[[start_col]]))
    peptide_df$End <- as.integer(to_numeric_safe(peptide_df[[end_col]]))

    if (!is.na(peptide_col)) {
      peptide_df$Peptide <- trimws(as.character(peptide_df[[peptide_col]]))
    } else {
      peptide_df$Peptide <- NA_character_
    }

    keep <- nzchar(peptide_df$Protein) & !is.na(peptide_df$Start) & !is.na(peptide_df$End)
    peptide_df <- peptide_df[keep, , drop = FALSE]

    non_contaminant <- !grepl("^contam|^CON__|contaminant",
                              peptide_df$Protein,
                              ignore.case = TRUE)
    peptide_df[non_contaminant, , drop = FALSE]
  }

  read_fasta_sequences <- function(path) {
    lines <- readLines(path, warn = FALSE)
    header_idx <- grep("^>", lines)

    if (length(header_idx) == 0) {
      stop("FASTA file has no headers beginning with '>'.")
    }

    protein <- character(length(header_idx))
    sequence <- character(length(header_idx))

    for (i in seq_along(header_idx)) {
      start <- header_idx[i]
      end <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
      header <- sub("^>", "", lines[start])
      protein[i] <- strsplit(header, "\\s+")[[1]][1]
      sequence_lines <- if ((start + 1) <= end) lines[(start + 1):end] else character(0)
      sequence[i] <- toupper(gsub("\\s+", "", paste(sequence_lines, collapse = "")))
    }

    data.frame(Protein = protein, Sequence = sequence, stringsAsFactors = FALSE)
  }

  extract_motifs_from_merged <- function(merged_df) {
    nterm <- character(nrow(merged_df))
    cterm <- character(nrow(merged_df))
    nterm[] <- NA_character_
    cterm[] <- NA_character_

    for (i in seq_len(nrow(merged_df))) {
      protein_sequence <- merged_df$Sequence[i]
      start_pos <- merged_df$Start[i]
      end_pos <- merged_df$End[i]

      if (is.na(protein_sequence) || is.na(start_pos) || is.na(end_pos)) {
        next
      }

      nterm[i] <- substr(protein_sequence, max(1, start_pos - 4), min(nchar(protein_sequence), start_pos + 3))
      cterm[i] <- substr(protein_sequence, max(1, end_pos - 3), min(nchar(protein_sequence), end_pos + 4))
    }

    list(nterm = nterm, cterm = cterm)
  }

  filter_motifs_8mer <- function(merged_df, motifs) {
    out <- data.frame(
      Peptide = merged_df$Peptide,
      motif = motifs,
      stringsAsFactors = FALSE
    )

    out[!is.na(out$motif) & nchar(out$motif) == 8 & !grepl("X", out$motif), , drop = FALSE]
  }

  expand_motif_positions <- function(df) {
    if (nrow(df) == 0) {
      for (pos in POS8) {
        df[[pos]] <- character(0)
      }
      return(df)
    }

    df$P4 <- substr(df$motif, 1, 1)
    df$P3 <- substr(df$motif, 2, 2)
    df$P2 <- substr(df$motif, 3, 3)
    df$P1 <- substr(df$motif, 4, 4)
    df$p1 <- substr(df$motif, 5, 5)
    df$p2 <- substr(df$motif, 6, 6)
    df$p3 <- substr(df$motif, 7, 7)
    df$p4 <- substr(df$motif, 8, 8)
    df
  }

  make_dummy_position_table <- function(df8) {
    out <- df8[, c("motif", POS8), drop = FALSE]

    for (pos in POS8) {
      for (aa in AA20) {
        out[[paste0(pos, "_", aa)]] <- as.integer(out[[pos]] == aa)
      }
    }

    out
  }

  protein_background_and_aa_frequency <- function(
    protein_file_path,
    fasta_file_path,
    background_selection = "top_n",
    number_top_proteins = 100
  ) {
    protein_data <- read_input_table(protein_file_path)
    if (!"Combined Total Peptides" %in% names(protein_data) || !"Protein" %in% names(protein_data)) {
      stop("Protein file must contain Protein and Combined Total Peptides columns.")
    }

    protein_data$`Combined Total Peptides` <- to_numeric_safe(protein_data$`Combined Total Peptides`)
    protein_data <- protein_data[order(protein_data$`Combined Total Peptides`, decreasing = TRUE), , drop = FALSE]

    if (identical(background_selection, "all")) {
      selected_proteins <- protein_data
    } else {
      number_top_proteins <- max(1L, as.integer(number_top_proteins %||% 100))
      selected_proteins <- head(protein_data, number_top_proteins)
    }

    seq_tbl <- read_fasta_sequences(fasta_file_path)
    selected_sequences <- seq_tbl[
      match(selected_proteins$Protein, seq_tbl$Protein),
      ,
      drop = FALSE
    ]
    selected_sequences <- selected_sequences[
      !is.na(selected_sequences$Sequence),
      ,
      drop = FALSE
    ]

    all_seq <- paste(selected_sequences$Sequence, collapse = "")
    counts <- vapply(AA20, function(aa) {
      lengths(regmatches(all_seq, gregexpr(aa, all_seq, fixed = TRUE)))
    }, numeric(1))

    if (sum(counts) == 0) {
      stop("No selected protein sequences could be matched to the FASTA file.")
    }

    frequency <- data.frame(
      AminoAcid = names(counts),
      Count = as.numeric(counts),
      Frequency = as.numeric(counts / sum(counts)),
      stringsAsFactors = FALSE
    )

    list(
      selected_proteins = selected_proteins,
      amino_acid_frequency = frequency,
      probability = setNames(frequency$Frequency, frequency$AminoAcid)
    )
  }

  calculate_position_frequencies <- function(dummy_df) {
    count_matrix <- matrix(
      0,
      nrow = length(AA20),
      ncol = length(POS8),
      dimnames = list(AA20, POS8)
    )

    if (nrow(dummy_df) > 0) {
      for (pos in POS8) {
        for (aa in AA20) {
          column_name <- paste0(pos, "_", aa)
          if (column_name %in% names(dummy_df)) {
            count_matrix[aa, pos] <- sum(dummy_df[[column_name]] == 1, na.rm = TRUE)
          }
        }
      }
    }

    totals <- colSums(count_matrix, na.rm = TRUE)
    frequency_matrix <- sweep(
      count_matrix,
      2,
      ifelse(totals > 0, totals, 1),
      "/"
    )

    list(counts = count_matrix, frequencies = frequency_matrix)
  }

  sum_success_by_position <- function(dummy_df) {
    out <- setNames(as.list(rep(0, length(POS8) * length(AA20))),
                    as.vector(outer(POS8, AA20, paste, sep = "_")))

    if (nrow(dummy_df) > 0) {
      for (nm in names(out)) {
        if (nm %in% names(dummy_df)) {
          out[[nm]] <- sum(dummy_df[[nm]] == 1, na.rm = TRUE)
        }
      }
    }

    as.data.frame(out, check.names = FALSE)
  }

  create_position_frames <- function(sum_row, sample_size) {
    frames <- list()

    for (pos in POS8) {
      frames[[pos]] <- data.frame(
        aa = AA20,
        success = vapply(AA20, function(aa) {
          col <- paste0(pos, "_", aa)
          if (col %in% names(sum_row)) as.numeric(sum_row[[col]][1]) else 0
        }, numeric(1)),
        sample_size = sample_size,
        stringsAsFactors = FALSE
      )
    }

    frames
  }

  calculate_binomial_logodds <- function(frames, prob_aa) {
    mat <- matrix(NA_real_, nrow = length(AA20), ncol = length(POS8),
                  dimnames = list(AA20, POS8))

    for (pos in POS8) {
      df <- frames[[pos]]
      df$prob_aa <- as.numeric(prob_aa[df$aa])

      for (i in seq_len(nrow(df))) {
        s <- as.numeric(df$success[i])
        n <- as.numeric(df$sample_size[i])
        p <- as.numeric(df$prob_aa[i])

        if (is.na(s) || is.na(n) || is.na(p) || n <= 0) {
          next
        }

        db1 <- sum(stats::dbinom(s:n, n, p))
        db2 <- if (s == 0) {
          stats::dbinom(0, n, p)
        } else {
          sum(stats::dbinom(1:s, n, p))
        }

        if (db1 > 0 && db2 > 0) {
          mat[i, pos] <- -log(db1 / db2)
        }
      }
    }

    mat
  }

  make_logo_plot <- function(score_mat, motif_source = "protein_background", show_y_title = TRUE) {
    m <- as.matrix(score_mat)
    mode(m) <- "numeric"
    m[!is.finite(m)] <- NA_real_
    colnames(m) <- DISPLAY_X

    if (identical(motif_source, "frequency_counts")) {
      m[is.na(m)] <- 0

      return(
        ggseqlogo(m, method = "custom", seq_type = "aa") +
          theme(axis.line = element_line(colour = "black", linewidth = 1)) +
          scale_x_discrete(name = NULL, limits = colnames(m)) +
          scale_y_continuous(
            name = if (show_y_title) "Amino acid frequency" else NULL,
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.2),
            expand = expansion(mult = c(0, 0.02))
          ) +
          theme(
            axis.title.x = element_blank(),
            axis.title.y = if (show_y_title) element_text(face = "bold", size = 18) else element_blank(),
            axis.text.x = element_text(face = "bold", size = 16),
            axis.text.y = element_text(face = "bold", size = 16)
          )
      )
    }

    pos_stack <- colSums(pmax(m, 0), na.rm = TRUE)
    neg_stack <- colSums(pmax(-m, 0), na.rm = TRUE)
    max_abs_stack <- max(c(pos_stack, neg_stack), na.rm = TRUE)

    if (!is.finite(max_abs_stack) || max_abs_stack == 0) {
      max_abs_stack <- 50
    }

    step <- if (max_abs_stack > 400) 200 else 50
    lim <- max(step, ceiling(max_abs_stack / step) * step)
    x_limits <- colnames(m)

    ggseqlogo(m, method = "custom", seq_type = "aa") +
      theme(axis.line = element_line(colour = "black", linewidth = 1)) +
      scale_x_discrete(name = NULL, limits = x_limits) +
      scale_y_continuous(
        name = if (show_y_title) "Binomial log-odds" else NULL,
        limits = c(-lim, lim),
        breaks = seq(-lim, lim, by = step)
      ) +
      geom_hline(yintercept = 0, colour = "red", linewidth = 1) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = if (show_y_title) element_text(face = "bold", size = 18) else element_blank(),
        axis.text.x = element_text(face = "bold", size = 16),
        axis.text.y = element_text(face = "bold", size = 16)
      )
  }

  run_one_dataset <- function(
    peptide_df,
    label,
    seq_tbl,
    motif_source,
    outdir,
    prob_aa = NULL
  ) {
    safe_label <- safe_name(label)
    merged <- merge(
      peptide_df[, c("Protein", "Start", "End", "Peptide"), drop = FALSE],
      seq_tbl,
      by = "Protein",
      all.x = TRUE,
      sort = FALSE
    )

    write.csv(merged, file.path(outdir, paste0("merged_", safe_label, ".csv")), row.names = FALSE)

    motifs <- extract_motifs_from_merged(merged)
    n8 <- filter_motifs_8mer(merged, motifs$nterm)
    c8 <- filter_motifs_8mer(merged, motifs$cterm)

    n_dummy <- make_dummy_position_table(expand_motif_positions(n8))
    c_dummy <- make_dummy_position_table(expand_motif_positions(c8))

    write.csv(n8, file.path(outdir, paste0("ntermseq_", safe_label, ".csv")), row.names = FALSE)
    write.csv(c8, file.path(outdir, paste0("ctermseq_", safe_label, ".csv")), row.names = FALSE)
    write.csv(n_dummy, file.path(outdir, paste0("nterm_dummy_", safe_label, ".csv")), row.names = FALSE)
    write.csv(c_dummy, file.path(outdir, paste0("cterm_dummy_", safe_label, ".csv")), row.names = FALSE)

    if (identical(motif_source, "frequency_counts")) {
      n_frequency <- calculate_position_frequencies(n_dummy)
      c_frequency <- calculate_position_frequencies(c_dummy)

      score_n <- n_frequency$frequencies
      score_c <- c_frequency$frequencies

      write.csv(
        n_frequency$counts,
        file.path(outdir, paste0("amino_acid_counts_n_", safe_label, ".csv"))
      )
      write.csv(
        c_frequency$counts,
        file.path(outdir, paste0("amino_acid_counts_c_", safe_label, ".csv"))
      )
      write.csv(
        score_n,
        file.path(outdir, paste0("amino_acid_frequencies_n_", safe_label, ".csv"))
      )
      write.csv(
        score_c,
        file.path(outdir, paste0("amino_acid_frequencies_c_", safe_label, ".csv"))
      )
    } else {
      n_dist <- n_dummy[!duplicated(n_dummy$motif), , drop = FALSE]
      c_dist <- c_dummy[!duplicated(c_dummy$motif), , drop = FALSE]

      score_n <- calculate_binomial_logodds(
        create_position_frames(sum_success_by_position(n_dist), nrow(n_dist)),
        prob_aa
      )
      score_c <- calculate_binomial_logodds(
        create_position_frames(sum_success_by_position(c_dist), nrow(c_dist)),
        prob_aa
      )

      write.csv(score_n, file.path(outdir, paste0("ratio_n_", safe_label, ".csv")))
      write.csv(score_c, file.path(outdir, paste0("ratio_c_", safe_label, ".csv")))
    }

    p_n <- make_logo_plot(score_n, motif_source = motif_source, show_y_title = TRUE)
    p_c <- make_logo_plot(score_c, motif_source = motif_source, show_y_title = TRUE)

    ggsave(file.path(outdir, paste0("Nterm_", safe_label, ".tiff")), p_n,
           width = 6, height = 5, units = "in", dpi = 500)
    ggsave(file.path(outdir, paste0("Cterm_", safe_label, ".tiff")), p_c,
           width = 6, height = 5, units = "in", dpi = 500)

    list(
      label = label,
      score_n = score_n,
      score_c = score_c,
      motif_source = motif_source,
      plot_n = p_n,
      plot_c = p_c
    )
  }

  results_dir <- reactiveVal(NULL)
  plot_results <- reactiveVal(list())
  analysis_message <- reactiveVal("Awaiting analysis.")

  session$onSessionEnded(function() {
    directory <- results_dir()
    if (!is.null(directory) && dir.exists(directory)) {
      unlink(directory, recursive = TRUE, force = TRUE)
    }
  })

  peptide_data <- reactive({
    req(input$peptide_file)
    read_input_table(input$peptide_file$datapath)
  })

  selected_motif_source <- reactive({
    input$motif_source %||% "protein_background"
  })

  output$motif_source_ui <- renderUI({
    req(input$peptide_file)

    radioButtons(
      inputId = ns("motif_source"),
      label = "Extract proteolytic motif from:",
      choices = c(
        "Using top proteins as background" = "protein_background",
        "Frequency of counts" = "frequency_counts"
      ),
      selected = "protein_background"
    )
  })

  output$protein_background_ui <- renderUI({
    req(input$peptide_file)

    if (!identical(selected_motif_source(), "protein_background")) {
      return(NULL)
    }

    tagList(
      radioButtons(
        inputId = ns("background_selection"),
        label = "Protein background:",
        choices = c(
          "Select no. of top proteins" = "top_n",
          "Select all proteins" = "all"
        ),
        selected = "top_n"
      ),
      uiOutput(ns("top_protein_number_ui")),
      fileInput(
        inputId = ns("protein_file"),
        label = "Choose Protein Summary File (.csv/.tsv/.txt):",
        accept = c(".csv", ".tsv", ".txt")
      )
    )
  })

  output$top_protein_number_ui <- renderUI({
    req(input$peptide_file)

    if (!identical(selected_motif_source(), "protein_background") ||
        !identical(input$background_selection %||% "top_n", "top_n")) {
      return(NULL)
    }

    numericInput(
      inputId = ns("number_top_proteins"),
      label = "Number of top proteins:",
      value = 100,
      min = 1,
      step = 1
    )
  })

  output$fasta_file_ui <- renderUI({
    req(input$peptide_file)

    fileInput(
      inputId = ns("fasta_file"),
      label = "Choose FASTA File (.fasta):",
      accept = c(".fasta", ".fa", ".faa", ".txt")
    )
  })

  selected_quantity_type <- reactive({
    input$quantity_column_type %||% "MaxLFQ.Intensity"
  })

  quantitative_columns <- reactive({
    detect_quantitative_columns(peptide_data(), selected_quantity_type())
  })

  output$experiment_setup_ui <- renderUI({
    n_experiments <- input$n_experiments
    req(n_experiments)

    experiment_controls <- lapply(seq_len(n_experiments), function(i) {
      div(
        class = "p-2 mb-2 border rounded bg-light",
        textInput(
          inputId = ns(paste0("experiment_name_", i)),
          label = paste0("Experiment ", i, " Name:"),
          value = paste("Experiment", i)
        ),
        numericInput(
          inputId = ns(paste0("experiment_replicates_", i)),
          label = paste0("Experiment ", i, " Replicates:"),
          value = 2,
          min = 1,
          step = 1
        )
      )
    })

    tagList(experiment_controls)
  })

  experiment_info <- reactive({
    n_experiments <- input$n_experiments
    req(n_experiments)

    data.frame(
      ExperimentIndex = seq_len(n_experiments),
      ExperimentName = vapply(
        seq_len(n_experiments),
        function(i) input[[paste0("experiment_name_", i)]] %||% paste("Experiment", i),
        character(1)
      ),
      Replicates = vapply(
        seq_len(n_experiments),
        function(i) input[[paste0("experiment_replicates_", i)]] %||% 1,
        numeric(1)
      ),
      stringsAsFactors = FALSE
    )
  })

  output$input_check <- renderPrint({
    if (is.null(input$peptide_file)) {
     
      return(invisible(NULL))
    }

    cols <- quantitative_columns()
    info <- experiment_info()

    cat("Peptide file:", input$peptide_file$name, "\n")
    cat(
      "Motif calculation:",
      if (identical(selected_motif_source(), "frequency_counts")) {
        "Frequency of counts"
      } else {
        "Protein background with binomial log-odds"
      },
      "\n"
    )
    cat("Selected quantitative column type:", selected_quantity_type(), "\n")
    cat("Detected selected quantitative columns:", length(cols), "\n")
    if (length(cols) > 0) {
      cat(paste0("  ", seq_along(cols), ". ", cols, collapse = "\n"), "\n")
    }
    cat("\nTotal replicates requested:", sum(info$Replicates), "\n")
  })

  output$experiment_mapping_header <- renderText({
    paste("Experiment to", selected_quantity_type(), "Column Mapping")
  })

  output$valid_value_counts_header <- renderText({
    paste("Rows with Valid Numerical", selected_quantity_type(), "Values")
  })

  output$experiment_mapping <- renderTable({
    req(input$peptide_file)

    cols <- quantitative_columns()
    info <- experiment_info()

    start_col <- 1
    rows <- lapply(seq_len(nrow(info)), function(i) {
      n_reps <- info$Replicates[i]
      end_col <- start_col + n_reps - 1
      assigned_cols <- cols[start_col:end_col]

      out <- data.frame(
        Experiment = info$ExperimentName[i],
        Replicates = n_reps,
        SelectedColumnsByOrder = paste(assigned_cols[!is.na(assigned_cols)], collapse = "; "),
        stringsAsFactors = FALSE
      )

      start_col <<- end_col + 1
      out
    })

    do.call(rbind, rows)
  })

  valid_value_counts <- reactive({
    req(input$peptide_file)

    peptide_df <- prepare_peptide_table(peptide_data())
    cols <- quantitative_columns()
    info <- experiment_info()
    exp_cols <- assign_experiment_columns(cols, info)

    rows <- lapply(seq_len(nrow(info)), function(i) {
      exp_name <- info$ExperimentName[i]
      assigned_cols <- exp_cols[[i]]
      valid_mask <- row_has_any_valid_value(peptide_df, assigned_cols)

      data.frame(
        Experiment = exp_name,
        Replicates = info$Replicates[i],
        SelectedColumns = paste(assigned_cols, collapse = "; "),
        RowsWithAnyValidNumericValue = sum(valid_mask),
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)
  })

  output$valid_value_counts <- renderTable({
    valid_value_counts()
  })

  observeEvent(input$run_analysis, {
    req(input$peptide_file)

    motif_source <- selected_motif_source()

    if (is.null(input$fasta_file)) {
      analysis_message("Cannot run yet. Please upload the FASTA file.")
      return(invisible(NULL))
    }

    if (identical(motif_source, "protein_background") && is.null(input$protein_file)) {
      analysis_message("Cannot run yet. Please upload the protein summary file for the selected protein-background method.")
      return(invisible(NULL))
    }

    cols <- quantitative_columns()
    info <- experiment_info()
    total_replicates <- sum(info$Replicates)

    if (total_replicates != length(cols)) {
      analysis_message(paste(
        "Cannot run yet.",
        paste("Total replicates entered:", total_replicates),
        paste("Detected selected quantitative columns:", length(cols)),
        "These numbers should match so experiments can be assigned by selected-column order.",
        sep = "\n"
      ))
      return(invisible(NULL))
    }

    required_packages <- c("ggplot2", "ggseqlogo", "patchwork")
    missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_packages) > 0) {
      analysis_message(paste(
        "Cannot run motif plotting because these R packages are missing:",
        paste(missing_packages, collapse = ", "),
        "Install them first, then run again.",
        sep = "\n"
      ))
      return(invisible(NULL))
    }

    previous_outdir <- results_dir()
    if (!is.null(previous_outdir) && dir.exists(previous_outdir)) {
      unlink(previous_outdir, recursive = TRUE, force = TRUE)
    }
    results_dir(NULL)
    plot_results(list())

    outdir <- tempfile("proteolytic_motif_results_")
    dir.create(outdir, recursive = TRUE)

    tryCatch({
      peptide_df <- prepare_peptide_table(peptide_data())
      exp_cols <- assign_experiment_columns(cols, info)
      counts <- valid_value_counts()

      write.csv(info, file.path(outdir, "experiment_setup.csv"), row.names = FALSE)
      write.csv(data.frame(SelectedQuantitativeColumnType = selected_quantity_type()),
                file.path(outdir, "selected_quantitative_column_type.csv"),
                row.names = FALSE)
      write.csv(data.frame(SelectedColumnsInInputOrder = cols),
                file.path(outdir, "selected_columns_input_order.csv"),
                row.names = FALSE)
      write.csv(counts, file.path(outdir, "valid_value_counts.csv"), row.names = FALSE)

      seq_tbl <- read_fasta_sequences(input$fasta_file$datapath)
      prob_aa <- NULL

      if (identical(motif_source, "protein_background")) {
        background_selection <- input$background_selection %||% "top_n"
        number_top_proteins <- input$number_top_proteins %||% 100
        background_info <- protein_background_and_aa_frequency(
          protein_file_path = input$protein_file$datapath,
          fasta_file_path = input$fasta_file$datapath,
          background_selection = background_selection,
          number_top_proteins = number_top_proteins
        )
        prob_aa <- background_info$probability

        selected_background_name <- if (identical(background_selection, "all")) {
          "all_proteins_background.csv"
        } else {
          paste0("top_", as.integer(number_top_proteins), "_proteins_background.csv")
        }

        write.csv(
          background_info$selected_proteins,
          file.path(outdir, selected_background_name),
          row.names = FALSE
        )
        write.csv(
          background_info$amino_acid_frequency,
          file.path(outdir, "amino_acid_frequencies_protein_background.csv"),
          row.names = FALSE
        )
      }

      write.csv(
        data.frame(
          MotifCalculation = if (identical(motif_source, "frequency_counts")) {
            "Position-wise amino acid frequencies"
          } else {
            "Binomial log-odds using protein background"
          },
          ProteinBackgroundSelection = if (identical(motif_source, "protein_background")) {
            input$background_selection %||% "top_n"
          } else {
            "Not required"
          },
          NumberTopProteins = if (
            identical(motif_source, "protein_background") &&
            identical(input$background_selection %||% "top_n", "top_n")
          ) {
            input$number_top_proteins %||% 100
          } else {
            NA
          },
          stringsAsFactors = FALSE
        ),
        file.path(outdir, "motif_calculation_method.csv"),
        row.names = FALSE
      )

      plotted <- character(0)
      skipped <- character(0)
      results <- list()

      for (i in seq_len(nrow(info))) {
        exp_name <- info$ExperimentName[i]
        valid_mask <- row_has_any_valid_value(peptide_df, exp_cols[[i]])
        exp_peptides <- peptide_df[valid_mask, , drop = FALSE]
        write.csv(exp_peptides,
                  file.path(outdir, paste0("peptides_valid_", safe_name(exp_name), ".csv")),
                  row.names = FALSE)

        if (nrow(exp_peptides) == 0) {
          skipped <- c(skipped, exp_name)
          next
        }

        results[[exp_name]] <- run_one_dataset(
          peptide_df = exp_peptides,
          label = exp_name,
          seq_tbl = seq_tbl,
          motif_source = motif_source,
          outdir = outdir,
          prob_aa = prob_aa
        )
        plotted <- c(plotted, exp_name)
      }

      if (length(results) > 1) {
        p_n <- lapply(seq_along(results), function(i) {
          make_logo_plot(
            results[[i]]$score_n,
            motif_source = results[[i]]$motif_source,
            show_y_title = i == 1
          ) +
            ggtitle(results[[i]]$label) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
        })
        p_c <- lapply(seq_along(results), function(i) {
          make_logo_plot(
            results[[i]]$score_c,
            motif_source = results[[i]]$motif_source,
            show_y_title = i == 1
          ) +
            ggtitle(results[[i]]$label) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
        })

        ggsave(file.path(outdir, "Nterm_all_side_by_side.tiff"),
               wrap_plots(p_n, nrow = 1),
               width = 6 * length(results), height = 5, units = "in", dpi = 500)
        ggsave(file.path(outdir, "Cterm_all_side_by_side.tiff"),
               wrap_plots(p_c, nrow = 1),
               width = 6 * length(results), height = 5, units = "in", dpi = 500)
      }

      results_dir(outdir)
      plot_results(results)
      analysis_message(paste(
        "Analysis finished.",
        paste(
          "Motif calculation:",
          if (identical(motif_source, "frequency_counts")) {
            "Frequency of counts (protein summary file not used)"
          } else {
            "Binomial log-odds using selected protein background"
          }
        ),
        "Rows with valid numerical values:",
        paste(capture.output(print(counts[, c("Experiment", "RowsWithAnyValidNumericValue"), drop = FALSE],
                                   row.names = FALSE)), collapse = "\n"),
        if (length(plotted) > 0) paste("Motif plots created for:", paste(plotted, collapse = ", ")) else "No motif plots were created.",
        if (length(skipped) > 0) paste("Skipped experiments with no valid rows:", paste(skipped, collapse = ", ")) else NULL,
        paste("Output folder:", normalizePath(outdir)),
        sep = "\n"
      ))
    }, error = function(e) {
      unlink(outdir, recursive = TRUE, force = TRUE)
      analysis_message(paste("Analysis failed:", conditionMessage(e)))
    })
  })

  output$analysis_status <- renderPrint({
    cat(analysis_message(), "\n")
  })

  output$motif_plot_ui <- renderUI({
    results <- plot_results()
    if (length(results) == 0) {
      return(tags$p("Plots will appear here after analysis finishes."))
    }

    tagList(lapply(seq_along(results), function(i) {
      safe_label <- safe_name(results[[i]]$label)
      tagList(
        h5(results[[i]]$label),
        fluidRow(
          column(
            6,
            h6("N-terminal motif", class = "fw-bold text-center"),
            plotOutput(ns(paste0("plot_n_", safe_label)), height = "420px")
          ),
          column(
            6,
            h6("C-terminal motif", class = "fw-bold text-center"),
            plotOutput(ns(paste0("plot_c_", safe_label)), height = "420px")
          )
        )
      )
    }))
  })

  observe({
    results <- plot_results()
    if (length(results) == 0) {
      return(invisible(NULL))
    }

    for (result in results) {
      local({
        current <- result
        safe_label <- safe_name(current$label)
        output[[paste0("plot_n_", safe_label)]] <- renderPlot({
          current$plot_n
        })
        output[[paste0("plot_c_", safe_label)]] <- renderPlot({
          current$plot_c
        })
      })
    }
  })

  output$download_results_zip <- downloadHandler(
    filename = function() {
      paste0("proteolytic_motif_results_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(results_dir())
      zip_contents <- list.files(
        results_dir(),
        recursive = TRUE,
        full.names = FALSE
      )
      zip::zipr(zipfile = file, files = zip_contents, root = results_dir())
    }
  )

  })
}
