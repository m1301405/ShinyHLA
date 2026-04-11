## Fastq Upload & Warning
observeEvent(input$fastq_upload_button, {
  tmp_dir <- tempdir()
  bam_dir <- file.path(tmp_dir, "BAM")
  input_dir <- file.path(tmp_dir, "Input")
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8)
  )
  w$show()
  
  if (is.null(input$fastq_gz_upload)) {
    output$file_count <- renderText("No files selected")
    shinyalert("Error", "No files selected.", type = "error")
    w$hide()
    return()
  }
  
  required_files <- 2
  num_files <- nrow(input$fastq_gz_upload)
  
  if (num_files != required_files) {
    shinyalert("Error", paste("Please select exactly", required_files, "files. You selected", num_files, "."), type = "error")
  } else {
    if (dir.exists(bam_dir)) {
      unlink(bam_dir, recursive = TRUE)
    }
    dir.create(bam_dir, showWarnings = FALSE)
    
    if (dir.exists(input_dir)) {
      unlink(input_dir, recursive = TRUE)
    }
    dir.create(input_dir, showWarnings = FALSE)
    
    # Copy uploaded files to the Input directory
    file.copy(input$fastq_gz_upload$datapath, input_dir, overwrite = TRUE)
    
    # Rename uploaded files with meaningful names
    system(paste("mv", file.path(input_dir, "0.gz"), file.path(input_dir, "sample_1.fastq.gz")))
    system(paste("mv", file.path(input_dir, "1.gz"), file.path(input_dir, "sample_2.fastq.gz")))
    
    shinyalert("Success", paste("FASTQ Files Upload Successful."), type = "success")
  }
  w$hide()
})

## BAM Upload & Warning & Convert
observeEvent(input$bam_upload_button, {
  tmp_dir <- tempdir()
  bam_dir <- file.path(tmp_dir, "BAM")
  input_dir <- file.path(tmp_dir, "Input")
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) # Transparent background
  )
  w$show() # Show the loading animation
  
  if (is.null(input$bam_upload)) {
    shinyalert("Error", "No files selected.", type = "error")
    w$hide()
    return()
  } else {
    if (dir.exists(bam_dir)) {
      unlink(bam_dir, recursive = TRUE)
    }
    dir.create(bam_dir, showWarnings = FALSE)
    
    # Convert BAM to Fastq
    # file.copy(input$bam_upload$datapath, bam_dir, overwrite = TRUE)
    # bam_file <- file.path(bam_dir, "sample.bam")
    bam_file <- file.path(bam_dir, "sample.bam")
    file.copy(input$bam_upload$datapath, bam_file, overwrite = TRUE)

    system(paste("/arcasHLA/arcasHLA extract", bam_file, "-o", bam_dir, "-t", num_cores, "-v"))
    
    # Remove existing Input directory and recreate it
    if (dir.exists(input_dir)) {
      unlink(input_dir, recursive = TRUE)
    }
    dir.create(input_dir, showWarnings = FALSE)
    
    # Move sorted files to the Input directory
    system(paste("mv", file.path(bam_dir, "sample.sorted.1.fq.gz"), file.path(input_dir, "sample_1.fastq.gz")))
    system(paste("mv", file.path(bam_dir, "sample.sorted.2.fq.gz"), file.path(input_dir, "sample_2.fastq.gz")))
    
    shinyalert("Success", paste("Bam File Upload Successful."), type = "success")
  }
  w$hide()
})
# ------------------------------------------------------------
# Common reload core
# source_type: "history" or "demo"
# ------------------------------------------------------------
load_result_from_dir <- function(
    target_dir,
    display_id,
    source_type = c("history", "demo"),
    show_success_alert = TRUE,
    id_box_label = "JobID",
    mode_value = "reload",
    use_waiter = TRUE
) {
  source_type <- match.arg(source_type)
  
  target_rds <- file.path(target_dir, "result_data.rds")
  
  if (!dir.exists(target_dir)) {
    shinyalert("Error", paste("Target directory not found:", target_dir), type = "error")
    return(FALSE)
  }
  
  if (!file.exists(target_rds)) {
    shinyalert("Error", paste("result_data.rds not found in:", target_dir), type = "error")
    return(FALSE)
  }
  
  w <- NULL
  if (isTRUE(use_waiter)) {
    w <- Waiter$new(html = spin_loader(), color = transparent(0.8))
    w$show()
    on.exit(w$hide(), add = TRUE)
  }
  
  loaded_data <- readRDS(target_rds)
  
  # ----------------------------------------------------------
  # effective source
  # ----------------------------------------------------------
  effective_source <- source_type
  if (
    identical(source_type, "history") &&
    !is.null(loaded_data$source_type) &&
    nzchar(loaded_data$source_type)
  ) {
    effective_source <- loaded_data$source_type
  }
  
  current_result_source(effective_source)
  
  # ----------------------------------------------------------
  # Demo result -> auto write merged csv into demo pivot pool
  # ----------------------------------------------------------
  if (identical(effective_source, "demo")) {
    tryCatch(
      {
        write_demo_pivot_csv_from_result(
          target_dir = target_dir,
          loaded_data = loaded_data,
          display_id = display_id
        )
      },
      error = function(e) {
        warning(sprintf("Failed to prepare demo pivot csv for %s: %s", display_id, e$message))
      }
    )
  }
  
  # ----------------------------------------------------------
  # optional cleanup
  # ----------------------------------------------------------
  if (!is.null(loaded_data$table)) {
    drop_cols <- intersect(c("PharmGKB", "pharmgkb_allele"), colnames(loaded_data$table))
    if (length(drop_cols) > 0) {
      loaded_data$table <- loaded_data$table[, setdiff(colnames(loaded_data$table), drop_cols), drop = FALSE]
    }
  }
  
  # ----------------------------------------------------------
  # Resolve ZIP
  # ----------------------------------------------------------
  resolved_zip <- resolve_zip_from_loaded_data(target_dir, loaded_data)
  
  if (is.null(resolved_zip) || !file.exists(resolved_zip)) {
    shinyalert("Error", "This result exists but the archived ZIP file is missing.", type = "error")
    return(FALSE)
  }
  
  # ----------------------------------------------------------
  # Update reactive state
  # ----------------------------------------------------------
  current_job_id(display_id)
  package_value(loaded_data$tool)
  
  current_mode(mode_value)
  current_tmp_zip(NULL)
  current_history_zip(normalizePath(resolved_zip, winslash = "/", mustWork = FALSE))
  current_prefix(basename(resolved_zip))
  
  # ----------------------------------------------------------
  # UI: table
  # ----------------------------------------------------------
  output$hla_typing_table <- renderDataTable({
    render_hla_typing_datatable(loaded_data$table)
  })
  
  # ----------------------------------------------------------
  # UI: summary
  # ----------------------------------------------------------
  output$summary_boxes <- renderUI({
    fluidRow(
      column(
        width = 3,
        summaryBox("NGS Type", loaded_data$ngs_type, width = 12, icon = icon("dna"), style = "info")
      ),
      column(
        width = 3,
        summaryBox("IPD-IMGT/HLA Version", loaded_data$version, width = 12, icon = icon("code-branch"), style = "success")
      ),
      column(
        width = 3,
        summaryBox("Typing Tool", loaded_data$tool, width = 12, icon = icon("box"), style = "danger")
      ),
      column(
        width = 3,
        summaryBox(
          id_box_label,
          display_id,
          width = 12,
          icon = icon("fingerprint"),
          style = "primary"
        )
      )
    )
  })
  
  updateTabItems(session, "my_tabs", "hla")
  updateTabsetPanel(session, "hla_panel", selected = "Result")
  
  if (isTRUE(show_success_alert)) {
    shinyalert(
      "Success",
      paste(id_box_label, "loaded successfully:", display_id),
      type = "success"
    )
  }
  
  TRUE
}

# ============================================================
# Manual JobID reload
# ============================================================
observeEvent(input$jobid_load_button, {
  if (is.null(input$jobid_input) || !nzchar(input$jobid_input)) {
    shinyalert("Error", "Please enter a valid JobID.", type = "error")
    return()
  }
  
  jid <- trimws(input$jobid_input)
  target_dir <- file.path(HISTORY_ROOT, jid)
  
  load_result_from_dir(
    target_dir  = target_dir,
    display_id  = jid,
    source_type = "history"
  )
})
## Demo QC option
wes_demo_uploaded_once <- reactiveVal(FALSE)

is_current_demo_wes <- reactive({
  input$upload_choice == "demo" && input$demo_choice == "WES"
})

observeEvent(input$demo_upload_button, {
  if (input$upload_choice == "demo") {
    if (input$demo_choice == "WES") {
      wes_demo_uploaded_once(TRUE)
    } else if (input$demo_choice == "RNA-seq") {
      wes_demo_uploaded_once(FALSE)
    }
  }
})

output$show_qc_button <- reactive({
  (input$upload_choice == "bam" && input$bam_upload_button >= 1) ||
    (input$upload_choice == "demo" && is_current_demo_wes() && wes_demo_uploaded_once())
})
outputOptions(output, "show_qc_button", suspendWhenHidden = FALSE)
