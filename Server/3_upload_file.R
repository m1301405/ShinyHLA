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
    
    file.copy(input$bam_upload$datapath, bam_dir, overwrite = TRUE)
    
    # Convert BAM to Fastq
    bam_file <- file.path(bam_dir, "sample.bam")
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

#-----------------------------------------------------------------------
observeEvent(input$jobid_load_button, {
  if (is.null(input$jobid_input) || !nzchar(input$jobid_input)) {
    shinyalert("Error", "Please enter a valid JobID.", type = "error")
    return()
  }
  
  jid <- input$jobid_input
  HISTORY_ROOT <- "/root/shiny/History"
  
  target_dir <- file.path(HISTORY_ROOT, jid)
  target_rds <- file.path(target_dir, "result_data.rds")
  
  if (!file.exists(target_rds)) {
    shinyalert("Error", paste("JobID", jid, "not found or incomplete."), type = "error")
    return()
  }
  
  w <- Waiter$new(html = spin_loader(), color = transparent(0.8))
  w$show()
  
  loaded_data <- readRDS(target_rds)
  
  # ------------------------------------------------------------
  # Reload ZIP resolution (canonical first)
  # ------------------------------------------------------------
  canonical_zip <- file.path(target_dir, "output.zip")
  
  fallback_zip <- if (!is.null(loaded_data$zip_path) && nzchar(loaded_data$zip_path)) loaded_data$zip_path else NULL
  
  zip_candidates <- list.files(target_dir, pattern = "\\.zip$", full.names = TRUE)
  first_zip <- if (length(zip_candidates) > 0) zip_candidates[1] else NULL
  
  resolved_zip <- NULL
  if (file.exists(canonical_zip)) {
    resolved_zip <- canonical_zip
  } else if (!is.null(fallback_zip) && file.exists(fallback_zip)) {
    resolved_zip <- fallback_zip
  } else if (!is.null(first_zip) && file.exists(first_zip)) {
    resolved_zip <- first_zip
  }
  
  if (is.null(resolved_zip) || !file.exists(resolved_zip)) {
    w$hide()
    shinyalert("Error", "This JobID exists but the archived ZIP file is missing.", type = "error")
    return()
  }
  
  # ------------------------------------------------------------
  # Update reactive state for Reload mode
  # ------------------------------------------------------------
  current_job_id(jid)
  package_value(loaded_data$tool)
  
  current_mode("reload")
  current_tmp_zip(NULL) # reload 不用 tmp
  
  current_history_zip(normalizePath(resolved_zip, winslash = "/", mustWork = FALSE))
  current_prefix(basename(resolved_zip))
  
  # Reset IGV button to force user to click again to load/update IGV
  shinyjs::reset("igv_reference")
  
  # ------------------------------------------------------------
  # UI: table
  # ------------------------------------------------------------
  output$hla_typing_table <- renderDataTable({
    datatable(
      loaded_data$table,
      filter = "top",
      rownames = FALSE,
      selection = "none",
      class = "nowrap",
      options = list(
        pageLength = 8,
        autoWidth = TRUE,
        scrollX = TRUE,
        columnDefs = list(list(
          width = "100px",
          targets = "_all",
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data != null && data.length > 40 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 40) + '...</span>' : data;",
            "}"
          )
        ))
      )
    )
  })
  
  # UI: summary
  output$summary_boxes <- renderUI({
    fluidRow(
      column(width = 3, summaryBox("NGS Type", loaded_data$ngs_type, width = 12, icon = icon("dna"), style = "info")),
      column(width = 3, summaryBox("IPD-IMGT/HLA Version", loaded_data$version, width = 12, icon = icon("code-branch"), style = "success")),
      column(width = 3, summaryBox("Typing Tool", loaded_data$tool, width = 12, icon = icon("box"), style = "danger")),
      column(width = 3, summaryBox("JobID", jid, width = 12, icon = icon("fingerprint"), style = "primary"))
    )
  })
  
  updateTabItems(session, "my_tabs", "hla")
  
  w$hide()
  shinyalert("Success", paste("JobID:", jid, "loaded successfully."), type = "success")
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
