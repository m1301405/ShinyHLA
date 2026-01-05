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
# [新增] JobID 重新載入功能邏輯
#-----------------------------------------------------------------------
observeEvent(input$jobid_load_button, {
  # 1. 檢查輸入
  if (input$jobid_input == "") {
    shinyalert("Error", "Please enter a valid JobID.", type = "error")
    return()
  }
  
  # 2. 定義 RDS 檔案的路徑
  target_rds <- file.path(tmp_dir, "History", input$jobid_input, "result_data.rds")
  
  if (file.exists(target_rds)) {
    w <- Waiter$new(html = spin_loader(), color = transparent(0.8))
    w$show()
    
    # 3. 讀取持久化數據
    loaded_data <- readRDS(target_rds)
    
    # 4. 更新全域變數狀態
    current_job_id(input$jobid_input) # 設定後 IGV 會自動導向歷史目錄
    package_value(loaded_data$tool)
    
    # 5. 重新渲染 HLA Typing Table (包含序列截斷邏輯)
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
          columnDefs = list(list(
            targets = "_all",
            render = JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data != null && data.length > 20 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
              "}"
            )
          ))
        )
      )
    })
    
    # 6. 重新渲染等寬排版的 Summary Boxes (3-3-3-3 佈局)
    output$summary_boxes <- renderUI({
      fluidRow(
        column(width = 3, summaryBox("NGS Type", loaded_data$ngs_type, width = 12, icon = icon("dna"), style = "info")),
        column(width = 3, summaryBox("IPD-IMGT/HLA Version", loaded_data$version, width = 12, icon = icon("code-branch"), style = "success")),
        column(width = 3, summaryBox("Typing Tool", loaded_data$tool, width = 12, icon = icon("box"), style = "danger")),
        column(width = 3, summaryBox("JobID", input$jobid_input, width = 12, icon = icon("fingerprint"), style = "primary"))
      )
    })
    
    # 7. 恢復下載功能 (從歸檔路徑抓取檔案)
    output$hla_typing_download <- downloadHandler(
      filename = function() { basename(loaded_data$zip_path) },
      content = function(file) { file.copy(loaded_data$zip_path, file) }
    )
    
    # 8. 自動切換到 HLA Typing 頁面 (修正 app.R 中的 my_tabs ID)
    updateTabItems(session, "my_tabs", "hla") 
    
    w$hide()
    shinyalert("Success", paste("JobID:", input$jobid_input, "loaded successfully."), type = "success")
    
  } else {
    shinyalert("Error", paste("JobID", input$jobid_input, "not found or incomplete."), type = "error")
  }
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

