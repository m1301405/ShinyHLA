## pivotable
observeEvent(input$mhci_pivottable_button,{
  ori_dir <- getwd()
  tmp_dir <- tempdir()
  setwd(tmp_dir)
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) # transparent background
  )
  w$show() # show the waiter
  
  # 建立資料夾
  dir.create("./Pivotable", showWarnings = FALSE)
  # 定義函數-1
  extract_allele <- function(file_path) {
    hla_output <- read.table(file_path, quote="\"", comment.char="")
    hla_output <- subset(hla_output, grepl("^>", V1))
    hla_output <- subset(hla_output, !grepl("_protein", V1))
    hla_output <- gsub("[>_nucleotide]", "", hla_output)
    hla_output <- gsub("[\"\\(\\)]", "", hla_output)
    hla_output <- unlist(strsplit(hla_output, ", "))
    return(hla_output)
  }
  # 定義函數-2
  allele_merge_hla_class <- function(data) {
    # 初始化結果數據框
    result <- data.frame(MHC_class = character(), Allele = character(), stringsAsFactors = FALSE)
    
    # 遍歷數據
    for (i in 1:length(data)) {
      if (grepl("^A", data[i])) {
        class_name <- "HLA-A"
      } else if (grepl("^B", data[i])) {
        class_name <- "HLA-B"
      } else if (grepl("^C", data[i])) {
        class_name <- "HLA-C"
      } else if (grepl("^DPB1", data[i])) {
        class_name <- "HLA-DPB1"
      } else if (grepl("^DQA1", data[i])) {
        class_name <- "HLA-DQA1"
      } else if (grepl("^DQB1", data[i])) {
        class_name <- "HLA-DQB1"
      } else if (grepl("^DRB1", data[i])) {
        class_name <- "HLA-DRB1"
      } else if (grepl("^DPA1", data[i])) {
        class_name <- "HLA-DPA1"
      }
      
      if (class_name %in% result$MHC_class) {
        result[result$MHC_class == class_name, "Allele"] <- paste(result[result$MHC_class == class_name, "Allele"], data[i], sep = " ")
      } else {
        result <- rbind(result, data.frame(MHC_class = class_name, Allele = data[i], stringsAsFactors = FALSE))
      }
    }
    
    return(result)
    
  }
  
  # 獲取目錄下的所有文件名
  file_list <- list.files(path = "./Output", pattern = "*0.txt", full.names = TRUE)
  
  # 設定輸出目錄
  output_directory <- "./Pivotable"
  
  # 迭代所有文件並讀取內容
  for (file in file_list) {
    # 讀取每個文件的內容到一個dataframe
    file_data <-extract_allele(file)
    result <- allele_merge_hla_class(file_data)
    
    row_number <- nrow(result)
    # 建立Tool, IMGTHLA, Sequencing_type表格
    df <- data.frame(
      Tool = character(row_number),
      NGS_type = character(row_number),
      IMGTHLA_version = character(row_number),
      stringsAsFactors = FALSE)
    # 添入Sequence類型
    if (grepl("_WES",file)) {
      df[1:row_number, "NGS_type"] <- "WES"
    } else {
      df[1:row_number, "NGS_type"] <- "RNA-seq"
    }
    # 填入IMGTHLA版本
    version_match <- regmatches(file, regexpr("v[0-9]+\\.[0-9]+\\.[0-9]+", file))
    if (length(version_match) > 0) {
      df[1:row_number, "IMGTHLA_version"] <- version_match
    }
    # 添入工具類型
    df[1:row_number, "Tool"] <- case_when(
      grepl("optitype", file) ~ "OptiType",
      grepl("arcashla", file) ~ "arcasHLA",
      grepl("hlahd", file) ~ "HLA-HD",
      TRUE ~ "SpecHLA"
    )
    # 合併表格
    merged_output <- cbind(result,df)
    # 輸出表格
    output_file <- file.path(output_directory, paste0(sub("\\.txt$", "", basename(file)),"_merged",".csv"))
    write.csv(merged_output, file = output_file, row.names = FALSE)
  }
  
  # 獲取目錄下所有的_merged.csv文件的文件路徑
  file_paths <- list.files(path = "./Pivotable", pattern = "_merged\\.csv$", full.names = TRUE)
  
  # 初始化合併後的數據框
  merged_pivottable <- data.frame()
  
  # 對每個文件進行操作
  for (file_path in file_paths) {
    # 讀取文件
    df <- read.csv(file_path)
    
    # 將當前文件的數據框與之前的數據框進行合併
    merged_pivottable <- rbind(merged_pivottable, df)
  }
  
  # 將所有文件合併後的數據框保存到一個新文件中
  write.csv(merged_pivottable, file = "./Pivotable/merged_all.csv", row.names = FALSE)
  
  # 將MHC_class列分成兩組
  group1 <- c("HLA-A", "HLA-B", "HLA-C")
  group2 <- setdiff(unique(merged_pivottable$MHC_class), group1)
  
  # 分割數據
  data_group1 <- merged_pivottable[merged_pivottable$MHC_class %in% group1,]
  data_group2 <- merged_pivottable[merged_pivottable$MHC_class %in% group2,]
  
  # 儲存為CSV文件
  write.csv(data_group1, file = "./Pivotable/MHC_classI.csv", row.names = FALSE)
  write.csv(data_group2, file = "./Pivotable/MHC_classII.csv", row.names = FALSE)
  
  
  # 執行MHC class I pivottabler
  if (length(input$mhci_col_list) > 2 || length(input$mhci_row_list) > 2 || length(input$mhci_col_list) == 0 || length(input$mhci_row_list) == 0) {
    shinyalert(
      title = "Warning",
      text = "One or more rank lists do not meet the requirements!",
      type = "error"
    )
  } else {
    output$mhci_pivot_table <- renderPivottabler({
      pt <- PivotTable$new()
      merged_all <- read_csv(file.path(tmp_dir,"Pivotable/MHC_classI.csv"))
      # 添加數據
      pt$addData(merged_all)
      
      if (length(input$mhci_col_list) >= 1) pt$addColumnDataGroups(input$mhci_col_list[1], addTotal = FALSE)
      if (length(input$mhci_col_list) >= 2) pt$addColumnDataGroups(input$mhci_col_list[2], addTotal = FALSE)
      
      if (length(input$mhci_row_list) >= 1) pt$addRowDataGroups(input$mhci_row_list[1], outlineBefore = list(isEmpty = TRUE, groupStyleDeclarations = list(color = "blue", width = "200px")), addTotal = FALSE)
      if (length(input$mhci_row_list) >= 2) pt$addRowDataGroups(input$mhci_row_list[2], addTotal = FALSE)
      
      # 定義計算
      pt$defineCalculation(
        calculationName = "Unique",
        summariseExpression = "paste(unique(Allele), collapse = ', ')"
      )
      
      # 評估樞紐分析表
      pt$evaluatePivot()
      
      # 建立excel
      wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
      addWorksheet(wb, "Data")
      pt$writeToExcelWorksheet(wb=wb, wsName="Data",
                               topRowNumber=2, leftMostColumnNumber=2, applyStyles=TRUE)
      saveWorkbook(wb, file = file.path(tmp_dir,"Pivotable/MHC_class_I_pivottable.xlsx"), overwrite = TRUE)
      
      # 建立Pivot table
      pivottabler(pt)
    })
  }
  # 隱藏等待動畫
  w$hide() 
  
  setwd(ori_dir)
})

observeEvent(input$mhcii_pivottable_button,{
  ori_dir <- getwd()
  tmp_dir <- tempdir()
  setwd(tmp_dir)
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) # transparent background
  )
  w$show() # show the waiter
  # 執行MHC class II pivottabler
  if (length(input$mhcii_col_list) > 2 || length(input$mhcii_row_list) > 2 || length(input$mhcii_col_list) == 0 || length(input$mhcii_row_list) == 0) {
    shinyalert(
      title = "Warning",
      text = "One or more rank lists do not meet the requirements!",
      type = "error"
    )
  } else {
    output$mhcii_pivot_table <- renderPivottabler({
      pt <- PivotTable$new()
      merged_all <- read_csv(file.path(tmp_dir,"Pivotable/MHC_classII.csv"))
      # 添加數據
      pt$addData(merged_all)
      
      if (length(input$mhcii_col_list) >= 1) pt$addColumnDataGroups(input$mhcii_col_list[1], addTotal = FALSE)
      if (length(input$mhcii_col_list) >= 2) pt$addColumnDataGroups(input$mhcii_col_list[2], addTotal = FALSE)
      
      if (length(input$mhcii_row_list) >= 1) pt$addRowDataGroups(input$mhcii_row_list[1], outlineBefore = list(isEmpty = TRUE, groupStyleDeclarations = list(color = "blue", width = "200px")), addTotal = FALSE)
      if (length(input$mhcii_row_list) >= 2) pt$addRowDataGroups(input$mhcii_row_list[2], addTotal = FALSE)
      
      # 定義計算
      pt$defineCalculation(
        calculationName = "Unique",
        summariseExpression = "paste(unique(Allele), collapse = ', ')"
      )
      
      # 評估樞紐分析表
      pt$evaluatePivot()
      
      # 建立excel
      wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
      addWorksheet(wb, "Data")
      pt$writeToExcelWorksheet(wb=wb, wsName="Data",
                               topRowNumber=2, leftMostColumnNumber=2, applyStyles=TRUE)
      saveWorkbook(wb, file = file.path(tmp_dir,"Pivotable/MHC_class_II_pivottable.xlsx"), overwrite = TRUE)
      
      # 建立Pivot table
      pivottabler(pt)
    })
  }
  # 隱藏等待動畫
  w$hide() 
  
  setwd(ori_dir)
})


## MHC class I pivottable download
output$mhci_pivot_download <- downloadHandler(
  tmp_dir <- tempdir(),
  filename = function() {
    "MHC_class_I_pivottable.xlsx"
  },
  content = function(file) {
    # 指定要下載的文件的路徑
    file.copy(file.path(tmp_dir, "Pivotable/MHC_class_I_pivottable.xlsx"), file)
  }
)

## MHC class II pivottable download
output$mhcii_pivot_download <- downloadHandler(
  tmp_dir <- tempdir(),
  filename = function() {
    "MHC_class_II_pivottable.xlsx"
  },
  content = function(file) {
    # 指定要下載的文件的路徑
    file.copy(file.path(tmp_dir, "Pivotable/MHC_class_II_pivottable.xlsx"), file)
  }
)

## file numbers warning of MHC class I
observe({
  mhci_warning_message <- NULL
  
  if (length(input$mhci_col_list) > 2) {
    mhci_warning_message <- "Warning: You can only have up to two items in the Columns list!"
  } else if (is.null(input$mhci_col_list) || length(input$mhci_col_list) == 0) {
    mhci_warning_message <- "Warning: The Columns list is empty!"
  }
  
  if (length(input$mhci_row_list) > 2) {
    if (is.null(mhci_warning_message)) {
      mhci_warning_message <- "Warning: You can only have up to two items in the Rows list!"
    } else {
      mhci_warning_message <- paste(mhci_warning_message, "Additionally, you can only have up to two items in the Rows list!", sep = " ")
    }
  } else if (is.null(input$mhci_row_list) || length(input$mhci_row_list) == 0) {
    if (is.null(mhci_warning_message)) {
      mhci_warning_message <- "Warning: The Rows list is empty!"
    } else {
      mhci_warning_message <- paste(mhci_warning_message, "Additionally, the Rows list is empty!", sep = " ")
    }
  }
  
  output$mhci_warning_message <- renderUI({
    if (!is.null(mhci_warning_message)) {
      tags$p(style = "color: red;", mhci_warning_message)
    } else {
      NULL
    }
  })
})


## file numbers warning of MHC class II
observe({
  mhcii_warning_message <- NULL
  
  if (length(input$mhcii_col_list) > 2) {
    mhcii_warning_message <- "Warning: You can only have up to two items in the Columns list!"
  } else if (is.null(input$mhcii_col_list) || length(input$mhcii_col_list) == 0) {
    mhcii_warning_message <- "Warning: The Columns list is empty!"
  }
  
  if (length(input$mhcii_row_list) > 2) {
    if (is.null(mhcii_warning_message)) {
      mhcii_warning_message <- "Warning: You can only have up to two items in the Rows list!"
    } else {
      mhcii_warning_message <- paste(mhcii_warning_message, "Additionally, you can only have up to two items in the Rows list!", sep = " ")
    }
  } else if (is.null(input$mhcii_row_list) || length(input$mhcii_row_list) == 0) {
    if (is.null(mhcii_warning_message)) {
      mhcii_warning_message <- "Warning: The Rows list is empty!"
    } else {
      mhcii_warning_message <- paste(mhcii_warning_message, "Additionally, the Rows list is empty!", sep = " ")
    }
  }
  
  output$mhcii_warning_message <- renderUI({
    if (!is.null(mhcii_warning_message)) {
      tags$p(style = "color: red;", mhcii_warning_message)
    } else {
      NULL
    }
  })
})

## MHC class i & ii tabpanel
observeEvent(input$mhci_pivottable_button, {
  updateTabsetPanel(session, "hla_panel", "MHC class I")
})
observeEvent(input$mhcii_pivottable_button, {
  updateTabsetPanel(session, "hla_panel", "MHC class II")
})
## typing tabpanel
observeEvent(input$hla_typing_button, {
  updateTabsetPanel(session, "hla_panel", "Result")
})
