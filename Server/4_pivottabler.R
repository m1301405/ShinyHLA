## pivotable
observeEvent(input$mhci_pivottable_button,{
  ori_dir <- getwd()
  tmp_dir <- tempdir()
  setwd(tmp_dir)
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8)
  )
  w$show()
  
  dir.create("./Pivotable", showWarnings = FALSE)
  # Define Function-1
  extract_allele <- function(file_path) {
    hla_output <- read.table(file_path, quote="\"", comment.char="")
    hla_output <- subset(hla_output, grepl("^>", V1))
    hla_output <- subset(hla_output, !grepl("_protein", V1))
    hla_output <- gsub("[>_nucleotide]", "", hla_output)
    hla_output <- gsub("[\"\\(\\)]", "", hla_output)
    hla_output <- unlist(strsplit(hla_output, ", "))
    return(hla_output)
  }
  # Define Function-2
  allele_merge_hla_class <- function(data) {
    result <- data.frame(MHC_class = character(), Allele = character(), stringsAsFactors = FALSE)
    
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
  
  # Get all file names in the directory
  file_list <- list.files(path = "./Output", pattern = "\\.txt$", full.names = TRUE)
  output_directory <- "./Pivotable"
  
  # Iterate through all files and read their contents
  for (file in file_list) {
    file_data <-extract_allele(file)
    result <- allele_merge_hla_class(file_data)
    row_number <- nrow(result)
    
    # Create Tool, IMGTHLA, and Sequencing_type tables
    df <- data.frame(
      Tool = character(row_number),
      NGS_type = character(row_number),
      IMGTHLA_version = character(row_number),
      stringsAsFactors = FALSE)
    # Insert Sequence types
    if (grepl("_WES",file)) {
      df[1:row_number, "NGS_type"] <- "WES"
    } else {
      df[1:row_number, "NGS_type"] <- "RNA-seq"
    }
    # Insert IMGTHLA version
    version_match <- regmatches(file, regexpr("v[0-9]+\\.[0-9]+\\.[0-9]+", file))
    if (length(version_match) > 0) {
      df[1:row_number, "IMGTHLA_version"] <- version_match
    }
    # Insert Tool type
    df[1:row_number, "Tool"] <- case_when(
      grepl("optitype", file) ~ "OptiType",
      grepl("arcashla", file) ~ "arcasHLA",
      grepl("hlahd", file) ~ "HLA-HD",
      TRUE ~ "T1K"
    )
    merged_output <- cbind(result,df)
    # export
    output_file <- file.path(output_directory, paste0(sub("\\.txt$", "", basename(file)),"_merged",".csv"))
    write.csv(merged_output, file = output_file, row.names = FALSE)
  }
  
  # Get file paths of all _merged.csv files in the directory
  file_paths <- list.files(path = "./Pivotable", pattern = "_merged\\.csv$", full.names = TRUE)
  
  # Initialize the merged dataframe
  merged_pivottable <- data.frame()
  
  # Process each file
  for (file_path in file_paths) {
    df <- read.csv(file_path)
    # Merge the current file's dataframe with the previous dataframe
    merged_pivottable <- rbind(merged_pivottable, df)
  }
  
  # Save the merged dataframe to a new file
  write.csv(merged_pivottable, file = "./Pivotable/merged_all.csv", row.names = FALSE)
  
  # Split the MHC_class column into two groups
  group1 <- c("HLA-A", "HLA-B", "HLA-C")
  group2 <- setdiff(unique(merged_pivottable$MHC_class), group1)
  
  # Split the data
  data_group1 <- merged_pivottable[merged_pivottable$MHC_class %in% group1,]
  data_group2 <- merged_pivottable[merged_pivottable$MHC_class %in% group2,]
  
  # Save as a CSV file
  write.csv(data_group1, file = "./Pivotable/MHC_classI.csv", row.names = FALSE)
  write.csv(data_group2, file = "./Pivotable/MHC_classII.csv", row.names = FALSE)
  
  # Check if the file is empty
  if (file.size("./Pivotable/MHC_classI.csv") == 0) {
    shinyalert(
      title = "Warning",
      text = "No genotyping results available for integration.",
      type = "error"
    )
    w$hide()
    setwd(ori_dir)
    return()
  }
  
  # Use tryCatch to prevent `read.csv()` from failing
  mhc_classI_data <- tryCatch(
    read.csv("./Pivotable/MHC_classI.csv", stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  
  # Check if `read.csv()` still results in an empty dataset
  if (nrow(mhc_classI_data) == 0) {
    shinyalert(
      title = "Warning",
      text = "No genotyping results available for integration.",
      type = "error"
    )
    w$hide()
    setwd(ori_dir)
    return()
  }
  
  # MHC class I pivottabler
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
      pt$addData(merged_all)

      if (length(input$mhci_col_list) >= 1) pt$addColumnDataGroups(input$mhci_col_list[1], addTotal = FALSE)
      if (length(input$mhci_col_list) >= 2) pt$addColumnDataGroups(input$mhci_col_list[2], addTotal = FALSE)

      if (length(input$mhci_row_list) >= 1) pt$addRowDataGroups(input$mhci_row_list[1], outlineBefore = list(isEmpty = TRUE, groupStyleDeclarations = list(color = "blue", width = "200px")), addTotal = FALSE)
      if (length(input$mhci_row_list) >= 2) pt$addRowDataGroups(input$mhci_row_list[2], addTotal = FALSE)

      pt$defineCalculation(
        calculationName = "Unique",
        summariseExpression = "paste(unique(Allele), collapse = ', ')"
      )

      pt$evaluatePivot()

      # excel
      wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
      addWorksheet(wb, "Data")
      pt$writeToExcelWorksheet(wb=wb, wsName="Data",
                               topRowNumber=2, leftMostColumnNumber=2, applyStyles=TRUE)
      saveWorkbook(wb, file = file.path(tmp_dir,"Pivotable/MHC_class_I_pivottable.xlsx"), overwrite = TRUE)

      # Pivot table
      pivottabler(pt)
    })
  }
  w$hide() 
  
  setwd(ori_dir)
})

observeEvent(input$mhcii_pivottable_button,{
  ori_dir <- getwd()
  tmp_dir <- tempdir()
  setwd(tmp_dir)
  
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8)
  )
  w$show()
  
  dir.create("./Pivotable", showWarnings = FALSE)
  # Define Function-1
  extract_allele <- function(file_path) {
    hla_output <- read.table(file_path, quote="\"", comment.char="")
    hla_output <- subset(hla_output, grepl("^>", V1))
    hla_output <- subset(hla_output, !grepl("_protein", V1))
    hla_output <- gsub("[>_nucleotide]", "", hla_output)
    hla_output <- gsub("[\"\\(\\)]", "", hla_output)
    hla_output <- unlist(strsplit(hla_output, ", "))
    return(hla_output)
  }
  # Define Function-2
  allele_merge_hla_class <- function(data) {
    result <- data.frame(MHC_class = character(), Allele = character(), stringsAsFactors = FALSE)
    
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
  
  file_list <- list.files(path = "./Output", pattern = "*0.txt", full.names = TRUE)
  
  output_directory <- "./Pivotable"
  
  for (file in file_list) {
    file_data <-extract_allele(file)
    result <- allele_merge_hla_class(file_data)
    
    row_number <- nrow(result)
    
    df <- data.frame(
      Tool = character(row_number),
      NGS_type = character(row_number),
      IMGTHLA_version = character(row_number),
      stringsAsFactors = FALSE)
    
    if (grepl("_WES",file)) {
      df[1:row_number, "NGS_type"] <- "WES"
    } else {
      df[1:row_number, "NGS_type"] <- "RNA-seq"
    }

    version_match <- regmatches(file, regexpr("v[0-9]+\\.[0-9]+\\.[0-9]+", file))
    if (length(version_match) > 0) {
      df[1:row_number, "IMGTHLA_version"] <- version_match
    }
    
    df[1:row_number, "Tool"] <- case_when(
      grepl("optitype", file) ~ "OptiType",
      grepl("arcashla", file) ~ "arcasHLA",
      grepl("hlahd", file) ~ "HLA-HD",
      TRUE ~ "T1K"
    )
    
    merged_output <- cbind(result,df)
    
    output_file <- file.path(output_directory, paste0(sub("\\.txt$", "", basename(file)),"_merged",".csv"))
    write.csv(merged_output, file = output_file, row.names = FALSE)
  }
  
  file_paths <- list.files(path = "./Pivotable", pattern = "_merged\\.csv$", full.names = TRUE)
  
  merged_pivottable <- data.frame()
  
  for (file_path in file_paths) {
    df <- read.csv(file_path)
    
    merged_pivottable <- rbind(merged_pivottable, df)
  }
  
  write.csv(merged_pivottable, file = "./Pivotable/merged_all.csv", row.names = FALSE)
  
  group1 <- c("HLA-A", "HLA-B", "HLA-C")
  group2 <- setdiff(unique(merged_pivottable$MHC_class), group1)
  
  data_group1 <- merged_pivottable[merged_pivottable$MHC_class %in% group1,]
  data_group2 <- merged_pivottable[merged_pivottable$MHC_class %in% group2,]
  
  write.csv(data_group1, file = "./Pivotable/MHC_classI.csv", row.names = FALSE)
  write.csv(data_group2, file = "./Pivotable/MHC_classII.csv", row.names = FALSE)
  
  # Check if the file is empty
  if (file.size("./Pivotable/MHC_classII.csv") == 0) {
    shinyalert(
      title = "Warning",
      text = "No genotyping results available for integration.",
      type = "error"
    )
    w$hide()
    setwd(ori_dir)
    return()
  }
  
  # Use tryCatch to prevent `read.csv()` from failing
  mhc_classII_data <- tryCatch(
    read.csv("./Pivotable/MHC_classII.csv", stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  
  # Check if `read.csv()` still results in an empty dataset
  if (nrow(mhc_classII_data) == 0) {
    shinyalert(
      title = "Warning",
      text = "No genotyping results available for integration.",
      type = "error"
    )
    w$hide()
    setwd(ori_dir)
    return()
  }
  
  # MHC class II pivottabler
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
      
      pt$addData(merged_all)
      
      if (length(input$mhcii_col_list) >= 1) pt$addColumnDataGroups(input$mhcii_col_list[1], addTotal = FALSE)
      if (length(input$mhcii_col_list) >= 2) pt$addColumnDataGroups(input$mhcii_col_list[2], addTotal = FALSE)
      
      if (length(input$mhcii_row_list) >= 1) pt$addRowDataGroups(input$mhcii_row_list[1], outlineBefore = list(isEmpty = TRUE, groupStyleDeclarations = list(color = "blue", width = "200px")), addTotal = FALSE)
      if (length(input$mhcii_row_list) >= 2) pt$addRowDataGroups(input$mhcii_row_list[2], addTotal = FALSE)
      
      pt$defineCalculation(
        calculationName = "Unique",
        summariseExpression = "paste(unique(Allele), collapse = ', ')"
      )
      
      pt$evaluatePivot()
      
      wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
      addWorksheet(wb, "Data")
      pt$writeToExcelWorksheet(wb=wb, wsName="Data",
                               topRowNumber=2, leftMostColumnNumber=2, applyStyles=TRUE)
      saveWorkbook(wb, file = file.path(tmp_dir,"Pivotable/MHC_class_II_pivottable.xlsx"), overwrite = TRUE)
      
      pivottabler(pt)
    })
  }
  w$hide() 
  
  setwd(ori_dir)
})

## MHC class I pivottable download
output$mhci_pivot_download <- downloadHandler(
  filename = function() {
    "MHC_class_I_pivottable.xlsx"
  },
  content = function(file) {
    pivot_xlsx_path <- file.path(tempdir(), "Pivotable/MHC_class_I_pivottable.xlsx")
    
    # Check if the PivotTable file exists
    if (!file.exists(pivot_xlsx_path)) {
      shinyalert(
        title = "Warning",
        text = "No genotyping results available for download.",
        type = "error"
      )
      return()
    }
    # download
    file.copy(pivot_xlsx_path, file)
  }
)

## MHC class II pivottable download
output$mhcii_pivot_download <- downloadHandler(
  tmp_dir <- tempdir(),
  filename = function() {
    "MHC_class_II_pivottable.xlsx"
  },
  content = function(file) {
    pivot_xlsx_path <- file.path(tempdir(), "Pivotable/MHC_class_II_pivottable.xlsx")
    
    # Check if the PivotTable file exists
    if (!file.exists(pivot_xlsx_path)) {
      shinyalert(
        title = "Warning",
        text = "No genotyping results available for download.",
        type = "error"
      )
      return()
    }
    # download
    file.copy(pivot_xlsx_path, file)
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
