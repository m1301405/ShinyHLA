# ============================================================
# Demo Pivot settings
# ============================================================
DEMO_PIVOT_POOL <- "/root/shiny/Demo/pivottable"
dir.create(DEMO_PIVOT_POOL, recursive = TRUE, showWarnings = FALSE)

current_result_source <- reactiveVal("normal")

# ------------------------------------------------------------
# Infer HLA class from allele
# ------------------------------------------------------------
infer_mhc_class <- function(allele) {
  allele <- trimws(as.character(allele))
  
  if (grepl("^A\\*", allele)) {
    "HLA-A"
  } else if (grepl("^B\\*", allele)) {
    "HLA-B"
  } else if (grepl("^C\\*", allele)) {
    "HLA-C"
  } else if (grepl("^DPA1\\*", allele)) {
    "HLA-DPA1"
  } else if (grepl("^DPB1\\*", allele)) {
    "HLA-DPB1"
  } else if (grepl("^DQA1\\*", allele)) {
    "HLA-DQA1"
  } else if (grepl("^DQB1\\*", allele)) {
    "HLA-DQB1"
  } else if (grepl("^DRB1\\*", allele)) {
    "HLA-DRB1"
  } else {
    NA_character_
  }
}

# ------------------------------------------------------------
# Merge alleles by class
# Output aligns with normal *_merged.csv
# ------------------------------------------------------------
allele_merge_hla_class_common <- function(data) {
  result <- data.frame(MHC_class = character(), Allele = character(), stringsAsFactors = FALSE)
  
  if (length(data) == 0) return(result)
  
  for (i in seq_along(data)) {
    class_name <- infer_mhc_class(data[i])
    if (is.na(class_name)) next
    
    if (class_name %in% result$MHC_class) {
      result[result$MHC_class == class_name, "Allele"] <-
        paste(result[result$MHC_class == class_name, "Allele"], data[i], sep = " ")
    } else {
      result <- rbind(
        result,
        data.frame(MHC_class = class_name, Allele = data[i], stringsAsFactors = FALSE)
      )
    }
  }
  
  result
}

# ------------------------------------------------------------
# Extract allele from normal txt
# ------------------------------------------------------------
extract_allele_from_txt_common <- function(file_path) {
  hla_output <- read.table(file_path, quote = "\"", comment.char = "", stringsAsFactors = FALSE)
  hla_output <- subset(hla_output, grepl("^>", V1))
  hla_output <- subset(hla_output, !grepl("_protein", V1))
  hla_output <- gsub("^>", "", hla_output$V1)
  hla_output <- gsub("_nucleotide$", "", hla_output)
  hla_output <- gsub("[\"\\(\\)]", "", hla_output)
  hla_output <- unlist(strsplit(hla_output, ", "))
  trimws(hla_output)
}

# ------------------------------------------------------------
# Extract allele from fasta (demo)
# >A*01:01:01_nucleotide -> A*01:01:01
# ------------------------------------------------------------
extract_allele_from_fasta_common <- function(fasta_path) {
  lines <- readLines(fasta_path, warn = FALSE)
  headers <- grep("^>", lines, value = TRUE)
  
  if (length(headers) == 0) return(character(0))
  
  headers <- gsub("^>", "", headers)
  headers <- headers[!grepl("_protein$", headers)]
  headers <- gsub("_nucleotide$", "", headers)
  headers <- gsub("[\"\\(\\)]", "", headers)
  trimws(headers)
}

# ------------------------------------------------------------
# Find nucleotide.fasta in demo result folder
# unzip if needed
# ------------------------------------------------------------
locate_demo_fasta_common <- function(target_dir, display_id) {
  direct_nuc <- file.path(target_dir, "nucleotide.fasta")
  direct_pro <- file.path(target_dir, "protein.fasta")
  
  if (file.exists(direct_nuc)) return(direct_nuc)
  if (file.exists(direct_pro)) return(direct_pro)
  
  zip_candidates <- list.files(target_dir, pattern = "\\.zip$", full.names = TRUE)
  if (length(zip_candidates) == 0) {
    stop(sprintf("No zip file found in demo folder: %s", target_dir))
  }
  
  extract_dir <- file.path(tempdir(), "demo_pivot_extract", display_id)
  unlink(extract_dir, recursive = TRUE, force = TRUE)
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  
  unzip(zip_candidates[1], exdir = extract_dir)
  
  nuc <- list.files(extract_dir, pattern = "^nucleotide\\.fasta$", recursive = TRUE, full.names = TRUE)
  if (length(nuc) > 0) return(nuc[1])
  
  pro <- list.files(extract_dir, pattern = "^protein\\.fasta$", recursive = TRUE, full.names = TRUE)
  if (length(pro) > 0) return(pro[1])
  
  stop(sprintf("No nucleotide.fasta or protein.fasta found in %s", target_dir))
}

# ------------------------------------------------------------
# Write one demo *_merged.csv into demo pivot pool
# ------------------------------------------------------------
write_demo_pivot_csv_from_result <- function(target_dir, loaded_data, display_id) {
  dir.create(DEMO_PIVOT_POOL, recursive = TRUE, showWarnings = FALSE)
  
  fasta_path <- locate_demo_fasta_common(target_dir, display_id)
  alleles <- extract_allele_from_fasta_common(fasta_path)
  
  if (length(alleles) == 0) {
    stop(sprintf("No alleles parsed from fasta: %s", fasta_path))
  }
  
  result <- allele_merge_hla_class_common(alleles)
  if (nrow(result) == 0) {
    stop(sprintf("No valid class rows generated for %s", display_id))
  }
  
  row_number <- nrow(result)
  
  df <- data.frame(
    Tool = rep(loaded_data$tool, row_number),
    NGS_type = rep(loaded_data$ngs_type, row_number),
    IMGTHLA_version = rep(loaded_data$version, row_number),
    stringsAsFactors = FALSE
  )
  
  merged_output <- cbind(result, df)
  
  tool_stub <- dplyr::case_when(
    loaded_data$tool == "OptiType" ~ "optitype",
    loaded_data$tool == "arcasHLA" ~ "arcashla",
    loaded_data$tool == "HLA-HD" ~ "hlahd",
    TRUE ~ "t1k"
  )
  
  ngs_stub <- ifelse(loaded_data$ngs_type == "RNA-seq", "RNAseq", loaded_data$ngs_type)
  
  output_file <- file.path(
    DEMO_PIVOT_POOL,
    paste0(tool_stub, "_", ngs_stub, "_", loaded_data$version, "_", display_id, "_merged.csv")
  )
  
  write.csv(merged_output, file = output_file, row.names = FALSE)
  output_file
}

# ------------------------------------------------------------
# Write merged_all / classI / classII
# ------------------------------------------------------------
write_split_pivot_files_common <- function(merged_pivottable, output_directory) {
  write.csv(merged_pivottable, file = file.path(output_directory, "merged_all.csv"), row.names = FALSE)
  
  group1 <- c("HLA-A", "HLA-B", "HLA-C")
  group2 <- setdiff(unique(merged_pivottable$MHC_class), group1)
  
  data_group1 <- merged_pivottable[merged_pivottable$MHC_class %in% group1, , drop = FALSE]
  data_group2 <- merged_pivottable[merged_pivottable$MHC_class %in% group2, , drop = FALSE]
  
  write.csv(data_group1, file = file.path(output_directory, "MHC_classI.csv"), row.names = FALSE)
  write.csv(data_group2, file = file.path(output_directory, "MHC_classII.csv"), row.names = FALSE)
}

# ------------------------------------------------------------
# Build Pivotable workspace from normal Output/*.txt
# ------------------------------------------------------------
prepare_normal_pivot_workspace <- function(tmp_dir) {
  output_directory <- file.path(tmp_dir, "Pivotable")
  unlink(output_directory, recursive = TRUE, force = TRUE)
  dir.create(output_directory, showWarnings = FALSE)
  
  file_list <- list.files(path = file.path(tmp_dir, "Output"), pattern = "\\.txt$", full.names = TRUE)
  if (length(file_list) == 0) return(FALSE)
  
  for (file in file_list) {
    file_data <- extract_allele_from_txt_common(file)
    result <- allele_merge_hla_class_common(file_data)
    row_number <- nrow(result)
    if (row_number == 0) next
    
    df <- data.frame(
      Tool = character(row_number),
      NGS_type = character(row_number),
      IMGTHLA_version = character(row_number),
      stringsAsFactors = FALSE
    )
    
    if (grepl("_WES", file, ignore.case = TRUE)) {
      df[1:row_number, "NGS_type"] <- "WES"
    } else {
      df[1:row_number, "NGS_type"] <- "RNA-seq"
    }
    
    version_match <- regmatches(file, regexpr("v[0-9]+\\.[0-9]+\\.[0-9]+", file))
    if (length(version_match) > 0) {
      df[1:row_number, "IMGTHLA_version"] <- version_match
    }
    
    df[1:row_number, "Tool"] <- dplyr::case_when(
      grepl("optitype", file, ignore.case = TRUE) ~ "OptiType",
      grepl("arcashla", file, ignore.case = TRUE) ~ "arcasHLA",
      grepl("hlahd", file, ignore.case = TRUE) ~ "HLA-HD",
      TRUE ~ "T1K"
    )
    
    merged_output <- cbind(result, df)
    output_file <- file.path(output_directory, paste0(sub("\\.txt$", "", basename(file)), "_merged.csv"))
    write.csv(merged_output, file = output_file, row.names = FALSE)
  }
  
  file_paths <- list.files(path = output_directory, pattern = "_merged\\.csv$", full.names = TRUE)
  if (length(file_paths) == 0) return(FALSE)
  
  merged_pivottable <- do.call(rbind, lapply(file_paths, function(fp) {
    read.csv(fp, stringsAsFactors = FALSE)
  }))
  
  if (nrow(merged_pivottable) == 0) return(FALSE)
  
  write_split_pivot_files_common(merged_pivottable, output_directory)
  TRUE
}

# ------------------------------------------------------------
# Build Pivotable workspace from demo pivot pool
# ------------------------------------------------------------
prepare_demo_pivot_workspace <- function(tmp_dir) {
  output_directory <- file.path(tmp_dir, "Pivotable")
  unlink(output_directory, recursive = TRUE, force = TRUE)
  dir.create(output_directory, showWarnings = FALSE)
  
  file_paths <- list.files(path = DEMO_PIVOT_POOL, pattern = "_merged\\.csv$", full.names = TRUE)
  if (length(file_paths) == 0) return(FALSE)
  
  merged_pivottable <- do.call(rbind, lapply(file_paths, function(fp) {
    read.csv(fp, stringsAsFactors = FALSE)
  }))
  
  if (nrow(merged_pivottable) == 0) return(FALSE)
  
  write_split_pivot_files_common(merged_pivottable, output_directory)
  TRUE
}