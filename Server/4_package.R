#----------------------------------------------------------
execution_times <- reactiveValues(times = character())
#----------------------------------------------------------
getCurrentTime <- function() {
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M")
  return(formatted_time)
}
#----------------------------------------------------------
append_execution_log <- function(job_id, ngs_type, version, tool, execution_seconds) {
  log_entry <- paste(
    getCurrentTime(),
    "| JobID=", job_id,
    "|", ngs_type,
    "|", version,
    "|", tool,
    "| Execution:", execution_seconds, "s"
  )
  
  execution_times$times <- c(log_entry, execution_times$times)
  
  output$hla_log <- renderUI({
    HTML(paste(
      "<span style='font-size:12px;'>",
      paste(execution_times$times, collapse = "<br>"),
      "</span>"
    ))
  })
}

#----------------------------------------------------------
resolve_zip_from_loaded_data <- function(target_dir, loaded_data) {
  canonical_zip <- file.path(target_dir, "output.zip")
  
  fallback_zip <- if (!is.null(loaded_data$zip_path) && nzchar(loaded_data$zip_path)) {
    loaded_data$zip_path
  } else {
    NULL
  }
  
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
  
  resolved_zip
}

#----------------------------------------------------------
archive_demo_result_to_history <- function(target_dir, demo_job_id) {
  history_dir <- file.path(HISTORY_ROOT, demo_job_id)
  dir.create(history_dir, recursive = TRUE, showWarnings = FALSE)
  
  target_rds <- file.path(target_dir, "result_data.rds")
  if (!file.exists(target_rds)) {
    stop(sprintf("result_data.rds not found in demo target_dir: %s", target_dir))
  }
  
  loaded_data <- readRDS(target_rds)
  resolved_zip <- resolve_zip_from_loaded_data(target_dir, loaded_data)
  
  if (is.null(resolved_zip) || !file.exists(resolved_zip)) {
    stop(sprintf("No zip file found in demo target_dir: %s", target_dir))
  }
  
  history_zip <- file.path(history_dir, basename(resolved_zip))
  ok <- file.copy(resolved_zip, history_zip, overwrite = TRUE)
  
  if (!ok || !file.exists(history_zip)) {
    stop("Failed to archive demo ZIP into /root/shiny/History.")
  }
  
  loaded_data$job_id <- demo_job_id
  loaded_data$zip_path <- normalizePath(history_zip, winslash = "/", mustWork = FALSE)
  loaded_data$source_type <- "demo"
  loaded_data$demo_file_folder <- demo_job_id
  
  saveRDS(
    loaded_data,
    file = file.path(history_dir, "result_data.rds")
  )
  
  list(
    history_dir = history_dir,
    history_zip = normalizePath(history_zip, winslash = "/", mustWork = FALSE),
    loaded_data = loaded_data
  )
}
#----------------------------------------------------------
observeEvent(input$hla_typing_button, {
  # ----------------------------------------------------------
  # Demo mode: simulate analysis, then write to History and load
  # ----------------------------------------------------------
  if (identical(source_mode(), "demo")) {
    req(input$sequence, input$imgthla, input$package)
    
    ngs_type <- trimws(as.character(input$sequence))
    version  <- trimws(as.character(input$imgthla))
    tool     <- trimws(as.character(input$package))
    
    demo_file_folder <- tryCatch(
      resolve_demo_folder(
        ngs_type = ngs_type,
        version  = version,
        tool     = tool
      ),
      error = function(e) {
        shinyalert("Error", e$message, type = "error")
        return(NULL)
      }
    )
    
    if (is.null(demo_file_folder)) return()
    
    target_dir <- tryCatch(
      build_demo_target_dir(
        ngs_type = ngs_type,
        demo_file_folder = demo_file_folder
      ),
      error = function(e) {
        shinyalert("Error", e$message, type = "error")
        return(NULL)
      }
    )
    
    if (is.null(target_dir)) return()
    
    start_time <- Sys.time()
    
    w <- Waiter$new(
      html = spin_loader(),
      color = transparent(0.8)
    )
    w$show()
    on.exit(w$hide(), add = TRUE)
    
    Sys.sleep(5)
    
    demo_archive <- tryCatch(
      archive_demo_result_to_history(
        target_dir = target_dir,
        demo_job_id = demo_file_folder
      ),
      error = function(e) {
        shinyalert("Error", e$message, type = "error")
        return(NULL)
      }
    )
    
    if (is.null(demo_archive)) return()
    
    current_job_id(demo_file_folder)
    
    ok <- load_result_from_dir(
      target_dir = demo_archive$history_dir,
      display_id = demo_file_folder,
      source_type = "demo",
      show_success_alert = FALSE,
      id_box_label = "JobID",
      mode_value = "run",
      use_waiter = FALSE
    )
    
    if (!ok) return()
    
    end_time <- Sys.time()
    execution_time <- difftime(end_time, start_time, units = "secs")
    execution_time <- round(as.numeric(execution_time), 0)
    
    append_execution_log(
      job_id = demo_file_folder,
      ngs_type = ngs_type,
      version = version,
      tool = tool,
      execution_seconds = execution_time
    )
    
    return()
  }
  # ----------------------------------------------------------
  # Normal mode:
  # ----------------------------------------------------------
  HISTORY_ROOT <- "/root/shiny/History"
  dir.create(HISTORY_ROOT, recursive = TRUE, showWarnings = FALSE)
  
  ori_dir <- getwd()
  WORK_ROOT <- tempdir()
  
  setwd(WORK_ROOT)
  
  new_id <- paste0(format(Sys.time(), "%Y%m%d"), "-", sample(1000:9999, 1))
  current_job_id(new_id)
  
  history_dir <- file.path(HISTORY_ROOT, new_id)
  dir.create(history_dir, recursive = TRUE, showWarnings = FALSE)
  
  #----------------------------------------------------------
  input_dir <- file.path(WORK_ROOT, "Input")
  input_files <- list.files(input_dir, full.names = TRUE)
  
  if (length(input_files) == 0) {
    shinyalert(
      title = "Error",
      text = "No input files found. Please upload files before submitting.",
      type = "error"
    )
    setwd(ori_dir)
    return()
  }
  #----------------------------------------------------------
  start_time <- Sys.time()
  #----------------------------------------------------------
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8)
  )
  w$show()
  #----------------------------------------------------------
  dir.create("./Output", showWarnings = FALSE)
  dir.create("./Download", showWarnings = FALSE)
  #----------------------------------------------------------
  if (input$package == "OptiType") {
    if (dir.exists("./Output/optitype")) {
      unlink("./Output/optitype", recursive = TRUE)
    }
    dir.create("./Output/optitype", showWarnings = FALSE)
    #----------------------------------------------------------
    if (input$sequence == "WES") {
      system("OptiTypePipeline.py -i ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -d -p sample -o ./Output/optitype")
    } else if (input$sequence == "RNA-seq") {
      system("OptiTypePipeline.py -i ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -r -p sample -o ./Output/optitype")
    }
    #----------------------------------------------------------
    optitype <- read.table(file = "./Output/optitype/sample_result.tsv", sep = "\t", header = TRUE)
    # whether the dataframe is empty
    if (nrow(optitype) == 0) {
      output$hla_typing_table <- renderDataTable({
        datatable(
          data.frame(Message = "HLA typing unsuccessful: Insufficient sequencing reads detected."),
          class = "nowrap"
        )
      })
      
      # Define file and directory paths
      output_file <- paste("./Output/optitype_", input$sequence, "_", input$imgthla, ".txt", sep = "")
      zip_folder_pathway <- paste("./Download/optitype_", input$sequence, "_", input$imgthla, "_zip", sep = "")
      output_file_zip_pathway <- paste("Download/optitype_", input$sequence, "_", input$imgthla, "_output.zip", sep = "")
      merged_file <- paste("./Pivotable/optitype_", input$sequence, "_", input$imgthla, "_merged.csv", sep = "")
      # igv_folder <- "./IGV/optitype"
      
      # Remove the previous successful typing file (output_file) if it exists
      if (file.exists(output_file)) {
        file.remove(output_file)
      }
      
      # Remove the previous typing folder (zip_folder_pathway) if it exists
      if (dir.exists(zip_folder_pathway)) {
        unlink(zip_folder_pathway, recursive = TRUE)
      }
      
      # Remove the previous zip file (output_file_zip_pathway) if it exists
      if (file.exists(output_file_zip_pathway)) {
        file.remove(output_file_zip_pathway)
      }
      
      # Remove the previous successful typing file (merged_file) if it exists
      if (file.exists(merged_file)) {
        file.remove(merged_file)
      }
      
      # # Remove the previous successful typing folder (igv_folder) if it exists
      # if (dir.exists(igv_folder)) {
      #   unlink(igv_folder, recursive = TRUE)
      # }
    } else {
      optitype_mhc_table <- data.frame(
        Allele = character(6),
        nucleotide = character(6),
        protein = character(6),
        stringsAsFactors = FALSE
      )
      
      optitype_mhc_table[1, 1] <- if (is.na(optitype[1, 2]) || optitype[1, 2] == "") "" else optitype[1, 2]
      optitype_mhc_table[2, 1] <- if (is.na(optitype[1, 3]) || optitype[1, 3] == "") "" else optitype[1, 3]
      optitype_mhc_table[3, 1] <- if (is.na(optitype[1, 4]) || optitype[1, 4] == "") "" else optitype[1, 4]
      optitype_mhc_table[4, 1] <- if (is.na(optitype[1, 5]) || optitype[1, 5] == "") "" else optitype[1, 5]
      optitype_mhc_table[5, 1] <- if (is.na(optitype[1, 6]) || optitype[1, 6] == "") "" else optitype[1, 6]
      optitype_mhc_table[6, 1] <- if (is.na(optitype[1, 7]) || optitype[1, 7] == "") "" else optitype[1, 7]

      optitype_mhc_table <- optitype_mhc_table[order(optitype_mhc_table$Allele), , drop = FALSE]
      optitype_mhc_table <- optitype_mhc_table %>% filter(Allele != "")
      
      #----------------------------------------------------------
      for (i in 1:nrow(optitype_mhc_table)) {
        imgthla_version <- input$imgthla
        seq <- optitype_mhc_table[i, 1]
        allele_type <- substr(seq, 1, 1)
        allele_number <- substr(seq, 3, nchar(seq))
        IMGTHLA_nuc <- switch(allele_type,
                              "A" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/A_nuc.fasta", sep = ""),
                              "B" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/B_nuc.fasta", sep = ""),
                              "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_nuc.fasta", sep = "")
        )
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        }
        optitype_mhc_table[i, 2] <- system(cmd, intern = TRUE)
      }
      #----------------------------------------------------------
      for (i in 1:nrow(optitype_mhc_table)) {
        imgthla_version <- input$imgthla
        seq <- optitype_mhc_table[i, 1]
        allele_type <- substr(seq, 1, 1)
        allele_number <- substr(seq, 3, nchar(seq))
        IMGTHLA_prot <- switch(allele_type,
                               "A" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/A_prot.fasta", sep = ""),
                               "B" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/B_prot.fasta", sep = ""),
                               "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_prot.fasta", sep = "")
        )
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_prot, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_prot, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_prot, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_prot, " - | sed -n '2p' -")
        }
        optitype_mhc_table[i, 3] <- system(cmd, intern = TRUE)
      }
      #----------------------------------------------------------
      optitype_mhc_table <- optitype_mhc_table[order(optitype_mhc_table$Allele), , drop = FALSE]
      optitype_mhc_table <- optitype_mhc_table %>% filter(Allele != "")
      #----------------------------------------------------------
      output$hla_typing_table <- renderDataTable({
        render_hla_typing_datatable(optitype_mhc_table)
      })
      #----------------------------------------------------------
      sequence_file <- NULL
      for (i in seq_len(nrow(optitype_mhc_table))) {
        allele_i <- optitype_mhc_table$Allele[i]
        nucleotide_i <- optitype_mhc_table$nucleotide[i]
        protein_i <- optitype_mhc_table$protein[i]
        
        nucleotide_header <- paste0(">", allele_i, "_nucleotide")
        protein_header <- paste0(">", allele_i, "_protein")
        
        seqs <- rbind(
          c(nucleotide_header),
          c(nucleotide_i),
          c(protein_header),
          c(protein_i)
        )
        
        sequence_file <- rbind(sequence_file, seqs)
      }
      output_file <- paste("./Output/optitype_", input$sequence, "_", input$imgthla, ".txt", sep = "")
      write(sequence_file, file = output_file)
      #----------------------------------------------------------
      zip_folder_pathway <- file.path(WORK_ROOT, "Download", paste0("optitype_", input$sequence, "_", input$imgthla, "_zip"))
      dir.create(zip_folder_pathway, recursive = TRUE, showWarnings = FALSE)
      #----------------------------------------------------------
      lines <- readLines(output_file)
      dna_lines <- c()
      protein_lines <- c()
      #----------------------------------------------------------
      for (i in seq(1, length(lines), by = 2)) {
        header <- lines[i]
        sequence <- lines[i + 1]
        
        if (grepl("_nucleotide", header)) {
          dna_lines <- c(dna_lines, header, sequence)
        } else if (grepl("_protein", header)) {
          protein_lines <- c(protein_lines, header, sequence)
        }
      }
      #----------------------------------------------------------
      dna_file <- file.path(zip_folder_pathway, "nucleotide.fasta")
      protein_file <- file.path(zip_folder_pathway, "protein.fasta")
      writeLines(dna_lines, dna_file)
      writeLines(protein_lines, protein_file)
      #----------------------------------------------------------
      setwd(zip_folder_pathway)
      output_file_zip <- zip(file.path(WORK_ROOT, paste("Download/optitype_", input$sequence, "_", input$imgthla, "_output.zip", sep = "")), files = ".")
      output_file_zip_pathway <- paste("Download/optitype_", input$sequence, "_", input$imgthla, "_output.zip", sep = "")
      setwd(WORK_ROOT)
      #----------------------------------------------------------
      zip_name <- paste0("optitype_", input$sequence, "_", input$imgthla, "_output.zip")
      zip_src <- file.path(WORK_ROOT, "Download", zip_name)
      
      if (!file.exists(zip_src)) {
        w$hide()
        shinyalert("Error", paste("ZIP not generated:", zip_src), type = "error")
        setwd(ori_dir)
        return()
      }
      
      current_tmp_zip(normalizePath(zip_src, winslash = "/", mustWork = FALSE))
      current_mode("run")
      
      history_zip <- file.path(history_dir, basename(zip_src))
      ok <- file.copy(zip_src, history_zip, overwrite = TRUE)
      
      if (!ok || !file.exists(history_zip)) {
        w$hide()
        shinyalert("Error", "Failed to archive ZIP into /root/shiny/History.", type = "error")
        setwd(ori_dir)
        return()
      }
      
      current_history_zip(normalizePath(history_zip, winslash = "/", mustWork = FALSE))
      current_prefix(basename(zip_src))
      
      saveRDS(list(
        table = optitype_mhc_table,
        ngs_type = input$sequence,
        version = input$imgthla,
        tool = input$package,
        job_id = new_id,
        zip_path = normalizePath(history_zip, winslash = "/", mustWork = FALSE)
      ), file = file.path(history_dir, "result_data.rds"))
    }
    #----------------------------------------------------------
  } else if (input$package == "arcasHLA") {
    if (dir.exists("./Output/arcashla")) {
      unlink("./Output/arcashla", recursive = TRUE)
    }
    dir.create("./Output/arcashla", showWarnings = FALSE)
    #----------------------------------------------------------
    if (input$imgthla == "v3.14.0") {
      system(paste("/arcasHLA/arcasHLA genotype14 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.49.0") {
      system(paste("/arcasHLA/arcasHLA genotype49 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.50.0") {
      system(paste("/arcasHLA/arcasHLA genotype50 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.51.0") {
      system(paste("/arcasHLA/arcasHLA genotype51 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.52.0") {
      system(paste("/arcasHLA/arcasHLA genotype52 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.53.0") {
      system(paste("/arcasHLA/arcasHLA genotype53 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.54.0") {
      system(paste("/arcasHLA/arcasHLA genotype54 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.55.0") {
      system(paste("/arcasHLA/arcasHLA genotype55 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.56.0") {
      system(paste("/arcasHLA/arcasHLA genotype56 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.57.0") {
      system(paste("/arcasHLA/arcasHLA genotype57 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.58.0") {
      system(paste("/arcasHLA/arcasHLA genotype58 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.59.0") {
      system(paste("/arcasHLA/arcasHLA genotype59 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    } else if (input$imgthla == "v3.60.0") {
      system(paste("/arcasHLA/arcasHLA genotype60 ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v", sep = ""))
    }
    #----------------------------------------------------------
    arcashla_genotype_path <- paste("./Output/arcashla/sample_1.genotype.json", sep = "")
    arcashla_genotype <- fromJSON(file = arcashla_genotype_path)
    #----------------------------------------------------------
    # whether the dataframe is empty
    if (length(arcashla_genotype) == 0) {
      output$hla_typing_table <- renderDataTable({
        datatable(
          data.frame(Message = "HLA typing unsuccessful: Insufficient sequencing reads detected."),
          class = "nowrap"
        )
      })
      
      # Define file and directory paths
      output_file <- paste("./Output/arcashla_", input$sequence, "_", input$imgthla, ".txt", sep = "")
      zip_folder_pathway <- paste("./Download/arcashla_", input$sequence, "_", input$imgthla, "_zip", sep = "")
      output_file_zip_pathway <- paste("Download/arcashla_", input$sequence, "_", input$imgthla, "_output.zip", sep = "")
      merged_file <- paste("./Pivotable/arcashla_", input$sequence, "_", input$imgthla, "_merged.csv", sep = "")
      
      # Remove the previous successful typing file (output_file) if it exists
      if (file.exists(output_file)) {
        file.remove(output_file)
      }
      
      # Remove the previous typing folder (zip_folder_pathway) if it exists
      if (dir.exists(zip_folder_pathway)) {
        unlink(zip_folder_pathway, recursive = TRUE)
      }
      
      # Remove the previous zip file (output_file_zip_pathway) if it exists
      if (file.exists(output_file_zip_pathway)) {
        file.remove(output_file_zip_pathway)
      }
      
      # Remove the previous successful typing file (merged_file) if it exists
      if (file.exists(merged_file)) {
        file.remove(merged_file)
      }
    } else {
      arcashla_mhc_A <- arcashla_genotype[["A"]]
      arcashla_mhc_B <- arcashla_genotype[["B"]]
      arcashla_mhc_C <- arcashla_genotype[["C"]]
      arcashla_mhc_DPB1 <- arcashla_genotype[["DPB1"]]
      arcashla_mhc_DQA1 <- arcashla_genotype[["DQA1"]]
      arcashla_mhc_DQB1 <- arcashla_genotype[["DQB1"]]
      arcashla_mhc_DRB1 <- arcashla_genotype[["DRB1"]]
      arcashla_mhc_all <- rbind(arcashla_mhc_A, arcashla_mhc_B, arcashla_mhc_C, arcashla_mhc_DPB1, arcashla_mhc_DQA1, arcashla_mhc_DQB1, arcashla_mhc_DRB1)
      arcashla_mhc_all <- as.data.frame(arcashla_mhc_all)
      #----------------------------------------------------------
      arcashla_stacked <- stack(arcashla_mhc_all)
      arcashla_stacked <- data.frame(allele = arcashla_stacked$values)
      arcashla_stacked <- arcashla_stacked[order(arcashla_stacked$allele), , drop = FALSE]
      if ("-" %in% arcashla_stacked$allele) {
        arcashla_stacked <- subset(arcashla_stacked, allele != "-")
      }
      arcashla_mhc_table <- data.frame(
        Allele = character(nrow(arcashla_stacked)),
        nucleotide = character(nrow(arcashla_stacked)),
        protein = character(nrow(arcashla_stacked))
      )
      for (i in 1:nrow(arcashla_stacked)) {
        arcashla_mhc_table[i, 1] <- arcashla_stacked[i, ]
      }
      #---------------------------------------------------------- sorting
      arcashla_mhc_table <- arcashla_mhc_table[order(arcashla_mhc_table$Allele), , drop = FALSE]
      #----------------------------------------------------------
      imgthla_version <- input$imgthla
      for (i in 1:nrow(arcashla_mhc_table)) {
        seq <- arcashla_mhc_table[i, 1]
        # Determine whether it is a type A, B, or C gene
        if (grepl("^[ABC]\\*", seq)) {
          allele_type <- substr(seq, 1, 1)
          allele_number <- sub("^[ABC]\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_nuc <- switch(allele_type,
                                "A" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/A_nuc.fasta", sep = ""),
                                "B" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/B_nuc.fasta", sep = ""),
                                "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_nuc.fasta", sep = "")
          )
        } else {
          # Other types of genes
          allele_type <- substr(seq, 1, 4)
          allele_number <- sub("^.{4}\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_nuc <- switch(allele_type,
                                "DPB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPB1_nuc.fasta", sep = ""),
                                "DQA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQA1_nuc.fasta", sep = ""),
                                "DQB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQB1_nuc.fasta", sep = ""),
                                "DRB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DRB1_nuc.fasta", sep = "")
          )
        }
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        }
        result <- system(cmd, intern = TRUE)
        if (length(result) == 0) {
          warning(paste("No result found for allele_number:", allele_number, "in sequence:", seq))
          arcashla_mhc_table[i, 2] <- NA # or some other appropriate default value or handling
        } else {
          arcashla_mhc_table[i, 2] <- result
        }
      }
      #----------------------------------------------------------
      imgthla_version <- input$imgthla
      for (i in 1:nrow(arcashla_mhc_table)) {
        seq <- arcashla_mhc_table[i, 1]
        # Determine whether it is a type A, B, or C gene
        if (grepl("^[ABC]\\*", seq)) {
          allele_type <- substr(seq, 1, 1)
          allele_number <- sub("^[ABC]\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_prot <- switch(allele_type,
                                 "A" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/A_prot.fasta", sep = ""),
                                 "B" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/B_prot.fasta", sep = ""),
                                 "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_prot.fasta", sep = "")
          )
        } else {
          # Other types of genes
          allele_type <- substr(seq, 1, 4)
          allele_number <- sub("^.{4}\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_prot <- switch(allele_type,
                                 "DPB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPB1_prot.fasta", sep = ""),
                                 "DQA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQA1_prot.fasta", sep = ""),
                                 "DQB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQB1_prot.fasta", sep = ""),
                                 "DRB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DRB1_prot.fasta", sep = "")
          )
        }
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_prot, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_prot, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_prot, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_prot, " - | sed -n '2p' -")
        }
        result <- system(cmd, intern = TRUE)
        if (length(result) == 0) {
          warning(paste("No result found for allele_number:", allele_number, "in sequence:", seq))
          arcashla_mhc_table[i, 3] <- NA
        } else {
          arcashla_mhc_table[i, 3] <- result
        }
      }
      #---------------------------------------------------------- sorting
      arcashla_mhc_table <- arcashla_mhc_table[order(arcashla_mhc_table$Allele), , drop = FALSE]
      arcashla_mhc_table <- arcashla_mhc_table[, c("Allele", "nucleotide", "protein"), drop = FALSE]
      #----------------------------------------------------------
      output$hla_typing_table <- renderDataTable({
        render_hla_typing_datatable(arcashla_mhc_table)
      })
      #----------------------------------------------------------
      sequence_file <- NULL
      for (i in seq_len(nrow(arcashla_mhc_table))) {
        allele_i <- arcashla_mhc_table$Allele[i]
        nucleotide_i <- arcashla_mhc_table$nucleotide[i]
        protein_i <- arcashla_mhc_table$protein[i]
        
        nucleotide_header <- paste0(">", allele_i, "_nucleotide")
        protein_header <- paste0(">", allele_i, "_protein")
        
        seqs <- rbind(
          c(nucleotide_header),
          c(nucleotide_i),
          c(protein_header),
          c(protein_i)
        )
        
        sequence_file <- rbind(sequence_file, seqs)
      }
      output_file <- paste("./Output/arcashla_", input$sequence, "_", input$imgthla, ".txt", sep = "")
      write(sequence_file, file = output_file)
      #----------------------------------------------------------
      zip_folder_pathway <- file.path(WORK_ROOT, "Download", paste0("arcashla_", input$sequence, "_", input$imgthla, "_zip"))
      dir.create(zip_folder_pathway, recursive = TRUE, showWarnings = FALSE)
      #----------------------------------------------------------
      lines <- readLines(output_file)
      dna_lines <- c()
      protein_lines <- c()
      #----------------------------------------------------------
      for (i in seq(1, length(lines), by = 2)) {
        header <- lines[i]
        sequence <- lines[i + 1]
        if (grepl("_nucleotide", header)) {
          dna_lines <- c(dna_lines, header, sequence)
        } else if (grepl("_protein", header)) {
          protein_lines <- c(protein_lines, header, sequence)
        }
      }
      dna_file <- file.path(zip_folder_pathway, "nucleotide.fasta")
      protein_file <- file.path(zip_folder_pathway, "protein.fasta")
      writeLines(dna_lines, dna_file)
      writeLines(protein_lines, protein_file)
      #----------------------------------------------------------
      setwd(zip_folder_pathway)
      output_file_zip <- zip(file.path(WORK_ROOT, paste("Download/arcashla_", input$sequence, "_", input$imgthla, "_output.zip", sep = "")), files = ".")
      output_file_zip_pathway <- paste("Download/arcashla_", input$sequence, "_", input$imgthla, "_output.zip", sep = "")
      setwd(WORK_ROOT)
      #----------------------------------------------------------
      zip_name <- paste0("arcashla_", input$sequence, "_", input$imgthla, "_output.zip")
      zip_src <- file.path(WORK_ROOT, "Download", zip_name)
      
      if (!file.exists(zip_src)) {
        w$hide()
        shinyalert("Error", paste("ZIP not generated:", zip_src), type = "error")
        setwd(ori_dir)
        return()
      }
      
      current_tmp_zip(normalizePath(zip_src, winslash = "/", mustWork = FALSE))
      current_mode("run")
      
      history_zip <- file.path(history_dir, basename(zip_src))
      ok <- file.copy(zip_src, history_zip, overwrite = TRUE)
      
      if (!ok || !file.exists(history_zip)) {
        w$hide()
        shinyalert("Error", "Failed to archive ZIP into /root/shiny/History.", type = "error")
        setwd(ori_dir)
        return()
      }
      
      current_history_zip(normalizePath(history_zip, winslash = "/", mustWork = FALSE))
      current_prefix(basename(history_zip))
      
      saveRDS(list(
        table = arcashla_mhc_table,
        ngs_type = input$sequence,
        version = input$imgthla,
        tool = input$package,
        job_id = new_id,
        zip_path = normalizePath(history_zip, winslash = "/", mustWork = FALSE)
      ), file = file.path(history_dir, "result_data.rds"))
    }
    #----------------------------------------------------------
  } else if (input$package == "HLA-HD") {
    if (dir.exists("./Output/hlahd")) {
      unlink("./Output/hlahd", recursive = TRUE)
    }
    dir.create("./Output/hlahd", showWarnings = FALSE)
    #----------------------------------------------------------
    if (input$imgthla == "v3.14.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary14 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.49.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary49 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.50.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary50 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.51.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary51 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.52.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary52 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.53.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary53 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.54.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary54 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.55.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary55 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.56.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary56 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.57.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary57 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.58.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary58 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.59.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary59 sample ./Output/hlahd", sep = ""))
    } else if (input$imgthla == "v3.60.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample_1.fastq.gz ./Input/sample_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary60 sample ./Output/hlahd", sep = ""))
    }
    #----------------------------------------------------------
    if (input$imgthla == "v3.14.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary14/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.49.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary49/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.50.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary50/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.51.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary51/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.52.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary52/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.53.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary53/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.54.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary54/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.55.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary55/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.56.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary56/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.57.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary57/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.58.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary58/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.59.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary59/all_exon_N150.fasta"
    } else if (input$imgthla == "v3.60.0") {
      reference_fasta <- "/hlahd.1.7.0/IMGTHLA_version/dictionary60/all_exon_N150.fasta"
    } else {
      stop("Unknown IMGTHLA version")
    }
    #-----------------------------------------------------------
    system("head -n 8 ./Output/hlahd/sample/result/sample_final.result.txt > ./Output/hlahd/sample_final_mhc.result.txt")
    
    # Read the file and remove lines containing 'Couldn't read result file.'
    lines <- readLines("./Output/hlahd/sample_final_mhc.result.txt")
    
    # Filter out lines containing 'Couldn't read result file.'
    valid_lines <- lines[!grepl("Couldn't read result file.", lines)]
    
    # Write the filtered content to a new temporary file
    temp_file <- "./Output/hlahd/sample_final_mhc_clean.result.txt"
    writeLines(valid_lines, temp_file)
    
    # Read the cleaned file
    hlahd <- read.table(temp_file, header = FALSE)
    hlahd <- hlahd[, -1]
    hlahd$V2 <- sub("HLA-", "", hlahd$V2)
    hlahd$V3 <- sub("HLA-", "", hlahd$V3)
    
    # Remove intermediate files
    file.remove(temp_file)
    #-----------------------------------------------------------
    hlahd_stacked <- stack(hlahd)
    hlahd_stacked <- data.frame(allele = hlahd_stacked$values)
    hlahd_stacked <- hlahd_stacked[order(hlahd_stacked$allele), , drop = FALSE]
    # remove "-", "Not" and "typed"
    values_to_remove <- c("-", "Not", "typed")
    hlahd_stacked <- subset(hlahd_stacked, !allele %in% values_to_remove)
    #----------------------------------------------------------
    # whether the dataframe is empty
    if (nrow(hlahd_stacked) == 0) {
      output$hla_typing_table <- renderDataTable({
        datatable(
          data.frame(Message = "HLA typing unsuccessful: Insufficient sequencing reads detected."),
          class = "nowrap"
        )
      })
      
      # Define file and directory paths
      output_file <- paste("./Output/hlahd_", input$sequence, "_", input$imgthla, ".txt", sep = "")
      zip_folder_pathway <- paste("./Download/hlahd_", input$sequence, "_", input$imgthla, "_zip", sep = "")
      output_file_zip_pathway <- paste("Download/hlahd_", input$sequence, "_", input$imgthla, "_output.zip", sep = "")
      merged_file <- paste("./Pivotable/hlahd_", input$sequence, "_", input$imgthla, "_merged.csv", sep = "")
      # igv_folder <- "./IGV/hlahd"
      
      # Remove the previous successful typing file (output_file) if it exists
      if (file.exists(output_file)) {
        file.remove(output_file)
      }
      
      # Remove the previous typing folder (zip_folder_pathway) if it exists
      if (dir.exists(zip_folder_pathway)) {
        unlink(zip_folder_pathway, recursive = TRUE)
      }
      
      # Remove the previous zip file (output_file_zip_pathway) if it exists
      if (file.exists(output_file_zip_pathway)) {
        file.remove(output_file_zip_pathway)
      }
      
      # Remove the previous successful typing file (merged_file) if it exists
      if (file.exists(merged_file)) {
        file.remove(merged_file)
      }
      
    } else {
      hlahd_mhc_table <- data.frame(
        Allele = character(nrow(hlahd_stacked)),
        nucleotide = character(nrow(hlahd_stacked)),
        protein = character(nrow(hlahd_stacked))
      )
      for (i in 1:nrow(hlahd_stacked)) {
        hlahd_mhc_table[i, 1] <- hlahd_stacked[i, ]
      }
      # distinct remove multiple DQB1
      hlahd_mhc_table_unique <- hlahd_mhc_table %>%
        filter(grepl("^DQB1", Allele)) %>%
        distinct(Allele)
      
      # unique DQB1 add to table
      hlahd_mhc_table <- hlahd_mhc_table %>%
        filter(!grepl("^DQB1", Allele)) %>%
        bind_rows(hlahd_mhc_table_unique)
      # sorting
      hlahd_mhc_table <- hlahd_mhc_table[order(hlahd_mhc_table$Allele), , drop = FALSE]
      #----------------------------------------------------------
      imgthla_version <- input$imgthla
      for (i in 1:nrow(hlahd_mhc_table)) {
        seq <- hlahd_mhc_table[i, 1]
        # Determine whether it is a type A, B, or C gene
        if (grepl("^[ABC]\\*", seq)) {
          allele_type <- substr(seq, 1, 1)
          allele_number <- sub("^[ABC]\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_nuc <- switch(allele_type,
                                "A" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/A_nuc.fasta", sep = ""),
                                "B" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/B_nuc.fasta", sep = ""),
                                "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_nuc.fasta", sep = "")
          )
        } else {
          # Other types of genes
          allele_type <- substr(seq, 1, 4)
          allele_number <- sub("^.{4}\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_nuc <- switch(allele_type,
                                "DPA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPA1_nuc.fasta", sep = ""),
                                "DPB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPB1_nuc.fasta", sep = ""),
                                "DQA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQA1_nuc.fasta", sep = ""),
                                "DQB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQB1_nuc.fasta", sep = ""),
                                "DRB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DRB1_nuc.fasta", sep = "")
          )
        }
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        }
        result <- system(cmd, intern = TRUE)
        if (length(result) == 0) {
          warning(paste("No result found for allele_number:", allele_number, "in sequence:", seq))
          hlahd_mhc_table[i, 2] <- NA # or some other appropriate default value or handling
        } else {
          hlahd_mhc_table[i, 2] <- result
        }
      }
      #----------------------------------------------------------
      imgthla_version <- input$imgthla
      for (i in 1:nrow(hlahd_mhc_table)) {
        seq <- hlahd_mhc_table[i, 1]
        # Determine whether it is a type A, B, or C gene
        if (grepl("^[ABC]\\*", seq)) {
          allele_type <- substr(seq, 1, 1)
          allele_number <- sub("^[ABC]\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_prot <- switch(allele_type,
                                 "A" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/A_prot.fasta", sep = ""),
                                 "B" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/B_prot.fasta", sep = ""),
                                 "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_prot.fasta", sep = "")
          )
        } else {
          # Other types of genes
          allele_type <- substr(seq, 1, 4)
          allele_number <- sub("^.{4}\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_prot <- switch(allele_type,
                                 "DPA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPA1_prot.fasta", sep = ""),
                                 "DPB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPB1_prot.fasta", sep = ""),
                                 "DQA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQA1_prot.fasta", sep = ""),
                                 "DQB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQB1_prot.fasta", sep = ""),
                                 "DRB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DRB1_prot.fasta", sep = "")
          )
        }
        
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_prot, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_prot, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_prot, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_prot, " - | sed -n '2p' -")
        }
        
        result <- system(cmd, intern = TRUE)
        
        if (length(result) == 0) {
          warning(paste("No result found for allele_number:", allele_number, "in sequence:", seq))
          hlahd_mhc_table[i, 3] <- NA # or some other appropriate default value or handling
        } else {
          hlahd_mhc_table[i, 3] <- result
        }
      }
      #----------------------------------------------------------
      hlahd_mhc_table <- hlahd_mhc_table[, c("Allele", "nucleotide", "protein"), drop = FALSE]
      #----------------------------------------------------------
      output$hla_typing_table <- renderDataTable({
        render_hla_typing_datatable(hlahd_mhc_table)
      })
      #----------------------------------------------------------
      sequence_file <- NULL
      for (i in seq_len(nrow(hlahd_mhc_table))) {
        allele_i <- hlahd_mhc_table$Allele[i]
        nucleotide_i <- hlahd_mhc_table$nucleotide[i]
        protein_i <- hlahd_mhc_table$protein[i]
        
        nucleotide_header <- paste0(">", allele_i, "_nucleotide")
        protein_header <- paste0(">", allele_i, "_protein")
        
        seqs <- rbind(
          c(nucleotide_header),
          c(nucleotide_i),
          c(protein_header),
          c(protein_i)
        )
        
        sequence_file <- rbind(sequence_file, seqs)
      }
      output_file <- paste("./Output/hlahd_", input$sequence, "_", input$imgthla, ".txt", sep = "")
      write(sequence_file, file = output_file)
      #----------------------------------------------------------
      zip_folder_pathway <- file.path(WORK_ROOT, "Download", paste0("hlahd_", input$sequence, "_", input$imgthla, "_zip"))
      dir.create(zip_folder_pathway, recursive = TRUE, showWarnings = FALSE)
      #----------------------------------------------------------
      lines <- readLines(output_file)
      dna_lines <- c()
      protein_lines <- c()
      for (i in seq(1, length(lines), by = 2)) {
        header <- lines[i]
        sequence <- lines[i + 1]
        if (grepl("_nucleotide", header)) {
          dna_lines <- c(dna_lines, header, sequence)
        } else if (grepl("_protein", header)) {
          protein_lines <- c(protein_lines, header, sequence)
        }
      }
      #----------------------------------------------------------
      dna_file <- file.path(zip_folder_pathway, "nucleotide.fasta")
      protein_file <- file.path(zip_folder_pathway, "protein.fasta")
      writeLines(dna_lines, dna_file)
      writeLines(protein_lines, protein_file)
      #----------------------------------------------------------
      setwd(zip_folder_pathway)
      output_file_zip <- zip(file.path(WORK_ROOT, paste("Download/hlahd_", input$sequence, "_", input$imgthla, "_output.zip", sep = "")), files = ".")
      output_file_zip_pathway <- paste("Download/hlahd_", input$sequence, "_", input$imgthla, "_output.zip", sep = "")
      setwd(WORK_ROOT)
      #----------------------------------------------------------
      zip_name <- paste0("hlahd_", input$sequence, "_", input$imgthla, "_output.zip")
      zip_src <- file.path(WORK_ROOT, "Download", zip_name)
      
      if (!file.exists(zip_src)) {
        w$hide()
        shinyalert("Error", paste("ZIP not generated:", zip_src), type = "error")
        setwd(ori_dir)
        return()
      }
      
      current_tmp_zip(normalizePath(zip_src, winslash = "/", mustWork = FALSE))
      current_mode("run")
      
      history_zip <- file.path(history_dir, basename(zip_src))
      ok <- file.copy(zip_src, history_zip, overwrite = TRUE)
      
      if (!ok || !file.exists(history_zip)) {
        w$hide()
        shinyalert("Error", "Failed to archive ZIP into /root/shiny/History.", type = "error")
        setwd(ori_dir)
        return()
      }
      
      current_history_zip(normalizePath(history_zip, winslash = "/", mustWork = FALSE))
      current_prefix(basename(history_zip))
      
      saveRDS(list(
        table = hlahd_mhc_table,
        ngs_type = input$sequence,
        version = input$imgthla,
        tool = input$package,
        job_id = new_id,
        zip_path = normalizePath(history_zip, winslash = "/", mustWork = FALSE)
      ), file = file.path(history_dir, "result_data.rds"))
    }
  } else if (input$package == "T1K") {
    convert_fa_to_fq <- function(fa_path, fq_path) {
      awk_cmd <- sprintf(
        "awk 'BEGIN{OFS=\"\\n\"} /^>/ {h=$0; next} {qual=\"\"; for(i=1;i<=length($0);i++) qual=qual\"I\"; print h, $0, \"+\", qual}' %s > %s",
        shQuote(fa_path), shQuote(fq_path)
      )
      system(awk_cmd)
    }
    #----------------------------------------------------------
    seq_suffix <- ifelse(input$sequence == "WES", "dna", "rna")
    imgt_version_number <- gsub("\\.", "", sub("^v", "", input$imgthla))
    imgt_dir <- file.path("/T1K/IMGTHLA", input$imgthla)
    origin_fa <- file.path(imgt_dir, paste0("IMGTHLA_", imgt_version_number, "_", seq_suffix, "_seq.fa"))
    #----------------------------------------------------------
    base_output_dir <- file.path(WORK_ROOT, "Output", "t1k")
    # igv_out_dir <- file.path(WORK_ROOT, "IGV", "t1k")
    
    if (dir.exists(base_output_dir)) unlink(base_output_dir, recursive = TRUE)
    dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # if (dir.exists(igv_out_dir)) unlink(igv_out_dir, recursive = TRUE)
    # dir.create(igv_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    output_dir <- file.path(base_output_dir, "t1k")
    #----------------------------------------------------------
    cmd <- paste(
      "perl /T1K/run-t1k",
      "-t", num_cores,
      "--preset hla",
      "-1 ./Input/sample_1.fastq.gz",
      "-2 ./Input/sample_2.fastq.gz",
      "-f", origin_fa,
      "-o", output_dir
    )
    system(cmd)
    #----------------------------------------------------------
    t1k_raw <- read.delim("./Output/t1k/t1k_genotype.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    keep_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1")
    
    extract_alleles <- function(row) {
      gene <- row[[1]]
      if (!gene %in% keep_genes) {
        return(NULL)
      }
      
      alleles <- c()
      if (row[[3]] != "." && row[[3]] != "0" && row[[3]] != "-") alleles <- c(alleles, row[[3]])
      if (row[[6]] != "." && row[[6]] != "0" && row[[6]] != "-") alleles <- c(alleles, row[[6]])
      
      if (length(row) >= 9 && grepl("\\*", row[[9]])) {
        third <- unlist(strsplit(row[[9]], ";"))
        if (length(third) >= 1 && third[1] != "." && third[1] != "-") {
          alleles <- c(alleles, third[1])
        }
      }
      
      if (length(alleles) == 0) {
        return(NULL)
      }
      
      data.frame(
        Allele = alleles,
        nucleotide = "",
        protein = "",
        stringsAsFactors = FALSE
      )
    }
    #----------------------------------------------------------
    t1k_mhc_table <- do.call(rbind, lapply(1:nrow(t1k_raw), function(i) extract_alleles(t1k_raw[i, ])))
    t1k_mhc_table <- subset(t1k_mhc_table, !(Allele %in% c("-", ".", "0")))
    t1k_mhc_table <- t1k_mhc_table[order(t1k_mhc_table$Allele), , drop = FALSE]
    #----------------------------------------------------------
    #---------------nucleotide
    imgthla_version <- input$imgthla
    for (i in 1:nrow(t1k_mhc_table)) {
      seq <- t1k_mhc_table[i, 1]
      
      gene_name <- sub("^HLA-([A-Z0-9]+)\\*.*", "\\1", seq)
      allele_number <- sub("^HLA-[A-Z0-9]+\\*([0-9:]+)", "\\1", seq)
      
      fasta_path <- file.path("/root/shiny/IMGTHLA", imgthla_version, "fasta", paste0(gene_name, "_nuc.fasta"))
      
      if (grepl("N", allele_number)) {
        cmd <- paste(
          "grep", allele_number, fasta_path,
          "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq",
          fasta_path, "-", "| sed -n '2p'"
        )
      } else {
        cmd <- paste(
          "grep", allele_number, fasta_path,
          "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq",
          fasta_path, "-", "| sed -n '2p'"
        )
      }
      
      result <- system(cmd, intern = TRUE)
      
      if (length(result) == 0) {
        warning(paste("No nucleotide found for:", seq))
        t1k_mhc_table[i, 2] <- NA
      } else {
        t1k_mhc_table[i, 2] <- result
      }
    }
    #---------------protein
    for (i in 1:nrow(t1k_mhc_table)) {
      seq <- t1k_mhc_table[i, 1]
      
      gene_name <- sub("^HLA-([A-Z0-9]+)\\*.*", "\\1", seq)
      allele_number <- sub("^HLA-[A-Z0-9]+\\*([0-9:]+)", "\\1", seq)
      
      fasta_path <- file.path("/root/shiny/IMGTHLA", imgthla_version, "fasta", paste0(gene_name, "_prot.fasta"))
      
      if (grepl("N", allele_number)) {
        cmd <- paste(
          "grep", allele_number, fasta_path,
          "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq",
          fasta_path, "-", "| sed -n '2p'"
        )
      } else {
        cmd <- paste(
          "grep", allele_number, fasta_path,
          "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq",
          fasta_path, "-", "| sed -n '2p'"
        )
      }
      
      result <- system(cmd, intern = TRUE)
      
      if (length(result) == 0) {
        warning(paste("No protein found for:", seq))
        t1k_mhc_table[i, 3] <- NA
      } else {
        t1k_mhc_table[i, 3] <- result
      }
    }
    #----------------------------------------------------------
    t1k_mhc_table <- t1k_mhc_table[, c("Allele", "nucleotide", "protein"), drop = FALSE]
    #----------------------------------------------------------
    output$hla_typing_table <- renderDataTable({
      render_hla_typing_datatable(t1k_mhc_table)
    })
    #----------------------------------------------------------
    sequence_file <- NULL
    
    for (i in seq_len(nrow(t1k_mhc_table))) {
      allele_i <- t1k_mhc_table$Allele[i]
      nucleotide_i <- t1k_mhc_table$nucleotide[i]
      protein_i <- t1k_mhc_table$protein[i]
      
      short_allele_i <- sub("^HLA-", "", allele_i)
      
      nucleotide_header <- paste0(">", short_allele_i, "_nucleotide")
      protein_header <- paste0(">", short_allele_i, "_protein")
      
      seqs <- rbind(
        c(nucleotide_header),
        c(nucleotide_i),
        c(protein_header),
        c(protein_i)
      )
      
      sequence_file <- rbind(sequence_file, seqs)
    }
    
    output_file <- file.path(WORK_ROOT, "Output", paste0("t1k_", input$sequence, "_", input$imgthla, ".txt"))
    write(sequence_file, file = output_file)
    #----------------------------------------------------------
    zip_folder_pathway <- file.path(WORK_ROOT, paste0("Download/t1k_", input$sequence, "_", input$imgthla, "_zip"))
    dir.create(zip_folder_pathway, recursive = TRUE, showWarnings = FALSE)
    
    output_file <- file.path(WORK_ROOT, paste0("Output/t1k_", input$sequence, "_", input$imgthla, ".txt"))
    lines <- readLines(output_file)
    
    dna_lines <- c()
    protein_lines <- c()
    
    for (i in seq(1, length(lines), by = 2)) {
      header <- lines[i]
      sequence <- lines[i + 1]
      
      if (grepl("_nucleotide", header)) {
        dna_lines <- c(dna_lines, header, sequence)
      } else if (grepl("_protein", header)) {
        protein_lines <- c(protein_lines, header, sequence)
      }
    }
    
    dna_file <- file.path(zip_folder_pathway, "nucleotide.fasta")
    protein_file <- file.path(zip_folder_pathway, "protein.fasta")
    writeLines(dna_lines, dna_file)
    writeLines(protein_lines, protein_file)
    #----------------------------------------------------------
    setwd(zip_folder_pathway)
    output_file_zip_pathway <- file.path(WORK_ROOT, paste0("Download/t1k_", input$sequence, "_", input$imgthla, "_output.zip"))
    zip(zipfile = output_file_zip_pathway, files = list.files(".", full.names = FALSE))
    setwd(WORK_ROOT)
    #----------------------------------------------------------
    zip_src <- output_file_zip_pathway
    history_zip <- file.path(history_dir, basename(zip_src))
    ok <- file.copy(zip_src, history_zip, overwrite = TRUE)
    
    if (!ok || !file.exists(history_zip)) {
      w$hide()
      shinyalert("Error", "Failed to archive ZIP into /root/shiny/History.", type = "error")
      setwd(ori_dir)
      return()
    }
    
    current_history_zip(normalizePath(history_zip, winslash = "/", mustWork = FALSE))
    current_prefix(basename(history_zip))
    
    saveRDS(list(
      table = t1k_mhc_table,
      ngs_type = input$sequence,
      version = input$imgthla,
      tool = input$package,
      job_id = new_id,
      zip_path = normalizePath(history_zip, winslash = "/", mustWork = FALSE)
    ), file = file.path(history_dir, "result_data.rds"))
  }
  #----------------------------------------------------------
  w$hide()
  #----------------------------------------------------------
  # Calculate execution time
  end_time <- Sys.time()
  execution_time <- difftime(end_time, start_time, units = "secs")
  execution_time <- round(as.numeric(execution_time), 0)
  
  append_execution_log(
    job_id = current_job_id(),
    ngs_type = input$sequence,
    version = input$imgthla,
    tool = input$package,
    execution_seconds = execution_time
  )
  #----------------------------------------------------------
  current_result_source("normal")
  #----------------------------------------------------------
  setwd(ori_dir)
})

#----------------------------------------------------------
# Download Handler (Global Scope)
#----------------------------------------------------------
resolve_download_zip <- function() {
  tmp <- current_tmp_zip()
  if (!is.null(tmp) && nzchar(tmp) && file.exists(tmp)) {
    return(tmp)
  }
  
  hz <- current_history_zip()
  if (!is.null(hz) && nzchar(hz) && file.exists(hz)) {
    return(hz)
  }
  
  jid <- current_job_id()
  if (!is.null(jid) && nzchar(jid)) {
    rds <- file.path("/root/shiny/History", jid, "result_data.rds")
    if (file.exists(rds)) {
      ld <- readRDS(rds)
      if (!is.null(ld$zip_path) && file.exists(ld$zip_path)) {
        return(ld$zip_path)
      }
    }
  }
  NULL
}

output$hla_typing_download <- downloadHandler(
  filename = function() {
    nm <- current_prefix()
    if (is.null(nm) || !nzchar(nm)) {
      paste0("HLA_typing_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    } else {
      nm
    }
  },
  content = function(file) {
    src <- resolve_download_zip()
    shiny::validate(need(
      !is.null(src) && file.exists(src),
      "No downloadable result is currently loaded. Please run or reload a job."
    ))
    file.copy(src, file, overwrite = TRUE)
  }
)
