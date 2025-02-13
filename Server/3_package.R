#----------------------------------------------------------
execution_times <- reactiveValues(times = character())
#----------------------------------------------------------
getCurrentTime <- function() {
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M")
  return(formatted_time)
}
#----------------------------------------------------------
observeEvent(input$hla_typing_button,{
  ori_dir <- getwd()
  tmp_dir <- tempdir()
  setwd(tmp_dir)
  
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
  dir.create("./IGV", showWarnings = FALSE)
  #----------------------------------------------------------
  if (input$package == "OptiType") {
    if (dir.exists("./Output/optitype")) {
      unlink("./Output/optitype", recursive = TRUE)
    }
    dir.create("./Output/optitype", showWarnings = FALSE)
    #----------------------------------------------------------
    if (input$sequence == "WES") {
      system("OptiTypePipeline.py -i ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -d -p sample1 -o ./Output/optitype")
    }else if (input$sequence == "RNA-seq") {
      system("OptiTypePipeline.py -i ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -r -p sample1 -o ./Output/optitype")
    }
    #----------------------------------------------------------
    if (dir.exists("./IGV/optitype")) {
      unlink("./IGV/optitype", recursive = TRUE)
    }
    dir.create("./IGV/optitype", showWarnings = FALSE)
    #----------------------------------------------------------
    # --- Merge alignment files ---
    # Sort _1.bam files if they exist
    bam1_files <- Sys.glob("./Output/optitype/*_1.bam")
    if (length(bam1_files) > 0) {
      sort_cmd1 <- "samtools sort -T ./IGV/optitype/temp_1 -o ./IGV/optitype/optitype_1_sorted.bam ./Output/optitype/*_1.bam"
      system(sort_cmd1)
    } else {
      message("[Warning] No *_1.bam files found. Skipping sort for *_1.bam.")
    }
    
    # Sort _2.bam files if they exist
    bam2_files <- Sys.glob("./Output/optitype/*_2.bam")
    if (length(bam2_files) > 0) {
      sort_cmd2 <- "samtools sort -T ./IGV/optitype/temp_2 -o ./IGV/optitype/optitype_2_sorted.bam ./Output/optitype/*_2.bam"
      system(sort_cmd2)
    } else {
      message("[Warning] No *_2.bam files found. Skipping sort for *_2.bam.")
    }
    
    # Merge sorted BAM files if both exist
    sorted_bam1 <- "./IGV/optitype/optitype_1_sorted.bam"
    sorted_bam2 <- "./IGV/optitype/optitype_2_sorted.bam"
    if (file.exists(sorted_bam1) && file.exists(sorted_bam2)) {
      merge_cmd <- paste("samtools merge ./IGV/optitype/optitype_merge.bam", sorted_bam1, sorted_bam2)
      system(merge_cmd)
    } else {
      message("[Warning] One or both sorted BAM files are missing. Merge not executed.")
    }
    
    # Remove intermediate files (only remove if they exist)
    if (length(bam1_files) > 0) file.remove(bam1_files)
    if (length(bam2_files) > 0) file.remove(bam2_files)
    if (file.exists(sorted_bam1)) file.remove(sorted_bam1)
    if (file.exists(sorted_bam2)) file.remove(sorted_bam2)
    
    # --- Extract the "main sequence" from the alignment file ---
    
    merged_bam <- "./IGV/optitype/optitype_merge.bam"
    merged_sam <- "./IGV/optitype/optitype_merge.sam"
    if (file.exists(merged_bam)) {
      view_cmd <- paste("samtools view -h", merged_bam, "> ./IGV/optitype/optitype_merge.sam")
      system(view_cmd)
    } else {
      message("[Warning] Merged BAM file does not exist. Cannot convert to SAM.")
    }
    
    filtered_sam <- "./IGV/optitype/optitype_merge_filtered.sam"
    if (file.exists(merged_sam)) {
      awk_cmd <- paste("awk '$2 == 0 || $1 ~ /^@/'", merged_sam, "> ./IGV/optitype/optitype_merge_filtered.sam")
      system(awk_cmd)
    } else {
      message("[Warning] SAM file not found. Skipping filtering step.")
    }
    
    filtered_bam <- "./IGV/optitype/optitype_merge_filtered.bam"
    if (file.exists(filtered_sam)) {
      conv_cmd <- paste("samtools view -Sb", filtered_sam, "> ./IGV/optitype/optitype_merge_filtered.bam")
      system(conv_cmd)
    } else {
      message("[Warning] Filtered SAM file not found. Cannot convert to BAM.")
    }
    
    if (file.exists(merged_bam)) file.remove(merged_bam)
    if (file.exists(merged_sam)) file.remove(merged_sam)
    if (file.exists(filtered_sam)) file.remove(filtered_sam)
    
    # --- Extract data based on the matched reference ---
    
    # Choose the reference FASTA file based on input$sequence
    if (input$sequence == "WES") {
      reference_fasta <- "/root/miniconda3/envs/optitype/share/optitype-1.3.2-3/data/hla_reference_dna.fasta"
    } else {
      reference_fasta <- "/root/miniconda3/envs/optitype/share/optitype-1.3.2-3/data/hla_reference_rna.fasta"
    }
    
    input_bam <- file.path(getwd(), "IGV/optitype/optitype_merge_filtered.bam")
    igv_fasta <- file.path(getwd(), "IGV/optitype/optitype_hla_reference_igv.fasta")
    
    # Execute the Python script if the input BAM exists
    if (file.exists(input_bam)) {
      py_cmd <- paste("python /root/shiny/Server/optitype_igv_genome.py", 
                      input_bam, reference_fasta, igv_fasta)
      system(py_cmd)
    } else {
      message("[Warning] Input BAM file does not exist. Python script not executed.")
    }
    
    # Create an index of the reference FASTA if it exists
    if (file.exists(igv_fasta)) {
      index_cmd <- paste("samtools faidx", igv_fasta)
      system(index_cmd)
    } else {
      message("[Warning] IGV FASTA file does not exist. Cannot create index.")
    }
    #----------------------------------------------------------
    optitype = read.table(file = "./Output/optitype/sample1_result.tsv", sep = '\t', header = TRUE)
    # whether the dataframe is empty
    if (nrow(optitype) == 0) {
      output$hla_typing_table <- renderDataTable({
        datatable(
          data.frame(Message = "HLA typing unsuccessful: Insufficient sequencing reads detected."),
          class = 'nowrap'
        )
      })
      
      # Define file and directory paths
      output_file = paste("./Output/optitype_",input$sequence,"_",input$imgthla,".txt",sep = "")
      zip_folder_pathway = paste("./Download/optitype_",input$sequence,"_",input$imgthla,"_zip",sep = "")
      output_file_zip_pathway <- paste("Download/optitype_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")
      merged_file = paste("./Pivotable/optitype_",input$sequence,"_",input$imgthla,"_merged.csv",sep = "")
      
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
      optitype_mhc_table <- data.frame(
        Allele = character(6),
        nucleotide = character(6),
        protein = character(6),
        stringsAsFactors = FALSE
      )
      
      optitype_mhc_table[1,1] <- if (is.na(optitype[1,2]) || optitype[1,2] == "") "" else optitype[1,2]
      optitype_mhc_table[2,1] <- if (is.na(optitype[1,3]) || optitype[1,3] == "") "" else optitype[1,3]
      optitype_mhc_table[3,1] <- if (is.na(optitype[1,4]) || optitype[1,4] == "") "" else optitype[1,4]
      optitype_mhc_table[4,1] <- if (is.na(optitype[1,5]) || optitype[1,5] == "") "" else optitype[1,5]
      optitype_mhc_table[5,1] <- if (is.na(optitype[1,6]) || optitype[1,6] == "") "" else optitype[1,6]
      optitype_mhc_table[6,1] <- if (is.na(optitype[1,7]) || optitype[1,7] == "") "" else optitype[1,7]
      
      optitype_mhc_table <- optitype_mhc_table[order(optitype_mhc_table$Allele), , drop = FALSE]
      optitype_mhc_table <- optitype_mhc_table %>% filter(Allele != "")
      #----------------------------------------------------------
      for (i in 1:nrow(optitype_mhc_table)) {
        imgthla_version = input$imgthla
        seq = optitype_mhc_table[i,1]
        allele_type = substr(seq,1,1)
        allele_number = substr(seq,3,nchar(seq))
        IMGTHLA_nuc <- switch(allele_type, "A"=paste("/root/shiny/IMGTHLA/",imgthla_version,"/fasta/A_nuc.fasta",sep = ""), "B"=paste("/root/shiny/IMGTHLA/",imgthla_version,"/fasta/B_nuc.fasta",sep = ""), "C"=paste("/root/shiny/IMGTHLA/",imgthla_version,"/fasta/C_nuc.fasta",sep = ""))
        if (grepl("N", allele_number)) {
          cmd <- paste("grep",allele_number,IMGTHLA_nuc,"| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq",IMGTHLA_nuc," - | sed -n '2p' -")
        } else {
          cmd <- paste("grep",allele_number,IMGTHLA_nuc,"| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq",IMGTHLA_nuc," - | sed -n '2p' -")
        }
        optitype_mhc_table[i,2] <- system(cmd, intern = TRUE)
      }
      #----------------------------------------------------------
      for (i in 1:nrow(optitype_mhc_table)) {
        imgthla_version = input$imgthla
        seq = optitype_mhc_table[i,1]
        allele_type = substr(seq,1,1)
        allele_number = substr(seq,3,nchar(seq))
        IMGTHLA_prot <- switch(allele_type, "A"=paste("/root/shiny/IMGTHLA/",imgthla_version,"/fasta/A_prot.fasta",sep = ""), "B"=paste("/root/shiny/IMGTHLA/",imgthla_version,"/fasta/B_prot.fasta",sep = ""), "C"=paste("/root/shiny/IMGTHLA/",imgthla_version,"/fasta/C_prot.fasta",sep = ""))
        if (grepl("N", allele_number)) {
          cmd <- paste("grep",allele_number,IMGTHLA_prot,"| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq",IMGTHLA_prot," - | sed -n '2p' -")
        } else {
          cmd <- paste("grep",allele_number,IMGTHLA_prot,"| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq",IMGTHLA_prot," - | sed -n '2p' -")
        }
        optitype_mhc_table[i,3] <- system(cmd, intern = TRUE)
      }
      #----------------------------------------------------------
      output$hla_typing_table <- renderDataTable({
        datatable(optitype_mhc_table,
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
                        "}")
                    )))
        )
      })
      #----------------------------------------------------------
      sequence_file <- NULL
      for (i in 1:nrow(optitype_mhc_table)) {
        nucleotide <- paste(">", optitype_mhc_table[i, 1], "_nucleotide", sep = "")
        protein <- paste(">", optitype_mhc_table[i, 1], "_protein", sep = "")
        seqs <- rbind(rbind(nucleotide, optitype_mhc_table[i, 2]), rbind(protein, optitype_mhc_table[i, 3]))
        sequence_file <- rbind(sequence_file, seqs)
      }
      output_file = paste("./Output/optitype_",input$sequence,"_",input$imgthla,".txt",sep = "")
      write(sequence_file,file = output_file)
      #----------------------------------------------------------
      dir.create(paste("./Download/optitype_",input$sequence,"_",input$imgthla,"_zip",sep = ""), showWarnings = FALSE)
      zip_folder_pathway = paste("./Download/optitype_",input$sequence,"_",input$imgthla,"_zip",sep = "")
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
      output_file_zip <- zip(file.path(tmp_dir,paste("Download/optitype_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")), files = ".")
      output_file_zip_pathway <- paste("Download/optitype_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")
      setwd(tmp_dir)
      #----------------------------------------------------------
      output$hla_typing_download <- downloadHandler(
        filename = function() {
          paste("optitype_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")},
        content = function(file) {
          file.copy(file.path(tmp_dir,output_file_zip_pathway), file)
        }
      )
    }
    #----------------------------------------------------------
  }else if (input$package == "arcasHLA") {
    if (dir.exists("./Output/arcashla")) {
      unlink("./Output/arcashla", recursive = TRUE)
    }
    dir.create("./Output/arcashla", showWarnings = FALSE)
    #----------------------------------------------------------
    if (input$imgthla == "v3.14.0") {
      system(paste("/arcasHLA/arcasHLA genotype14 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    } else if (input$imgthla == "v3.49.0") {
      system(paste("/arcasHLA/arcasHLA genotype49 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    } else if (input$imgthla == "v3.50.0") {
      system(paste("/arcasHLA/arcasHLA genotype50 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    } else if (input$imgthla == "v3.51.0") {
      system(paste("/arcasHLA/arcasHLA genotype51 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    } else if (input$imgthla == "v3.52.0") {
      system(paste("/arcasHLA/arcasHLA genotype52 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    } else if (input$imgthla == "v3.53.0") {
      system(paste("/arcasHLA/arcasHLA genotype53 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    } else if (input$imgthla == "v3.54.0") {
      system(paste("/arcasHLA/arcasHLA genotype54 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    } else if (input$imgthla == "v3.55.0") {
      system(paste("/arcasHLA/arcasHLA genotype55 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    } else if (input$imgthla == "v3.56.0") {
      system(paste("/arcasHLA/arcasHLA genotype56 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    } else if (input$imgthla == "v3.57.0") {
      system(paste("/arcasHLA/arcasHLA genotype57 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    } else if (input$imgthla == "v3.58.0") {
      system(paste("/arcasHLA/arcasHLA genotype58 ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o ./Output/arcashla -t ", num_cores, " -v",sep = ""))
    }
    #----------------------------------------------------------
    arcashla_genotype_path <- paste("./Output/arcashla/sample1_1.genotype.json",sep = "")
    arcashla_genotype <- fromJSON(file = arcashla_genotype_path)
    #----------------------------------------------------------
    # whether the dataframe is empty
    if (length(arcashla_genotype) == 0) {
      output$hla_typing_table <- renderDataTable({
        datatable(
          data.frame(Message = "HLA typing unsuccessful: Insufficient sequencing reads detected."),
          class = 'nowrap'
        )
      })
      
      # Define file and directory paths
      output_file = paste("./Output/arcashla_",input$sequence,"_",input$imgthla,".txt",sep = "")
      zip_folder_pathway = paste("./Download/arcashla_",input$sequence,"_",input$imgthla,"_zip",sep = "")
      output_file_zip_pathway <- paste("Download/arcashla_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")
      merged_file = paste("./Pivotable/arcashla_",input$sequence,"_",input$imgthla,"_merged.csv",sep = "")
      
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
      arcashla_mhc_all <- rbind(arcashla_mhc_A,arcashla_mhc_B,arcashla_mhc_C,arcashla_mhc_DPB1,arcashla_mhc_DQA1,arcashla_mhc_DQB1,arcashla_mhc_DRB1)
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
        protein = character(nrow(arcashla_stacked)))
      for (i in 1:nrow(arcashla_stacked)) {
        arcashla_mhc_table[i, 1] <- arcashla_stacked[i,]
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
                                "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_nuc.fasta", sep = ""))
        } else {
          # Other types of genes
          allele_type <- substr(seq, 1, 4)
          allele_number <- sub("^.{4}\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_nuc <- switch(allele_type, 
                                "DPB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPB1_nuc.fasta", sep = ""), 
                                "DQA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQA1_nuc.fasta", sep = ""), 
                                "DQB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQB1_nuc.fasta", sep = ""), 
                                "DRB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DRB1_nuc.fasta", sep = ""))
        }
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        }
        result <- system(cmd, intern = TRUE)
        if (length(result) == 0) {
          warning(paste("No result found for allele_number:", allele_number, "in sequence:", seq))
          arcashla_mhc_table[i, 2] <- NA  # or some other appropriate default value or handling
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
                                 "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_prot.fasta", sep = ""))
        } else {
          # Other types of genes
          allele_type <- substr(seq, 1, 4)
          allele_number <- sub("^.{4}\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_prot <- switch(allele_type, 
                                 "DPB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPB1_prot.fasta", sep = ""), 
                                 "DQA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQA1_prot.fasta", sep = ""), 
                                 "DQB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQB1_prot.fasta", sep = ""), 
                                 "DRB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DRB1_prot.fasta", sep = ""))
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
      #----------------------------------------------------------
      output$hla_typing_table <- renderDataTable({
        datatable(arcashla_mhc_table,
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
                        "}")
                    )))
        )
      })
      #----------------------------------------------------------
      sequence_file <- NULL
      for (i in 1:nrow(arcashla_mhc_table)) {
        nucleotide <- paste(">", arcashla_mhc_table[i, 1], "_nucleotide", sep = "")
        protein <- paste(">", arcashla_mhc_table[i, 1], "_protein", sep = "")
        seqs <- rbind(rbind(nucleotide, arcashla_mhc_table[i, 2]), rbind(protein, arcashla_mhc_table[i, 3]))
        sequence_file <- rbind(sequence_file, seqs)
      }
      output_file = paste("./Output/arcashla_",input$sequence,"_",input$imgthla,".txt",sep = "")
      write(sequence_file,file = output_file)
      #----------------------------------------------------------
      dir.create(paste("./Download/arcashla_",input$sequence,"_",input$imgthla,"_zip",sep = ""), showWarnings = FALSE)
      zip_folder_pathway = paste("./Download/arcashla_",input$sequence,"_",input$imgthla,"_zip",sep = "")
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
      output_file_zip <- zip(file.path(tmp_dir,paste("Download/arcashla_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")), files = ".")
      output_file_zip_pathway <- paste("Download/arcashla_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")
      setwd(tmp_dir)
      #----------------------------------------------------------
      output$hla_typing_download <- downloadHandler(
        filename = function() {
          paste("arcashla_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")},
        content = function(file) {
          file.copy(file.path(tmp_dir,output_file_zip_pathway), file)
        }
      )
    }
    #----------------------------------------------------------
  }else if (input$package == "HLA-HD") {
    if (dir.exists("./Output/hlahd")) {
      unlink("./Output/hlahd", recursive = TRUE)
    }
    dir.create("./Output/hlahd", showWarnings = FALSE)
    #----------------------------------------------------------
    if (input$imgthla == "v3.14.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary14 sample1 ./Output/hlahd",sep = ""))
    } else if (input$imgthla == "v3.49.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary49 sample1 ./Output/hlahd",sep = ""))
    } else if (input$imgthla == "v3.50.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary50 sample1 ./Output/hlahd",sep = ""))
    } else if (input$imgthla == "v3.51.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary51 sample1 ./Output/hlahd",sep = ""))
    } else if (input$imgthla == "v3.52.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary52 sample1 ./Output/hlahd",sep = ""))
    } else if (input$imgthla == "v3.53.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary53 sample1 ./Output/hlahd",sep = ""))
    } else if (input$imgthla == "v3.54.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary54 sample1 ./Output/hlahd",sep = ""))
    } else if (input$imgthla == "v3.55.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary55 sample1 ./Output/hlahd",sep = ""))
    } else if (input$imgthla == "v3.56.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary56 sample1 ./Output/hlahd",sep = ""))
    } else if (input$imgthla == "v3.57.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary57 sample1 ./Output/hlahd",sep = ""))
    } else if (input$imgthla == "v3.58.0") {
      system(paste("hlahd.sh -t ", num_cores, " -m 100 -c 0.95 -f /hlahd.1.7.0/freq_data ./Input/sample1_1.fastq.gz ./Input/sample1_2.fastq.gz /hlahd.1.7.0/HLA_gene.split.txt /hlahd.1.7.0/IMGTHLA_version/dictionary58 sample1 ./Output/hlahd",sep = ""))
    }
    #-----------------------------------------------------------
    if (dir.exists("./IGV/hlahd")) {
      unlink("./IGV/hlahd", recursive = TRUE)
    }
    dir.create("./IGV/hlahd", showWarnings = FALSE)
    #----------------------------------------------------------
    # Merge two SAM files
    system("cat ./Output/hlahd/sample1/mapfile/sample1.all.R1.pmap.NM.sam ./Output/hlahd/sample1/mapfile/sample1.all.R2.pmap.NM.sam > ./IGV/hlahd/sample1.all.merged.sam")
    
    # Extract headers and move them to the beginning of the file
    system("grep '^@' ./IGV/hlahd/sample1.all.merged.sam > ./IGV/hlahd/header.sam")
    system("grep -v '^@' ./IGV/hlahd/sample1.all.merged.sam > ./IGV/hlahd/body.sam")
    system("cat ./IGV/hlahd/header.sam ./IGV/hlahd/body.sam > ./IGV/hlahd/sample1.all.corrected.merged.sam")
    
    # Keep only the first occurrence of each sequence name and remove duplicates
    system("awk '!seen[$0]++' ./IGV/hlahd/sample1.all.corrected.merged.sam > ./IGV/hlahd/sample1.all.rmdup.sam")
    
    # Convert the merged SAM file to BAM format
    system("samtools view -Sb ./IGV/hlahd/sample1.all.rmdup.sam > ./IGV/hlahd/sample1.all.rmdup.bam")
    
    # Sort the merged BAM file
    system("samtools sort -o ./IGV/hlahd/sample1.all.sorted.bam ./IGV/hlahd/sample1.all.rmdup.bam")
    
    # Remove intermediate files
    file.remove("./IGV/hlahd/sample1.all.merged.sam")
    file.remove("./IGV/hlahd/sample1.all.corrected.merged.sam")
    file.remove("./IGV/hlahd/sample1.all.rmdup.sam")
    file.remove("./IGV/hlahd/sample1.all.rmdup.bam")
    file.remove("./IGV/hlahd/header.sam")
    file.remove("./IGV/hlahd/body.sam")
    
    # Extract exon names matched with the reference from the mapping results as a list
    system("samtools view ./IGV/hlahd/sample1.all.sorted.bam | cut -f 3 | grep -E '^[ABC]:Exon' > ./IGV/hlahd/hlahd_mhc_i_list.txt")
    system("samtools view ./IGV/hlahd/sample1.all.sorted.bam | cut -f 3 | grep -E '^D(PA|PB|QA|QB|RB)1:Exon' > ./IGV/hlahd/hlahd_mhc_ii_list.txt")
    
    # Create hla_reference.fasta for sample1.all.sorted.bam
    input_list_I <- file.path(getwd(), "IGV/hlahd/hlahd_mhc_i_list.txt")
    input_list_II <- file.path(getwd(), "IGV/hlahd/hlahd_mhc_ii_list.txt")
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
    } else {
      stop("Unknown IMGTHLA version")
    }
    #----------------------------------------------------------
    igv_fasta_I <- file.path(getwd(), "IGV/hlahd/hlahd_hla_reference_I_igv.fasta")
    igv_fasta_II <- file.path(getwd(), "IGV/hlahd/hlahd_hla_reference_II_igv.fasta")
    # Create IGV Genome reference
    system(paste("python /root/shiny/Server/hlahd_igv_genome.py", 
                 input_list_I, 
                 reference_fasta, 
                 igv_fasta_I, 
                 sep = " "))
    system(paste("python /root/shiny/Server/hlahd_igv_genome.py", 
                 input_list_II, 
                 reference_fasta, 
                 igv_fasta_II, 
                 sep = " "))
    
    # Step 1: Modify BAM file using samtools and sed
    system('samtools view -h ./IGV/hlahd/sample1.all.sorted.bam | 
       sed "s/A:Exon/A_Exon/g; s/B:Exon/B_Exon/g; s/C:Exon/C_Exon/g; s/DPA1:Exon/DPA1_Exon/g; 
            s/DPB1:Exon/DPB1_Exon/g; s/DQA1:Exon/DQA1_Exon/g; s/DQB1:Exon/DQB1_Exon/g; 
            s/DRB1:Exon/DRB1_Exon/g" > ./IGV/hlahd/sample1.modified.sam')
    
    # Step 2: Convert modified SAM file back to BAM format
    system('samtools view -Sb ./IGV/hlahd/sample1.modified.sam > ./IGV/hlahd/sample1.modified.bam')
    
    # Step 3: Modify HLA reference I FASTA file by replacing ':' with '_'
    system('mv ./IGV/hlahd/hlahd_hla_reference_I_igv.fasta ./IGV/hlahd/hlahd_hla_reference_I_igv_colon.fasta')
    system('sed "s/:Exon/_Exon/g" ./IGV/hlahd/hlahd_hla_reference_I_igv_colon.fasta > ./IGV/hlahd/hlahd_hla_reference_I_igv.fasta')
    
    # Step 4: Modify HLA reference II FASTA file by replacing ':' with '_'
    system('mv ./IGV/hlahd/hlahd_hla_reference_II_igv.fasta ./IGV/hlahd/hlahd_hla_reference_II_igv_colon.fasta')
    system('sed "s/:Exon/_Exon/g" ./IGV/hlahd/hlahd_hla_reference_II_igv_colon.fasta > ./IGV/hlahd/hlahd_hla_reference_II_igv.fasta')
    
    # Create reference.fasta index
    fasta_files <- c(
      "./IGV/hlahd/hlahd_hla_reference_I_igv.fasta",
      "./IGV/hlahd/hlahd_hla_reference_II_igv.fasta"
    )
    
    # Function to check and index FASTA
    index_fasta <- function(fasta_path) {
      if (file.exists(fasta_path) && file.size(fasta_path) > 0) {
        system(paste("samtools faidx", fasta_path))
      }
    }
    
    # Run indexing for each file
    lapply(fasta_files, index_fasta)
    
    # Remove intermediate files
    file.remove(input_list_I)
    file.remove(input_list_II)
    file.remove("./IGV/hlahd/sample1.all.sorted.bam")
    file.remove("./IGV/hlahd/sample1.modified.sam")
    file.remove("./IGV/hlahd/hlahd_hla_reference_I_igv_colon.fasta")
    file.remove("./IGV/hlahd/hlahd_hla_reference_II_igv_colon.fasta")
    #-----------------------------------------------------------
    system("head -n 8 ./Output/hlahd/sample1/result/sample1_final.result.txt > ./Output/hlahd/sample1_final_mhc.result.txt")
    
    # Read the file and remove lines containing 'Couldn't read result file.'
    lines <- readLines("./Output/hlahd/sample1_final_mhc.result.txt")
    
    # Filter out lines containing 'Couldn't read result file.'
    valid_lines <- lines[!grepl("Couldn't read result file.", lines)]
    
    # Write the filtered content to a new temporary file
    temp_file <- "./Output/hlahd/sample1_final_mhc_clean.result.txt"
    writeLines(valid_lines, temp_file)
    
    # Read the cleaned file
    hlahd <- read.table(temp_file, header = FALSE)
    hlahd <- hlahd[ ,-1]
    hlahd$V2 <- sub("HLA-", "",hlahd$V2)
    hlahd$V3 <- sub("HLA-", "",hlahd$V3)
    
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
          class = 'nowrap'
        )
      })
      
      # Define file and directory paths
      output_file = paste("./Output/hlahd_",input$sequence,"_",input$imgthla,".txt",sep = "")
      zip_folder_pathway = paste("./Download/hlahd_",input$sequence,"_",input$imgthla,"_zip",sep = "")
      output_file_zip_pathway <- paste("Download/hlahd_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")
      merged_file = paste("./Pivotable/hlahd_",input$sequence,"_",input$imgthla,"_merged.csv",sep = "")
      
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
        protein = character(nrow(hlahd_stacked)))
      for (i in 1:nrow(hlahd_stacked)) {
        hlahd_mhc_table[i, 1] <- hlahd_stacked[i,]
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
                                "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_nuc.fasta", sep = ""))
        } else {
          # Other types of genes
          allele_type <- substr(seq, 1, 4)
          allele_number <- sub("^.{4}\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_nuc <- switch(allele_type, 
                                "DPA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPA1_nuc.fasta", sep = ""), 
                                "DPB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPB1_nuc.fasta", sep = ""), 
                                "DQA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQA1_nuc.fasta", sep = ""), 
                                "DQB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQB1_nuc.fasta", sep = ""), 
                                "DRB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DRB1_nuc.fasta", sep = ""))
        }
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        }
        result <- system(cmd, intern = TRUE)
        if (length(result) == 0) {
          warning(paste("No result found for allele_number:", allele_number, "in sequence:", seq))
          hlahd_mhc_table[i, 2] <- NA  # or some other appropriate default value or handling
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
                                 "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_prot.fasta", sep = ""))
        } else {
          # Other types of genes
          allele_type <- substr(seq, 1, 4)
          allele_number <- sub("^.{4}\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_prot <- switch(allele_type, 
                                 "DPA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPA1_prot.fasta", sep = ""), 
                                 "DPB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPB1_prot.fasta", sep = ""), 
                                 "DQA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQA1_prot.fasta", sep = ""), 
                                 "DQB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQB1_prot.fasta", sep = ""), 
                                 "DRB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DRB1_prot.fasta", sep = ""))
        }
        
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_prot, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_prot, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_prot, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_prot, " - | sed -n '2p' -")
        }
        
        result <- system(cmd, intern = TRUE)
        
        if (length(result) == 0) {
          warning(paste("No result found for allele_number:", allele_number, "in sequence:", seq))
          hlahd_mhc_table[i, 3] <- NA  # or some other appropriate default value or handling
        } else {
          hlahd_mhc_table[i, 3] <- result
        }
      }
      #----------------------------------------------------------
      output$hla_typing_table <- renderDataTable({
        datatable(hlahd_mhc_table,
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
                        "}")
                    )))
        )
      })
      #----------------------------------------------------------
      sequence_file <- NULL
      for (i in 1:nrow(hlahd_mhc_table)) {
        nucleotide <- paste(">", hlahd_mhc_table[i, 1], "_nucleotide", sep = "")
        protein <- paste(">", hlahd_mhc_table[i, 1], "_protein", sep = "")
        seqs <- rbind(rbind(nucleotide, hlahd_mhc_table[i, 2]), rbind(protein, hlahd_mhc_table[i, 3]))
        sequence_file <- rbind(sequence_file, seqs)
      }
      output_file = paste("./Output/hlahd_",input$sequence,"_",input$imgthla,".txt",sep = "")
      write(sequence_file,file = output_file)
      #----------------------------------------------------------
      dir.create(paste("./Download/hlahd_",input$sequence,"_",input$imgthla,"_zip",sep = ""), showWarnings = FALSE)
      zip_folder_pathway = paste("./Download/hlahd_",input$sequence,"_",input$imgthla,"_zip",sep = "")
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
      output_file_zip <- zip(file.path(tmp_dir,paste("Download/hlahd_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")), files = ".")
      output_file_zip_pathway <- paste("Download/hlahd_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")
      setwd(tmp_dir)
      #----------------------------------------------------------
      output$hla_typing_download <- downloadHandler(
        filename = function() {
          paste("hlahd_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")},
        content = function(file) {
          file.copy(file.path(tmp_dir,output_file_zip_pathway),file)
        }
      )
    }
    #----------------------------------------------------------
  }else if (input$package == "SpecHLA") {
    if (dir.exists("./Output/spechla")) {
      unlink("./Output/spechla", recursive = TRUE)
    }
    dir.create("./Output/spechla", showWarnings = FALSE)
    system(paste("/SpecHLA/spechla.sh ", input$imgthla, " -u 1 -j ", num_cores, " -n sample1 -1 ./Input/sample1_1.fastq.gz -2 ./Input/sample1_2.fastq.gz -o ./Output/spechla",sep = ""))
    #----------------------------------------------------------
    if (dir.exists("./IGV/spechla")) {
      unlink("./IGV/spechla", recursive = TRUE)
    }
    dir.create("./IGV/spechla", showWarnings = FALSE)
    
    # Insert into the reference file
    # Copy the reference FASTA file if it exists
    source_ref <- "/SpecHLA/SpecHLA57/db/ref/hla.ref.extend.fa"
    dest_ref   <- file.path(tmp_dir, "IGV", "spechla", "specHLA_hla_reference.fasta")
    if (file.exists(source_ref)) {
      cp_cmd <- paste("cp", source_ref, dest_ref)
      system(cp_cmd)
    } else {
      message("[Warning] Source reference file '", source_ref, "' does not exist. Skipping copy.")
    }
    
    # Copy the reference index (.fai) file if it exists
    source_fai <- "/SpecHLA/SpecHLA57/db/ref/hla.ref.extend.fa.fai"
    dest_fai   <- file.path(tmp_dir, "IGV", "spechla", "specHLA_hla_reference.fasta.fai")
    if (file.exists(source_fai)) {
      cp_fai_cmd <- paste("cp", source_fai, dest_fai)
      system(cp_fai_cmd)
    } else {
      message("[Warning] Source reference index file '", source_fai, "' does not exist. Skipping copy.")
    }
    
    # Insert into the alignment results
    
    # Move the BAM file if it exists
    source_bam <- "/sample1.merge.bam"
    dest_bam   <- file.path(tmp_dir, "IGV", "spechla", "sample1.merge.bam")
    if (file.exists(source_bam)) {
      mv_cmd <- paste("mv", source_bam, dest_bam)
      system(mv_cmd)
    } else {
      message("[Warning] Source BAM file '", source_bam, "' does not exist. Skipping move.")
    }
    
    # Read the result file if it exists
    result_file <- "./Output/spechla/sample1/hla.result.txt"
    if (file.exists(result_file)) {
      spechla <- read.table(result_file)
    } else {
      message("[Warning] Result file '", result_file, "' does not exist. Skipping read.")
    }
    
    # whether the dataframe is empty
    if (nrow(spechla) == 0) {
      output$hla_typing_table <- renderDataTable({
        datatable(
          data.frame(Message = "HLA typing unsuccessful: Insufficient sequencing reads detected."),
          class = 'nowrap'
        )
      })
      
      # Define file and directory paths
      output_file = paste("./Output/spechla_",input$sequence,"_",input$imgthla,".txt",sep = "")
      zip_folder_pathway = paste("./Download/spechla_",input$sequence,"_",input$imgthla,"_zip",sep = "")
      output_file_zip_pathway <- paste("Download/spechla_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")
      merged_file = paste("./Pivotable/spechla_",input$sequence,"_",input$imgthla,"_merged.csv",sep = "")
      
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
      spechla_mhc_table <- data.frame(
        Allele = character(16),
        nucleotide = character(16),
        protein = character(16))
      for (i in 1:16) {
        spechla_mhc_table[i, 1] <- spechla[2, (i + 1)]
      }
      if ("-" %in% spechla_mhc_table$Allele) {
        spechla_mhc_table <- subset(spechla_mhc_table, Allele != "-")
      }
      # sorting
      spechla_mhc_table <- spechla_mhc_table[order(spechla_mhc_table$Allele), , drop = FALSE]
      imgthla_version <- input$imgthla
      for (i in 1:nrow(spechla_mhc_table)) {
        seq <- spechla_mhc_table[i, 1]
        # Determine whether it is a type A, B, or C gene
        if (grepl("^[ABC]\\*", seq)) {
          allele_type <- substr(seq, 1, 1)
          allele_number <- sub("^[ABC]\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_nuc <- switch(allele_type,
                                "A" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/A_nuc.fasta", sep = ""),
                                "B" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/B_nuc.fasta", sep = ""),
                                "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_nuc.fasta", sep = ""))
        } else {
          # Other types of genes
          allele_type <- substr(seq, 1, 4)
          allele_number <- sub("^.{4}\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_nuc <- switch(allele_type,
                                "DPA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPA1_nuc.fasta", sep = ""),
                                "DPB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPB1_nuc.fasta", sep = ""),
                                "DQA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQA1_nuc.fasta", sep = ""),
                                "DQB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQB1_nuc.fasta", sep = ""),
                                "DRB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DRB1_nuc.fasta", sep = ""))
        }
        #----------------------------------------------------------
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_nuc, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_nuc, " - | sed -n '2p' -")
        }
        #----------------------------------------------------------
        result <- system(cmd, intern = TRUE)
        #----------------------------------------------------------
        if (length(result) == 0) {
          warning(paste("No result found for allele_number:", allele_number, "in sequence:", seq))
          spechla_mhc_table[i, 2] <- NA
        } else {
          spechla_mhc_table[i, 2] <- result
        }
      }
      #----------------------------------------------------------
      imgthla_version <- input$imgthla
      for (i in 1:nrow(spechla_mhc_table)) {
        seq <- spechla_mhc_table[i, 1]
        # Determine whether it is a type A, B, or C gene
        if (grepl("^[ABC]\\*", seq)) {
          allele_type <- substr(seq, 1, 1)
          allele_number <- sub("^[ABC]\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_prot <- switch(allele_type,
                                 "A" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/A_prot.fasta", sep = ""),
                                 "B" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/B_prot.fasta", sep = ""),
                                 "C" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/C_prot.fasta", sep = ""))
        } else {
          # Other types of genes
          allele_type <- substr(seq, 1, 4)
          allele_number <- sub("^.{4}\\*([0-9:]+)", "\\1", seq)
          IMGTHLA_prot <- switch(allele_type,
                                 "DPA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPA1_prot.fasta", sep = ""),
                                 "DPB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DPB1_prot.fasta", sep = ""),
                                 "DQA1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQA1_prot.fasta", sep = ""),
                                 "DQB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DQB1_prot.fasta", sep = ""),
                                 "DRB1" = paste("/root/shiny/IMGTHLA/", imgthla_version, "/fasta/DRB1_prot.fasta", sep = ""))
        }
        #----------------------------------------------------------
        if (grepl("N", allele_number)) {
          cmd <- paste("grep", allele_number, IMGTHLA_prot, "| sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_prot, " - | sed -n '2p' -")
        } else {
          cmd <- paste("grep", allele_number, IMGTHLA_prot, "| grep -v 'N' | sort -t: -k2,2n -k3,3n -k4,4n | sed 's/^>//g' | seqtk subseq", IMGTHLA_prot, " - | sed -n '2p' -")
        }
        #----------------------------------------------------------
        result <- system(cmd, intern = TRUE)
        #----------------------------------------------------------
        if (length(result) == 0) {
          warning(paste("No result found for allele_number:", allele_number, "in sequence:", seq))
          spechla_mhc_table[i, 3] <- NA  # or some other appropriate default value or handling
        } else {
          spechla_mhc_table[i, 3] <- result
        }
      }
      #----------------------------------------------------------
      output$hla_typing_table <- renderDataTable({
        datatable(spechla_mhc_table,
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
                        "}")
                    )))
                  )
      })
      #----------------------------------------------------------
      # Assume sequence_file has been defined as an empty dataframe or matrix
      sequence_file <- NULL
      for (i in 1:nrow(spechla_mhc_table)) {
        seq = spechla_mhc_table[i, 1]
        nucleotide <- paste(">", seq, "_nucleotide", sep = "")
        protein <- paste(">", seq, "_protein", sep = "")
        seqs <- rbind(rbind(nucleotide, spechla_mhc_table[i, 2]), rbind(protein, spechla_mhc_table[i, 3]))
        sequence_file <- rbind(sequence_file, seqs)
      }
      #----------------------------------------------------------
      output_file = paste("./Output/spechla_",input$sequence,"_",input$imgthla,".txt",sep = "")
      write(sequence_file,file = output_file)
      #----------------------------------------------------------
      dir.create(paste("./Download/spechla_",input$sequence,"_",input$imgthla,"_zip",sep = ""), showWarnings = FALSE)
      zip_folder_pathway = paste("./Download/spechla_",input$sequence,"_",input$imgthla,"_zip",sep = "")
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
      output_file_zip <- zip(file.path(tmp_dir,paste("Download/spechla_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")), files = ".")
      output_file_zip_pathway <- paste("Download/spechla_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")
      setwd(tmp_dir)
      #----------------------------------------------------------
      output$hla_typing_download <- downloadHandler(
        filename = function() {
          paste("spechla_",input$sequence,"_",input$imgthla,"_output.zip",sep = "")},
        content = function(file) {
          file.copy(file.path(tmp_dir,output_file_zip_pathway), file)
        }
      )
    }
  }
  w$hide()
  #----------------------------------------------------------
  # Calculate execution time
  end_time <- Sys.time()
  execution_time <- difftime(end_time, start_time, units = "secs")
  execution_time <- round(as.numeric(execution_time), 0)
  
  # Record the current time and execution time
  log_entry <- paste(getCurrentTime(), "-",input$sequence,"|",input$imgthla,"|",input$package,"|","Execution Time:", execution_time, "seconds.")
  
  # Append the new execution time to the existing log, adding from top to bottom
  execution_times$times <- c(log_entry, execution_times$times)
  
  # Display all execution times line by line
  output$hla_log <- renderUI({
    HTML(paste(
      "<span style='font-size:12px;'>",
      paste(execution_times$times, collapse = "<br>"),
      "</span>"
    ))
  })
  #----------------------------------------------------------
  setwd(ori_dir)
})