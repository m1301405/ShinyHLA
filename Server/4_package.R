#----------------------------------------------------------
execution_times <- reactiveValues(times = character())
#----------------------------------------------------------
getCurrentTime <- function() {
  current_time <- Sys.time()
  formatted_time <- format(current_time, "%Y-%m-%d %H:%M")
  return(formatted_time)
}
#----------------------------------------------------------
observeEvent(input$hla_typing_button, {
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
  dir.create("./IGV", showWarnings = FALSE)
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
    
    # Extract the "main sequence" from the alignment file
    
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
    
    ## Clean BAM (remove 0x800 + 0x100)
    clean_unsorted <- "./IGV/optitype/optitype_clean.unsorted.bam"
    clean_bam <- "./IGV/optitype/optitype_clean.bam"
    
    if (file.exists(filtered_bam)) {
      
      # Remove supplementary + secondary (keep header)
      cmd_view <- paste(
        "samtools view -bh -F 0x900",
        shQuote(filtered_bam),
        "-o",
        shQuote(clean_unsorted)
      )
      system(cmd_view)
      
      # Sort for safety
      cmd_sort <- paste(
        "samtools sort -o",
        shQuote(clean_bam),
        shQuote(clean_unsorted)
      )
      system(cmd_sort)
      
      # Index
      system(paste("samtools index", shQuote(clean_bam)))
      
      # Remove intermediate file
      if (file.exists(clean_unsorted)) {
        file.remove(clean_unsorted)
      }
      
    } else {
      message("[Warning] filtered BAM not found. Clean step skipped.")
    }
    
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
      py_cmd <- paste(
        "python /root/shiny/Server/optitype_igv_genome.py",
        input_bam, reference_fasta, igv_fasta
      )
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
      igv_folder <- "./IGV/optitype"
      
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
      
      # Remove the previous successful typing folder (igv_folder) if it exists
      if (dir.exists(igv_folder)) {
        unlink(igv_folder, recursive = TRUE)
      }
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
      sequence_file <- NULL
      for (i in 1:nrow(optitype_mhc_table)) {
        nucleotide <- paste(">", optitype_mhc_table[i, 1], "_nucleotide", sep = "")
        protein <- paste(">", optitype_mhc_table[i, 1], "_protein", sep = "")
        seqs <- rbind(rbind(nucleotide, optitype_mhc_table[i, 2]), rbind(protein, optitype_mhc_table[i, 3]))
        sequence_file <- rbind(sequence_file, seqs)
      }
      output_file <- paste("./Output/optitype_", input$sequence, "_", input$imgthla, ".txt", sep = "")
      write(sequence_file, file = output_file)
      #----------------------------------------------------------
      # Extract typed allele list (safe version)
      allele_txt <- output_file
      allele_list_file <- "./IGV/optitype/typed_alleles.txt"
      
      if (file.exists(allele_txt)) {
        
        extract_cmd <- paste(
          "grep '^>'", shQuote(allele_txt),
          "| sed 's/^>//'",
          "| sed 's/_nucleotide$//'",
          "| sed 's/_protein$//'",
          "| sort -u",
          ">",
          shQuote(allele_list_file)
        )
        
        system(extract_cmd)
        
        if (!file.exists(allele_list_file)) {
          message("[Error] Failed to generate typed allele list.")
        }
        
      } else {
        message("[Warning] Allele txt not found.")
      }
      #----------------------------------------------------------
      # NEW: Build contig mapping (safe)
      allelelist_path <- file.path("/root/shiny/IMGTHLA", input$imgthla, "Allelelist.txt")
      contig_map_file <- "./IGV/optitype/contig2allele4.tsv"
      
      if (file.exists(allelelist_path)) {
        
        awk_cmd <- paste(
          "awk -v OFS='\\t' 'NF>=2 {print $1,$2}'",
          shQuote(allelelist_path),
          ">",
          shQuote(contig_map_file)
        )
        
        system(awk_cmd)
        
        if (!file.exists(contig_map_file)) {
          message("[Error] Failed to generate contig mapping file.")
        }
        
      } else {
        message("[Warning] Allelelist.txt not found.")
      }
      #----------------------------------------------------------
      # Map typed 2-field → contig IDs
      allele_list_file   <- "./IGV/optitype/typed_alleles.txt"
      contig_map_file    <- "./IGV/optitype/contig2allele4.tsv"
      target_contig_file <- "./IGV/optitype/target_contigs.txt"
      
      if (file.exists(allele_list_file) && file.exists(contig_map_file)) {
        
        join_cmd <- paste(
          "awk -F'\\t' 'NR==FNR{t[$1]=1;next}",
          "{split($2,x,\":\");key=x[1]\":\"x[2]; if(key in t) print $1}'",
          shQuote(allele_list_file),
          shQuote(contig_map_file),
          "| sort -u >",
          shQuote(target_contig_file)
        )
        
        system(join_cmd)
        
        if (file.exists(target_contig_file)) {
          n_contigs <- length(readLines(target_contig_file))
          message("[Info] Target contigs detected: ", n_contigs)
        }
        
        if (!file.exists(target_contig_file)) {
          message("[Error] Failed to generate target_contigs.txt")
        }
        
      } else {
        message("[Warning] Missing typed_alleles.txt or contig2allele4.tsv")
      }
      #----------------------------------------------------------
      # --- NEW: keep only contigs with mapped reads > 0, then intersect with target_contigs.txt ---
      target_contigs <- "./IGV/optitype/target_contigs.txt"
      
      bam_contigs_mapped <- "./IGV/optitype/bam_contigs_mapped.txt"
      target_contigs_mapped_present <- "./IGV/optitype/target_contigs_mapped_present.txt"
      
      if (file.exists(clean_bam) && file.exists(target_contigs)) {
        
        system(paste(
          "samtools idxstats", shQuote(clean_bam),
          "| awk -F'\\t' '$3>0 {print $1}' >",
          shQuote(bam_contigs_mapped)
        ))
        
        system(paste(
          "grep -Fxf", shQuote(bam_contigs_mapped), shQuote(target_contigs),
          ">", shQuote(target_contigs_mapped_present)
        ))
        
      } else {
        message("[Warning] clean_bam or target_contigs missing; cannot build mapped-present contigs list.")
      }
      #----------------------------------------------------------
      # NEW: Build target BAM (mapped contigs only)
      target_bam <- "./IGV/optitype/optitype_target.bam"
      
      if (file.exists(clean_bam) && file.exists(target_contigs_mapped_present)) {
        
        ctgs <- readLines(target_contigs_mapped_present)
        ctgs <- ctgs[nzchar(ctgs)]
        
        if (length(ctgs) > 0) {
          
          view_cmd <- paste(
            "samtools view -b",
            shQuote(clean_bam),
            paste(shQuote(ctgs), collapse = " ")
          )
          
          sort_cmd <- paste(
            "samtools sort -o",
            shQuote(target_bam),
            "-"
          )
          
          system(paste(view_cmd, "|", sort_cmd))
          system(paste("samtools index", shQuote(target_bam)))
          
        } else {
          message("[Warning] No mapped target contigs found.")
        }
        
      } else {
        message("[Warning] clean_bam or target_contigs_mapped_present missing.")
      }
      #----------------------------------------------------------
      # NEW: Build rename map (2-field)
      rename_map <- "./IGV/optitype/rename_map.tsv"
      target_contigs_mapped_present <- "./IGV/optitype/target_contigs_mapped_present.txt"
      contig_map_file <- "./IGV/optitype/contig2allele4.tsv"
      
      if (file.exists(target_contigs_mapped_present) && file.exists(contig_map_file)) {
        
        awk_cmd <- paste(
          "awk -F'\\t' '",
          "NR==FNR {keep[$1]=1; next}",
          "{",
          "ctg=$1; allele=$2;",
          "if (ctg in keep) {",
          "split(allele,x,\":\");",
          "key=x[1]\":\"x[2];",
          "gsub(\":\",\"-\",key);",
          "print ctg\"\\tHLA-\"key",
          "}",
          "}'",
          target_contigs_mapped_present,
          contig_map_file,
          ">",
          rename_map
        )
        
        system(awk_cmd)
        
      } else {
        message("[Warning] Cannot build rename_map.tsv")
      }
      #----------------------------------------------------------
      # NEW: Reheader BAM for IGV
      target_bam <- "./IGV/optitype/optitype_target.bam"
      final_bam  <- "./IGV/optitype/optitype.igv.bam"
      
      header_file <- "./IGV/optitype/header.sam"
      header_renamed <- "./IGV/optitype/header.renamed.sam"
      
      if (file.exists(target_bam) && file.exists(rename_map)) {
        
        # Extract header
        system(paste("samtools view -H", shQuote(target_bam), ">", shQuote(header_file)))
        
        # Rename SN: safely without breaking tabs
        awk_cmd <- paste(
          "awk '",
          "NR==FNR {map[$1]=$2; next}",
          "{",
          "for (old in map) {",
          "gsub(\"SN:\" old, \"SN:\" map[old])",
          "}",
          "print",
          "}'",
          rename_map,
          header_file,
          ">",
          header_renamed
        )
        
        system(awk_cmd)
        
        # Apply new header
        system(paste(
          "samtools reheader",
          shQuote(header_renamed),
          shQuote(target_bam),
          ">",
          shQuote(final_bam)
        ))
        
        system(paste("samtools index", shQuote(final_bam)))
        
      } else {
        message("[Warning] Cannot reheader BAM")
      }
      #----------------------------------------------------------
      # NEW: Build IGV FASTA (renamed)
      igv_fasta_original <- "./IGV/optitype/optitype_hla_reference_igv.fasta"
      final_fasta <- "./IGV/optitype/optitype.igv.fasta"
      tmp_fasta <- "./IGV/optitype/optitype.igv.tmp.fasta"
      
      rename_contigs <- "./IGV/optitype/rename_contigs.txt"
      
      if (file.exists(rename_map) && file.exists(igv_fasta_original)) {
        
        # Extract original contig names
        system(paste("cut -f1", shQuote(rename_map), ">", shQuote(rename_contigs)))
        
        # Extract sequences
        system(paste(
          "samtools faidx",
          shQuote(igv_fasta_original),
          "`cat", rename_contigs, "`",
          ">",
          shQuote(tmp_fasta)
        ))
        
        # Rename FASTA header
        awk_cmd <- paste(
          "awk '",
          "NR==FNR {map[$1]=$2; next}",
          "{",
          "if($0 ~ /^>/){",
          "name=substr($0,2);",
          "print \">\" map[name]",
          "} else {print}",
          "}'",
          rename_map,
          tmp_fasta,
          ">",
          final_fasta
        )
        
        system(awk_cmd)
        
        system(paste("samtools faidx", shQuote(final_fasta)))
        
      } else {
        message("[Warning] Cannot build IGV FASTA")
      }
      #----------------------------------------------------------
      # FINAL CLEANUP (single BAM architecture)
      files_to_remove <- c(
        header_file,
        header_renamed,
        tmp_fasta,
        rename_contigs,
        bam_contigs_mapped,
        target_contigs,
        target_contigs_mapped_present,
        rename_map,
        contig_map_file,
        allele_list_file,
        "./IGV/optitype/optitype_target.bam",
        "./IGV/optitype/optitype_target.bam.bai",
        "./IGV/optitype/optitype_clean.bam",
        "./IGV/optitype/optitype_clean.bam.bai"
      )
      
      for (f in files_to_remove) {
        if (file.exists(f)) file.remove(f)
      }
      #----------------------------------------------------------
      # Extract contig list and caculate coverage/meandepth
      bam   <- "./IGV/optitype/optitype.igv.bam"
      fasta <- "./IGV/optitype/optitype.igv.fasta"
      
      igv_contigs <- "./IGV/optitype/igv_contigs.txt"
      qc_tsv      <- "./IGV/optitype/coverage_tmp.tsv"
      
      # contig list
      system(paste(
        "grep '^>'", shQuote(fasta),
        "| sed 's/^>//' >",
        shQuote(igv_contigs)
      ))
      
      # QC tsv
      system(paste(
        "printf 'HLA\\tCoverage\\tMeanDepth\\n' >", shQuote(qc_tsv), "&&",
        "while IFS= read -r ctg; do",
        "[ -z \"$ctg\" ] && continue;",
        "samtools depth -a -r \"$ctg\"", shQuote(bam),
        "| awk -v C=\"$ctg\" '{sum+=$3; n++; if($3>0) cov++} END{",
        "if(n==0){printf \"%s\\tNA\\tNA\\n\", C}",
        "else{printf \"%s\\t%.6f\\t%.6f\\n\", C, cov/n, sum/n}",
        "}' >>", shQuote(qc_tsv), ";",
        "done <", shQuote(igv_contigs)
      ))
      #----------------------------------------------------------
      # merge to optitype_mhc_table
      qc <- read.table(qc_tsv, sep="\t", header=TRUE, stringsAsFactors=FALSE)
      
      to_igv_key <- function(a){
        parts <- strsplit(a, ":", fixed=TRUE)[[1]]
        a2 <- if (length(parts) >= 2) paste(parts[1], parts[2], sep=":") else a
        paste0("HLA-", gsub(":", "-", a2, fixed=TRUE))
      }
      
      optitype_mhc_table$IGV_KEY <- vapply(optitype_mhc_table$Allele, to_igv_key, character(1))
      
      qc <- read.table(qc_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      
      need_cols <- c("HLA", "Coverage", "MeanDepth")
      if (!all(need_cols %in% colnames(qc))) {
        stop("[Error] qc_tsv header invalid. Got: ", paste(colnames(qc), collapse = ", "))
      }
      
      qc <- qc[!duplicated(qc$HLA), ]
      
      to_igv_key <- function(a){
        parts <- strsplit(a, ":", fixed = TRUE)[[1]]
        a2 <- if (length(parts) >= 2) paste(parts[1], parts[2], sep=":") else a
        paste0("HLA-", gsub(":", "-", a2, fixed = TRUE))
      }
      optitype_mhc_table$IGV_KEY <- vapply(optitype_mhc_table$Allele, to_igv_key, character(1))
      
      m <- match(optitype_mhc_table$IGV_KEY, qc$HLA)
      optitype_mhc_table$Coverage  <- qc$Coverage[m]
      optitype_mhc_table$MeanDepth <- qc$MeanDepth[m]
      
      #----------------------------------------------------------
      # modify coverage (use raw numeric, then format for display)
      optitype_mhc_table$Coverage_raw  <- suppressWarnings(as.numeric(optitype_mhc_table$Coverage))
      optitype_mhc_table$MeanDepth_raw <- suppressWarnings(as.numeric(optitype_mhc_table$MeanDepth))
      
      optitype_mhc_table$Coverage <- ifelse(
        is.na(optitype_mhc_table$Coverage_raw),
        "NA",
        sprintf("%.2f%%", optitype_mhc_table$Coverage_raw * 100)
      )
      
      optitype_mhc_table$MeanDepth <- ifelse(
        is.na(optitype_mhc_table$MeanDepth_raw),
        "NA",
        sprintf("%.2f", optitype_mhc_table$MeanDepth_raw)
      )
      
      optitype_mhc_table$Coverage_color <- ifelse(
        is.na(optitype_mhc_table$Coverage_raw), "na",
        ifelse(optitype_mhc_table$Coverage_raw >= 0.8, "good",
               ifelse(optitype_mhc_table$Coverage_raw >= 0.5, "mid", "bad"))
      )
      

      if ("nucleotide" %in% colnames(optitype_mhc_table)) {
        colnames(optitype_mhc_table)[colnames(optitype_mhc_table)=="nucleotide"] <- "Nucleotide"
      }
      if ("protein" %in% colnames(optitype_mhc_table)) {
        colnames(optitype_mhc_table)[colnames(optitype_mhc_table)=="protein"] <- "Protein"
      }
      
      optitype_mhc_table <- optitype_mhc_table[, c(
        "Allele",
        "Coverage",
        "MeanDepth",
        "Nucleotide",
        "Protein",
        "Coverage_color"
      )]
      #----------------------------------------------------------
      output$hla_typing_table <- renderDataTable({
        datatable(
          optitype_mhc_table,
          rownames = FALSE,
          selection = "none",
          class = "nowrap",
          options = list(
            pageLength = 8,
            autoWidth = TRUE,
            scrollX = TRUE,
            columnDefs = list(
              list(targets = which(colnames(optitype_mhc_table)=="Coverage_color") - 1, visible = FALSE),
              
              list(
                targets = which(colnames(optitype_mhc_table)=="Coverage") - 1,
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if(type !== 'display') return data;",
                  "  var color = row[row.length - 1];",
                  "  var cls = '';",
                  "  if(color === 'good') cls='color:#1b7f3a;';",
                  "  else if(color === 'mid') cls='color:#b8860b;';",
                  "  else if(color === 'bad') cls='color:#b22222;';",
                  "  else cls='color:#666;';",
                  "  return '<span style=\"' + cls + '\">' + data + '</span>';",
                  "}"
                )
              ),
              
              list(
                width = "100px",
                targets = "_all",
                render = JS(
                  "function(data, type, row, meta) {",
                  "return type === 'display' && data != null && data.length > 30 ?",
                  "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                  "}"
                )
              )
            )
          )
        )
      })
      #----------------------------------------------------------
      if (file.exists(qc_tsv)) file.remove(qc_tsv)
      if (file.exists(igv_contigs)) file.remove(igv_contigs)
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
      backup_igv_files("OptiType", new_id)
      
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
      #----------------------------------------------------------
      sequence_file <- NULL
      for (i in 1:nrow(arcashla_mhc_table)) {
        nucleotide <- paste(">", arcashla_mhc_table[i, 1], "_nucleotide", sep = "")
        protein <- paste(">", arcashla_mhc_table[i, 1], "_protein", sep = "")
        seqs <- rbind(rbind(nucleotide, arcashla_mhc_table[i, 2]), rbind(protein, arcashla_mhc_table[i, 3]))
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
    }
    #-----------------------------------------------------------
    if (dir.exists("./IGV/hlahd")) {
      unlink("./IGV/hlahd", recursive = TRUE)
    }
    dir.create("./IGV/hlahd", showWarnings = FALSE)
    #----------------------------------------------------------
    # Merge two SAM files
    system("cat ./Output/hlahd/sample/mapfile/sample.all.R1.pmap.NM.sam ./Output/hlahd/sample/mapfile/sample.all.R2.pmap.NM.sam > ./IGV/hlahd/sample.all.merged.sam")
    
    # Extract headers and move them to the beginning of the file
    system("grep '^@' ./IGV/hlahd/sample.all.merged.sam > ./IGV/hlahd/header.sam")
    system("grep -v '^@' ./IGV/hlahd/sample.all.merged.sam > ./IGV/hlahd/body.sam")
    system("cat ./IGV/hlahd/header.sam ./IGV/hlahd/body.sam > ./IGV/hlahd/sample.all.corrected.merged.sam")
    
    # Keep only the first occurrence of each sequence name
    system("awk '!seen[$0]++' ./IGV/hlahd/sample.all.corrected.merged.sam > ./IGV/hlahd/sample.all.rmdup.sam")
    
    # Convert the merged SAM file to BAM format
    system("samtools view -Sb ./IGV/hlahd/sample.all.rmdup.sam > ./IGV/hlahd/sample.all.rmdup.bam")

    # Sort the merged BAM file
    system("samtools sort -o ./IGV/hlahd/sample.all.sorted.bam ./IGV/hlahd/sample.all.rmdup.bam")
    
    # Remove intermediate files
    file.remove("./IGV/hlahd/sample.all.merged.sam")
    file.remove("./IGV/hlahd/sample.all.corrected.merged.sam")
    file.remove("./IGV/hlahd/sample.all.rmdup.sam")
    file.remove("./IGV/hlahd/sample.all.rmdup.bam")
    file.remove("./IGV/hlahd/header.sam")
    file.remove("./IGV/hlahd/body.sam")
    
    # Extract exon names matched with the reference from the mapping results as a list
    system("samtools view ./IGV/hlahd/sample.all.sorted.bam | cut -f 3 | grep -E '^[ABC]:Exon' > ./IGV/hlahd/hlahd_mhc_i_list.txt")
    system("samtools view ./IGV/hlahd/sample.all.sorted.bam | cut -f 3 | grep -E '^D(PA|PB|QA|QB|RB)1:Exon' > ./IGV/hlahd/hlahd_mhc_ii_list.txt")

    # Create hla_reference.fasta for sample.all.sorted.bam
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
                 sep = " "
    ))
    system(paste("python /root/shiny/Server/hlahd_igv_genome.py",
                 input_list_II,
                 reference_fasta,
                 igv_fasta_II,
                 sep = " "
    ))
    
    # Step 1: Modify BAM file using samtools and sed
    system('samtools view -h ./IGV/hlahd/sample.all.sorted.bam |
       sed "s/A:Exon/A_Exon/g; s/B:Exon/B_Exon/g; s/C:Exon/C_Exon/g; s/DPA1:Exon/DPA1_Exon/g;
            s/DPB1:Exon/DPB1_Exon/g; s/DQA1:Exon/DQA1_Exon/g; s/DQB1:Exon/DQB1_Exon/g;
            s/DRB1:Exon/DRB1_Exon/g" > ./IGV/hlahd/sample.modified.sam')
    
    # Step 2: Convert modified SAM file back to BAM format
    system("samtools view -Sb ./IGV/hlahd/sample.modified.sam > ./IGV/hlahd/sample.modified.bam")
    
    # Step 3: Modify HLA reference I FASTA file by replacing ':' with '_'
    system("mv ./IGV/hlahd/hlahd_hla_reference_I_igv.fasta ./IGV/hlahd/hlahd_hla_reference_I_igv_colon.fasta")
    system('sed "s/:Exon/_Exon/g" ./IGV/hlahd/hlahd_hla_reference_I_igv_colon.fasta > ./IGV/hlahd/hlahd_hla_reference_I_igv.fasta')
    
    # Step 4: Modify HLA reference II FASTA file by replacing ':' with '_'
    system("mv ./IGV/hlahd/hlahd_hla_reference_II_igv.fasta ./IGV/hlahd/hlahd_hla_reference_II_igv_colon.fasta")
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
    file.remove("./IGV/hlahd/sample.all.sorted.bam")
    file.remove("./IGV/hlahd/sample.modified.sam")
    file.remove("./IGV/hlahd/hlahd_hla_reference_I_igv_colon.fasta")
    file.remove("./IGV/hlahd/hlahd_hla_reference_II_igv_colon.fasta")
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
      igv_folder <- "./IGV/hlahd"
      
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
      
      # Remove the previous successful typing folder (igv_folder) if it exists
      if (dir.exists(igv_folder)) {
        unlink(igv_folder, recursive = TRUE)
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
      sequence_file <- NULL
      for (i in 1:nrow(hlahd_mhc_table)) {
        nucleotide <- paste(">", hlahd_mhc_table[i, 1], "_nucleotide", sep = "")
        protein <- paste(">", hlahd_mhc_table[i, 1], "_protein", sep = "")
        seqs <- rbind(rbind(nucleotide, hlahd_mhc_table[i, 2]), rbind(protein, hlahd_mhc_table[i, 3]))
        sequence_file <- rbind(sequence_file, seqs)
      }
      output_file <- paste("./Output/hlahd_", input$sequence, "_", input$imgthla, ".txt", sep = "")
      write(sequence_file, file = output_file)
      #----------------------------------------------------------
      # Step0 (HLA-HD exon-level demo BAM + FASTQ for bowtie2 remap)
      igv_dir <- "./IGV/hlahd"
      dir.create(igv_dir, recursive = TRUE, showWarnings = FALSE)
      
      r1_sam <- "./Output/hlahd/sample/mapfile/sample.all.R1.pmap.NM.sam"
      r2_sam <- "./Output/hlahd/sample/mapfile/sample.all.R2.pmap.NM.sam"
      
      merged_sam <- file.path(igv_dir, "sample.all.singleheader.merged.sam")
      all_sorted_bam <- file.path(igv_dir, "sample.all.sorted.bam")
      
      mod_sam <- file.path(igv_dir, "sample.modified.sam")
      mod_bam <- file.path(igv_dir, "sample.modified.bam")
      
      uniq1_bam <- file.path(igv_dir, "sample.uniq1.bam")
      reads_fq_gz <- file.path(igv_dir, "all.fq.gz")
      
      # (0) sanity
      if (!file.exists(r1_sam) || file.size(r1_sam) == 0) stop("[Step0] Missing R1 sam: ", r1_sam)
      if (!file.exists(r2_sam) || file.size(r2_sam) == 0) stop("[Step0] Missing R2 sam: ", r2_sam)
      
      # (1) single-header merge (header from R1 only; body from R1+R2)
      cmd_merge <- paste0(
        "(grep '^@' ", shQuote(r1_sam), " ; ",
        "grep -v '^@' ", shQuote(r1_sam), " ; ",
        "grep -v '^@' ", shQuote(r2_sam), ") > ",
        shQuote(merged_sam)
      )
      system(cmd_merge)
      
      if (!file.exists(merged_sam) || file.size(merged_sam) == 0) {
        stop("[Step0] single-header merged sam not generated: ", merged_sam)
      }
      
      # (2) merged SAM -> sorted BAM (demo exon-level) + index
      cmd_all_bam <- paste0(
        "samtools view -bS ", shQuote(merged_sam), " | ",
        "samtools sort -o ", shQuote(all_sorted_bam), " -"
      )
      system(cmd_all_bam)
      system(paste("samtools index", shQuote(all_sorted_bam)))
      
      if (!file.exists(all_sorted_bam) || file.size(all_sorted_bam) == 0) {
        stop("[Step0] sample.all.sorted.bam not generated: ", all_sorted_bam)
      }
      
      # (3) modify contig names :Exon -> _Exon (keep your original sed logic)
      cmd_mod <- paste0(
        "samtools view -h ", shQuote(all_sorted_bam), " | ",
        "sed \"s/A:Exon/A_Exon/g; s/B:Exon/B_Exon/g; s/C:Exon/C_Exon/g; ",
        "s/DPA1:Exon/DPA1_Exon/g; s/DPB1:Exon/DPB1_Exon/g; ",
        "s/DQA1:Exon/DQA1_Exon/g; s/DQB1:Exon/DQB1_Exon/g; ",
        "s/DRB1:Exon/DRB1_Exon/g\" > ",
        shQuote(mod_sam)
      )
      system(cmd_mod)
      
      system(paste("samtools view -bS", shQuote(mod_sam), "-o", shQuote(mod_bam)))
      system(paste("samtools index", shQuote(mod_bam)))
      
      if (!file.exists(mod_bam) || file.size(mod_bam) == 0) {
        stop("[Step0] sample.modified.bam not generated: ", mod_bam)
      }
      
      # (4) uniq1 (mapped only + first record per QNAME) -> BAM + index
      cmd_uniq1 <- paste0(
        "samtools view -h -F 4 ", shQuote(mod_bam), " | ",
        "awk 'BEGIN{FS=OFS=\"\\t\"} /^@/ {print; next} !seen[$1]++ {print}' | ",
        "samtools view -b -o ", shQuote(uniq1_bam), " -"
      )
      system(cmd_uniq1)
      system(paste("samtools index", shQuote(uniq1_bam)))
      
      if (!file.exists(uniq1_bam) || file.size(uniq1_bam) == 0) {
        stop("[Step0] sample.uniq1.bam not generated: ", uniq1_bam)
      }
      
      # (5) uniq1 -> single-end FASTQ.gz for bowtie2 (-U)
      cmd_fastq <- paste0(
        "samtools fastq -n ", shQuote(uniq1_bam),
        " | gzip -c > ", shQuote(reads_fq_gz)
      )
      system(cmd_fastq)
      
      if (!file.exists(reads_fq_gz) || file.size(reads_fq_gz) == 0) {
        stop("[Step0] all.fq.gz not generated: ", reads_fq_gz)
      } else {
        message("[Step0] reads FASTQ ready: ", reads_fq_gz)
      }
      #----------------------------------------------------------
      # Step1 (HLA-HD IGV prep): typed_alleles.txt + nucleotide.raw.fasta
      # Paths (Shiny paths)
      hlahd_txt   <- paste0("./Output/hlahd_", input$sequence, "_", input$imgthla, ".txt")
      
      igv_dir     <- "./IGV/hlahd"
      alleles_txt <- file.path(igv_dir, "typed_alleles.txt")
      nuc_raw_fa  <- file.path(igv_dir, "hlahd.nucleotide.raw.fasta")
      
      rename_map  <- file.path(igv_dir, "hlahd.rename_map.tsv")
      fa_all      <- file.path(igv_dir, "hlahd.igv.all.fasta")
      fa_i        <- file.path(igv_dir, "hlahd.igv.classI.fasta")
      fa_ii       <- file.path(igv_dir, "hlahd.igv.classII.fasta")
      
      # [Step0 outputs sanity check] - Step0 should already generate these
      step0_mod_bam  <- file.path(igv_dir, "sample.modified.bam")
      step0_reads_fq <- file.path(igv_dir, "all.fq.gz")
      
      if (!file.exists(step0_mod_bam) || file.size(step0_mod_bam) == 0) {
        message("[Info] Step0 demo BAM not found yet: ", step0_mod_bam,
                " (OK if you haven't run Step0, but Step2 will need reads FASTQ.)")
      }
      if (!file.exists(step0_reads_fq) || file.size(step0_reads_fq) == 0) {
        message("[Warning] Step0 reads FASTQ missing: ", step0_reads_fq,
                " (Step2 bowtie2 will FAIL without this.)")
      }
      
      # Safety: ensure IGV folder exists
      dir.create(igv_dir, recursive = TRUE, showWarnings = FALSE)
      
      if (!file.exists(hlahd_txt) || file.size(hlahd_txt) == 0) {
        message("[Warning] HLA-HD output txt missing or empty: ", hlahd_txt)
      } else {
        
        # 1) typed_alleles.txt (unique allele list)
        cmd_alleles <- paste(
          "grep '^>'", shQuote(hlahd_txt),
          "| sed 's/^>//'",
          "| sed 's/_nucleotide$//'",
          "| sed 's/_protein$//'",
          "| sort -u >", shQuote(alleles_txt)
        )
        system(cmd_alleles)
        
        # 2) nucleotide-only FASTA (raw, headers are original: A*01:01:01 etc.)
        #    keep only _nucleotide records and strip suffix
        cmd_nuc <- paste(
          "awk '",
          "BEGIN{keep=0}",
          "(/^>/){",
          "  keep = ($0 ~ /_nucleotide$/);",
          "  if(keep){ h=$0; sub(/_nucleotide$/,\"\",h); print h }",
          "  next",
          "}",
          "{ if(keep) print }",
          "'",
          shQuote(hlahd_txt),
          ">",
          shQuote(nuc_raw_fa)
        )
        system(cmd_nuc)
        
        # 3) rename map: old allele -> IGV contig (HLA- + ':' -> '-')
        cmd_map <- paste(
          "awk 'BEGIN{FS=OFS=\"\\t\"}",
          "/^>/{old=substr($0,2); new=\"HLA-\" old; gsub(/:/,\"-\",new); print old, new}'",
          shQuote(nuc_raw_fa),
          ">",
          shQuote(rename_map)
        )
        system(cmd_map)
        
        # 4) build IGV FASTA (ALL) with renamed headers
        cmd_fa_all <- paste0(
          "awk -v MAP=", shQuote(rename_map),
          " 'BEGIN{FS=OFS=\"\\t\"; while((getline < MAP)>0){m[$1]=$2} close(MAP); keep=0}",
          " /^>/{old=substr($0,2); if(old in m){print \">\" m[old]; keep=1}else{keep=0}; next}",
          " {if(keep) print}' ",
          shQuote(nuc_raw_fa),
          " > ",
          shQuote(fa_all)
        )
        rc <- system(cmd_fa_all)
        if (rc != 0) stop("[HLA-HD IGV Step1] awk rename->fasta failed")
        
        # 5) split Class I / Class II FASTA
        cmd_fa_i <- paste(
          "awk '",
          "/^>/{name=substr($0,2); if(name ~ /^HLA-(A|B|C)\\*/) keep=1; else keep=0}",
          "{if(keep) print}",
          "'",
          shQuote(fa_all),
          ">",
          shQuote(fa_i)
        )
        system(cmd_fa_i)
        
        cmd_fa_ii <- paste(
          "awk '",
          "/^>/{name=substr($0,2); if(name ~ /^HLA-(DPA1|DPB1|DQA1|DQB1|DRB1)\\*/) keep=1; else keep=0}",
          "{if(keep) print}",
          "'",
          shQuote(fa_all),
          ">",
          shQuote(fa_ii)
        )
        system(cmd_fa_ii)
        
        # 6) faidx for IGV
        if (file.exists(fa_all) && file.size(fa_all) > 0) system(paste("samtools faidx", shQuote(fa_all)))
        if (file.exists(fa_i)   && file.size(fa_i)   > 0) system(paste("samtools faidx", shQuote(fa_i)))
        if (file.exists(fa_ii)  && file.size(fa_ii)  > 0) system(paste("samtools faidx", shQuote(fa_ii)))
        
        # Optional quick sanity logs
        if (file.exists(rename_map)) {
          message("[Info] rename_map rows: ", length(readLines(rename_map)))
        }
        if (file.exists(fa_i)) {
          message("[Info] Class I contigs: ", sum(grepl("^>", readLines(fa_i))))
        }
        if (file.exists(fa_ii)) {
          message("[Info] Class II contigs: ", sum(grepl("^>", readLines(fa_ii))))
        }
      }
      # [Hard checks before Step2]
      if (!file.exists(fa_all) || file.size(fa_all) == 0) {
        stop("[HLA-HD IGV Step1] Failed to build IGV FASTA: ", fa_all)
      }
      if (!file.exists(step0_reads_fq) || file.size(step0_reads_fq) == 0) {
        stop("[HLA-HD IGV Step1] Missing reads FASTQ for Step2: ", step0_reads_fq,
             " (Run Step0 to generate all.fq.gz)")
      }
      #----------------------------------------------------------
      # Step2: bowtie2 align reads to allele-centric genome FASTA
      igv_dir   <- "./IGV/hlahd"
      igv_fasta <- file.path(igv_dir, "hlahd.igv.all.fasta")
      reads_fq  <- file.path(igv_dir, "all.fq.gz")  # <-- if your Shiny pipeline stores reads elsewhere, change here
      
      idx_prefix <- file.path(igv_dir, "hlahd_WES_igv")
      raw_bam    <- file.path(igv_dir, "allelecentric.raw.bam")
      
      # Safety: ensure dir exists
      if (!dir.exists(igv_dir)) {
        dir.create(igv_dir, recursive = TRUE, showWarnings = FALSE)
      }
      
      # Basic checks
      if (!file.exists(igv_fasta) || file.size(igv_fasta) == 0) {
        stop("[HLA-HD IGV Step2] Missing IGV FASTA: ", igv_fasta)
      }
      if (!file.exists(reads_fq) || file.size(reads_fq) == 0) {
        stop("[HLA-HD IGV Step2] Missing reads FASTQ: ", reads_fq,
             " (You need to put/copy your all.fq.gz into ./IGV/hlahd/ before running Step2)")
      }
      
      # Clean old bowtie2 index files (avoid mismatch with new FASTA)
      old_idx <- Sys.glob(paste0(idx_prefix, ".*.bt2*"))
      if (length(old_idx) > 0) file.remove(old_idx)
      
      # 2-1) Build bowtie2 index (ok to overwrite/rebuild)
      build_cmd <- paste(
        "bowtie2-build",
        shQuote(igv_fasta),
        shQuote(idx_prefix)
      )
      
      rc_build <- system(build_cmd)
      if (rc_build != 0) {
        stop("[HLA-HD IGV Step2] bowtie2-build failed. cmd=", build_cmd)
      }
      
      # 2-2) Align (keep all alignments with -a), sort to BAM
      # Note: capture bowtie2 log into file to help debugging
      bowtie2_log <- file.path(igv_dir, "bowtie2.step2.log")
      bt2_threads <- min(as.integer(num_cores), 8)
      align_cmd <- paste(
        "bowtie2 --end-to-end --very-sensitive -a -p", bt2_threads,
        "-x", shQuote(idx_prefix),
        "-U", shQuote(reads_fq),
        "2>", shQuote(bowtie2_log),
        "| samtools sort -o", shQuote(raw_bam), "-"
      )
      
      rc_align <- system(align_cmd)
      if (rc_align != 0) {
        stop("[HLA-HD IGV Step2] bowtie2 alignment failed. Check log: ", bowtie2_log)
      }      
      # 2-3) Index BAM
      if (!file.exists(raw_bam) || file.size(raw_bam) == 0) {
        stop("[HLA-HD IGV Step2] raw bam not produced/empty: ", raw_bam,
             "\nCheck log: ", bowtie2_log)
      }
      rc_index <- system(paste("samtools index", shQuote(raw_bam)))
      if (rc_index != 0 || !file.exists(paste0(raw_bam, ".bai"))) {
        stop("[HLA-HD IGV Step2] samtools index failed: ", raw_bam)
      }
      
      # 2-4) Ensure FASTA index exists (for IGV genome track usage)
      if (!file.exists(paste0(igv_fasta, ".fai"))) {
        system(paste("samtools faidx", shQuote(igv_fasta)))
      }
      
      # Optional: quick sanity prints into Shiny log
      if (file.exists(raw_bam)) {
        message("[HLA-HD IGV Step2] raw bam created: ", raw_bam)
        # mapped contigs >0 (first few)
        idxstats_out <- tryCatch(system(paste("samtools idxstats", shQuote(raw_bam), "| awk '$3>0{print $1\"\\t\"$3}' | head"), intern = TRUE),
                                 error = function(e) character(0))
        if (length(idxstats_out) > 0) {
          message("[HLA-HD IGV Step2] mapped contigs head:\n", paste(idxstats_out, collapse = "\n"))
        }
      } else {
        stop("[HLA-HD IGV Step2] raw bam not produced: ", raw_bam,
             "\nCheck log: ", bowtie2_log)
      }
      #----------------------------------------------------------
      # Step3: Remove secondary + supplementary (keep primary only)
      igv_dir <- "./IGV/hlahd"
      raw_bam <- file.path(igv_dir, "allelecentric.raw.bam")
      primary_bam <- file.path(igv_dir, "allelecentric.primaryOnly.bam")
      
      if (!file.exists(raw_bam) || file.size(raw_bam) == 0) {
        stop("[HLA-HD IGV Step3] Missing raw bam: ", raw_bam)
      }
      
      sort_threads <- min(as.integer(num_cores), 8)
      cmd_step3 <- paste(
        "samtools view -h -F 0x900", shQuote(raw_bam),
        "| samtools sort -@", sort_threads, "-o", shQuote(primary_bam), "-"
      )
      
      rc3 <- system(cmd_step3)
      if (rc3 != 0) {
        stop("[HLA-HD IGV Step3] samtools view/sort failed. cmd=", cmd_step3)
      }
      
      rc3i <- system(paste("samtools index", shQuote(primary_bam)))
      if (rc3i != 0 || !file.exists(paste0(primary_bam, ".bai"))) {
        stop("[HLA-HD IGV Step3] samtools index failed: ", primary_bam)
      }
      
      if (!file.exists(primary_bam) || file.size(primary_bam) == 0) {
        stop("[HLA-HD IGV Step3] primaryOnly bam not produced: ", primary_bam)
      } else {
        message("[HLA-HD IGV Step3] primaryOnly bam created: ", primary_bam)
      }
      #----------------------------------------------------------
      # Step4: Keep 1 best alignment per read (best NM, tie-breaker MAPQ)
      primary_bam <- "./IGV/hlahd/allelecentric.primaryOnly.bam"
      bestnm_bam  <- "./IGV/hlahd/allelecentric.bestNM.bam"
      
      if (file.exists(primary_bam) && file.size(primary_bam) > 0) {
        awk_prog <- paste(
          "BEGIN{FS=OFS=\"\\t\"; cur=\"\"; best_line=\"\"; best_nm=999999; best_mapq=-1}",
          "/^@/ {print; next}",
          "{",
          "  q=$1;",
          "  nm=999999;",
          "  for(i=12;i<=NF;i++){ if($i ~ /^NM:i:/){ nm=substr($i,6)+0; break } }",
          "  mapq=$5+0;",
          "  if(cur==\"\") { cur=q }",
          "  if(q!=cur){",
          "    if(best_line!=\"\") print best_line;",
          "    cur=q; best_line=$0; best_nm=nm; best_mapq=mapq;",
          "  } else {",
          "    if(nm < best_nm || (nm==best_nm && mapq > best_mapq)){",
          "      best_line=$0; best_nm=nm; best_mapq=mapq;",
          "    }",
          "  }",
          "}",
          "END{ if(best_line!=\"\") print best_line }",
          sep="\n"
        )
        
        sort_threads <- min(as.integer(num_cores), 8)
        
        cmd <- paste0(
          "samtools sort -n -@", sort_threads, " -O SAM ", shQuote(primary_bam), " | ",
          "awk ", shQuote(awk_prog), " | ",
          "samtools view -b - | ",
          "samtools sort -@", sort_threads, " -o ", shQuote(bestnm_bam), " -"
        )
        rc4 <- system(cmd)
        if (rc4 != 0) stop("[HLA-HD IGV Step4] command failed: ", cmd)
        
        # index
        if (file.exists(bestnm_bam) && file.size(bestnm_bam) > 0) {
          rc4i <- system(paste("samtools index", shQuote(bestnm_bam)))
          if (rc4i != 0 || !file.exists(paste0(bestnm_bam, ".bai"))) {
            stop("[HLA-HD IGV Step4] samtools index failed: ", bestnm_bam)
          }
        } else {
          warning("[HLA-HD IGV][Step4] allelecentric.bestNM.bam not generated or empty.")
        }
      } else {
        warning("[HLA-HD IGV][Step4] primaryOnly BAM missing/empty: ", primary_bam)
      }
      #----------------------------------------------------------
      # Step5: CR filter (NM/readlen <= 0.05) on bestNM BAM
      bestnm_bam <- "./IGV/hlahd/allelecentric.bestNM.bam"
      cr_bam     <- "./IGV/hlahd/allelecentric.CR.bam"
      
      if (file.exists(bestnm_bam) && file.size(bestnm_bam) > 0) {
        
        awk_prog <- paste(
          "BEGIN{FS=OFS=\"\\t\"}",
          "/^@/ {print; next}",
          "{",
          "  nm=999999;",
          "  for(i=12;i<=NF;i++){ if($i ~ /^NM:i:/){ nm=substr($i,6); break } }",
          "  readlen=length($10);",
          "  if(readlen>0 && (nm/readlen) <= 0.05) print;",
          "}",
          sep = "\n"
        )
        
        
        sort_threads <- min(as.integer(num_cores), 8)
        
        cmd <- paste0(
          "samtools view -h ", shQuote(bestnm_bam), " | ",
          "awk ", shQuote(awk_prog), " | ",
          "samtools view -b - | ",
          "samtools sort -@", sort_threads, " -o ", shQuote(cr_bam), " -"
        )
        
        rc5 <- system(cmd)
        if (rc5 != 0) stop("[HLA-HD IGV Step5] command failed: ", cmd)
        
        if (file.exists(cr_bam) && file.size(cr_bam) > 0) {
          rc5i <- system(paste("samtools index", shQuote(cr_bam)))
          if (rc5i != 0 || !file.exists(paste0(cr_bam, ".bai"))) {
            stop("[HLA-HD IGV Step5] samtools index failed: ", cr_bam)
          }
          message("[HLA-HD IGV][Step5] CR BAM generated: ", cr_bam)
        } else {
          warning("[HLA-HD IGV][Step5] allelecentric.CR.bam not generated or empty.")
        }
        
      } else {
        warning("[HLA-HD IGV][Step5] bestNM BAM missing/empty: ", bestnm_bam)
      }
      #----------------------------------------------------------
      # Step6: QC (Coverage / MeanDepth) on allele-centric BAM
      bam   <- file.path(igv_dir, "allelecentric.CR.bam")
      fasta <- file.path(igv_dir, "hlahd.igv.all.fasta")
      
      contigs_all <- file.path(igv_dir, "igv_contigs_all.txt")
      bam_mapped  <- file.path(igv_dir, "bam_contigs_mapped.txt")
      contigs_qc  <- file.path(igv_dir, "igv_contigs_for_qc.txt")
      qc_tsv      <- file.path(igv_dir, "coverage_tmp.tsv")
      
      # Safety checks
      if (!file.exists(bam) || file.size(bam) == 0) stop("[HLA-HD QC] Missing BAM: ", bam)
      if (!file.exists(fasta) || file.size(fasta) == 0) stop("[HLA-HD QC] Missing FASTA: ", fasta)
      
      # 1) contigs from FASTA (truth set)
      system(paste(
        "grep '^>'", shQuote(fasta),
        "| sed 's/^>//' >",
        shQuote(contigs_all)
      ))
      
      # 2) contigs with mapped reads > 0 from BAM (safety filter)
      system(paste(
        "samtools idxstats", shQuote(bam),
        "| awk -F'\\t' '$3>0{print $1}' >",
        shQuote(bam_mapped)
      ))
      
      # 3) intersect => contigs used for QC
      system(paste(
        "grep -Fxf", shQuote(bam_mapped), shQuote(contigs_all),
        ">",
        shQuote(contigs_qc)
      ))
      
      # sanity logs (optional)
      message("[Info] contigs in FASTA: ", length(readLines(contigs_all)))
      message("[Info] contigs mapped>0 in BAM: ", length(readLines(bam_mapped)))
      message("[Info] contigs used for QC: ", length(readLines(contigs_qc)))
      
      # 4) compute Coverage/MeanDepth per contig
      system(paste(
        "printf 'Contig\\tCoverage\\tMeanDepth\\n' >", shQuote(qc_tsv), "&&",
        "while IFS= read -r ctg; do",
        "  [ -z \"$ctg\" ] && continue;",
        "  samtools depth -a -r \"$ctg\"", shQuote(bam),
        "  | awk -v C=\"$ctg\" '{sum+=$3; n++; if($3>0) cov++} END{",
        "      if(n==0){printf \"%s\\tNA\\tNA\\n\", C}",
        "      else{printf \"%s\\t%.6f\\t%.6f\\n\", C, cov/n, sum/n}",
        "    }' >>", shQuote(qc_tsv), ";",
        "done <", shQuote(contigs_qc)
      ))
      
      if (!file.exists(qc_tsv) || file.size(qc_tsv) == 0) {
        stop("[HLA-HD QC] qc_tsv not generated: ", qc_tsv)
      }
      
      # 5) merge QC back to hlahd_mhc_table
      qc <- read.table(qc_tsv, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
      
      need_cols <- c("Contig", "Coverage", "MeanDepth")
      if (!all(need_cols %in% colnames(qc))) {
        stop("[HLA-HD QC] qc_tsv header invalid. Got: ", paste(colnames(qc), collapse=","))
      }
      qc <- qc[!duplicated(qc$Contig), ]
      
      # table allele -> contig key：A*01:01:01 -> HLA-A*01-01-01
      to_igv_key_hlahd <- function(a){
        paste0("HLA-", gsub(":", "-", a, fixed=TRUE))
      }
      hlahd_mhc_table$IGV_KEY <- vapply(hlahd_mhc_table$Allele, to_igv_key_hlahd, character(1))
      
      m <- match(hlahd_mhc_table$IGV_KEY, qc$Contig)
      hlahd_mhc_table$Coverage  <- qc$Coverage[m]
      hlahd_mhc_table$MeanDepth <- qc$MeanDepth[m]
      
      # 6) formatting + color class
      hlahd_mhc_table$Coverage_raw  <- suppressWarnings(as.numeric(hlahd_mhc_table$Coverage))
      hlahd_mhc_table$MeanDepth_raw <- suppressWarnings(as.numeric(hlahd_mhc_table$MeanDepth))
      
      hlahd_mhc_table$Coverage <- ifelse(
        is.na(hlahd_mhc_table$Coverage_raw),
        "NA",
        sprintf("%.2f%%", hlahd_mhc_table$Coverage_raw * 100)
      )
      hlahd_mhc_table$MeanDepth <- ifelse(
        is.na(hlahd_mhc_table$MeanDepth_raw),
        "NA",
        sprintf("%.2f", hlahd_mhc_table$MeanDepth_raw)
      )
      
      hlahd_mhc_table$Coverage_color <- ifelse(
        is.na(hlahd_mhc_table$Coverage_raw), "na",
        ifelse(hlahd_mhc_table$Coverage_raw >= 0.8, "good",
               ifelse(hlahd_mhc_table$Coverage_raw >= 0.5, "mid", "bad"))
      )
      #----------------------------------------------------------
      # --- normalize column names for display (HLA-HD) ---
      if ("nucleotide" %in% colnames(hlahd_mhc_table) && !("Nucleotide" %in% colnames(hlahd_mhc_table))) {
        colnames(hlahd_mhc_table)[colnames(hlahd_mhc_table) == "nucleotide"] <- "Nucleotide"
      }
      if ("protein" %in% colnames(hlahd_mhc_table) && !("Protein" %in% colnames(hlahd_mhc_table))) {
        colnames(hlahd_mhc_table)[colnames(hlahd_mhc_table) == "protein"] <- "Protein"
      }
      #----------------------------------------------------------
      show_cols <- c("Allele", "Coverage", "MeanDepth", "Nucleotide", "Protein", "Coverage_color")
      show_cols <- show_cols[show_cols %in% colnames(hlahd_mhc_table)]
      
      if ("Coverage_color" %in% show_cols) {
        show_cols <- c(setdiff(show_cols, "Coverage_color"), "Coverage_color")
      }
      hlahd_mhc_table <- hlahd_mhc_table[, show_cols, drop = FALSE]
      #----------------------------------------------------------
      output$hla_typing_table <- renderDataTable({
        datatable(
          hlahd_mhc_table,
          rownames = FALSE,
          selection = "none",
          class = "nowrap",
          options = list(
            pageLength = 8,
            autoWidth = TRUE,
            scrollX = TRUE,
            columnDefs = list(
              list(
                targets = which(colnames(hlahd_mhc_table) == "Coverage_color") - 1,
                visible = FALSE
              ),
              list(
                targets = which(colnames(hlahd_mhc_table) == "Coverage") - 1,
                render = JS(
                  "function(data, type, row, meta) {",
                  "  if(type !== 'display') return data;",
                  "  var color = row[row.length - 1];",
                  "  var cls = '';",
                  "  if(color === 'good') cls='color:#1b7f3a;';",
                  "  else if(color === 'mid') cls='color:#b8860b;';",
                  "  else if(color === 'bad') cls='color:#b22222;';",
                  "  else cls='color:#666;';",
                  "  return '<span style=\"' + cls + '\">' + data + '</span>';",
                  "}"
                )
              ),
              
              list(
                width = "100px",
                targets = "_all",
                render = JS(
                  "function(data, type, row, meta) {",
                  "return type === 'display' && data != null && data.length > 30 ?",
                  "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                  "}"
                )
              )
            )
          )
        )
      })
      #----------------------------------------------------------
      # [CLEANUP] remove HLA-HD IGV intermediates in ./IGV/hlahd
      # --- keep list (final deliverables + your preferred originals) ---
      keep_files <- c(
        # final allele-centric references for IGV
        file.path(igv_dir, "hlahd.igv.classI.fasta"),
        file.path(igv_dir, "hlahd.igv.classI.fasta.fai"),
        file.path(igv_dir, "hlahd.igv.classII.fasta"),
        file.path(igv_dir, "hlahd.igv.classII.fasta.fai"),
        
        # final BAM for IGV
        file.path(igv_dir, "allelecentric.CR.bam"),
        file.path(igv_dir, "allelecentric.CR.bam.bai"),
        
        # keep Step0 "initial IGV processing成果" (as you listed)
        file.path(igv_dir, "sample.modified.bam"),
        file.path(igv_dir, "sample.modified.bam.bai"),
        
        # keep old/reference FASTA (if you still produce these filenames)
        file.path(igv_dir, "hlahd_hla_reference_I_igv.fasta"),
        file.path(igv_dir, "hlahd_hla_reference_I_igv.fasta.fai"),
        file.path(igv_dir, "hlahd_hla_reference_II_igv.fasta"),
        file.path(igv_dir, "hlahd_hla_reference_II_igv.fasta.fai"),
        
        # (recommended) keep for reload/debug
        file.path(igv_dir, "typed_alleles.txt"),
        file.path(igv_dir, "hlahd.rename_map.tsv")
      )
      
      # keep only those that exist (avoid setdiff issues)
      keep_files <- unique(keep_files[file.exists(keep_files)])
      
      # list all top-level files under ./IGV/hlahd (do not recurse)
      all_files <- list.files(igv_dir, full.names = TRUE, recursive = FALSE)
      
      # remove everything not in keep list
      to_remove <- setdiff(all_files, keep_files)
      
      # additionally remove bowtie2 index (always intermediates)
      bt2_idx <- Sys.glob(file.path(igv_dir, "hlahd_WES_igv*.bt2*"))
      to_remove <- unique(c(to_remove, bt2_idx))
      to_remove <- to_remove[file.exists(to_remove)]
      
      if (length(to_remove) > 0) {
        ok_rm <- file.remove(to_remove)
        if (any(!ok_rm)) {
          warning("[Cleanup] Failed to remove some files:\n", paste(to_remove[!ok_rm], collapse = "\n"))
        } else {
          message("[Cleanup] Removed intermediates: ", length(to_remove), " files.")
        }
      } else {
        message("[Cleanup] Nothing to remove.")
      }
      
      # optional sanity: ensure key outputs still exist before backup
      must_exist <- c(
        file.path(igv_dir, "hlahd.igv.classI.fasta"),
        file.path(igv_dir, "hlahd.igv.classI.fasta.fai"),
        file.path(igv_dir, "hlahd.igv.classII.fasta"),
        file.path(igv_dir, "hlahd.igv.classII.fasta.fai"),
        file.path(igv_dir, "allelecentric.CR.bam"),
        file.path(igv_dir, "allelecentric.CR.bam.bai")
      )
      missing <- must_exist[!file.exists(must_exist)]
      if (length(missing) > 0) {
        stop("[Cleanup] Missing required final IGV files after cleanup:\n", paste(missing, collapse = "\n"))
      }
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
      backup_igv_files("HLA-HD", new_id)
      
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
    igv_fa <- file.path(imgt_dir, paste0("IMGTHLA_", imgt_version_number, "_", seq_suffix, "_seq_igv.fa"))
    #----------------------------------------------------------
    base_output_dir <- file.path(WORK_ROOT, "Output", "t1k")
    igv_out_dir <- file.path(WORK_ROOT, "IGV", "t1k")
    
    if (dir.exists(base_output_dir)) unlink(base_output_dir, recursive = TRUE)
    dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (dir.exists(igv_out_dir)) unlink(igv_out_dir, recursive = TRUE)
    dir.create(igv_out_dir, recursive = TRUE, showWarnings = FALSE)
    
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
    fa1 <- file.path(base_output_dir, "t1k_aligned_1.fa")
    fa2 <- file.path(base_output_dir, "t1k_aligned_2.fa")
    
    system(paste("/T1K/replace_fasta_colon.sh", origin_fa, igv_fa))
    system(paste("bwa index", igv_fa))
    
    fq1 <- file.path(igv_out_dir, "t1k_aligned_1.fq")
    fq2 <- file.path(igv_out_dir, "t1k_aligned_2.fq")
    convert_fa_to_fq(fa1, fq1)
    convert_fa_to_fq(fa2, fq2)
    
    bam_path <- file.path(igv_out_dir, "t1k_mapping.bam")
    system(paste(
      "bwa mem -t", num_cores, igv_fa, fq1, fq2,
      "| samtools sort -@", num_cores, "-o", bam_path
    ))
    system(paste("samtools index", bam_path))
    
    system(paste("/T1K/split_hla_fasta.sh", igv_fa, igv_out_dir))
    classI_fa <- file.path(igv_out_dir, "hla_classI_ABC.fasta")
    classII_fa <- file.path(igv_out_dir, "hla_classII_D.fasta")
    system(paste("samtools faidx", classI_fa))
    system(paste("samtools faidx", classII_fa))
    #----------------------------------------------------------
    unlink(c(fq1, fq2))
    file.remove(paste0(igv_fa, c("", ".amb", ".ann", ".bwt", ".pac", ".sa")))
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
    output$hla_typing_table <- renderDataTable({
      datatable(t1k_mhc_table,
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
    #----------------------------------------------------------
    sequence_file <- NULL
    
    for (i in 1:nrow(t1k_mhc_table)) {
      seq <- t1k_mhc_table[i, "Allele"]
      nuc <- t1k_mhc_table[i, "nucleotide"]
      prot <- t1k_mhc_table[i, "protein"]
      
      short_seq <- sub("^HLA-", "", seq)
      
      nucleotide_header <- paste0(">", short_seq, "_nucleotide")
      protein_header <- paste0(">", short_seq, "_protein")
      
      seqs <- rbind(c(nucleotide_header), c(nuc), c(protein_header), c(prot))
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
    backup_igv_files("T1K", new_id)
    
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
  
  # Get current JobID
  jid <- current_job_id()
  
  # Record the current time, jobid and execution time
  log_entry <- paste(
    getCurrentTime(),
    "| JobID=", jid,
    "|", input$sequence,
    "|", input$imgthla,
    "|", input$package,
    "| Execution:", execution_time, "s"
  )
  
  # Append the new execution time to the existing log (newest on top)
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
