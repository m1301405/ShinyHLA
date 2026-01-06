ori_dir <- getwd()

HISTORY_ROOT <- "/root/shiny/History"
dir.create(HISTORY_ROOT, recursive = TRUE, showWarnings = FALSE)

tmp_dir <- tempdir()
IGV_TMP_ROOT <- file.path(tmp_dir, "IGV")

directory_o <- file.path(IGV_TMP_ROOT, "optitype")
directory_h <- file.path(IGV_TMP_ROOT, "hlahd")
directory_t <- file.path(IGV_TMP_ROOT, "t1k")

#----------------------------------------------------

get_active_igv_dir <- function(tool) {
  jid <- current_job_id()
  
  if (is.null(jid) || jid %in% c("No Job Active", "New Analysis", "")) {
    return(switch(tool,
                  "OptiType" = directory_o,
                  "HLA-HD"   = directory_h,
                  "T1K"      = directory_t
    ))
  }
  
  file.path(HISTORY_ROOT, jid, "IGV")
}


backup_igv_files <- function(tool, job_id) {
  source_dir <- switch(tool,
                       "OptiType" = directory_o,
                       "HLA-HD"   = directory_h,
                       "T1K"      = directory_t
  )
  
  target_dir <- file.path(HISTORY_ROOT, job_id, "IGV")
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  files_to_move <- list.files(source_dir, full.names = TRUE)
  if (length(files_to_move) == 0) {
    return(invisible(FALSE))
  }
  
  ok <- file.rename(files_to_move, file.path(target_dir, basename(files_to_move)))
  if (any(!ok)) {
    need_copy <- files_to_move[!ok]
    file.copy(need_copy, target_dir, overwrite = TRUE)
    file.remove(need_copy)
  }
  
  invisible(TRUE)
}

#----------------------------------------------------
options_values <- reactiveValues(
  optitype_igv_options = NULL,
  hlahd_igv_i_options = NULL,
  hlahd_igv_ii_options = NULL,
  t1k_igv_i_options = NULL,
  t1k_igv_ii_options = NULL
)

#----------------------------------------------------
hg38_igv_i_options <- parseAndValidateGenomeSpec(
  genomeName = "hg38",
  initialLocus = c("HLA-A", "HLA-B", "HLA-C")
)
hg38_igv_ii_options <- parseAndValidateGenomeSpec(
  genomeName = "hg38",
  initialLocus = c("HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1")
)

load_and_show_region_hg38_I <- function(igv_id) {
  bam_file <- file.path(tmp_dir, "BAM", "sample.bam")
  if (file.exists(bam_file)) {
    x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
    loadBamTrackFromLocalData(session, id = igv_id, trackName = "Bam file", data = x)
    showGenomicRegion(session, id = igv_id, region = hg38_igv_i_options)
  }
}
load_and_show_region_hg38_II <- function(igv_id) {
  bam_file <- file.path(tmp_dir, "BAM", "sample.bam")
  if (file.exists(bam_file)) {
    x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
    loadBamTrackFromLocalData(session, id = igv_id, trackName = "bam file", data = x)
    showGenomicRegion(session, id = igv_id, region = hg38_igv_ii_options)
  }
}

output$hg38_igv_i <- renderIgvShiny({
  igvShiny(hg38_igv_i_options)
})
output$hg38_igv_ii <- renderIgvShiny({
  igvShiny(hg38_igv_ii_options)
})

#----------------------------------------------------
observeEvent(
  {
    input$igv_reference
    package_value()
  },
  {
    req(input$igv_reference >= 1)
    
    w <- Waiter$new(html = spin_loader(), color = transparent(0.8))
    w$show()
    
    updateTabsetPanel(session, "igv_tabsetPanel", "Typing Package Results")
    
    current_tool <- package_value()
    if (is.null(current_tool) || current_tool == "") {
      w$hide()
      return()
    }
    
    active_dir <- get_active_igv_dir(current_tool)
    
    if (current_tool == "OptiType") {
      fasta_path_o <- file.path(active_dir, "optitype_hla_reference_igv.fasta")
      fastaindex_path_o <- file.path(active_dir, "optitype_hla_reference_igv.fasta.fai")
      
      if (!file.exists(fasta_path_o)) {
        shinyalert("Error", "Required FASTA file is missing for OptiType in current directory.", type = "error")
        w$hide()
        return()
      }
      
      options_values$optitype_igv_options <- parseAndValidateGenomeSpec(
        genomeName = "optitype", initialLocus = "All", stockGenome = FALSE, dataMode = "localFiles",
        fasta = fasta_path_o, fastaIndex = fastaindex_path_o
      )
      output$alignment_igv <- renderIgvShiny({
        igvShiny(options_values$optitype_igv_options)
      })
      shinyalert("Success", "OptiType reference uploaded successfully.", type = "success")
    } else if (current_tool == "HLA-HD") {
      # HLA-HD I & II
      fasta_path_h_i <- file.path(active_dir, "hlahd_hla_reference_I_igv.fasta")
      fastaindex_path_h_i <- file.path(active_dir, "hlahd_hla_reference_I_igv.fasta.fai")
      fasta_path_h_ii <- file.path(active_dir, "hlahd_hla_reference_II_igv.fasta")
      fastaindex_path_h_ii <- file.path(active_dir, "hlahd_hla_reference_II_igv.fasta.fai")
      
      if (!file.exists(fasta_path_h_i)) {
        shinyalert("Error", "Required FASTA file is missing for HLA-HD.", type = "error")
        w$hide()
        return()
      }
      
      options_values$hlahd_igv_i_options <- parseAndValidateGenomeSpec(
        genomeName = "hlahd", initialLocus = "All", stockGenome = FALSE, dataMode = "localFiles",
        fasta = fasta_path_h_i, fastaIndex = fastaindex_path_h_i
      )
      options_values$hlahd_igv_ii_options <- parseAndValidateGenomeSpec(
        genomeName = "hlahd", initialLocus = "All", stockGenome = FALSE, dataMode = "localFiles",
        fasta = fasta_path_h_ii, fastaIndex = fastaindex_path_h_ii
      )
      output$alignment_igv <- renderIgvShiny({
        igvShiny(options_values$hlahd_igv_i_options)
      })
      output$alignment_ii_igv <- renderIgvShiny({
        igvShiny(options_values$hlahd_igv_ii_options)
      })
      shinyalert("Success", "HLA-HD reference uploaded successfully.", type = "success")
    } else if (current_tool == "T1K") {
      fasta_path_t_i <- file.path(active_dir, "hla_classI_ABC.fasta")
      fastaindex_path_t_i <- file.path(active_dir, "hla_classI_ABC.fasta.fai")
      fasta_path_t_ii <- file.path(active_dir, "hla_classII_D.fasta")
      fastaindex_path_t_ii <- file.path(active_dir, "hla_classII_D.fasta.fai")
      
      if (!file.exists(fasta_path_t_i)) {
        shinyalert("Error", "Required FASTA file is missing for T1K.", type = "error")
        w$hide()
        return()
      }
      
      options_values$t1k_igv_i_options <- parseAndValidateGenomeSpec(
        genomeName = "t1k", initialLocus = "All", stockGenome = FALSE, dataMode = "localFiles",
        fasta = fasta_path_t_i, fastaIndex = fastaindex_path_t_i
      )
      options_values$t1k_igv_ii_options <- parseAndValidateGenomeSpec(
        genomeName = "t1k", initialLocus = "All", stockGenome = FALSE, dataMode = "localFiles",
        fasta = fasta_path_t_ii, fastaIndex = fastaindex_path_t_ii
      )
      output$alignment_igv <- renderIgvShiny({
        igvShiny(options_values$t1k_igv_i_options)
      })
      output$alignment_ii_igv <- renderIgvShiny({
        igvShiny(options_values$t1k_igv_ii_options)
      })
      shinyalert("Success", "T1K reference uploaded successfully.", type = "success")
    }
    
    w$hide()
  }
)

#-------------------------------------------------
### OptiType
load_and_show_region_optitype <- function(igv_id) {
  active_dir <- get_active_igv_dir("OptiType")
  bam_file <- file.path(active_dir, "optitype_merge_filtered.bam")
  if (file.exists(bam_file)) {
    x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
    loadBamTrackFromLocalData(session, id = igv_id, trackName = "OptiType Alignment", data = x)
    showGenomicRegion(session, id = igv_id, region = options_values$optitype_igv_options)
  }
}

### HLA-HD
load_and_show_region_hlahd_I <- function(igv_id) {
  active_dir <- get_active_igv_dir("HLA-HD")
  bam_file <- file.path(active_dir, "sample.modified.bam")
  if (file.exists(bam_file)) {
    x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
    loadBamTrackFromLocalData(session, id = igv_id, trackName = "HLA-HD Class I", data = x)
    showGenomicRegion(session, id = igv_id, region = options_values$hlahd_igv_i_options)
  }
}

load_and_show_region_hlahd_II <- function(igv_id) {
  active_dir <- get_active_igv_dir("HLA-HD")
  bam_file <- file.path(active_dir, "sample.modified.bam")
  if (file.exists(bam_file)) {
    x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
    loadBamTrackFromLocalData(session, id = igv_id, trackName = "HLA-HD Class II", data = x)
    showGenomicRegion(session, id = igv_id, region = options_values$hlahd_igv_ii_options)
  }
}

### T1K
load_and_show_region_t1k_I <- function(igv_id) {
  active_dir <- get_active_igv_dir("T1K")
  bam_file <- file.path(active_dir, "t1k_mapping.bam")
  if (file.exists(bam_file)) {
    x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
    loadBamTrackFromLocalData(session, id = igv_id, trackName = "T1K Class I", data = x)
    showGenomicRegion(session, id = igv_id, region = options_values$t1k_igv_i_options)
  }
}

load_and_show_region_t1k_II <- function(igv_id) {
  active_dir <- get_active_igv_dir("T1K")
  bam_file <- file.path(active_dir, "t1k_mapping.bam")
  if (file.exists(bam_file)) {
    x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
    loadBamTrackFromLocalData(session, id = igv_id, trackName = "T1K Class II", data = x)
    showGenomicRegion(session, id = igv_id, region = options_values$t1k_igv_ii_options)
  }
}

observeEvent(input$igv_bam_button, {
  w <- Waiter$new(html = spin_loader(), color = transparent(0.8))
  w$show()
  
  tool <- package_value()
  if (tool == "OptiType") {
    load_and_show_region_optitype("alignment_igv")
    shinyalert("Success", "BAM File uploaded successfully for OptiType.", type = "success")
  } else if (tool == "HLA-HD") {
    load_and_show_region_hlahd_I("alignment_igv")
    load_and_show_region_hlahd_II("alignment_ii_igv")
    shinyalert("Success", "BAM File uploaded successfully for HLA-HD.", type = "success")
  } else if (tool == "T1K") {
    load_and_show_region_t1k_I("alignment_igv")
    load_and_show_region_t1k_II("alignment_ii_igv")
    shinyalert("Success", "BAM File uploaded successfully for T1K.", type = "success")
  }
  
  w$hide()
})

#------------------------------------------
observeEvent(input$hg38_igv_i_bam_button, {
  w <- Waiter$new(html = spin_loader(), color = transparent(0.8))
  w$show()
  bam_file <- file.path(tmp_dir, "BAM", "sample.bam")
  if (file.exists(bam_file)) {
    load_and_show_region_hg38_I("hg38_igv_i")
    shinyalert("Success", "BAM File Uploaded Successfully.", type = "success")
  } else {
    shinyalert("Error", "Please Upload BAM File!", type = "error")
  }
  w$hide()
})

observeEvent(input$hg38_igv_ii_bam_button, {
  w <- Waiter$new(html = spin_loader(), color = transparent(0.8))
  w$show()
  bam_file <- file.path(tmp_dir, "BAM", "sample.bam")
  if (file.exists(bam_file)) {
    load_and_show_region_hg38_II("hg38_igv_ii")
    shinyalert("Success", "BAM File Uploaded Successfully.", type = "success")
  } else {
    shinyalert("Error", "Please Upload BAM File!", type = "error")
  }
  w$hide()
})
