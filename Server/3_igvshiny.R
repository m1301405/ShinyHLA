## igv
# Define three directory paths
ori_dir <- getwd()
tmp_dir <- tempdir()

directory_o <- file.path(tmp_dir, "IGV", "optitype")
directory_h <- file.path(tmp_dir, "IGV", "hlahd")
directory_t <- file.path(tmp_dir, "IGV", "t1k")
#----------------------------------------------------
# Use reactiveValues to store globally shared variables
options_values <- reactiveValues(
  optitype_igv_options = NULL,
  hlahd_igv_i_options = NULL,
  hlahd_igv_ii_options = NULL,
  t1k_igv_i_options = NULL,
  t1k_igv_ii_options = NULL
)

#----------------------------------------------------
# Local BAM file upload
hg38_igv_i_options <- parseAndValidateGenomeSpec(
  genomeName = "hg38",
  initialLocus = c("HLA-A", "HLA-B", "HLA-C"))
hg38_igv_ii_options <- parseAndValidateGenomeSpec(
  genomeName = "hg38",
  initialLocus = c("HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"))

load_and_show_region_hg38_I <- function(igv_id) {
  bam_file <- file.path(tmp_dir, "BAM", "sample.bam")  
  bam_name <- "Bam file"
  
  if (file.exists(bam_file)) {
    x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
    loadBamTrackFromLocalData(session, id = igv_id, trackName = bam_name, data = x)
    showGenomicRegion(session, id = igv_id, region = hg38_igv_i_options)
  }
}
load_and_show_region_hg38_II <- function(igv_id) {
  bam_file <- file.path(tmp_dir, "BAM", "sample.bam")  
  bam_name <- "bam file"
  
  if (file.exists(bam_file)) {
    x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
    loadBamTrackFromLocalData(session, id = igv_id, trackName = bam_name, data = x)
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
# Perform IGV visualization for individual packages
observeEvent(input$igv_reference, {
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) # Transparent background
  )
  w$show() # Show the waiter
  
  updateTabsetPanel(session, "igv_tabsetPanel", "Typing Package Results")
  
  if (package_value() == "OptiType") {
    # OptiType
    fasta_path_o <- file.path(directory_o, "optitype_hla_reference_igv.fasta")
    fastaindex_path_o <- file.path(directory_o, "optitype_hla_reference_igv.fasta.fai")
    
    if (!file.exists(fasta_path_o) || !file.exists(fastaindex_path_o)) {
      shinyalert("Error", "Required FASTA or index file is missing for OptiType.", type = "error")
      w$hide()
      return()
    }
    
    options_values$optitype_igv_options <- parseAndValidateGenomeSpec(
      genomeName = "optitype",
      initialLocus = "All",
      stockGenome = FALSE,
      dataMode = "localFiles",
      fasta = fasta_path_o,
      fastaIndex = fastaindex_path_o
    )
    
    if (is.null(options_values$optitype_igv_options)) {
      shinyalert("Error", "Failed to parse and validate genome specification for OptiType.", type = "error")
      w$hide()
      return()
    }
    
    output$alignment_igv <- renderIgvShiny({
      igvShiny(options_values$optitype_igv_options)
    })
    
    shinyalert("Success", "OptiType reference uploaded successfully.", type = "success")
    
  } else if (package_value() == "HLA-HD") {
    # HLA-HD Class I
    fasta_path_h_i <- file.path(directory_h, "hlahd_hla_reference_I_igv.fasta")
    fastaindex_path_h_i <- file.path(directory_h, "hlahd_hla_reference_I_igv.fasta.fai")
    
    if (!file.exists(fasta_path_h_i) || !file.exists(fastaindex_path_h_i)) {
      shinyalert("Error", "Required FASTA or index file is missing for HLA-HD Class I and II.", type = "error")
      w$hide()
      return()
    }
    
    options_values$hlahd_igv_i_options <- parseAndValidateGenomeSpec(
      genomeName = "hlahd",
      initialLocus = "All",
      stockGenome = FALSE,
      dataMode = "localFiles",
      fasta = fasta_path_h_i,
      fastaIndex = fastaindex_path_h_i
    )
    
    if (is.null(options_values$hlahd_igv_i_options)) {
      shinyalert("Error", "Failed to parse and validate genome specification for HLA-HD Class I.", type = "error")
      w$hide()
      return()
    }
    
    output$alignment_igv <- renderIgvShiny({
      igvShiny(options_values$hlahd_igv_i_options)
    })
    
    # HLA-HD Class II
    fasta_path_h_ii <- file.path(directory_h, "hlahd_hla_reference_II_igv.fasta")
    fastaindex_path_h_ii <- file.path(directory_h, "hlahd_hla_reference_II_igv.fasta.fai")
    
    if (!file.exists(fasta_path_h_ii) || !file.exists(fastaindex_path_h_ii)) {
      shinyalert("Error", "Required FASTA or index file is missing for HLA-HD Class II.", type = "error")
      w$hide()
      return()
    }
    
    options_values$hlahd_igv_ii_options <- parseAndValidateGenomeSpec(
      genomeName = "hlahd",
      initialLocus = "All",
      stockGenome = FALSE,
      dataMode = "localFiles",
      fasta = fasta_path_h_ii,
      fastaIndex = fastaindex_path_h_ii
    )
    
    if (is.null(options_values$hlahd_igv_ii_options)) {
      shinyalert("Error", "Failed to parse and validate genome specification for HLA-HD Class II.", type = "error")
      w$hide()
      return()
    }
    
    output$alignment_ii_igv <- renderIgvShiny({
      igvShiny(options_values$hlahd_igv_ii_options)
    })
    
    shinyalert("Success", "HLA-HD reference uploaded successfully.", type = "success")
    
  } else if (package_value() == "T1K") {
    # T1K Class I
    fasta_path_t_i <- file.path(directory_t, "hla_classI_ABC.fasta")
    fastaindex_path_t_i <- file.path(directory_t, "hla_classI_ABC.fasta.fai")
    
    if (!file.exists(fasta_path_t_i) || !file.exists(fastaindex_path_t_i)) {
      shinyalert("Error", "Required FASTA or index file is missing for T1K Class I and II.", type = "error")
      w$hide()
      return()
    }
    
    options_values$t1k_igv_i_options <- parseAndValidateGenomeSpec(
      genomeName = "t1k",
      initialLocus = "All",
      stockGenome = FALSE,
      dataMode = "localFiles",
      fasta = fasta_path_t_i,
      fastaIndex = fastaindex_path_t_i
    )
    
    if (is.null(options_values$t1k_igv_i_options)) {
      shinyalert("Error", "Failed to parse and validate genome specification for T1K Class I.", type = "error")
      w$hide()
      return()
    }
    
    output$alignment_igv <- renderIgvShiny({
      igvShiny(options_values$t1k_igv_i_options)
    })
    
    # T1K Class II
    fasta_path_t_ii <- file.path(directory_t, "hla_classII_D.fasta")
    fastaindex_path_t_ii <- file.path(directory_t, "hla_classII_D.fasta.fai")
    
    if (!file.exists(fasta_path_t_ii) || !file.exists(fastaindex_path_t_ii)) {
      shinyalert("Error", "Required FASTA or index file is missing for T1K Class II.", type = "error")
      w$hide()
      return()
    }
    
    options_values$t1k_igv_ii_options <- parseAndValidateGenomeSpec(
      genomeName = "t1k",
      initialLocus = "All",
      stockGenome = FALSE,
      dataMode = "localFiles",
      fasta = fasta_path_t_ii,
      fastaIndex = fastaindex_path_t_ii
    )
    
    if (is.null(options_values$t1k_igv_ii_options)) {
      shinyalert("Error", "Failed to parse and validate genome specification for T1K Class II.", type = "error")
      w$hide()
      return()
    }
    
    output$alignment_ii_igv <- renderIgvShiny({
      igvShiny(options_values$t1k_igv_ii_options)
    })
    
    shinyalert("Success", "T1K reference uploaded successfully.", type = "success")
  }
  
  w$hide() # Hide the waiter
})
#-------------------------------------------------
## BAM file upload for each tabPanel
### OptiType
load_and_show_region_optitype <- function(igv_id) {
  bam_file <- file.path(directory_o, "optitype_merge_filtered.bam")
  bam_name <- "OptiType Alignment File"
  
  if (!file.exists(bam_file)) {
    shinyalert("Error", paste("BAM File for OptiType not found!"), type = "error")
    return()
  }
  
  x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
  loadBamTrackFromLocalData(session, id = igv_id, trackName = bam_name, data = x)
  showGenomicRegion(session, id = igv_id, region = options_values$optitype_igv_options)
}

### HLA-HD Class I
load_and_show_region_hlahd_I <- function(igv_id) {
  bam_file <- file.path(directory_h, "sample.modified.bam")
  bam_name <- "HLA-HD Class I Alignment File"
  
  if (!file.exists(bam_file)) {
    shinyalert("Error", paste("BAM File for HLA-HD Class I not found!"), type = "error")
    return()
  }
  
  x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
  loadBamTrackFromLocalData(session, id = igv_id, trackName = bam_name, data = x)
  showGenomicRegion(session, id = igv_id, region = options_values$hlahd_igv_i_options)
}

### HLA-HD Class II
load_and_show_region_hlahd_II <- function(igv_id) {
  bam_file <- file.path(directory_h, "sample.modified.bam")
  bam_name <- "HLA-HD Class II Alignment File"
  
  if (!file.exists(bam_file)) {
    shinyalert("Error", paste("BAM File for HLA-HD Class II not found!"), type = "error")
    return()
  }
  
  x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
  loadBamTrackFromLocalData(session, id = igv_id, trackName = bam_name, data = x)
  showGenomicRegion(session, id = igv_id, region = options_values$hlahd_igv_ii_options)
}

### T1K Class I
load_and_show_region_t1k_I <- function(igv_id) {
  bam_file <- file.path(directory_t, "t1k_mapping.bam")
  bam_name <- "T1K Class I Alignment File"
  
  if (!file.exists(bam_file)) {
    shinyalert("Error", "BAM File for T1K Class I not found!", type = "error")
    return()
  }
  
  x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
  loadBamTrackFromLocalData(session, id = igv_id, trackName = bam_name, data = x)
  showGenomicRegion(session, id = igv_id, region = options_values$t1k_igv_i_options)
}

### T1K Class II
load_and_show_region_t1k_II <- function(igv_id) {
  bam_file <- file.path(directory_t, "t1k_mapping.bam")
  bam_name <- "T1K Class II Alignment File"
  
  if (!file.exists(bam_file)) {
    shinyalert("Error", "BAM File for T1K Class II not found!", type = "error")
    return()
  }
  
  x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
  loadBamTrackFromLocalData(session, id = igv_id, trackName = bam_name, data = x)
  showGenomicRegion(session, id = igv_id, region = options_values$t1k_igv_ii_options)
}


observeEvent(input$igv_bam_button, {
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8)
  )
  w$show()
  
  if (package_value() == "OptiType") {
    bam_file <- file.path(directory_o, "optitype_merge_filtered.bam")
    if (!file.exists(bam_file)) {
      shinyalert("Error", "BAM file for OptiType not found!", type = "error")
      w$hide()
      return()
    }
    
    load_and_show_region_optitype("alignment_igv")
    shinyalert("Success", "BAM File uploaded successfully for OptiType.", type = "success")
    
  } else if (package_value() == "HLA-HD") {
    bam_file <- file.path(directory_h, "sample.modified.bam")
    if (!file.exists(bam_file)) {
      shinyalert("Error", "BAM file for HLA-HD not found!", type = "error")
      w$hide()
      return()
    }
    
    load_and_show_region_hlahd_I("alignment_igv")
    load_and_show_region_hlahd_II("alignment_ii_igv")
    shinyalert("Success", "BAM File uploaded successfully for HLA-HD.", type = "success")
    
  } else if (package_value() == "T1K") {
    bam_file <- file.path(directory_t, "t1k_mapping.bam")
    if (!file.exists(bam_file)) {
      shinyalert("Error", "BAM file for T1K not found!", type = "error")
      w$hide()
      return()
    }
    
    load_and_show_region_t1k_I("alignment_igv")
    load_and_show_region_t1k_II("alignment_ii_igv")
    shinyalert("Success", "BAM File uploaded successfully for T1K.", type = "success")
  }
  
  w$hide()
})
#------------------------------------------
observeEvent(input$hg38_igv_i_bam_button, {
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) 
  )
  w$show()
  
  bam_file <- file.path(tmp_dir, "BAM", "sample.bam")  
  
  if (file.exists(bam_file)) {
    load_and_show_region_hg38_I("hg38_igv_i")
    shinyalert("Success", paste("BAM File Uploaded Successfully."), type = "success")
  } else {
    shinyalert("Error", paste("Please Upload BAM File!"), type = "error")
  }
  
  w$hide() 
})

observeEvent(input$hg38_igv_ii_bam_button, {
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) 
  )
  w$show()
  
  bam_file <- file.path(tmp_dir, "BAM", "sample.bam")  
  
  if (file.exists(bam_file)) {
    load_and_show_region_hg38_II("hg38_igv_ii")
    shinyalert("Success", paste("BAM File Uploaded Successfully."), type = "success")
  } else {
    shinyalert("Error", paste("Please Upload BAM File!"), type = "error")
  }
  
  w$hide() 
})

