## igv
# Define three directory paths
ori_dir <- getwd()
tmp_dir <- tempdir()

directory_o <- file.path(tmp_dir, "IGV", "optitype")
directory_h <- file.path(tmp_dir, "IGV", "hlahd")
directory_s <- file.path(tmp_dir, "IGV", "spechla")
#----------------------------------------------------
# Use reactiveValues to store globally shared variables
options_values <- reactiveValues(
  optitype_igv_options = NULL,
  hlahd_igv_i_options = NULL,
  hlahd_igv_ii_options = NULL,
  spechla_igv_options = NULL
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
  bam_file <- file.path(tmp_dir, "BAM", "0.bam")  
  bam_name <- "Bam file"
  
  if (file.exists(bam_file)) {
    x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
    loadBamTrackFromLocalData(session, id = igv_id, trackName = bam_name, data = x)
    showGenomicRegion(session, id = igv_id, region = hg38_igv_i_options)
  }
}
load_and_show_region_hg38_II <- function(igv_id) {
  bam_file <- file.path(tmp_dir, "BAM", "0.bam")  
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
  
  updateTabsetPanel(session, "igv_tabsetPanel", "Package")
  
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
    
  } else if (package_value() == "SpecHLA") {
    # SpecHLA
    fasta_path_s <- file.path(directory_s, "specHLA_hla_reference.fasta")
    fastaindex_path_s <- file.path(directory_s, "specHLA_hla_reference.fasta.fai")
    
    if (!file.exists(fasta_path_s) || !file.exists(fastaindex_path_s)) {
      shinyalert("Error", "Required FASTA or index file is missing for SpecHLA.", type = "error")
      w$hide()
      return()
    }
    
    options_values$spechla_igv_options <- parseAndValidateGenomeSpec(
      genomeName = "spechla",
      initialLocus = "All",
      stockGenome = FALSE,
      dataMode = "localFiles",
      fasta = fasta_path_s,
      fastaIndex = fastaindex_path_s
    )
    
    if (is.null(options_values$spechla_igv_options)) {
      shinyalert("Error", "Failed to parse and validate genome specification for SpecHLA.", type = "error")
      w$hide()
      return()
    }
    
    output$alignment_igv <- renderIgvShiny({
      igvShiny(options_values$spechla_igv_options)
    })
    
    shinyalert("Success", "SpecHLA reference uploaded successfully.", type = "success")
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
  bam_file <- file.path(directory_h, "sample1.modified.bam")
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
  bam_file <- file.path(directory_h, "sample1.modified.bam")
  bam_name <- "HLA-HD Class II Alignment File"
  
  if (!file.exists(bam_file)) {
    shinyalert("Error", paste("BAM File for HLA-HD Class II not found!"), type = "error")
    return()
  }
  
  x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
  loadBamTrackFromLocalData(session, id = igv_id, trackName = bam_name, data = x)
  showGenomicRegion(session, id = igv_id, region = options_values$hlahd_igv_ii_options)
}

### SpecHLA
load_and_show_region_spechla <- function(igv_id) {
  bam_file <- file.path(directory_s, "sample1.merge.bam")
  bam_name <- "SpecHLA Alignment File"
  
  if (!file.exists(bam_file)) {
    shinyalert("Error", paste("BAM File for SpecHLA not found!"), type = "error")
    return()
  }
  
  x <- readGAlignments(bam_file, param = Rsamtools::ScanBamParam(what = "seq"))
  loadBamTrackFromLocalData(session, id = igv_id, trackName = bam_name, data = x)
  showGenomicRegion(session, id = igv_id, region = options_values$spechla_igv_options)
}

observeEvent(input$igv_bam_button, {
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8)
  )
  w$show()
  
  if (package_value() == "OptiType") {
    success_optitype <- file.exists(file.path(directory_o, "optitype_merge_filtered.bam"))
    load_and_show_region_optitype("alignment_igv")
    
    if (success_optitype) {
      shinyalert("Success", paste("BAM File uploaded successfully for OptiType."), type = "success")
    }
  } else if (package_value() == "HLA-HD") {
    success_hla <- file.exists(file.path(directory_h, "sample1.modified.bam"))

    load_and_show_region_hlahd_I("alignment_igv")
    load_and_show_region_hlahd_II("alignment_ii_igv")
    
    if (success_hla) {
      shinyalert("Success", paste("BAM File uploaded successfully for HLA-HD."), type = "success")
    }
    
  } else if (package_value() == "SpecHLA") {
    success_spechla <- file.exists(file.path(directory_s, "sample1.merge.bam"))
    load_and_show_region_spechla("alignment_igv")
    
    if (success_spechla) {
      shinyalert("Success", paste("BAM File uploaded successfully for SpecHLA."), type = "success")
    }
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
  
  bam_file <- file.path(tmp_dir, "BAM", "0.bam")  
  
  if (file.exists(bam_file)) {
    load_and_show_region_hg38_I("hg38_igv_i")
    shinyalert("Success", paste("BAM File uploaded successfully."), type = "success")
  } else {
    shinyalert("Error", paste("Please upload the BAM file!"), type = "error")
  }
  
  w$hide() 
})

observeEvent(input$hg38_igv_ii_bam_button, {
  w <- Waiter$new(
    html = spin_loader(),
    color = transparent(0.8) 
  )
  w$show()
  
  bam_file <- file.path(tmp_dir, "BAM", "0.bam")  
  
  if (file.exists(bam_file)) {
    load_and_show_region_hg38_II("hg38_igv_ii")
    shinyalert("Success", paste("BAM File uploaded successfully."), type = "success")
  } else {
    shinyalert("Error", paste("Please upload the BAM file!"), type = "error")
  }
  
  w$hide() 
})

