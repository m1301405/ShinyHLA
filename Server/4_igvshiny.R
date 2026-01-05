## igv server logic

# 1. 定義目錄路徑
ori_dir <- getwd()
tmp_dir <- tempdir()

# 暫存區路徑 (用於新分析)
directory_o <- file.path(tmp_dir, "IGV", "optitype")
directory_h <- file.path(tmp_dir, "IGV", "hlahd")
directory_t <- file.path(tmp_dir, "IGV", "t1k")

#----------------------------------------------------
# [新增功能] 動態取得當前套件的檔案目錄
# 會根據 current_job_id() 狀態自動切換「暫存區」或「歷史歸檔區」
get_active_igv_dir <- function(tool) {
  # 如果目前是正在分析的狀態，使用 tmp_dir 下的暫存路徑
  if (current_job_id() == "No Job Active" || current_job_id() == "New Analysis") {
    return(switch(tool, 
                  "OptiType" = directory_o, 
                  "HLA-HD"   = directory_h, 
                  "T1K"      = directory_t))
  } else {
    # 如果載入了歷史 JobID，則指向該 ID 下歸檔的 IGV 資料夾
    return(file.path(tmp_dir, "History", current_job_id(), "IGV"))
  }
}

# [新增功能] 備份 IGV 檔案到歷史目錄的 Helper Function
# 此函數應在 4_srv_package.R 每個套件分析結尾呼叫，將檔案從 tmp 移至持久化目錄
backup_igv_files <- function(tool, job_id) {
  source_dir <- switch(tool, 
                       "OptiType" = directory_o, 
                       "HLA-HD"   = directory_h, 
                       "T1K"      = directory_t)
  
  target_dir <- file.path(tmp_dir, "History", job_id, "IGV")
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 複製該套件目錄下的所有 BAM, FASTA, FAI 檔案
  files_to_copy <- list.files(source_dir, full.names = TRUE)
  if (length(files_to_copy) > 0) {
    file.copy(files_to_copy, target_dir, overwrite = TRUE)
  }
}

#----------------------------------------------------
# 使用 reactiveValues 存储全局共享的变量
options_values <- reactiveValues(
  optitype_igv_options = NULL,
  hlahd_igv_i_options = NULL,
  hlahd_igv_ii_options = NULL,
  t1k_igv_i_options = NULL,
  t1k_igv_ii_options = NULL
)

#----------------------------------------------------
# Local BAM file upload (固定使用 tmp_dir)
hg38_igv_i_options <- parseAndValidateGenomeSpec(
  genomeName = "hg38",
  initialLocus = c("HLA-A", "HLA-B", "HLA-C"))
hg38_igv_ii_options <- parseAndValidateGenomeSpec(
  genomeName = "hg38",
  initialLocus = c("HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"))

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

output$hg38_igv_i <- renderIgvShiny({ igvShiny(hg38_igv_i_options) })
output$hg38_igv_ii <- renderIgvShiny({ igvShiny(hg38_igv_ii_options) })

#----------------------------------------------------
# 進行單獨套件的IGV展示
observeEvent(input$igv_reference, {
  w <- Waiter$new(html = spin_loader(), color = transparent(0.8))
  w$show()
  
  updateTabsetPanel(session, "igv_tabsetPanel", "Typing Package Results")
  
  # 獲取動態目錄
  current_tool <- package_value()
  active_dir <- get_active_igv_dir(current_tool)
  
  if (current_tool == "OptiType") {
    fasta_path_o <- file.path(active_dir, "optitype_hla_reference_igv.fasta")
    fastaindex_path_o <- file.path(active_dir, "optitype_hla_reference_igv.fasta.fai")
    
    if (!file.exists(fasta_path_o)) {
      shinyalert("Error", "Required FASTA file is missing for OptiType in current directory.", type = "error")
      w$hide(); return()
    }
    
    options_values$optitype_igv_options <- parseAndValidateGenomeSpec(
      genomeName = "optitype", initialLocus = "All", stockGenome = FALSE, dataMode = "localFiles",
      fasta = fasta_path_o, fastaIndex = fastaindex_path_o
    )
    output$alignment_igv <- renderIgvShiny({ igvShiny(options_values$optitype_igv_options) })
    shinyalert("Success", "OptiType reference uploaded successfully.", type = "success")
    
  } else if (current_tool == "HLA-HD") {
    # HLA-HD I & II
    fasta_path_h_i <- file.path(active_dir, "hlahd_hla_reference_I_igv.fasta")
    fastaindex_path_h_i <- file.path(active_dir, "hlahd_hla_reference_I_igv.fasta.fai")
    fasta_path_h_ii <- file.path(active_dir, "hlahd_hla_reference_II_igv.fasta")
    fastaindex_path_h_ii <- file.path(active_dir, "hlahd_hla_reference_II_igv.fasta.fai")
    
    if (!file.exists(fasta_path_h_i)) {
      shinyalert("Error", "Required FASTA file is missing for HLA-HD.", type = "error")
      w$hide(); return()
    }
    
    options_values$hlahd_igv_i_options <- parseAndValidateGenomeSpec(
      genomeName = "hlahd", initialLocus = "All", stockGenome = FALSE, dataMode = "localFiles",
      fasta = fasta_path_h_i, fastaIndex = fastaindex_path_h_i
    )
    options_values$hlahd_igv_ii_options <- parseAndValidateGenomeSpec(
      genomeName = "hlahd", initialLocus = "All", stockGenome = FALSE, dataMode = "localFiles",
      fasta = fasta_path_h_ii, fastaIndex = fastaindex_path_h_ii
    )
    output$alignment_igv <- renderIgvShiny({ igvShiny(options_values$hlahd_igv_i_options) })
    output$alignment_ii_igv <- renderIgvShiny({ igvShiny(options_values$hlahd_igv_ii_options) })
    shinyalert("Success", "HLA-HD reference uploaded successfully.", type = "success")
    
  } else if (current_tool == "T1K") {
    fasta_path_t_i <- file.path(active_dir, "hla_classI_ABC.fasta")
    fastaindex_path_t_i <- file.path(active_dir, "hla_classI_ABC.fasta.fai")
    fasta_path_t_ii <- file.path(active_dir, "hla_classII_D.fasta")
    fastaindex_path_t_ii <- file.path(active_dir, "hla_classII_D.fasta.fai")
    
    if (!file.exists(fasta_path_t_i)) {
      shinyalert("Error", "Required FASTA file is missing for T1K.", type = "error")
      w$hide(); return()
    }
    
    options_values$t1k_igv_i_options <- parseAndValidateGenomeSpec(
      genomeName = "t1k", initialLocus = "All", stockGenome = FALSE, dataMode = "localFiles",
      fasta = fasta_path_t_i, fastaIndex = fastaindex_path_t_i
    )
    options_values$t1k_igv_ii_options <- parseAndValidateGenomeSpec(
      genomeName = "t1k", initialLocus = "All", stockGenome = FALSE, dataMode = "localFiles",
      fasta = fasta_path_t_ii, fastaIndex = fastaindex_path_t_ii
    )
    output$alignment_igv <- renderIgvShiny({ igvShiny(options_values$t1k_igv_i_options) })
    output$alignment_ii_igv <- renderIgvShiny({ igvShiny(options_values$t1k_igv_ii_options) })
    shinyalert("Success", "T1K reference uploaded successfully.", type = "success")
  }
  w$hide()
})

#-------------------------------------------------
## 每一個tabpanel的BAM載入 (皆使用動態目錄)

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

# BAM 按鈕監聽邏輯
observeEvent(input$igv_bam_button, {
  w <- Waiter$new(html = spin_loader(), color = transparent(0.8)); w$show()
  
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
# hg38 BAM 顯示邏輯 (固定使用 tmp_dir)
observeEvent(input$hg38_igv_i_bam_button, {
  w <- Waiter$new(html = spin_loader(), color = transparent(0.8)); w$show()
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
  w <- Waiter$new(html = spin_loader(), color = transparent(0.8)); w$show()
  bam_file <- file.path(tmp_dir, "BAM", "sample.bam")  
  if (file.exists(bam_file)) {
    load_and_show_region_hg38_II("hg38_igv_ii")
    shinyalert("Success", "BAM File Uploaded Successfully.", type = "success")
  } else {
    shinyalert("Error", "Please Upload BAM File!", type = "error")
  }
  w$hide() 
})