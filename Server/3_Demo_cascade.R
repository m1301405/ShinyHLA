# ============================================================
# Common settings
# ============================================================
DEMO_LIST_CSV <- "/root/shiny/Demo_list.csv"
DEMO_ROOT     <- "/root/shiny/Demo"
HISTORY_ROOT  <- "/root/shiny/History"

# source_mode:
#   "fastq" = normal FASTQ upload workflow
#   "bam"   = normal BAM/CRAM upload workflow
#   "demo"  = precomputed demo reload workflow
source_mode <- reactiveVal("fastq")
# ------------------------------------------------------------
# Upload button -> mode switch
# ------------------------------------------------------------
demo_uploaded <- reactiveVal(FALSE)
demo_selected_ngs <- reactiveVal(NULL)

observeEvent(input$fastq_upload_button, {
  source_mode("fastq")
  demo_uploaded(FALSE)
  demo_selected_ngs(NULL)
})

observeEvent(input$bam_upload_button, {
  source_mode("bam")
  demo_uploaded(FALSE)
  demo_selected_ngs(NULL)
})

observeEvent(input$demo_upload_button, {
  source_mode("demo")
})
# ------------------------------------------------------------
# Helper: read Demo list
# ------------------------------------------------------------
read_demo_list <- function() {
  if (!file.exists(DEMO_LIST_CSV)) {
    stop(sprintf("Demo list file not found: %s", DEMO_LIST_CSV))
  }
  
  df <- read.csv(
    DEMO_LIST_CSV,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    strip.white = TRUE
  )
  
  colnames(df) <- trimws(colnames(df))
  
  required_cols <- c("NGS_type", "IMGTHLA_version", "Typing_package", "Demo_file_folder")
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Demo list is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  df
}

# ------------------------------------------------------------
# Helper: convert NGS_type -> demo subdir
# ------------------------------------------------------------
ngs_to_demo_subdir <- function(ngs_type) {
  ngs_type <- trimws(as.character(ngs_type))
  
  if (identical(ngs_type, "WES")) return("wes")
  if (identical(ngs_type, "RNA-seq")) return("rna-seq")
  
  stop(sprintf("Unsupported NGS_type for Demo path routing: %s", ngs_type))
}

# ------------------------------------------------------------
# Helper: resolve unique Demo_file_folder
# ------------------------------------------------------------
resolve_demo_folder <- function(ngs_type, version, tool) {
  demo_df <- read_demo_list()
  
  matched <- subset(
    demo_df,
    trimws(NGS_type) == trimws(ngs_type) &
      trimws(IMGTHLA_version) == trimws(version) &
      trimws(Typing_package) == trimws(tool)
  )
  
  if (nrow(matched) == 0) {
    stop(sprintf(
      "No Demo result found for NGS_type=%s, IMGTHLA_version=%s, Typing_package=%s",
      ngs_type, version, tool
    ))
  }
  
  if (nrow(matched) > 1) {
    stop(sprintf(
      "Multiple Demo results found for NGS_type=%s, IMGTHLA_version=%s, Typing_package=%s",
      ngs_type, version, tool
    ))
  }
  
  matched$Demo_file_folder[[1]]
}

# ------------------------------------------------------------
# Helper: build demo target_dir
# ------------------------------------------------------------
build_demo_target_dir <- function(ngs_type, demo_file_folder) {
  subdir <- ngs_to_demo_subdir(ngs_type)
  file.path(DEMO_ROOT, subdir, demo_file_folder)
}