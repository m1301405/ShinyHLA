## =====================================================
## app_globals.R
## Global options, shared datasets, and file paths
## =====================================================

## Upload file size limit
options(shiny.maxRequestSize = 5000 * 1024^2)

## System CPU
num_cores <- min(parallel::detectCores(), 8)

## Cascade filter
HLA_typing <- readxl::read_excel("/root/shiny/HLA_package_list.xlsx")

## PharmGKB root directory
PHARMGKB_DIR <- "/root/shiny/PharmGKB"

## =====================================================
## PharmGKB datasets
## =====================================================

summary_hla <- read.delim(
  file.path(PHARMGKB_DIR, "summaryAnnotationHLA.tsv"),
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

variant_hla <- read.delim(
  file.path(PHARMGKB_DIR, "variantAnnotationHLA.tsv"),
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

druglabel_hla <- read.delim(
  file.path(PHARMGKB_DIR, "druglabelHLA.tsv"),
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

drugs <- read.delim(
  file.path(PHARMGKB_DIR, "drugs.tsv"),
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

phenotypes <- read.delim(
  file.path(PHARMGKB_DIR, "phenotypes.tsv"),
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)