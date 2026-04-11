## app_helpers.R
# ============================================================
# Demo / Reload common settings
# ============================================================
DEMO_LIST_CSV <- "/root/shiny/Demo_list.csv"
DEMO_ROOT     <- "/root/shiny/Demo"
HISTORY_ROOT  <- "/root/shiny/History"

demo_mode <- reactiveVal(FALSE)

# ------------------------------------------------------------
# Helper: Read Demo list
# ------------------------------------------------------------
read_demo_list <- function() {
  if (!file.exists(DEMO_LIST_CSV)) {
    stop(sprintf("Demo list file not found: %s", DEMO_LIST_CSV))
  }
  
  df <- read.csv(DEMO_LIST_CSV, stringsAsFactors = FALSE, sep = ",")
  
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
# Helper: NGS_type -> Demo folder name
# WES -> wes, RNA-seq -> rna-seq
# ------------------------------------------------------------
ngs_to_demo_subdir <- function(ngs_type) {
  ngs_type <- trimws(as.character(ngs_type))
  
  if (identical(ngs_type, "WES")) {
    return("wes")
  }
  if (identical(ngs_type, "RNA-seq")) {
    return("rna-seq")
  }
  
  stop(sprintf("Unsupported NGS_type for Demo path routing: %s", ngs_type))
}

# ------------------------------------------------------------
# Helper: 根據條件查找唯一 Demo_file_folder
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
# Helper: Demo target_dir
# ------------------------------------------------------------
build_demo_target_dir <- function(ngs_type, demo_file_folder) {
  subdir <- ngs_to_demo_subdir(ngs_type)
  file.path(DEMO_ROOT, subdir, demo_file_folder)
}

## =====================================================
## Shared helper functions for ShinyHLA + PharmGKB
## =====================================================

safe_trim <- function(x) {
  x <- as.character(x)
  trimws(x)
}

clean_empty <- function(x) {
  x <- safe_trim(x)
  x[x == ""] <- NA
  x
}

collapse_one <- function(x, sep = "; ") {
  x <- unique(clean_empty(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  paste(x, collapse = sep)
}

split_clean_unique <- function(x, pattern = "[;,]") {
  x <- clean_empty(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(character(0))
  out <- unlist(strsplit(x, pattern))
  out <- trimws(out)
  out <- out[out != ""]
  unique(out)
}

collapse_pretty <- function(x, sep = "; ") {
  x <- split_clean_unique(x)
  if (length(x) == 0) return("N/A")
  x <- sort(x)
  paste(x, collapse = sep)
}

truncate_text <- function(x, width = 110) {
  x <- safe_trim(x)
  ifelse(nchar(x) > width, paste0(substr(x, 1, width), "..."), x)
}

clean_pvalue_text <- function(x) {
  x <- as.character(x)
  
  out <- sapply(x, function(v) {
    if (is.na(v) || trimws(v) == "") return("N/A")
    
    vals <- unlist(strsplit(v, ";"))
    vals <- trimws(vals)
    
    vals <- gsub("^=", "", vals)
    vals <- gsub("=$", "", vals)
    vals <- gsub("^= ", "", vals)
    vals <- gsub(" = ", "", vals)
    vals <- vals[vals != ""]
    
    if (length(vals) == 0) return("N/A")
    paste(vals, collapse = "; ")
  })
  
  return(out)
}

safe_allele_filename <- function(x) {
  x <- safe_trim(x)
  x <- gsub("\\*", "_", x)
  x <- gsub(":", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
  x
}

allele_to_gene <- function(x) {
  x <- safe_trim(x)
  x <- sub("\\*.*$", "", x)
  x
}

download_button_compact <- function(outputId, label = "Download") {
  div(
    style = "display:inline-block;",
    downloadButton(outputId, label, class = "btn btn-default btn-sm")
  )
}

safe_external_link <- function(href, label, class = NULL, style = NULL) {
  tags$a(
    href = href,
    target = "_blank",
    rel = "noopener noreferrer",
    onclick = "event.stopPropagation();",
    class = class,
    style = style,
    label
  )
}

## =====================================================
## URL / external resource helpers
## =====================================================

make_summary_search_url <- function(allele) {
  paste0(
    "https://www.clinpgx.org/search?query=",
    utils::URLencode(allele, reserved = TRUE)
  )
}

make_variant_record_url <- function(id) {
  paste0("https://www.clinpgx.org/variantAnnotation/", id)
}

make_druglabel_url <- function(id) {
  paste0("https://www.clinpgx.org/labelAnnotation/", id)
}

## =====================================================
## PMID / literature helpers
## =====================================================

make_pmid_links <- function(pmid_string) {
  if (is.na(pmid_string) || trimws(pmid_string) == "") {
    return("N/A")
  }
  
  pmids <- unlist(strsplit(as.character(pmid_string), "[;,]"))
  pmids <- trimws(pmids)
  pmids <- pmids[pmids != ""]
  pmids <- unique(pmids)
  
  if (length(pmids) == 0) {
    return("N/A")
  }
  
  link_tags <- lapply(seq_along(pmids), function(i) {
    pmid <- pmids[i]
    tagList(
      safe_external_link(
        href = paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/"),
        label = pmid
      ),
      if (i < length(pmids)) HTML("; ")
    )
  })
  
  do.call(tagList, link_tags)
}

make_literature_links <- function(x) {
  x <- safe_trim(x)
  if (length(x) == 0 || all(is.na(x)) || all(x == "")) return("N/A")
  
  vals <- unlist(strsplit(paste(x, collapse = "; "), "[;]"))
  vals <- trimws(vals)
  vals <- vals[vals != ""]
  vals <- unique(vals)
  
  if (length(vals) == 0) return("N/A")
  
  out <- lapply(seq_along(vals), function(i) {
    one <- vals[i]
    
    if (grepl("^PMCID:", one, ignore.case = TRUE)) {
      pmcid <- sub("^PMCID:", "", one, ignore.case = TRUE)
      href <- paste0("https://pmc.ncbi.nlm.nih.gov/articles/", pmcid, "/")
      tagList(
        safe_external_link(href = href, label = one),
        if (i < length(vals)) HTML("; ")
      )
    } else {
      pmid <- sub("^PMID:", "", one, ignore.case = TRUE)
      pmid <- trimws(pmid)
      href <- paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/")
      label <- if (grepl("^PMID:", one, ignore.case = TRUE)) one else paste0("PMID:", pmid)
      tagList(
        safe_external_link(href = href, label = label),
        if (i < length(vals)) HTML("; ")
      )
    }
  })
  
  do.call(tagList, out)
}

clean_quoted_text_block <- function(x) {
  x <- safe_trim(x)
  
  if (length(x) == 0) return(x)
  
  out <- vapply(x, function(one) {
    if (is.na(one) || one == "") return("N/A")
    
    one <- gsub("\\s+", " ", one)
    one <- trimws(one)
    
    one <- gsub('"[[:space:]]+"', '"; "', one)
    one <- gsub('^"(.*)"$', "\\1", one)
    one <- gsub('";[[:space:]]*"', "; ", one)
    one <- gsub('^"+|"+$', "", one)
    one <- trimws(one)
    
    one
  }, FUN.VALUE = character(1))
  
  out
}

## =====================================================
## Drug helpers
## =====================================================

normalize_drug_name <- function(x) {
  x <- safe_trim(x)
  x <- gsub('^"(.*)"$', "\\1", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x <- tolower(x)
  x
}

split_outside_quotes <- function(x, separators = c(",", ";")) {
  x <- safe_trim(x)
  x <- x[!is.na(x) & x != ""]
  
  if (length(x) == 0) return(character(0))
  
  out_all <- character(0)
  
  for (one_string in x) {
    chars <- strsplit(one_string, "", fixed = TRUE)[[1]]
    buf <- character(0)
    pieces <- character(0)
    in_quotes <- FALSE
    i <- 1
    
    while (i <= length(chars)) {
      ch <- chars[i]
      
      if (ch == '"') {
        if (in_quotes && i < length(chars) && chars[i + 1] == '"') {
          buf <- c(buf, '"')
          i <- i + 2
          next
        } else {
          in_quotes <- !in_quotes
          i <- i + 1
          next
        }
      }
      
      if (!in_quotes && ch %in% separators) {
        piece <- paste(buf, collapse = "")
        piece <- trimws(piece)
        if (nzchar(piece)) {
          pieces <- c(pieces, piece)
        }
        buf <- character(0)
        i <- i + 1
        next
      }
      
      buf <- c(buf, ch)
      i <- i + 1
    }
    
    piece <- paste(buf, collapse = "")
    piece <- trimws(piece)
    if (nzchar(piece)) {
      pieces <- c(pieces, piece)
    }
    
    out_all <- c(out_all, pieces)
  }
  
  out_all <- trimws(out_all)
  out_all <- gsub('^"(.*)"$', "\\1", out_all)
  out_all <- gsub("\\s+", " ", out_all)
  out_all <- out_all[out_all != ""]
  out_all
}

split_drug_names <- function(x) {
  split_outside_quotes(x, separators = c(",", ";"))
}

split_drug_names_unique <- function(x) {
  unique(split_drug_names(x))
}

expand_drug_vector <- function(x) {
  vals <- split_drug_names(x)
  vals <- safe_trim(vals)
  vals <- vals[!is.na(vals) & vals != ""]
  vals
}

find_drug_id <- function(drug_name) {
  drug_name_norm <- normalize_drug_name(drug_name)
  
  hit <- drugs[
    normalize_drug_name(drugs$Name) == drug_name_norm,
    ,
    drop = FALSE
  ]
  
  if (nrow(hit) == 0) {
    return(NA_character_)
  }
  
  as.character(hit$`PharmGKB Accession Id`[1])
}

make_drug_links <- function(drug_string) {
  if (is.na(drug_string) || trimws(drug_string) == "") {
    return("N/A")
  }
  
  vals <- split_drug_names_unique(drug_string)
  
  if (length(vals) == 0) {
    return("N/A")
  }
  
  out <- lapply(seq_along(vals), function(i) {
    one <- vals[i]
    drug_id <- find_drug_id(one)
    
    link_part <- if (!is.na(drug_id) && drug_id != "") {
      safe_external_link(
        href = paste0("https://www.clinpgx.org/drug/", drug_id),
        label = one
      )
    } else {
      one
    }
    
    tagList(
      link_part,
      if (i < length(vals)) HTML("; ")
    )
  })
  
  do.call(tagList, out)
}

make_top_drugs_df <- function(x, top_n = 5) {
  drug_vec <- expand_drug_vector(x)
  
  if (length(drug_vec) == 0) {
    return(data.frame(drug = character(0), n = integer(0), stringsAsFactors = FALSE))
  }
  
  data.frame(drug = drug_vec, stringsAsFactors = FALSE) %>%
    dplyr::filter(!is.na(drug), drug != "") %>%
    dplyr::count(drug, sort = TRUE) %>%
    dplyr::slice_head(n = top_n)
}

## =====================================================
## Variant helpers
## =====================================================

make_variant_links <- function(variant_string, evidence_id_string) {
  if (is.na(variant_string) || trimws(variant_string) == "") {
    return("N/A")
  }
  
  variants <- unlist(strsplit(as.character(variant_string), ","))
  variants <- trimws(variants)
  variants <- variants[variants != ""]
  
  evidence_ids <- unlist(strsplit(as.character(evidence_id_string), ";"))
  evidence_ids <- trimws(evidence_ids)
  evidence_ids <- evidence_ids[evidence_ids != ""]
  
  out <- lapply(seq_along(variants), function(i) {
    v <- variants[i]
    
    link_part <- if (i <= length(evidence_ids) && evidence_ids[i] != "") {
      safe_external_link(
        href = paste0("https://www.clinpgx.org/variantAnnotation/", evidence_ids[i]),
        label = v
      )
    } else {
      v
    }
    
    tagList(
      link_part,
      if (i < length(variants)) HTML(", ")
    )
  })
  
  do.call(tagList, out)
}

clean_variant_phenotype <- function(x) {
  x <- safe_trim(x)
  x <- gsub("^Side Effect:", "", x)
  x <- gsub("^Disease:", "", x)
  x <- gsub("^Other:", "", x)
  x <- gsub("^Efficacy:", "", x)
  x <- gsub("^PK:", "", x)
  trimws(x)
}

split_clean_variant_phenotypes <- function(x) {
  x <- safe_trim(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(character(0))
  
  x <- gsub(',[[:space:]]*"?[[:space:]]*Disease:', '; Disease:', x)
  x <- gsub(',[[:space:]]*"?[[:space:]]*Side Effect:', '; Side Effect:', x)
  x <- gsub(',[[:space:]]*"?[[:space:]]*Other:', '; Other:', x)
  x <- gsub(',[[:space:]]*"?[[:space:]]*Efficacy:', '; Efficacy:', x)
  x <- gsub(',[[:space:]]*"?[[:space:]]*PK:', '; PK:', x)
  
  vals <- unlist(strsplit(x, ";", fixed = FALSE))
  vals <- trimws(vals)
  vals <- vals[vals != ""]
  
  vals <- gsub('^"+|"+$', "", vals)
  vals <- gsub("^'+|'+$", "", vals)
  vals <- trimws(vals)
  
  vals <- gsub("^Side Effect:", "", vals)
  vals <- gsub("^Disease:", "", vals)
  vals <- gsub("^Other:", "", vals)
  vals <- gsub("^Efficacy:", "", vals)
  vals <- gsub("^PK:", "", vals)
  
  vals <- gsub('^"+|"+$', "", vals)
  vals <- gsub("^'+|'+$", "", vals)
  vals <- trimws(vals)
  
  vals <- vals[vals != ""]
  vals <- unique(vals)
  
  vals
}

## =====================================================
## Phenotype helpers
## =====================================================

normalize_phenotype_name <- function(x) {
  x <- safe_trim(x)
  x <- tolower(x)
  x
}

find_phenotype_id <- function(phenotype_name) {
  phenotype_name_norm <- normalize_phenotype_name(phenotype_name)
  
  hit <- phenotypes[
    normalize_phenotype_name(phenotypes$Name) == phenotype_name_norm,
    ,
    drop = FALSE
  ]
  
  if (nrow(hit) == 0) {
    return(NA_character_)
  }
  
  as.character(hit$`PharmGKB Accession Id`[1])
}

make_phenotype_links <- function(phenotype_string, use_breaks = TRUE) {
  if (is.na(phenotype_string) || trimws(phenotype_string) == "") {
    return("N/A")
  }
  
  vals <- unlist(strsplit(as.character(phenotype_string), ";"))
  vals <- trimws(vals)
  vals <- vals[vals != ""]
  vals <- unique(vals)
  
  if (length(vals) == 0) {
    return("N/A")
  }
  
  out <- lapply(seq_along(vals), function(i) {
    one <- vals[i]
    pid <- find_phenotype_id(one)
    
    link_part <- if (!is.na(pid) && pid != "") {
      safe_external_link(
        href = paste0("https://www.clinpgx.org/disease/", pid),
        label = one
      )
    } else {
      one
    }
    
    sep_part <- NULL
    if (i < length(vals)) {
      sep_part <- if (use_breaks) tags$br() else HTML("; ")
    }
    
    tagList(link_part, sep_part)
  })
  
  do.call(tagList, out)
}

## =====================================================
## Gene / drug label helpers
## =====================================================

make_gene_links <- function(x) {
  vals <- split_clean_unique(x)
  if (length(vals) == 0) return("N/A")
  
  out <- lapply(seq_along(vals), function(i) {
    g <- vals[i]
    tagList(
      safe_external_link(
        href = paste0("https://www.clinpgx.org/gene/", utils::URLencode(g, reserved = TRUE)),
        label = g
      ),
      if (i < length(vals)) HTML("; ")
    )
  })
  
  do.call(tagList, out)
}

make_druglabel_title_cell <- function(title, url,
                                      has_prescribing = "",
                                      has_dosing = "",
                                      has_alternate = "",
                                      has_other = "",
                                      biomarker_flag = "",
                                      match_scope = "Allele-level") {
  
  badge <- function(label, bg = "#f2c94c", fg = "#222") {
    tags$span(
      style = paste0(
        "display:inline-block; margin:4px 5px 0 0; padding:2px 8px;",
        "border-radius:7px; font-size:11px; font-weight:600;",
        "line-height:1.4; background:", bg, "; color:", fg, ";"
      ),
      label
    )
  }
  
  tagList(
    safe_external_link(href = url, label = title),
    tags$div(
      style = "margin-top:6px;",
      if (!is.na(has_alternate)   && safe_trim(has_alternate)   != "") badge("Alternate Drug"),
      if (!is.na(has_prescribing) && safe_trim(has_prescribing) != "") badge("Prescribing Info"),
      if (!is.na(has_dosing)      && safe_trim(has_dosing)      != "") badge("Dosing Info"),
      if (!is.na(has_other)       && safe_trim(has_other)       != "") badge("Other Guidance"),
      if (!is.na(biomarker_flag)  && safe_trim(biomarker_flag)  != "") badge("FDA Biomarker", bg = "#b7d7ea"),
      if (!is.na(match_scope)     && safe_trim(match_scope)     != "") {
        if (match_scope == "Gene-level") {
          badge("Gene-level", bg = "#e9ecef", fg = "#444")
        } else {
          badge("Allele-level", bg = "#fdecea", fg = "#b71c1c")
        }
      }
    )
  )
}

make_pgx_level_badge <- function(x) {
  x <- safe_trim(x)
  if (is.na(x) || x == "") return("N/A")
  
  bg <- dplyr::case_when(
    x == "Testing Required"    ~ "#c0392b",
    x == "Testing Recommended" ~ "#e67e22",
    x == "Actionable PGx"      ~ "#f2c94c",
    x == "Informative PGx"     ~ "#d9d9d9",
    TRUE                       ~ "#d9d9d9"
  )
  
  fg <- ifelse(x %in% c("Testing Required", "Testing Recommended"), "white", "#222")
  
  tags$span(
    style = paste0(
      "display:inline-block;",
      "padding:2px 7px;",
      "margin:2px 4px 0 0;",
      "border-radius:6px;",
      "font-size:10.5px;",
      "font-weight:600;",
      "line-height:1.3;",
      "background:", bg, ";",
      "color:", fg, ";"
    ),
    x
  )
}

normalize_typing_allele_to_pharmgkb <- function(x) {
  x <- safe_trim(x)
  
  if (is.na(x) || x == "") return(NA_character_)
  
  x <- sub("^HLA-", "", x, ignore.case = TRUE)
  
  if (!grepl("\\*", x)) return(NA_character_)
  
  gene_part <- sub("\\*.*$", "", x)
  field_part <- sub("^.*\\*", "", x)
  
  fields <- unlist(strsplit(field_part, ":", fixed = TRUE))
  fields <- trimws(fields)
  fields <- fields[fields != ""]
  
  if (length(fields) < 2) return(NA_character_)
  
  two_field <- paste(fields[1:2], collapse = ":")
  
  paste0("HLA-", gene_part, "*", two_field)
}

make_pharmgkb_row_button <- function(allele_query) {
  if (is.na(allele_query) || allele_query == "") return("")
  
  allele_escaped <- gsub("'", "\\\\'", allele_query)
  
  onclick_js <- sprintf(
    paste0(
      "this.disabled = true;",
      "this.innerHTML = '<i class=\"fas fa-spinner fa-spin\"></i> Searching...';",
      "Shiny.setInputValue(\"pharmgkb_from_table\", \"%s\", {priority: \"event\"});"
    ),
    allele_escaped
  )
  
  as.character(
    tags$button(
      type = "button",
      class = "btn btn-xs btn-danger pharmgkb-search-btn",
      title = "Click to load PharmGKB annotation",
      onclick = onclick_js,
      "Search"
    )
  )
}


prepare_hla_table_for_display <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  drop_cols <- intersect(c("PharmGKB", "pharmgkb_allele"), colnames(df))
  if (length(drop_cols) > 0) {
    df <- df[, setdiff(colnames(df), drop_cols), drop = FALSE]
  }
  
  core_cols <- c("Allele", "nucleotide", "protein")
  df <- df[, intersect(core_cols, colnames(df)), drop = FALSE]
  
  df$pharmgkb_allele <- vapply(
    df$Allele,
    normalize_typing_allele_to_pharmgkb,
    FUN.VALUE = character(1)
  )
  
  df$PharmGKB <- vapply(
    df$pharmgkb_allele,
    make_pharmgkb_row_button,
    FUN.VALUE = character(1)
  )
  
  df <- df[, c("PharmGKB", "Allele", "nucleotide", "protein", "pharmgkb_allele"), drop = FALSE]
  
  df
}


render_hla_typing_datatable <- function(df) {
  df_display <- prepare_hla_table_for_display(df)
  
  DT::datatable(
    df_display,
    filter = "none",
    rownames = FALSE,
    selection = "none",
    escape = FALSE,
    class = "nowrap",
    options = list(
      pageLength = 8,
      autoWidth = TRUE,
      scrollX = TRUE,
      columnDefs = list(
        list(
          targets = 0,
          orderable = FALSE,
          searchable = FALSE,
          width = "110px"
        ),
        list(
          targets = 4,
          visible = FALSE,
          searchable = FALSE
        ),
        list(
          width = "100px",
          targets = c(1, 2, 3),
          render = JS(
            "function(data, type, row, meta) {",
            "  if (type === 'display' && data != null && data.length > 25) {",
            "    return '<span title=\"' + data + '\">' + data.substr(0, 25) + '...</span>';",
            "  }",
            "  return data;",
            "}"
          )
        )
      )
    )
  )
}
