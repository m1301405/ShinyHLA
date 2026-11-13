pharmgkb_panel_visible <- reactiveVal(FALSE)
selected_allele <- reactiveVal(NULL)

output$summary_annotation_overview <- renderUI({
  allele <- selected_allele()
  req(allele)
  
  sdat <- summary_hla[
    !is.na(summary_hla$search_tag) &
      trimws(summary_hla$search_tag) == trimws(allele),
    ,
    drop = FALSE
  ]
  
  if (nrow(sdat) > 0) {
    sdat <- unique(sdat)
  }
  
  if (nrow(sdat) == 0) {
    return(
      div(
        class = "summary-card",
        style = "width:100%; min-width:100%; box-sizing:border-box;",
        p("No summary annotation found for this allele in the local SummaryAnnotationHLA dataset.")
      )
    )
  }
  
  ## -----------------------------
  ## Evidence level count summary
  ## -----------------------------
  level_order <- c("1A", "1B", "2A", "2B", "3", "4")
  
  level_dat <- unique(
    sdat[, c("summary_annotation_id", "level_of_evidence"), drop = FALSE]
  )
  level_dat$level_of_evidence <- safe_trim(level_dat$level_of_evidence)
  level_dat <- level_dat[
    !is.na(level_dat$level_of_evidence) & level_dat$level_of_evidence != "",
    ,
    drop = FALSE
  ]
  
  if (nrow(level_dat) == 0) {
    evidence_levels <- "N/A"
  } else {
    level_tab <- table(factor(level_dat$level_of_evidence, levels = level_order))
    level_tab <- level_tab[level_tab > 0]
    
    if (length(level_tab) == 0) {
      evidence_levels <- "N/A"
    } else {
      evidence_levels <- paste0(
        "Level ", names(level_tab), " (", as.integer(level_tab), ")",
        collapse = "; "
      )
    }
  }
  
  modifiers <- collapse_pretty(sdat$level_modifiers)
  drugs_txt <- collapse_pretty(sdat$drugs)
  pheno_txt <- collapse_pretty(sdat$phenotypes)
  
  total_summary <- length(unique(sdat$summary_annotation_id))
  total_pmids <- length(unique(split_clean_unique(sdat$aggregated_evidence_pmids)))
  drug_vec <- expand_drug_vector(sdat$drugs)
  total_drugs <- length(unique(drug_vec))
  
  ext_url <- make_summary_search_url(allele)
  
  tagList(
    div(
      class = "summary-card",
      
      div(
        style = "display:flex; justify-content:flex-end; align-items:flex-start; margin-bottom:10px;",
        span(class = "summary-badge", "Summary Annotation Available")
      ),
      
      div(
        class = "overview-flex",
        
        ## LEFT
        div(
          class = "overview-left",
          
          tags$p(
            style = "font-size:15px; color:#444; margin:0 0 18px 0; line-height:1.7;",
            "Summary annotations integrate evidence from multiple variant annotations and clinical prescribing guidelines to describe genotype-based drug response risks. Each record represents a curated interpretation of the phenotypic impact supported by literature evidence levels and annotation scores."
          ),
          
          div(
            class = "overview-metrics",
            div(class = "mini-card", div(class = "mini-title", "Summary Records"), div(class = "mini-value", total_summary)),
            div(class = "mini-card", div(class = "mini-title", "Unique Drugs"), div(class = "mini-value", total_drugs)),
            div(class = "mini-card", div(class = "mini-title", "Unique Literature"), div(class = "mini-value", total_pmids))
          ),
          
          h4("Clinical Overview"),
          p(tags$b("Associated drug(s): "), drugs_txt),
          p(tags$b("Phenotype / risk: "), pheno_txt),
          p(tags$b("Evidence level(s): "), evidence_levels),
          p(tags$b("Level modifier(s): "), ifelse(modifiers == "", "N/A", modifiers)),
          
          div(
            style = "display:flex; flex-wrap:wrap; gap:10px; align-items:center; margin-top:18px;",
            safe_external_link(
              href = ext_url,
              label = "Search full record in ClinPGX",
              class = "btn btn-danger btn-sm"
            ),
            download_button_compact("download_summary_csv", "Download")
          )
        ),
        
        ## RIGHT
        div(
          class = "overview-right",
          plotOutput("summary_annotation_plot", height = "320px")
        )
      )
    )
  )
})

output$summary_annotation_tbl <- reactable::renderReactable({
  allele <- selected_allele()
  req(allele)
  
  sdat <- summary_hla[
    !is.na(summary_hla$search_tag) &
      trimws(summary_hla$search_tag) == trimws(allele),
    ,
    drop = FALSE
  ]
  
  if (nrow(sdat) > 0) {
    sdat <- unique(sdat)
  }
  
  if (nrow(sdat) == 0) {
    return(
      reactable::reactable(
        data.frame(Message = "No summary annotation found for this allele."),
        searchable = FALSE,
        pagination = FALSE,
        bordered = TRUE,
        striped = TRUE
      )
    )
  }
  
  keep_df <- data.frame(
    summary_annotation_id     = sdat$summary_annotation_id,
    level_of_evidence         = sdat$level_of_evidence,
    variant                   = sdat$variant_haplotypes,
    gene                      = sdat$gene,
    drugs                     = sdat$drugs,
    phenotype_category        = sdat$phenotype_category,
    phenotypes                = sdat$phenotypes,
    specialty_population      = sdat$specialty_population,
    latest_history_date       = sdat$latest_history_date,
    summary_url               = sdat$summary_url,
    annotation_text           = sdat$annotation_text,
    allele_function           = sdat$allele_function,
    aggregated_evidence_ids   = sdat$aggregated_evidence_ids,
    aggregated_evidence_types = sdat$aggregated_evidence_types,
    aggregated_evidence_pmids = sdat$aggregated_evidence_pmids,
    stringsAsFactors = FALSE
  )
  
  keep_df <- unique(keep_df)
  
  parent_df <- keep_df %>%
    dplyr::group_by(summary_annotation_id, variant) %>%
    dplyr::summarise(
      level_of_evidence         = collapse_one(level_of_evidence),
      gene                      = collapse_one(gene),
      drugs                     = collapse_one(drugs),
      phenotype_category        = collapse_one(phenotype_category),
      phenotypes                = collapse_one(phenotypes),
      specialty_population      = collapse_one(specialty_population),
      latest_history_date       = collapse_one(latest_history_date),
      annotation_text           = collapse_one(annotation_text),
      allele_function           = collapse_one(allele_function),
      aggregated_evidence_ids   = collapse_one(aggregated_evidence_ids),
      aggregated_evidence_types = collapse_one(aggregated_evidence_types),
      aggregated_evidence_pmids = collapse_one(aggregated_evidence_pmids),
      summary_url               = collapse_one(summary_url),
      pediatric                 = ifelse(
        any(grepl("Pediatric", specialty_population, ignore.case = TRUE), na.rm = TRUE),
        "TRUE", "FALSE"
      ),
      .groups = "drop"
    )
  
  parent_df <- parent_df %>%
    dplyr::select(
      summary_url,
      level_of_evidence,
      variant,
      gene,
      drugs,
      phenotype_category,
      phenotypes,
      pediatric,
      everything()
    )
  
  reactable::reactable(
    parent_df,
    searchable = FALSE,   # ⭐ turn off top search
    filterable = FALSE,
    striped = TRUE,
    bordered = TRUE,
    highlight = TRUE,
    compact = TRUE,
    defaultPageSize = 10,
    columns = list(
      summary_url = reactable::colDef(
        name = "",
        minWidth = 90,
        sortable = FALSE,
        filterable = FALSE,
        cell = function(value) {
          safe_external_link(
            href = value,
            label = "Details",
            class = "btn btn-primary btn-sm"
          )
        }
      ),
      
      level_of_evidence = reactable::colDef(
        name = "Level",
        minWidth = 80
      ),
      
      variant = reactable::colDef(
        name = "Variant",
        minWidth = 200
      ),
      
      gene = reactable::colDef(
        name = "Gene",
        minWidth = 100
      ),
      
      drugs = reactable::colDef(
        name = "Drugs",
        minWidth = 160,
        cell = function(value) {
          make_drug_links(value)
        }
      ),
      
      phenotype_category = reactable::colDef(
        name = "Phenotype Categories",
        minWidth = 140
      ),
      
      phenotypes = reactable::colDef(
        name = "Phenotype",
        minWidth = 260,
        cell = function(value) {
          make_phenotype_links(value, use_breaks = TRUE)
        }
      ),
      
      pediatric = reactable::colDef(
        name = "Pediatric",
        minWidth = 90
      ),
      
      summary_annotation_id = reactable::colDef(show = FALSE),
      specialty_population = reactable::colDef(show = FALSE),
      latest_history_date = reactable::colDef(show = FALSE),
      annotation_text = reactable::colDef(show = FALSE),
      allele_function = reactable::colDef(show = FALSE),
      aggregated_evidence_ids = reactable::colDef(show = FALSE),
      aggregated_evidence_types = reactable::colDef(show = FALSE),
      aggregated_evidence_pmids = reactable::colDef(show = FALSE)
    ),
    
    details = function(index) {
      row <- parent_df[index, ]
      
      tags$div(
        style = "padding: 16px 20px; background:#fafafa; border-left:4px solid #dc3545;",
        
        h5("Interpretation"),
        p(ifelse(row$annotation_text == "", "N/A", row$annotation_text)),
        
        tags$hr(),
        
        p(tags$b("Last updated: "), ifelse(row$latest_history_date == "", "N/A", row$latest_history_date)),
        
        p(tags$b("Specialty population: "), ifelse(row$specialty_population == "", "N/A", row$specialty_population)),
        
        p(tags$b("Allele function: "), ifelse(row$allele_function == "", "N/A", row$allele_function)),
        
        p(tags$b("Evidence type: "), ifelse(row$aggregated_evidence_types == "", "N/A", row$aggregated_evidence_types)),
        
        tags$p(tags$b("PMID(s): "), make_pmid_links(row$aggregated_evidence_pmids))
      )
    }
  )
})

output$summary_annotation_plot <- renderPlot({
  allele <- selected_allele()
  req(allele)
  
  sdat <- summary_hla[
    !is.na(summary_hla$search_tag) &
      trimws(summary_hla$search_tag) == trimws(allele),
    ,
    drop = FALSE
  ]
  
  if (nrow(sdat) == 0) return(NULL)
  
  sdat <- unique(sdat[, c("summary_annotation_id", "level_of_evidence")])
  sdat$level_of_evidence <- safe_trim(sdat$level_of_evidence)
  sdat <- sdat[sdat$level_of_evidence != "", , drop = FALSE]
  
  if (nrow(sdat) == 0) return(NULL)
  
  level_order <- c("1A", "1B", "2A", "2B", "3", "4")
  
  level_tab <- table(factor(sdat$level_of_evidence, levels = level_order))
  
  plot_df <- data.frame(
    level = factor(level_order, levels = rev(level_order)),
    count = as.integer(level_tab[level_order]),
    stringsAsFactors = FALSE
  )
  
  level_colors <- c(
    "1A" = "#2fb36d",
    "1B" = "#7fc97f",
    "2A" = "#2f73d0",
    "2B" = "#6baed6",
    "3"  = "#f4b63a",
    "4"  = "#c53a3a"
  )
  
  max_count <- max(plot_df$count, na.rm = TRUE)
  y_breaks <- 0:max(1, max_count)
  
  ggplot(plot_df, aes(x = level, y = count, fill = level)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = level_colors, guide = "none", drop = FALSE) +
    scale_x_discrete(labels = function(x) paste("Level", x)) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = as.character(y_breaks),
      limits = c(0, max(y_breaks)),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "Evidence level",
      y = "Number of summary annotations"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_blank(),
      axis.text.y = element_text(size = 13),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )
})

output$variant_annotation_overview <- renderUI({
  
  allele <- selected_allele()
  req(allele)
  
  vdat <- variant_hla[
    !is.na(variant_hla$search_tag) &
      trimws(variant_hla$search_tag) == trimws(allele),
    ,
    drop = FALSE
  ]
  
  if (nrow(vdat) > 0) {
    vdat <- unique(vdat)
  }
  
  if (nrow(vdat) == 0) {
    return(
      div(
        class = "summary-card",
        style = "width:100%; min-width:100%; box-sizing:border-box;",
        p("No variant annotation found for this allele.")
      )
    )
  }
  
  drug_long <- make_variant_drug_long(
    drug_strings = vdat$drugs,
    record_ids   = vdat$variant_annotation_id
  )
  
  lit_chr <- safe_trim(vdat$pmid)
  pheno_vec <- split_clean_variant_phenotypes(vdat$phenotype)
  pop_chr <- safe_trim(vdat$biogeographical_groups)
  
  total_var <- length(unique(vdat$variant_annotation_id))
  total_drugs <- length(unique(drug_long$drug))
  total_pmids <- length(unique(split_clean_unique(lit_chr)))
  
  top_drugs_df <- make_top_drugs_from_long(drug_long, top_n = 5)
  
  top_drugs_text <- if (nrow(top_drugs_df) == 0) {
    "N/A"
  } else {
    paste0(top_drugs_df$drug, " (", top_drugs_df$n, ")", collapse = "; ")
  }
  
  top_pheno_df <- data.frame(pheno = pheno_vec, stringsAsFactors = FALSE) %>%
    dplyr::filter(!is.na(pheno), pheno != "") %>%
    dplyr::count(pheno, sort = TRUE) %>%
    dplyr::slice_head(n = 6)
  
  top_pheno_text <- if (nrow(top_pheno_df) == 0) {
    "N/A"
  } else {
    paste(top_pheno_df$pheno, collapse = "; ")
  }
  
  pop_vec <- split_clean_unique(pop_chr)
  pop_vec <- pop_vec[!tolower(pop_vec) %in% c("unknown", "unk")]
  pop_text <- if (length(pop_vec) == 0) {
    "N/A"
  } else {
    paste(head(pop_vec, 4), collapse = "; ")
  }
  
  ext_url <- make_summary_search_url(allele)
  
  tagList(
    div(
      class = "summary-card",
      
      div(
        style = "display:flex; justify-content:flex-end; align-items:flex-start; margin-bottom:10px;",
        span(class = "summary-badge", "Variant Annotation Available")
      ),
      
      div(
        class = "overview-flex",
        ## LEFT
        div(
          class = "overview-left",
          
          tags$p(
            style = "font-size:15px; color:#444; margin:0 0 18px 0; line-height:1.7;",
            "Variant annotations capture literature-level associations between genetic variants and drug response phenotypes, including toxicity, efficacy, and pharmacokinetic outcomes. Each entry reflects curated evidence from individual studies and provides detailed context for interpreting variant-specific pharmacogenomic signals."
          ),
          
          div(
            class = "overview-metrics",
            div(class = "mini-card", div(class = "mini-title", "Variant Records"), div(class = "mini-value", total_var)),
            div(class = "mini-card", div(class = "mini-title", "Unique Drugs"), div(class = "mini-value", total_drugs)),
            div(class = "mini-card", div(class = "mini-title", "Unique Literature"), div(class = "mini-value", total_pmids))
          ),
          
          h4("Evidence Landscape"),
          p(tags$b("Top drug signals: "), top_drugs_text),
          p(tags$b("Major phenotype domains: "), top_pheno_text),
          p(tags$b("Population diversity: "), pop_text),
          
          div(
            style = "display:flex; flex-wrap:wrap; gap:10px; align-items:center; margin-top:18px;",
            safe_external_link(
              href = ext_url,
              label = "Search full record in ClinPGX",
              class = "btn btn-danger btn-sm"
            ),
            download_button_compact("download_variant_csv", "Download")
          )
        ),
        
        ## RIGHT
        div(
          class = "overview-right",
          plotOutput("variant_annotation_plot", height = "340px")
        )
      )
    )
  )
})

output$variant_annotation_tbl <- reactable::renderReactable({
  
  allele <- selected_allele()
  req(allele)
  
  vdat <- variant_hla[
    !is.na(variant_hla$search_tag) &
      trimws(variant_hla$search_tag) == trimws(allele),
    ,
    drop = FALSE
  ]
  
  if (nrow(vdat) > 0) {
    vdat <- unique(vdat)
  }
  
  if (nrow(vdat) == 0) {
    return(
      reactable::reactable(
        data.frame(Message = "No variant annotation found for this allele."),
        searchable = FALSE,
        pagination = FALSE,
        bordered = TRUE,
        striped = TRUE
      )
    )
  }
  
  keep_df <- data.frame(
    variant_annotation_id   = as.character(vdat$variant_annotation_id),
    detail_url              = make_variant_record_url(as.character(vdat$variant_annotation_id)),
    variant                 = as.character(vdat$variant_haplotypes),
    literature              = as.character(vdat$pmid),
    drugs                   = as.character(vdat$drugs),
    significance            = as.character(vdat$significance),
    p_value                 = as.character(vdat$study_p_values),
    cases                   = as.character(vdat$study_cases),
    controls                = as.character(vdat$study_controls),
    phenotype_category      = as.character(vdat$phenotype_category),
    pediatric               = ifelse(
      grepl("Pediatric", safe_trim(vdat$specialty_population), ignore.case = TRUE),
      "TRUE", "FALSE"
    ),
    association             = as.character(vdat$sentence),
    notes                   = clean_quoted_text_block(as.character(vdat$notes)),
    specialty_population    = as.character(vdat$specialty_population),
    biogeographical_groups  = as.character(vdat$biogeographical_groups),
    study_types             = as.character(vdat$study_types),
    study_ratio_stat_types  = as.character(vdat$study_ratio_stat_types),
    study_ratio_stats       = as.character(vdat$study_ratio_stats),
    study_ci_start          = as.character(vdat$study_ci_start),
    study_ci_stop           = as.character(vdat$study_ci_stop),
    comparison_genotype     = as.character(vdat$comparison_genotype),
    stringsAsFactors = FALSE
  )
  
  keep_df <- unique(keep_df)
  
  keep_df$significance_simple <- dplyr::case_when(
    trimws(tolower(keep_df$significance)) == "yes" ~ "Yes",
    trimws(tolower(keep_df$significance)) == "no" ~ "No",
    TRUE ~ "Not stated"
  )
  
  keep_df$p_value <- clean_pvalue_text(keep_df$p_value)
  
  display_df <- keep_df[, c(
    "detail_url",
    "variant",
    "literature",
    "drugs",
    "significance_simple",
    "p_value",
    "cases",
    "controls",
    "phenotype_category",
    "pediatric"
  ), drop = FALSE]
  
  reactable::reactable(
    data = display_df,
    rownames = FALSE,
    searchable = FALSE,
    filterable = FALSE,
    striped = TRUE,
    bordered = TRUE,
    highlight = TRUE,
    compact = TRUE,
    defaultPageSize = 10,
    columns = list(
      detail_url = reactable::colDef(
        name = "",
        minWidth = 90,
        sortable = FALSE,
        cell = function(value) {
          safe_external_link(
            href = value,
            label = "Details",
            class = "btn btn-primary btn-sm"
          )
        }
      ),
      
      variant = reactable::colDef(
        name = "Variant",
        minWidth = 190
      ),
      
      literature = reactable::colDef(
        name = "Literature",
        minWidth = 130,
        cell = function(value) make_literature_links(value)
      ),
      
      drugs = reactable::colDef(
        name = "Drugs",
        minWidth = 150,
        cell = function(value) {
          make_drug_links(value)
        }
      ),
      
      significance_simple = reactable::colDef(
        name = "Significance",
        minWidth = 100
      ),
      
      p_value = reactable::colDef(
        name = "P-Value",
        minWidth = 100
      ),
      
      cases = reactable::colDef(
        name = "# of Cases",
        minWidth = 95
      ),
      
      controls = reactable::colDef(
        name = "# of Controls",
        minWidth = 105
      ),
      
      phenotype_category = reactable::colDef(
        name = "Phenotype Categories",
        minWidth = 150
      ),
      
      pediatric = reactable::colDef(
        name = "Pediatric",
        minWidth = 90
      )
    ),
    
    details = function(index) {
      row <- keep_df[index, ]
      
      tags$div(
        style = "padding: 16px 20px; background:#fafafa; border-left:4px solid #0d6efd;",
        
        h5("Association"),
        p(ifelse(is.na(row$association) || row$association == "", "N/A", row$association)),
        
        tags$hr(),
        
        h5("Additional Context"),
        p(
          tags$b("More details: "),
          ifelse(is.na(row$notes) || row$notes == "", "N/A", row$notes)
        ),
        p(
          tags$b("Biogeographical groups: "),
          ifelse(is.na(row$biogeographical_groups) || row$biogeographical_groups == "", "N/A", row$biogeographical_groups)
        ),
        p(
          tags$b("Specialty population: "),
          ifelse(is.na(row$specialty_population) || row$specialty_population == "", "N/A", row$specialty_population)
        ),
        p(
          tags$b("Comparison allele(s) / genotype(s): "),
          ifelse(is.na(row$comparison_genotype) || row$comparison_genotype == "", "N/A", row$comparison_genotype)
        ),
        
        tags$hr(),
        
        h5("Study Information"),
        p(
          tags$b("Study type: "),
          ifelse(is.na(row$study_types) || row$study_types == "", "N/A", row$study_types)
        ),
        p(
          tags$b("Ratio statistic: "),
          paste0(
            ifelse(is.na(row$study_ratio_stat_types) || row$study_ratio_stat_types == "", "N/A", row$study_ratio_stat_types),
            " / ",
            ifelse(is.na(row$study_ratio_stats) || row$study_ratio_stats == "", "N/A", row$study_ratio_stats)
          )
        ),
        p(
          tags$b("Confidence interval: "),
          paste0(
            ifelse(is.na(row$study_ci_start) || row$study_ci_start == "", "N/A", row$study_ci_start),
            " - ",
            ifelse(is.na(row$study_ci_stop) || row$study_ci_stop == "", "N/A", row$study_ci_stop)
          )
        )
      )
    }
  )
  
})

output$variant_annotation_plot <- renderPlot({
  allele <- selected_allele()
  req(allele)
  
  vdat <- variant_hla[
    !is.na(variant_hla$search_tag) &
      trimws(variant_hla$search_tag) == trimws(allele),
    ,
    drop = FALSE
  ]
  
  if (nrow(vdat) == 0) return(NULL)
  
  vdat <- unique(vdat)
  
  drug_long <- make_variant_drug_long(
    drug_strings = vdat$drugs,
    record_ids   = vdat$variant_annotation_id
  )
  
  plot_df <- make_top_drugs_from_long(drug_long, top_n = 5)
  
  if (nrow(plot_df) == 0) return(NULL)
  
  plot_df$drug <- factor(plot_df$drug, levels = rev(plot_df$drug))
  
  drug_colors <- c("#4db3d9", "#f46a2a", "#e83e8c", "#4f79bd", "#a6b93f")
  
  max_count <- max(plot_df$n, na.rm = TRUE)
  x_breaks <- 0:max(1, max_count)
  
  ggplot(plot_df, aes(x = drug, y = n, fill = drug)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(values = drug_colors[seq_len(nrow(plot_df))], guide = "none") +
    scale_y_continuous(
      breaks = x_breaks,
      labels = as.integer(x_breaks),
      limits = c(0, max(x_breaks)),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "Drug",
      y = "Number of variant annotations"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_blank(),
      axis.text.y = element_text(size = 13),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )
})

output$druglabel_annotation_overview <- renderUI({
  
  allele <- selected_allele()
  req(allele)
  
  allele_gene <- allele_to_gene(allele)
  
  ddat <- druglabel_hla[
    !is.na(druglabel_hla$search_tag) &
      trimws(druglabel_hla$search_tag) %in% c(trimws(allele), trimws(allele_gene)),
    ,
    drop = FALSE
  ]
  
  if (nrow(ddat) > 0) {
    ddat <- unique(ddat)
  }
  
  if (nrow(ddat) == 0) {
    return(
      div(
        class = "summary-card",
        style = "width:100%; min-width:100%; box-sizing:border-box;",
        p("No drug label annotation found for this allele.")
      )
    )
  }
  
  ddat$match_scope <- ifelse(trimws(ddat$search_tag) == trimws(allele), "Allele-level", "Gene-level")
  
  drug_vec <- expand_drug_vector(ddat$Chemicals)
  
  total_records <- length(unique(ddat$`PharmGKB ID`))
  total_drugs <- length(unique(drug_vec))
  total_allele_level <- sum(ddat$match_scope == "Allele-level", na.rm = TRUE)
  
  top_drugs_df <- make_top_drugs_df(ddat$Chemicals, top_n = 5)
  
  top_drugs_text <- if (nrow(top_drugs_df) == 0) {
    "N/A"
  } else {
    paste0(top_drugs_df$drug, " (", top_drugs_df$n, ")", collapse = "; ")
  }
  
  source_text <- collapse_pretty(ddat$Source)
  pgx_text    <- collapse_pretty(ddat$`Testing Level`)
  
  tagList(
    div(
      class = "summary-card",
      
      div(
        style = "display:flex; justify-content:flex-end; align-items:flex-start;",
        span(class = "summary-badge", "Drug Label Annotation Available")
      ),
      
      tags$p(
        style = "font-size:15px; color:#444; margin:6px 0 18px 0; line-height:1.7;",
        "Drug label annotations summarize regulatory pharmacogenomic information extracted from approved drug labels across major agencies. These records highlight clinically relevant recommendations, testing requirements, and actionable PGx guidance linked to HLA variation."
      ),
      
      div(
        style = "display:grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap:14px; margin-bottom:18px;",
        div(class = "mini-card", div(class = "mini-title", "Drug Label Records"), div(class = "mini-value", total_records)),
        div(class = "mini-card", div(class = "mini-title", "Unique Drugs"), div(class = "mini-value", total_drugs)),
        div(class = "mini-card", div(class = "mini-title", "Allele-level Labels"), div(class = "mini-value", total_allele_level))
      ),
      
      div(
        style = "display:flex; flex-wrap:wrap; gap:24px;",
        div(
          style = "flex:1; min-width:380px;",
          h4("Regulatory Landscape"),
          p(tags$b("Major regulatory agencies: "), source_text),
          p(tags$b("Top labeled drugs: "), top_drugs_text),
          p(tags$b("PGx action levels: "), pgx_text)
        ),
        div(
          style = "flex:0.9; min-width:280px;",
          h4("External Resource"),
          div(
            style = "display:flex; flex-wrap:wrap; gap:10px; align-items:center;",
            safe_external_link(
              href = make_summary_search_url(allele),
              label = "Search full record in ClinPGX",
              class = "btn btn-danger btn-sm"
            ),
            download_button_compact("download_druglabel_csv", "Download")
          )
        )
      )
    )
  )
})

output$druglabel_annotation_tbl <- reactable::renderReactable({
  
  allele <- selected_allele()
  req(allele)
  
  allele_gene <- allele_to_gene(allele)
  
  ddat <- druglabel_hla[
    !is.na(druglabel_hla$search_tag) &
      trimws(druglabel_hla$search_tag) %in% c(trimws(allele), trimws(allele_gene)),
    ,
    drop = FALSE
  ]
  
  if (nrow(ddat) > 0) {
    ddat <- unique(ddat)
  }
  
  if (nrow(ddat) == 0) {
    return(
      reactable::reactable(
        data.frame(Message = "No drug label annotation found for this allele."),
        searchable = FALSE,
        pagination = FALSE,
        bordered = TRUE,
        striped = TRUE
      )
    )
  }
  
  keep_df <- data.frame(
    pharmgkb_id                  = as.character(ddat$`PharmGKB ID`),
    detail_url                   = make_druglabel_url(as.character(ddat$`PharmGKB ID`)),
    pgx_level                    = as.character(ddat$`Testing Level`),
    source                       = as.character(ddat$Source),
    title                        = as.character(ddat$Name),
    genes                        = as.character(ddat$Genes),
    drugs                        = as.character(ddat$Chemicals),
    biomarker_flag               = as.character(ddat$`Biomarker Flag`),
    has_prescribing_info         = as.character(ddat$`Has Prescribing Info`),
    has_dosing_info              = as.character(ddat$`Has Dosing Info`),
    has_alternate_drug           = as.character(ddat$`Has Alternate Drug`),
    has_other_prescribing        = as.character(ddat$`Has Other Prescribing Guidance`),
    cancer_genome                = as.character(ddat$`Cancer Genome`),
    prescribing                  = as.character(ddat$Prescribing),
    variants_haplotypes          = as.character(ddat$`Variants/Haplotypes`),
    latest_history_date          = as.character(ddat$`Latest History Date (YYYY-MM-DD)`),
    search_tag                   = as.character(ddat$search_tag),
    match_scope                  = ifelse(trimws(ddat$search_tag) == trimws(allele), "Allele-level", "Gene-level"),
    stringsAsFactors = FALSE
  )
  
  keep_df <- unique(keep_df)
  
  keep_df$match_scope <- factor(
    keep_df$match_scope,
    levels = c("Allele-level", "Gene-level")
  )
  
  keep_df <- keep_df %>%
    arrange(match_scope, source, drugs, title)
  
  reactable::reactable(
    keep_df,
    searchable = FALSE,
    filterable = FALSE,
    striped = TRUE,
    bordered = TRUE,
    highlight = TRUE,
    compact = TRUE,
    defaultPageSize = 10,
    columns = list(
      detail_url = reactable::colDef(
        name = "",
        minWidth = 90,
        sortable = FALSE,
        cell = function(value) {
          safe_external_link(
            href = value,
            label = "Details",
            class = "btn btn-primary btn-sm"
          )
        }
      ),
      
      pgx_level = reactable::colDef(
        name = "PGx Level",
        minWidth = 150,
        cell = function(value) make_pgx_level_badge(value)
      ),
      
      source = reactable::colDef(
        name = "Source",
        minWidth = 110
      ),
      
      title = reactable::colDef(
        name = "Title",
        minWidth = 340,
        cell = function(value, index) {
          make_druglabel_title_cell(
            title            = value,
            url              = keep_df$detail_url[index],
            has_prescribing  = keep_df$has_prescribing_info[index],
            has_dosing       = keep_df$has_dosing_info[index],
            has_alternate    = keep_df$has_alternate_drug[index],
            has_other        = keep_df$has_other_prescribing[index],
            biomarker_flag   = keep_df$biomarker_flag[index],
            match_scope      = as.character(keep_df$match_scope[index])
          )
        }
      ),
      
      drugs = reactable::colDef(
        name = "Drugs",
        minWidth = 160,
        cell = function(value) {
          make_drug_links(value)
        }
      ),
      
      pharmgkb_id = reactable::colDef(show = FALSE),
      genes = reactable::colDef(show = FALSE),
      biomarker_flag = reactable::colDef(show = FALSE),
      has_prescribing_info = reactable::colDef(show = FALSE),
      has_dosing_info = reactable::colDef(show = FALSE),
      has_alternate_drug = reactable::colDef(show = FALSE),
      has_other_prescribing = reactable::colDef(show = FALSE),
      cancer_genome = reactable::colDef(show = FALSE),
      prescribing = reactable::colDef(show = FALSE),
      variants_haplotypes = reactable::colDef(show = FALSE),
      latest_history_date = reactable::colDef(show = FALSE),
      search_tag = reactable::colDef(show = FALSE),
      match_scope = reactable::colDef(show = FALSE)
    ),
    
    details = function(index) {
      row <- keep_df[index, ]
      
      tags$div(
        style = "padding: 16px 20px; background:#fafafa; border-left:4px solid #6c757d;",
        
        h5("Annotation Metadata"),
        p(tags$b("Match scope: "), row$match_scope),
        p(tags$b("Last updated: "), ifelse(is.na(row$latest_history_date) || row$latest_history_date == "", "N/A", row$latest_history_date)),
        p(tags$b("Biomarker flag: "), ifelse(is.na(row$biomarker_flag) || row$biomarker_flag == "", "N/A", row$biomarker_flag)),
        
        tags$hr(),
        
        h5("Label Content Summary"),
        p(tags$b("Prescribing info available: "), ifelse(is.na(row$has_prescribing_info) || row$has_prescribing_info == "", "FALSE", "TRUE")),
        p(tags$b("Dosing info available: "), ifelse(is.na(row$has_dosing_info) || row$has_dosing_info == "", "FALSE", "TRUE")),
        p(tags$b("Alternate drug available: "), ifelse(is.na(row$has_alternate_drug) || row$has_alternate_drug == "", "FALSE", "TRUE")),
        p(tags$b("Other prescribing guidance: "), ifelse(is.na(row$has_other_prescribing) || row$has_other_prescribing == "", "FALSE", "TRUE")),
        p(tags$b("Cancer genome flag: "), ifelse(is.na(row$cancer_genome) || row$cancer_genome == "", "FALSE", "TRUE")),
        
        tags$hr(),
        
        h5("HLA / Label Context"),
        p(tags$b("Genes: "), ifelse(is.na(row$genes) || row$genes == "", "N/A", row$genes)),
        p(tags$b("Variants / Haplotypes: "), ifelse(is.na(row$variants_haplotypes) || row$variants_haplotypes == "", "N/A", row$variants_haplotypes)),
        p(tags$b("Prescribing category: "), ifelse(is.na(row$prescribing) || row$prescribing == "", "N/A", row$prescribing)),
      )
    }
  )
})

output$download_summary_csv <- downloadHandler(
  filename = function() {
    allele <- selected_allele()
    paste0(safe_allele_filename(allele), "_summary_annotations.csv")
  },
  content = function(file) {
    allele <- selected_allele()
    
    sdat <- summary_hla[
      !is.na(summary_hla$search_tag) &
        trimws(summary_hla$search_tag) == trimws(allele),
      ,
      drop = FALSE
    ]
    
    if (nrow(sdat) > 0) {
      sdat <- unique(sdat)
    }
    
    utils::write.csv(sdat, file, row.names = FALSE, na = "")
  }
)

output$download_variant_csv <- downloadHandler(
  filename = function() {
    allele <- selected_allele()
    paste0(safe_allele_filename(allele), "_variant_annotations.csv")
  },
  content = function(file) {
    allele <- selected_allele()
    
    vdat <- variant_hla[
      !is.na(variant_hla$search_tag) &
        trimws(variant_hla$search_tag) == trimws(allele),
      ,
      drop = FALSE
    ]
    
    if (nrow(vdat) > 0) {
      vdat <- unique(vdat)
    }
    
    utils::write.csv(vdat, file, row.names = FALSE, na = "")
  }
)

output$download_druglabel_csv <- downloadHandler(
  filename = function() {
    allele <- selected_allele()
    paste0(safe_allele_filename(allele), "_drug_label_annotations.csv")
  },
  content = function(file) {
    allele <- selected_allele()
    allele_gene <- allele_to_gene(allele)
    
    ddat <- druglabel_hla[
      !is.na(druglabel_hla$search_tag) &
        trimws(druglabel_hla$search_tag) %in% c(trimws(allele), trimws(allele_gene)),
      ,
      drop = FALSE
    ]
    
    if (nrow(ddat) > 0) {
      ddat <- unique(ddat)
    }
    
    utils::write.csv(ddat, file, row.names = FALSE, na = "")
  }
)

observeEvent(input$pharmgkb_from_table, {
  req(input$pharmgkb_from_table)
  
  allele_clicked <- trimws(as.character(input$pharmgkb_from_table))
  
  pharmgkb_panel_visible(TRUE)
  
  if ("allele_input" %in% names(input)) {
    updateTextInput(
      session,
      inputId = "allele_input",
      value = allele_clicked
    )
  }
  
  selected_allele(allele_clicked)
  
  session$onFlushed(function() {
    updateTabsetPanel(session, "pharmgkb_tabs", selected = "Summary Annotations")
    session$sendCustomMessage("resetPharmgkbButtons", list())
  }, once = TRUE)
})

output$pharmgkb_box_ui <- renderUI({

  if (!pharmgkb_panel_visible()) {
    return(NULL)
  }

  column(
    width = 12,
    box(
      id = "pharmgkb_box",
      height = "auto",
      width = 12,
      elevation = 4,
      title = strong("PharmGKB Database"),
      collapsible = TRUE,
      solidHeader = FALSE,
      maximizable = TRUE,
      status = "warning",

      div(
        style = "width:100%; display:block;",
        tabsetPanel(
          id = "pharmgkb_tabs",
          type = "tabs",

          tabPanel(
            "Summary Annotations",
            tags$br(),
            uiOutput("summary_annotation_overview"),
            reactable::reactableOutput("summary_annotation_tbl")
          ),

          tabPanel(
            "Variant Annotations",
            tags$br(),
            uiOutput("variant_annotation_overview"),
            reactable::reactableOutput("variant_annotation_tbl")
          ),

          tabPanel(
            "Drug Label Annotations",
            tags$br(),
            uiOutput("druglabel_annotation_overview"),
            reactable::reactableOutput("druglabel_annotation_tbl")
          )
        )
      )
    )
  )
})

outputOptions(output, "summary_annotation_overview", suspendWhenHidden = FALSE)
outputOptions(output, "summary_annotation_tbl", suspendWhenHidden = FALSE)
outputOptions(output, "summary_annotation_plot", suspendWhenHidden = FALSE)

outputOptions(output, "variant_annotation_overview", suspendWhenHidden = FALSE)
outputOptions(output, "variant_annotation_tbl", suspendWhenHidden = FALSE)
outputOptions(output, "variant_annotation_plot", suspendWhenHidden = FALSE)

outputOptions(output, "druglabel_annotation_overview", suspendWhenHidden = FALSE)
outputOptions(output, "druglabel_annotation_tbl", suspendWhenHidden = FALSE)
