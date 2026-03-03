# 04_panel_optimization.R
# Statistically determine optimal panel sizes for African populations
# Based on GLOBOCAN 2022 Top 15 cancer types

library(tidyverse)
library(here)

# Ensure dplyr functions take precedence
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Load comprehensive gene database
gene_scores <- read_csv(
  here("data", "processed", "comprehensive_gene_database_top15.csv"),
  show_col_types = FALSE
)

# Load TCGA mutations
tcga_combined <- read_csv(
  here("data", "processed", "tcga_mutations_top15_africa.csv"),
  show_col_types = FALSE
)

cat("Data loaded successfully\n")
cat("  Genes:", nrow(gene_scores), "\n")
cat("  Mutations:", nrow(tcga_combined), "\n\n")

# Helper function: Quick coverage calculation ----

calculate_quick_coverage <- function(panel_genes, mutations_data) {
  
  panel_gene_list <- panel_genes$gene
  
  # Identify column names
  gene_col <- if ("gene" %in% names(mutations_data)) {
    "gene"
  } else if ("Hugo_Symbol" %in% names(mutations_data)) {
    "Hugo_Symbol"
  } else {
    stop("Cannot find gene column")
  }
  
  sample_col <- if ("sample_id" %in% names(mutations_data)) {
    "sample_id"
  } else if ("Tumor_Sample_Barcode" %in% names(mutations_data)) {
    "Tumor_Sample_Barcode"
  } else {
    stop("Cannot find sample ID column")
  }
  
  # Calculate coverage
  coverage_data <- mutations_data %>%
    mutate(in_panel = .data[[gene_col]] %in% panel_gene_list)
  
  total_samples <- n_distinct(coverage_data[[sample_col]])
  samples_with_panel <- n_distinct(coverage_data[[sample_col]][coverage_data$in_panel])
  
  mean_sample_coverage <- (samples_with_panel / total_samples) * 100
  
  return(list(
    mean_sample_coverage = mean_sample_coverage,
    total_samples = total_samples,
    samples_covered = samples_with_panel
  ))
}


# SECTION 1: Determine Optimal Panel Sizes Statistically ----


determine_optimal_panel_sizes <- function(gene_scores, mutations_data) {
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("STATISTICAL DETERMINATION OF OPTIMAL PANEL SIZES\n")
  cat(strrep("=", 80), "\n\n")
  
  # Prepare gene data with African prioritization
  gene_data <- gene_scores %>%
    mutate(
      # Calculate top 3 cancer relevance
      top3_cancer_freq = pmax(
        if_else(!is.na(freq_Breast), freq_Breast, 0),
        if_else(!is.na(freq_Prostate), freq_Prostate, 0),
        if_else(!is.na(freq_Cervical), freq_Cervical, 0),
        na.rm = TRUE
      ),
      
      # Actionability
      is_actionable = !is.na(oncokb_level) & oncokb_level %in% c("1", "2", "3"),
      
      # African priority score
      african_priority = case_when(
        is_actionable == TRUE & pancancer_frequency > 0.15 ~ clinical_relevance * 1.5,
        gene_category == "DNA Repair" ~ clinical_relevance * 1.4,
        top3_cancer_freq > 0.1 ~ clinical_relevance * 1.3,
        is_actionable == TRUE ~ clinical_relevance * 1.2,
        TRUE ~ clinical_relevance
      )
    ) %>%
    arrange(desc(african_priority))
  
  # METHOD 1: Elbow method
  cat("METHOD 1: Elbow Analysis\n")
  cat("Testing panel sizes from 20 to 200 genes...\n")
  
  panel_sizes_test <- seq(20, 200, by = 10)
  
  elbow_analysis <- map_dfr(panel_sizes_test, function(size) {
    
    if (size %% 50 == 0) cat("  Testing", size, "genes...\n")
    
    # Select top N genes
    panel_genes <- gene_data %>%
      slice_head(n = size)
    
    # Calculate coverage
    coverage_results <- calculate_quick_coverage(panel_genes, mutations_data)
    
    tibble(
      panel_size = size,
      mean_coverage = coverage_results$mean_sample_coverage,
      n_actionable = sum(panel_genes$is_actionable),
      n_dna_repair = sum(panel_genes$gene_category == "DNA Repair")
    )
  })
  
  # Calculate marginal benefit
  elbow_analysis <- elbow_analysis %>%
    mutate(
      marginal_coverage = c(NA, diff(mean_coverage)),
      marginal_coverage_per_gene = marginal_coverage / 10
    )
  
  # Find elbow point
  elbow_threshold <- max(elbow_analysis$marginal_coverage_per_gene, na.rm = TRUE) * 0.2
  elbow_point <- elbow_analysis %>%
    filter(marginal_coverage_per_gene < elbow_threshold) %>%
    slice_head(n = 1) %>%
    pull(panel_size)
  
  if (length(elbow_point) == 0) elbow_point <- 80
  
  cat("\n  Elbow point:", elbow_point, "genes\n\n")
  
  # METHOD 2: Coverage thresholds
  cat("METHOD 2: Coverage Threshold Analysis\n")
  
  threshold_70 <- elbow_analysis %>%
    filter(mean_coverage >= 70) %>%
    slice_head(n = 1) %>%
    pull(panel_size)
  
  threshold_80 <- elbow_analysis %>%
    filter(mean_coverage >= 80) %>%
    slice_head(n = 1) %>%
    pull(panel_size)
  
  threshold_90 <- elbow_analysis %>%
    filter(mean_coverage >= 90) %>%
    slice_head(n = 1) %>%
    pull(panel_size)
  
  if (length(threshold_70) == 0) threshold_70 <- 50
  if (length(threshold_80) == 0) threshold_80 <- 100
  if (length(threshold_90) == 0) threshold_90 <- 150
  
  cat("  70% coverage:", threshold_70, "genes\n")
  cat("  80% coverage:", threshold_80, "genes\n")
  cat("  90% coverage:", threshold_90, "genes\n\n")
  
  # METHOD 3: Actionable gene saturation
  cat("METHOD 3: Actionable Gene Saturation\n")
  
  total_actionable <- sum(gene_data$is_actionable, na.rm = TRUE)
  
  actionable_80 <- elbow_analysis %>%
    mutate(actionable_percent = (n_actionable / total_actionable) * 100) %>%
    filter(actionable_percent >= 80) %>%
    slice_head(n = 1) %>%
    pull(panel_size)
  
  if (length(actionable_80) == 0) actionable_80 <- 60
  
  cat("  Total actionable:", total_actionable, "\n")
  cat("  For 80% actionable:", actionable_80, "genes\n\n")
  
  # METHOD 4: Cost-effectiveness
  cat("METHOD 4: Cost-Effectiveness Optimization\n")
  
  cost_effectiveness <- elbow_analysis %>%
    mutate(
      total_cost = (panel_size * 5) + 150,
      efficiency_ratio = (n_actionable * mean_coverage) / total_cost
    )
  
  optimal_efficiency <- cost_effectiveness %>%
    slice_max(efficiency_ratio, n = 1) %>%
    pull(panel_size)
  
  cat("  Most efficient:", optimal_efficiency, "genes\n\n")
  
  # FINAL RECOMMENDATIONS
  cat(strrep("=", 80), "\n")
  cat("RECOMMENDED PANEL SIZES\n")
  cat(strrep("=", 80), "\n\n")
  
  # Small panel: minimum for essential coverage
  small_panel <- min(threshold_70, actionable_80, elbow_point, na.rm = TRUE)
  small_panel <- as.integer(round(small_panel / 5) * 5)  # Round to nearest 5
  small_panel <- max(30, min(60, small_panel))  # Bound 30-60
  
  # Medium panel: optimal cost-benefit
  medium_panel <- min(threshold_80, optimal_efficiency, na.rm = TRUE)
  medium_panel <- as.integer(round(medium_panel / 5) * 5)
  medium_panel <- max(70, min(110, medium_panel))  # Bound 70-110
  
  # Large panel: comprehensive
  large_panel <- min(threshold_90, 150, na.rm = TRUE)
  large_panel <- as.integer(round(large_panel / 5) * 5)
  large_panel <- max(130, min(180, large_panel))  # Bound 130-180
  
  # Ensure proper ordering and spacing
  if (medium_panel - small_panel < 20) {
    medium_panel <- small_panel + 30
  }
  if (large_panel - medium_panel < 30) {
    large_panel <- medium_panel + 40
  }
  
  cat("SMALL PANEL:", small_panel, "genes\n")
  cat("  Based on: 70% coverage threshold + actionable saturation\n")
  cat("  Target: Essential genes, resource-limited settings\n\n")
  
  cat("MEDIUM PANEL:", medium_panel, "genes\n")
  cat("  Based on: Optimal cost-effectiveness + 80% coverage\n")
  cat("  Target: Standard clinical use (RECOMMENDED)\n\n")
  
  cat("LARGE PANEL:", large_panel, "genes\n")
  cat("  Based on: 90% coverage + comprehensive profiling\n")
  cat("  Target: Research and comprehensive profiling\n\n")
  
  # Save analysis
  write_csv(
    elbow_analysis,
    here("results", "tables", "panel_size_elbow_analysis.csv")
  )
  
  write_csv(
    cost_effectiveness,
    here("results", "tables", "panel_size_cost_effectiveness.csv")
  )
  
  recommendation_summary <- tibble(
    panel_type = c("Small", "Medium", "Large"),
    recommended_size = c(small_panel, medium_panel, large_panel),
    rationale = c(
      "Essential actionable genes, minimum coverage",
      "Optimal cost-benefit, 80% coverage target",
      "Comprehensive profiling, 90% coverage target"
    ),
    statistical_basis = c(
      "Elbow + Actionable saturation",
      "Cost-effectiveness maximum",
      "Coverage threshold + Breakpoint"
    ),
    target_coverage = c(70, 80, 90),
    use_case = c(
      "Resource-limited settings",
      "Standard clinical use",
      "Research/comprehensive profiling"
    )
  )
  
  write_csv(
    recommendation_summary,
    here("results", "tables", "recommended_panel_sizes.csv")
  )
  
  return(list(
    recommended_sizes = c(small_panel, medium_panel, large_panel),
    elbow_analysis = elbow_analysis,
    cost_effectiveness = cost_effectiveness,
    summary = recommendation_summary
  ))
}

# Run statistical analysis
optimal_sizes <- determine_optimal_panel_sizes(gene_scores, tcga_combined)

# Extract recommended sizes as INTEGERS
small_panel_size <- as.integer(optimal_sizes$recommended_sizes[1])
medium_panel_size <- as.integer(optimal_sizes$recommended_sizes[2])
large_panel_size <- as.integer(optimal_sizes$recommended_sizes[3])

cat("Statistically-determined panel sizes:\n")
cat("  Small:", small_panel_size, "genes\n")
cat("  Medium:", medium_panel_size, "genes\n")
cat("  Large:", large_panel_size, "genes\n\n")

# SECTION 2: Select genes with African prioritization ----


select_panel_genes <- function(
    gene_scores,
    target_size = 50,
    min_priority_score = 0.3,
    african_weighted = TRUE
) {
  
  # ENSURE target_size is an integer
  target_size <- as.integer(round(target_size))
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat(sprintf("SELECTING %d-GENE PANEL FOR AFRICAN POPULATIONS\n", target_size))
  cat(strrep("=", 80), "\n\n")
  
  # Apply African-specific weighting
  if (african_weighted) {
    gene_scores <- gene_scores %>%
      mutate(
        # Calculate top 3 cancer relevance
        top3_cancer_freq = pmax(
          if_else(!is.na(freq_Breast), freq_Breast, 0),
          if_else(!is.na(freq_Prostate), freq_Prostate, 0),
          if_else(!is.na(freq_Cervical), freq_Cervical, 0),
          na.rm = TRUE
        ),
        
        # Actionability
        is_actionable = !is.na(oncokb_level) & oncokb_level %in% c("1", "2", "3"),
        
        african_priority = case_when(
          # Highest priority: Actionable + high frequency
          is_actionable == TRUE & pancancer_frequency > 0.15 ~ clinical_relevance * 1.5,
          
          # High priority: DNA repair
          gene_category == "DNA Repair" ~ clinical_relevance * 1.4,
          
          # High priority: Relevant to top 3 African cancers
          top3_cancer_freq > 0.1 ~ clinical_relevance * 1.3,
          
          # Moderate priority: Actionable genes
          is_actionable == TRUE ~ clinical_relevance * 1.2,
          
          # Standard priority
          TRUE ~ clinical_relevance
        )
      )
  } else {
    gene_scores <- gene_scores %>%
      mutate(
        is_actionable = !is.na(oncokb_level) & oncokb_level %in% c("1", "2", "3"),
        african_priority = clinical_relevance
      )
  }
  
  # TIER 1: Essential actionable genes
  tier1_genes <- gene_scores %>%
    filter(
      is_actionable == TRUE,
      african_priority >= min_priority_score * 4,
      gene_category %in% c("Oncogene", "Tumor Suppressor", "DNA Repair")
    ) %>%
    arrange(desc(african_priority))
  
  cat("Tier 1 (Essential Actionable):", nrow(tier1_genes), "genes\n")
  
  # TIER 2: DNA repair genes
  tier2_genes <- gene_scores %>%
    filter(
      gene_category == "DNA Repair",
      !gene %in% tier1_genes$gene,
      african_priority >= min_priority_score * 2
    ) %>%
    arrange(desc(african_priority))
  
  cat("Tier 2 (DNA Repair):", nrow(tier2_genes), "genes\n")
  
  # TIER 3: High-frequency drivers
  tier3_genes <- gene_scores %>%
    filter(
      !gene %in% c(tier1_genes$gene, tier2_genes$gene),
      gene_category %in% c("Oncogene", "Tumor Suppressor"),
      african_priority >= min_priority_score * 1.5,
      pancancer_frequency > 0.05
    ) %>%
    arrange(desc(african_priority))
  
  cat("Tier 3 (High-frequency drivers):", nrow(tier3_genes), "genes\n")
  
  # TIER 4: Chromatin remodeling
  tier4_genes <- gene_scores %>%
    filter(
      !gene %in% c(tier1_genes$gene, tier2_genes$gene, tier3_genes$gene),
      gene_category == "Chromatin Remodeling",
      african_priority >= min_priority_score
    ) %>%
    arrange(desc(african_priority))
  
  cat("Tier 4 (Chromatin Remodeling):", nrow(tier4_genes), "genes\n")
  
  # TIER 5: Additional priority genes
  tier5_genes <- gene_scores %>%
    filter(
      !gene %in% c(tier1_genes$gene, tier2_genes$gene, 
                   tier3_genes$gene, tier4_genes$gene),
      african_priority >= min_priority_score * 0.5
    ) %>%
    arrange(desc(african_priority))
  
  cat("Tier 5 (Additional priority):", nrow(tier5_genes), "genes\n")
  
  # Combine tiers to reach target size
  panel_genes <- bind_rows(
    tier1_genes %>% mutate(tier = "Tier 1: Essential Actionable"),
    tier2_genes %>% mutate(tier = "Tier 2: DNA Repair"),
    tier3_genes %>% mutate(tier = "Tier 3: High-Frequency Drivers"),
    tier4_genes %>% mutate(tier = "Tier 4: Chromatin Remodeling"),
    tier5_genes %>% mutate(tier = "Tier 5: Additional Priority")
  ) %>%
    distinct(gene, .keep_all = TRUE) %>%
    slice_head(n = target_size) %>%  # target_size is now guaranteed integer
    mutate(
      selection_rank = row_number(),
      african_optimized = TRUE
    )
  
  cat("\n")
  cat(strrep("-", 80), "\n")
  cat("FINAL PANEL COMPOSITION\n")
  cat(strrep("-", 80), "\n\n")
  
  # Summary by tier
  tier_summary <- panel_genes %>%
    count(tier) %>%
    mutate(percentage = round(n / target_size * 100, 1))
  
  print(tier_summary)
  
  cat("\nGene categories:\n")
  print(panel_genes %>% count(gene_category, sort = TRUE))
  
  cat("\nActionability:\n")
  cat("  Actionable:", sum(panel_genes$is_actionable, na.rm = TRUE), 
      sprintf("(%.1f%%)\n", sum(panel_genes$is_actionable, na.rm = TRUE) / target_size * 100))
  
  cat("\n")
  
  return(panel_genes)
}


# SECTION 3: Calculate panel coverage ----


calculate_panel_coverage <- function(panel_genes, all_mutations) {
  
  panel_gene_list <- panel_genes$gene
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("CALCULATING PANEL COVERAGE\n")
  cat(strrep("=", 80), "\n\n")
  
  cat("Panel size:", length(panel_gene_list), "genes\n")
  cat("Total mutations:", nrow(all_mutations), "\n\n")
  
  # Identify column names
  gene_col <- if ("gene" %in% names(all_mutations)) "gene" else "Hugo_Symbol"
  sample_col <- if ("sample_id" %in% names(all_mutations)) "sample_id" else "Tumor_Sample_Barcode"
  project_col <- if ("study_id" %in% names(all_mutations)) "study_id" else "cancer_type"
  globocan_col <- if ("globocan_rank" %in% names(all_mutations)) "globocan_rank" else NULL
  
  # Calculate coverage
  all_mutations <- all_mutations %>%
    mutate(in_panel = .data[[gene_col]] %in% panel_gene_list)
  
  # By project with GLOBOCAN ranking
  if (!is.null(globocan_col)) {
    coverage_stats <- all_mutations %>%
      group_by(
        project = .data[[project_col]],
        globocan_rank = .data[[globocan_col]]
      ) %>%
      summarise(
        total_mutations = n(),
        panel_mutations = sum(in_panel),
        coverage_percent = (panel_mutations / total_mutations) * 100,
        total_samples = n_distinct(.data[[sample_col]]),
        samples_with_panel_mutation = n_distinct(.data[[sample_col]][in_panel]),
        sample_coverage_percent = (samples_with_panel_mutation / total_samples) * 100,
        .groups = "drop"
      ) %>%
      arrange(globocan_rank)
  } else {
    coverage_stats <- all_mutations %>%
      group_by(project = .data[[project_col]]) %>%
      summarise(
        total_mutations = n(),
        panel_mutations = sum(in_panel),
        coverage_percent = (panel_mutations / total_mutations) * 100,
        total_samples = n_distinct(.data[[sample_col]]),
        samples_with_panel_mutation = n_distinct(.data[[sample_col]][in_panel]),
        sample_coverage_percent = (samples_with_panel_mutation / total_samples) * 100,
        .groups = "drop"
      )
  }
  
  cat("Coverage by cancer type:\n")
  print(coverage_stats %>%
          mutate(
            coverage_percent = round(coverage_percent, 2),
            sample_coverage_percent = round(sample_coverage_percent, 2)
          ))
  
  cat("\nSummary:\n")
  cat("  Mean sample coverage:", sprintf("%.2f%%\n", mean(coverage_stats$sample_coverage_percent)))
  
  if (!is.null(globocan_col)) {
    top3_avg <- coverage_stats %>%
      filter(globocan_rank <= 3) %>%
      summarise(avg = mean(sample_coverage_percent)) %>%
      pull(avg)
    cat("  Top 3 cancer coverage:", sprintf("%.2f%%\n", top3_avg))
  }
  
  cat("\n")
  
  return(coverage_stats)
}


# SECTION 4: Add pharmacogenomic variants ----


add_pharmacogenomic_variants <- function(panel_genes) {
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("ADDING PHARMACOGENOMIC VARIANTS\n")
  cat(strrep("=", 80), "\n\n")
  
  pgx_variants <- tribble(
    ~gene, ~variant, ~drug, ~effect, ~african_relevance,
    "CYP2D6", "*17, *29", "Tamoxifen", "Reduced metabolism", "High",
    "DPYD", "Y186C (rs115232898)", "5-Fluorouracil/Capecitabine", "Reduced DPD, increased toxicity", "High",
    "UGT1A1", "*28", "Irinotecan", "Reduced glucuronidation, increased toxicity", "High",
    "TPMT", "*8", "6-Mercaptopurine", "Reduced metabolism, increased toxicity", "High",
    "CYP2C8", "*2", "Paclitaxel", "Reduced metabolism, increased toxicity", "High",
    "CYP2B6", "*6", "Cyclophosphamide", "Altered metabolism, increased toxicity", "High",
    "G6PD", "African variants", "Multiple chemotherapies", "Hemolytic risk", "High",
    "CYP3A5", "*3", "Multiple", "Drug metabolism", "High",
    "CYP2C19", "*2, *3", "Clopidogrel", "Reduced activation", "Moderate",
    "SLCO1B1", "rs4149056", "Statins", "Myopathy risk", "Moderate"
  )
  
  # Add missing PGx genes
  pgx_genes <- unique(pgx_variants$gene)
  missing_pgx <- pgx_genes[!pgx_genes %in% panel_genes$gene]
  
  if (length(missing_pgx) > 0) {
    cat("Adding", length(missing_pgx), "PGx genes:", paste(missing_pgx, collapse = ", "), "\n\n")
    
    new_pgx_genes <- tibble(
      gene = missing_pgx,
      tier = "Tier 6: Pharmacogenomics",
      gene_category = "Pharmacogenomics",
      is_actionable = TRUE,
      african_priority = 0.9,
      selection_rank = max(panel_genes$selection_rank) + seq_along(missing_pgx),
      african_optimized = TRUE
    )
    
    panel_genes <- bind_rows(panel_genes, new_pgx_genes)
    cat("New panel size:", nrow(panel_genes), "genes\n\n")
  } else {
    cat("All PGx genes already in panel\n\n")
  }
  
  return(list(panel = panel_genes, pgx_variants = pgx_variants))
}


# SECTION 5: Generate panels with optimal sizes ----


cat("\n")
cat(strrep("=", 80), "\n")
cat("GENERATING PANELS\n")
cat(strrep("=", 80), "\n")

# Generate panels
panel_small <- select_panel_genes(gene_scores, target_size = small_panel_size)
write_csv(panel_small, here("results", "tables", 
                            sprintf("panel_%d_genes.csv", small_panel_size)))

panel_medium <- select_panel_genes(gene_scores, target_size = medium_panel_size)
write_csv(panel_medium, here("results", "tables", 
                             sprintf("panel_%d_genes.csv", medium_panel_size)))

panel_large <- select_panel_genes(gene_scores, target_size = large_panel_size)
write_csv(panel_large, here("results", "tables", 
                            sprintf("panel_%d_genes.csv", large_panel_size)))


# SECTION 6: Coverage analysis ----


coverage_small <- calculate_panel_coverage(panel_small, tcga_combined) %>%
  mutate(panel_size = small_panel_size)

coverage_medium <- calculate_panel_coverage(panel_medium, tcga_combined) %>%
  mutate(panel_size = medium_panel_size)

coverage_large <- calculate_panel_coverage(panel_large, tcga_combined) %>%
  mutate(panel_size = large_panel_size)

all_coverage <- bind_rows(coverage_small, coverage_medium, coverage_large)
write_csv(all_coverage, here("results", "tables", "panel_coverage_stats.csv"))


# SECTION 7: Final cost-benefit ----


cost_benefit_final <- tibble(
  panel_size = c(small_panel_size, medium_panel_size, large_panel_size),
  panel_type = c("Small", "Medium", "Large"),
  total_cost = (panel_size * 5) + 150,
  n_actionable = c(
    sum(panel_small$is_actionable, na.rm = TRUE),
    sum(panel_medium$is_actionable, na.rm = TRUE),
    sum(panel_large$is_actionable, na.rm = TRUE)
  ),
  n_dna_repair = c(
    sum(panel_small$gene_category == "DNA Repair", na.rm = TRUE),
    sum(panel_medium$gene_category == "DNA Repair", na.rm = TRUE),
    sum(panel_large$gene_category == "DNA Repair", na.rm = TRUE)
  ),
  mean_sample_coverage = c(
    mean(coverage_small$sample_coverage_percent),
    mean(coverage_medium$sample_coverage_percent),
    mean(coverage_large$sample_coverage_percent)
  ),
  top3_coverage = c(
    mean(coverage_small$sample_coverage_percent[coverage_small$globocan_rank <= 3]),
    mean(coverage_medium$sample_coverage_percent[coverage_medium$globocan_rank <= 3]),
    mean(coverage_large$sample_coverage_percent[coverage_large$globocan_rank <= 3])
  )
) %>%
  mutate(
    cost_per_actionable = round(total_cost / n_actionable, 2),
    coverage_per_dollar = round(mean_sample_coverage / total_cost, 4),
    benefit_score = round((n_actionable * mean_sample_coverage) / total_cost, 2)
  )

write_csv(cost_benefit_final, here("results", "tables", "cost_benefit_analysis.csv"))


# SECTION 8: Add PGx to medium panel ----


panel_medium_with_pgx <- add_pharmacogenomic_variants(panel_medium)

write_csv(
  panel_medium_with_pgx$panel,
  here("results", "tables", 
       sprintf("panel_%d_with_pgx.csv", medium_panel_size))
)

write_csv(
  panel_medium_with_pgx$pgx_variants,
  here("results", "tables", "pharmacogenomic_variants.csv")
)


# SECTION 9: Final summary ----


cat("\n")
cat(strrep("=", 80), "\n")
cat("PANEL OPTIMIZATION COMPLETE\n")
cat(strrep("=", 80), "\n\n")

cat("STATISTICALLY-DETERMINED PANELS:\n\n")

cat("Small (", small_panel_size, " genes) - Budget-Conscious:\n", sep = "")
cat("  Actionable:", sum(panel_small$is_actionable, na.rm = TRUE), "\n")
cat("  Coverage:", sprintf("%.1f%%\n", mean(coverage_small$sample_coverage_percent)))
cat("  Cost: $", (small_panel_size * 5) + 150, "\n\n")

cat("Medium (", medium_panel_size, " genes) - RECOMMENDED:\n", sep = "")
cat("  Actionable:", sum(panel_medium$is_actionable, na.rm = TRUE), "\n")
cat("  Coverage:", sprintf("%.1f%%\n", mean(coverage_medium$sample_coverage_percent)))
cat("  Cost: $", (medium_panel_size * 5) + 150, "\n\n")

cat("Large (", large_panel_size, " genes) - Comprehensive:\n", sep = "")
cat("  Actionable:", sum(panel_large$is_actionable, na.rm = TRUE), "\n")
cat("  Coverage:", sprintf("%.1f%%\n", mean(coverage_large$sample_coverage_percent)))
cat("  Cost: $", (large_panel_size * 5) + 150, "\n\n")


