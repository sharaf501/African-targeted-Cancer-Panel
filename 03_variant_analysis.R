# 03_variant_analysis.R
# Analyze mutation frequencies and prioritize genes

library(tidyverse)
library(maftools)
library(here)

# Load data
tcga_mutations <- readRDS(here("data", "processed", "tcga_mutations_top15_africa.rds"))
gene_lit_freq <- read_csv(here("data", "processed", "gene_literature_frequency.csv"))
cancer_incidence <- read_csv(here("data", "processed", "cancer_incidence_africa_globocan2022.csv"))
curated_genes <- read_csv(here("data", "processed", "curated_cancer_genes_top15.csv"))
comprehensive_db <- read_csv(here("data", "processed", "comprehensive_gene_database_top15.csv"))


# Function: Calculate mutation frequencies from TCGA ----


calculate_mutation_frequencies_from_df <- function(mutations_df) {
  
  cat("Calculating mutation frequencies from data frame...\n")
  cat("Input:", nrow(mutations_df), "mutations\n\n")
  
  # Check columns
  cat("Available columns:", paste(names(mutations_df), collapse = ", "), "\n\n")
  
  # Calculate frequencies by study and gene
  freq_by_study <- mutations_df %>%
    group_by(study_id, cancer_type, gene) %>%
    summarise(
      n_mutations = n(),
      n_samples = n_distinct(sample_id),
      .groups = "drop"
    )
  
  # Get total samples per study
  total_samples_per_study <- mutations_df %>%
    group_by(study_id, cancer_type) %>%
    summarise(
      total_samples = n_distinct(sample_id),
      .groups = "drop"
    )
  
  # Calculate mutation frequency
  freq_df <- freq_by_study %>%
    left_join(total_samples_per_study, by = c("study_id", "cancer_type")) %>%
    mutate(
      mutation_frequency = n_samples / total_samples
    ) %>%
    select(gene, cancer_type, study_id, n_samples, total_samples, mutation_frequency) %>%
    rename(MutatedSamples = n_samples) %>%
    arrange(cancer_type, desc(mutation_frequency))
  
  cat("Calculated frequencies for", n_distinct(freq_df$gene), "genes\n")
  cat("Across", n_distinct(freq_df$cancer_type), "cancer types\n\n")
  
  return(freq_df)
}

# Calculate frequencies
tcga_frequencies <- calculate_mutation_frequencies_from_df(tcga_mutations)

# Save
write_csv(tcga_frequencies, here("data", "processed", "tcga_mutation_frequencies.csv"))


# Function: Calculate pan-cancer mutation scores ----


calculate_pancancer_scores <- function(freq_df) {
  
  pancancer_scores <- freq_df %>%
    group_by(gene) %>%
    summarise(
      n_cancer_types = n_distinct(cancer_type),
      mean_mutation_freq = mean(mutation_frequency, na.rm = TRUE),
      max_mutation_freq = max(mutation_frequency, na.rm = TRUE),
      total_mutations = sum(MutatedSamples, na.rm = TRUE),
      total_samples = sum(total_samples, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pancancer_frequency = total_mutations / total_samples,
      breadth_score = n_cancer_types / max(n_cancer_types),
      frequency_score = pancancer_frequency / max(pancancer_frequency),
      combined_score = (breadth_score + frequency_score) / 2
    ) %>%
    arrange(desc(combined_score))
  
  return(pancancer_scores)
}

pancancer_scores <- calculate_pancancer_scores(tcga_frequencies)
write_csv(pancancer_scores, here("data", "processed", "pancancer_mutation_scores.csv"))


# Function: Integrate African-specific data ----


integrate_african_data <- function(
    pancancer_scores,
    literature_freq,
    cancer_incidence,
    curated_genes
) {
  
  # Map TCGA cancer types to incidence data
  cancer_mapping <- tibble(
    cancer_type = c("BRCA", "PRAD", "COAD", "READ", "LIHC", "ESCA", 
                    "CESC", "LUAD", "LUSC", "STAD", "OV", "BLCA", "PAAD", "UCEC"),
    incidence_cancer = c("breast cancer", "prostate cancer", "colorectal cancer",
                         "colorectal cancer", "liver cancer", "esophageal cancer",
                         "cervical cancer", "lung cancer", "lung cancer",
                         "gastric cancer", "ovarian cancer", "bladder cancer",
                         "pancreatic cancer", "endometrial cancer")
  )
  
  # Join with incidence weights
  weighted_scores <- tcga_frequencies %>%
    left_join(cancer_mapping, by = "cancer_type") %>%
    left_join(cancer_incidence, by = c("incidence_cancer" = "cancer_type")) %>%
    mutate(
      incidence_weight = coalesce(priority_score, 5) / 10,
      weighted_mutation_freq = mutation_frequency * incidence_weight
    ) %>%
    group_by(gene) %>%
    summarise(
      african_weighted_score = sum(weighted_mutation_freq, na.rm = TRUE),
      n_relevant_cancers = n(),
      .groups = "drop"
    )
  
  # Combine all scores
  integrated_scores <- pancancer_scores %>%
    left_join(weighted_scores, by = "gene") %>%
    left_join(literature_freq, by = c("gene" = "genes")) %>%
    left_join(curated_genes %>% 
                mutate(in_curated_list = TRUE) %>%
                select(gene, in_curated_list),
              by = "gene") %>%
    mutate(
      literature_count = coalesce(literature_count, 0),
      in_curated_list = coalesce(in_curated_list, FALSE),
      african_weighted_score = coalesce(african_weighted_score, 0),
      
      # Normalize scores
      literature_score = literature_count / max(literature_count, na.rm = TRUE),
      curated_score = as.numeric(in_curated_list),
      #african_score = african_weighted_score / max(african_weighted_score, na.rm = TRUE),
      african_score = percent_rank(african_weighted_score),  # Better distribution
      
      # Combined priority score
      priority_score = (
        combined_score * 0.3 +              # Pan-cancer importance
          african_score * 0.3 +             # African cancer relevance
          literature_score * 0.2 +          # Literature support
          curated_score * 0.2               # Expert curation
      )
    ) %>%
    arrange(desc(priority_score))
  
  return(integrated_scores)
}

integrated_scores <- integrate_african_data(
  pancancer_scores,
  gene_lit_freq,
  cancer_incidence,
  curated_genes
)

# Check the new distribution
cat("\nAfrican score after percent_rank normalization:\n")
print(summary(integrated_scores$african_score))

write_csv(integrated_scores, here("data", "processed", "integrated_gene_scores.csv"))


# Function: Add gene annotations ----


# Gene annotation function

annotate_genes <- function(gene_scores) {
  
  library(biomaRt)
  library(org.Hs.eg.db)
  
  cat("Annotating", nrow(gene_scores), "genes...\n")
  
  # Get Ensembl annotation
  cat("Connecting to Ensembl...\n")
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  # Query gene information in batches (Ensembl has limits)
  genes_to_query <- unique(gene_scores$gene)
  cat("Querying", length(genes_to_query), "genes from Ensembl...\n")
  
  gene_info <- tryCatch({
    getBM(
      attributes = c(
        'hgnc_symbol',
        'chromosome_name',
        'start_position',
        'end_position',
        'gene_biotype',
        'description'
      ),
      filters = 'hgnc_symbol',
      values = genes_to_query,
      mart = ensembl
    )
  }, error = function(e) {
    cat("Warning: Ensembl query failed:", e$message, "\n")
    return(tibble(
      hgnc_symbol = character(),
      chromosome_name = character(),
      start_position = integer(),
      end_position = integer(),
      gene_biotype = character(),
      description = character()
    ))
  })
  
  cat("Retrieved info for", nrow(gene_info), "genes from Ensembl\n")
  
  # Gene Ontology terms
  cat("Querying Gene Ontology terms...\n")
  
  go_terms <- tryCatch({
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = genes_to_query,
      columns = c("SYMBOL", "GENENAME", "GO"),
      keytype = "SYMBOL"
    ) %>%
      as_tibble() %>%
      filter(!is.na(SYMBOL), !is.na(GENENAME)) %>%
      group_by(SYMBOL) %>%
      summarise(
        gene_name = dplyr::first(GENENAME),  
        go_terms = paste(unique(na.omit(GO)), collapse = "; "),
        n_go_terms = n_distinct(GO),
        .groups = "drop"
      )
  }, error = function(e) {
    cat("Warning: GO annotation failed:", e$message, "\n")
    return(tibble(
      SYMBOL = character(),
      gene_name = character(),
      go_terms = character(),
      n_go_terms = integer()
    ))
  })
  
  cat("Retrieved GO terms for", nrow(go_terms), "genes\n\n")
  
  # Combine annotations
  cat("Combining annotations...\n")
  annotated_scores <- gene_scores %>%
    left_join(gene_info, by = c("gene" = "hgnc_symbol")) %>%
    left_join(go_terms, by = c("gene" = "SYMBOL")) %>%
    mutate(
      chromosome_name = if_else(!is.na(chromosome_name), 
                                paste0("chr", chromosome_name), 
                                NA_character_),
      gene_length = if_else(!is.na(start_position) & !is.na(end_position),
                            end_position - start_position,
                            NA_integer_)
    )
  
  cat("Annotation complete!\n")
  cat("  Genes with Ensembl data:", sum(!is.na(annotated_scores$chromosome_name)), "\n")
  cat("  Genes with GO terms:", sum(!is.na(annotated_scores$gene_name)), "\n\n")
  
  return(annotated_scores)
}

# Run annotation
cat("Starting gene annotation...\n")
annotated_scores <- annotate_genes(integrated_scores)

# Save
write_csv(annotated_scores, here("data", "processed", "annotated_gene_scores.csv"))
saveRDS(annotated_scores, here("data", "processed", "annotated_gene_scores.rds"))

cat("✓ Annotated scores saved!\n\n")

# Show sample
cat("Sample of annotated data:\n")
print(annotated_scores %>% 
        dplyr::select(gene, gene_name, chromosome_name, gene_length, priority_score) %>%
        head(10))


# Function: Categorize genes by function ----


categorize_genes <- function(annotated_scores) {
  
  # Load the comprehensive gene database
  gene_db <- read_csv(
    here("data", "processed", "comprehensive_gene_database_top15.csv"),
    show_col_types = FALSE
  )
  
  # Define additional gene categories for genes marked as "Other" in database
  oncogenes <- c(
    "KRAS", "NRAS", "HRAS", "BRAF", "EGFR", "ERBB2", "ERBB3", "MET", "ALK",
    "ROS1", "RET", "PIK3CA", "PIK3CB", "MYC", "MYCN", "CCND1", "CCND2", "MDM2",
    "FGFR1", "FGFR2", "FGFR3", "KIT", "PDGFRA", "JAK2", "AKT1", "MTOR",
    "IDH1", "IDH2", "ESR1", "AR", "NOTCH1", "FLT3", "ABL1", "NTRK1"
  )
  
  tumor_suppressors <- c(
    "TP53", "PTEN", "RB1", "APC", "VHL", "NF1", "NF2", "CDKN2A", "CDKN2B",
    "STK11", "SMAD4", "BRCA1", "BRCA2", "ATM", "CHEK2", "PALB2", "TSC1",
    "TSC2", "FBXW7", "ARID1A", "BAP1", "SETD2", "PBRM1"
  )
  
  dna_repair <- c(
    "BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2", "PALB2", "RAD51",
    "RAD51B", "RAD51C", "RAD51D", "BRIP1", "BARD1", "MLH1", "MSH2", "MSH3",
    "MSH6", "PMS2", "POLE", "POLD1", "MUTYH", "ERCC2", "ERCC4", "FANCA",
    "FANCC", "NBN"
  )
  
  chromatin_remodeling <- c(
    "ARID1A", "ARID1B", "ARID2", "SMARCA4", "SMARCB1", "PBRM1", "SETD2",
    "KMT2D", "KMT2C", "KMT2A", "KDM6A", "KDM5C", "CREBBP", "EP300", "BAP1",
    "ASXL1", "EZH2", "DNMT3A", "TET2"
  )
  
  signaling <- c(
    "MAP3K1", "MAP2K1", "MAP2K2", "MAPK1", "PTPN11", "SOS1", "IRS1", "IRS2"
  )
  
  cell_cycle <- c(
    "CCND1", "CCND2", "CCND3", "CCNE1", "CDK4", "CDK6", "CDK12", "CDKN2A",
    "CDKN2B", "RB1", "TP53", "MDM2", "MDM4"
  )
  
  # Join with database and categorize
  categorized <- annotated_scores %>%
    left_join(
      gene_db %>% dplyr::select(gene, gene_category, oncokb_level, is_oncogene, is_tumor_suppressor),
      by = "gene"
    ) %>%
    mutate(
      # Use database category, or classify "Other" genes using knowledge
      final_category = case_when(
        !is.na(gene_category) & gene_category != "Other" ~ gene_category,
        gene %in% dna_repair ~ "DNA Repair",
        gene %in% chromatin_remodeling ~ "Chromatin Remodeling",
        gene %in% oncogenes ~ "Oncogene",
        gene %in% tumor_suppressors ~ "Tumor Suppressor",
        gene %in% signaling ~ "Signaling Pathway",
        gene %in% cell_cycle ~ "Cell Cycle",
        TRUE ~ "Other"
      ),
      
      # Actionable if oncokb_level is 1 or 2
      is_actionable = !is.na(oncokb_level) & oncokb_level <= 3,
      
      # Add actionability detail
      actionability_level = case_when(
        oncokb_level == 1 ~ "FDA-recognized",
        oncokb_level == 2 ~ "Standard care",
        oncokb_level == 3 ~ "Investigational",
        TRUE ~ "Not actionable"
      )
    ) %>%
    dplyr::select(-gene_category) %>%
    dplyr::rename(gene_category = final_category)
  
  return(categorized)
}

# Usage
categorized_scores <- categorize_genes(annotated_scores)


write_csv(categorized_scores, here("data", "processed", "categorized_gene_scores.csv"))

cat("Total genes scored:", nrow(categorized_scores), "\n")
cat("Top 10 priority genes:\n")
print(categorized_scores %>% 
        dplyr::select(gene, priority_score, gene_category) %>% head(10))
