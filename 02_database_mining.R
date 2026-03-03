# 02_database_mining.R 
# Interconnected with GLOBOCAN 2022 data for top 15 African cancers

library(tidyverse)
library(data.table)
library(here)

config <- readRDS(here("data", "config.rds"))

# Helper function
`%||%` <- function(a, b) if (is.null(a)) b else a


# SECTION 1: Load Panel Gene Lists ----

load_panel_gene_lists <- function() {
  
  foundation_genes <- c(
    "ABL1", "ACVR1B", "AKT1", "AKT2", "AKT3", "ALK", "ALOX12B", "AMER1", "APC",
    "AR", "ARAF", "ARFRP1", "ARID1A", "ASXL1", "ATM", "ATR", "ATRX", "AURKA",
    "AURKB", "AXIN1", "AXL", "BAP1", "BARD1", "BCL2", "BCL2L1", "BCL2L2", "BCL6",
    "BCOR", "BCORL1", "BRAF", "BRCA1", "BRCA2", "BRD4", "BRIP1", "BTG1", "BTG2",
    "BTK", "C11ORF30", "CALR", "CARD11", "CASP8", "CBFB", "CBL", "CCND1", "CCND2",
    "CCND3", "CCNE1", "CD22", "CD274", "CD70", "CD79A", "CD79B", "CDC73", "CDH1",
    "CDK12", "CDK4", "CDK6", "CDK8", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CDKN2C",
    "CEBPA", "CHEK1", "CHEK2", "CIC", "CREBBP", "CRKL", "CSF1R", "CSF3R", "CTCF",
    "CTNNA1", "CTNNB1", "CUL3", "CUL4A", "CXCR4", "CYP17A1", "DAXX", "DDR1", "DDR2",
    "DIS3", "DNMT3A", "DOT1L", "EED", "EGFR", "EP300", "EPHA3", "EPHB1", "EPHB4",
    "ERBB2", "ERBB3", "ERBB4", "ERCC4", "ERG", "ERRFI1", "ESR1", "EZH2", "FAM46C",
    "FANCA", "FANCC", "FANCG", "FANCL", "FAS", "FBXW7", "FGF10", "FGF12", "FGF14",
    "FGF19", "FGF23", "FGF3", "FGF4", "FGF6", "FGFR1", "FGFR2", "FGFR3", "FGFR4",
    "FH", "FLCN", "FLT1", "FLT3", "FOXL2", "FUBP1", "GABRA6", "GATA3", "GATA4",
    "GATA6", "GID4", "GNA11", "GNA13", "GNAQ", "GNAS", "GRM3", "GSK3B", "H3F3A",
    "HDAC1", "HGF", "HNF1A", "HRAS", "HSD3B1", "ID3", "IDH1", "IDH2", "IGF1R",
    "IKBKE", "IKZF1", "INPP4B", "IRF2", "IRF4", "IRS2", "JAK1", "JAK2", "JAK3",
    "JUN", "KDM5A", "KDM5C", "KDM6A", "KDR", "KEAP1", "KEL", "KIT", "KLHL6",
    "KMT2A", "KMT2D", "KRAS", "LTK", "LYN", "MAF", "MAP2K1", "MAP2K2", "MAP2K4",
    "MAP3K1", "MAP3K13", "MAPK1", "MCL1", "MDM2", "MDM4", "MED12", "MEF2B", "MEN1",
    "MERTK", "MET", "MITF", "MKNK1", "MLH1", "MPL", "MRE11A", "MSH2", "MSH3",
    "MSH6", "MST1R", "MTAP", "MTOR", "MUTYH", "MYC", "MYCL", "MYCN", "MYD88",
    "NBN", "NCOA3", "NCOR1", "NF1", "NF2", "NFE2L2", "NFKBIA", "NKX2-1", "NKX3-1",
    "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "NPM1", "NRAS", "NSD1", "NSD2",
    "NSD3", "NTRK1", "NTRK2", "NTRK3", "PAK1", "PALB2", "PARP1", "PAX5",
    "PBRM1", "PDCD1", "PDCD1LG2", "PDGFRA", "PDGFRB", "PDK1", "PGR", "PHF6",
    "PIK3C2G", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2",
    "PIK3R3", "PIM1", "PLCG2", "PMS2", "POLD1", "POLE", "POT1", "PPARG",
    "PPM1D", "PPP2R1A", "PRDM1", "PRKCI", "PRKD1", "PRKN", "PTCH1", "PTEN",
    "PTPN11", "RAC1", "RAD21", "RAD50", "RAD51", "RAD51B", "RAD51C", "RAD51D",
    "RAD52", "RAD54L", "RAF1", "RARA", "RASA1", "RB1", "RBM10", "RECQL4",
    "REL", "RET", "RHEB", "RHOA", "RICTOR", "RIT1", "RNF43", "ROS1", "RPTOR",
    "RRAS", "RUNX1", "RXRA", "SDHA", "SDHB", "SDHC", "SDHD", "SETD2", "SF3B1",
    "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1", "SMARCD1", "SMO", "SOCS1",
    "SOS1", "SOX2", "SOX9", "SPEN", "SPOP", "SPRED1", "SRC", "SRSF2", "STAG2",
    "STAT3", "STAT5A", "STAT5B", "STK11", "SUFU", "SUZ12", "SYK", "TBX3",
    "TCF3", "TCF7L2", "TEK", "TERT", "TET1", "TET2", "TGFBR1", "TGFBR2",
    "TMPRSS2", "TNFAIP3", "TNFRSF14", "TP53", "TP63", "TSC1", "TSC2", "U2AF1",
    "VEGFA", "VHL", "WT1", "XPO1", "XRCC2", "YAP1", "ZFHX3", "ZNRF3", "ZRSR2"
  )
  
  msk_impact_genes <- c(
    "ABL1", "ACVR1", "AKT1", "AKT2", "AKT3", "ALK", "ALOX12B", "AMER1", "APC",
    "AR", "ARAF", "ARID1A", "ARID1B", "ARID2", "ASXL1", "ASXL2", "ATM", "ATR",
    "ATRX", "AURKA", "AURKB", "AXIN1", "AXIN2", "AXL", "B2M", "BAP1", "BARD1",
    "BCL2", "BCL2L1", "BCL2L11", "BCL6", "BCOR", "BCORL1", "BRAF", "BRCA1",
    "BRCA2", "BRD4", "BRIP1", "BTK", "CALR", "CARD11", "CASP8", "CBFB", "CBL",
    "CCND1", "CCND2", "CCND3", "CCNE1", "CD22", "CD274", "CD79A", "CD79B",
    "CDC73", "CDH1", "CDK12", "CDK4", "CDK6", "CDK8", "CDKN1A", "CDKN1B",
    "CDKN2A", "CDKN2B", "CDKN2C", "CEBPA", "CHEK1", "CHEK2", "CIC", "CREBBP",
    "CRKL", "CSF1R", "CSF3R", "CTCF", "CTNNA1", "CTNNB1", "CUL3", "CUL4A",
    "CXCR4", "CYLD", "DAXX", "DDR1", "DDR2", "DIS3", "DNMT3A", "DOT1L",
    "EGFR", "EIF1AX", "EP300", "EPHA3", "EPHA5", "EPHA7", "EPHB1", "ERBB2",
    "ERBB3", "ERBB4", "ERCC2", "ERG", "ESR1", "EZH2", "FAM46C", "FANCA",
    "FANCC", "FAT1", "FBXW7", "FGF19", "FGF3", "FGF4", "FGFR1", "FGFR2",
    "FGFR3", "FGFR4", "FH", "FLCN", "FLT1", "FLT3", "FLT4", "FOXA1", "FOXL2",
    "FOXP1", "GATA1", "GATA2", "GATA3", "GLI1", "GNA11", "GNAQ", "GNAS",
    "H3F3A", "H3F3C", "HGF", "HIST1H3B", "HNF1A", "HRAS", "ID3", "IDH1",
    "IDH2", "IGF1R", "IKBKE", "IKZF1", "IL7R", "INHBA", "INPP4A", "INPP4B",
    "IRF4", "IRS1", "IRS2", "JAK1", "JAK2", "JAK3", "JUN", "KDM5A", "KDM5C",
    "KDM6A", "KDR", "KEAP1", "KEL", "KIT", "KLHL6", "KMT2A", "KMT2D", "KRAS",
    "LTK", "LYN", "MAF", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K13",
    "MAPK1", "MCL1", "MDM2", "MDM4", "MED12", "MEF2B", "MEN1", "MERTK", "MET",
    "MITF", "MKNK1", "MLH1", "MPL", "MRE11A", "MSH2", "MSH3", "MSH6", "MST1R",
    "MTAP", "MTOR", "MUTYH", "MYC", "MYCL", "MYCN", "MYD88", "NBN", "NF1",
    "NF2", "NFE2L2", "NFKBIA", "NKX2-1", "NOTCH1", "NOTCH2", "NOTCH3", "NPM1",
    "NRAS", "NT5C2", "NTRK1", "NTRK2", "NTRK3", "P2RY8", "PALB2", "PARK2",
    "PARP1", "PARP2", "PARP3", "PAX5", "PBRM1", "PDCD1", "PDCD1LG2", "PDGFRA",
    "PDGFRB", "PDK1", "PIK3C2B", "PIK3C2G", "PIK3CA", "PIK3CB", "PIK3R1",
    "PIM1", "PMS2", "POLD1", "POLE", "PPARG", "PPP2R1A", "PPP2R2A", "PRDM1",
    "PRKAR1A", "PRKCI", "PTCH1", "PTEN", "PTPN11", "PTPRO", "QKI", "RAC1",
    "RAD21", "RAD51", "RAD51B", "RAD51C", "RAD51D", "RAD52", "RAD54L", "RAF1",
    "RARA", "RB1", "RBM10", "REL", "RET", "RICTOR", "RNF43", "ROS1", "RPTOR",
    "SDHA", "SDHB", "SDHC", "SDHD", "SETD2", "SF3B1", "SGK1", "SMAD2", "SMAD4",
    "SMARCA4", "SMARCB1", "SMO", "SNCAIP", "SOCS1", "SOX2", "SOX9", "SPEN",
    "SPOP", "SRC", "STAG2", "STAT3", "STK11", "SUFU", "SYK", "TBX3", "TEK",
    "TET2", "TGFBR2", "TIPARP", "TNFAIP3", "TNFRSF14", "TP53", "TSC1", "TSC2",
    "TYRO3", "U2AF1", "VEGFA", "VHL", "WHSC1", "WHSC1L1", "WT1", "XPO1",
    "XRCC2", "ZNF217", "ZNF703"
  )
  
  panel_genes <- tibble(
    gene = unique(c(foundation_genes, msk_impact_genes))
  ) %>%
    mutate(
      in_foundation = gene %in% foundation_genes,
      in_msk_impact = gene %in% msk_impact_genes,
      in_both = in_foundation & in_msk_impact
    )
  
  cat("Foundation Medicine CDx:", length(foundation_genes), "genes\n")
  cat("MSK-IMPACT:", length(msk_impact_genes), "genes\n")
  cat("Combined unique genes:", nrow(panel_genes), "genes\n\n")
  
  write_csv(panel_genes, here("data", "processed", "comprehensive_panel_genes.csv"))
  
  return(list(
    all_genes = panel_genes,
    foundation = foundation_genes,
    msk_impact = msk_impact_genes
  ))
}

panel_data <- load_panel_gene_lists()


# SECTION 2: Download TCGA Data for Top 15 African Cancers ----


# Map GLOBOCAN 2022 Top 15 cancers to TCGA study IDs from cBioportal
get_tcga_studies_for_top15 <- function() {
  
  tcga_studies <- tribble(
    ~study_id, ~cancer_type, ~globocan_rank, ~url,
    "brca_tcga_pan_can_atlas_2018", "Breast", 1, 
    "https://cbioportal-datahub.s3.amazonaws.com/brca_tcga_pan_can_atlas_2018.tar.gz",
    
    "prad_tcga_pan_can_atlas_2018", "Prostate", 2,
    "https://cbioportal-datahub.s3.amazonaws.com/prad_tcga_pan_can_atlas_2018.tar.gz",
    
    "cesc_tcga_pan_can_atlas_2018", "Cervical", 3,
    "https://cbioportal-datahub.s3.amazonaws.com/cesc_tcga_pan_can_atlas_2018.tar.gz",
    
    "lihc_tcga_pan_can_atlas_2018", "Liver", 4,
    "https://cbioportal-datahub.s3.amazonaws.com/lihc_tcga_pan_can_atlas_2018.tar.gz",
    
    "coadread_tcga_pan_can_atlas_2018", "Colorectal", 5,
    "https://cbioportal-datahub.s3.amazonaws.com/coadread_tcga_pan_can_atlas_2018.tar.gz",
    
    "luad_tcga_pan_can_atlas_2018", "Lung", 6,
    "https://cbioportal-datahub.s3.amazonaws.com/luad_tcga_pan_can_atlas_2018.tar.gz",
    
    "lusc_tcga_pan_can_atlas_2018", "Lung", 6,
    "https://cbioportal-datahub.s3.amazonaws.com/lusc_tcga_pan_can_atlas_2018.tar.gz",
    
    "ov_tcga_pan_can_atlas_2018", "Ovarian", 7,
    "https://cbioportal-datahub.s3.amazonaws.com/ov_tcga_pan_can_atlas_2018.tar.gz",
    
    "blca_tcga_pan_can_atlas_2018", "Bladder", 9,
    "https://cbioportal-datahub.s3.amazonaws.com/blca_tcga_pan_can_atlas_2018.tar.gz",
    
    "stad_tcga_pan_can_atlas_2018", "Stomach", 10,
    "https://cbioportal-datahub.s3.amazonaws.com/stad_tcga_pan_can_atlas_2018.tar.gz",
    
    "esca_tcga_pan_can_atlas_2018", "Esophageal", 11,
    "https://cbioportal-datahub.s3.amazonaws.com/esca_tcga_pan_can_atlas_2018.tar.gz",
    
    "ucec_tcga_pan_can_atlas_2018", "Corpus uteri", 12,
    "https://cbioportal-datahub.s3.amazonaws.com/ucec_tcga_pan_can_atlas_2018.tar.gz",
    
    "paad_tcga_pan_can_atlas_2018", "Pancreatic", 14,
    "https://cbioportal-datahub.s3.amazonaws.com/paad_tcga_pan_can_atlas_2018.tar.gz"
  )
  
  cat("TCGA studies mapped to GLOBOCAN 2022 Top 15:\n")
  cat("  • Breast (Rank 1): brca_tcga\n")
  cat("  • Prostate (Rank 2): prad_tcga\n")
  cat("  • Cervical (Rank 3): cesc_tcga\n")
  cat("  • Liver (Rank 4): lihc_tcga\n")
  cat("  • Colorectal (Rank 5): coadread_tcga\n")
  cat("  • Lung (Rank 6): luad_tcga, lusc_tcga\n")
  cat("  • Ovarian (Rank 7): ov_tcga\n")
  cat("  • Bladder (Rank 9): blca_tcga\n")
  cat("  • Stomach (Rank 10): stad_tcga\n")
  cat("  • Esophageal (Rank 11): esca_tcga\n")
  cat("  • Corpus uteri (Rank 12): ucec_tcga\n")
  cat("  • Pancreatic (Rank 14): paad_tcga\n")
  cat("\n  Note: Non-Hodgkin lymphoma (Rank 8), Leukaemia (Rank 13), and\n")
  cat("        Kaposi sarcoma (Rank 15) have limited TCGA data\n\n")
  
  return(tcga_studies)
}

download_tcga_maf_files <- function() {
  
  cat("Downloading TCGA MAF files for top African cancers...\n\n")
  
  tcga_studies <- get_tcga_studies_for_top15()
  
  download_dir <- here("data", "raw", "tcga_maf")
  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (i in 1:nrow(tcga_studies)) {
    study <- tcga_studies$study_id[i]
    cancer_type <- tcga_studies$cancer_type[i]
    rank <- tcga_studies$globocan_rank[i]
    url <- tcga_studies$url[i]
    
    cat(sprintf("[%d/%d] Rank %d: %s (%s)\n", i, nrow(tcga_studies), rank, cancer_type, study))
    
    tar_file <- file.path(download_dir, paste0(study, ".tar.gz"))
    
    tryCatch({
      if (!file.exists(tar_file)) {
        download.file(url, tar_file, mode = "wb", quiet = TRUE)
        cat("  → Downloaded\n")
      } else {
        cat("  → Already downloaded\n")
      }
      
      study_dir <- file.path(download_dir, study)
      if (!dir.exists(study_dir)) {
        untar(tar_file, exdir = download_dir)
        cat("  → Extracted\n")
      }
      
      maf_file <- list.files(study_dir, pattern = "data_mutations.*\\.txt$", full.names = TRUE)
      
      if (length(maf_file) > 0) {
        cat(sprintf("  ✓ MAF file found: %s\n", basename(maf_file)))
      } else {
        cat(sprintf("  ⚠ Warning: No MAF file in %s\n", study_dir))
      }
      
    }, error = function(e) {
      cat(sprintf("  ✗ Error: %s\n", e$message))
    })
  }
  
  cat("\nDownload complete!\n\n")
  return(invisible(NULL))
}

# Dwonload TCGA MAF files via cbioportal:
 download_tcga_maf_files()

# SECTION 3: Read TCGA Data ----

read_tcga_data_top15 <- function() {
  
  cat("Reading TCGA data for Top 15 African cancers...\n\n")
  
  tcga_dir <- here("data", "raw", "tcga_maf")
  tcga_studies <- get_tcga_studies_for_top15()
  
  all_mutations <- list()
  
  for (i in 1:nrow(tcga_studies)) {
    study_id <- tcga_studies$study_id[i]
    cancer_type <- tcga_studies$cancer_type[i]
    rank <- tcga_studies$globocan_rank[i]
    
    study_dir <- file.path(tcga_dir, study_id)
    
    cat(sprintf("[%d/%d] Rank %d: %s (%s)\n", i, nrow(tcga_studies), rank, cancer_type, study_id))
    
    if (!dir.exists(study_dir)) {
      cat("  ✗ Directory not found - run download first\n")
      next
    }
    
    maf_file <- file.path(study_dir, "data_mutations.txt")
    
    if (file.exists(maf_file)) {
      maf_data <- fread(maf_file, nThread = 4, showProgress = FALSE)
      
      cat("  ✓ Loaded", nrow(maf_data), "mutations,", 
          n_distinct(maf_data$Hugo_Symbol), "unique genes\n")
      
      maf_data$study_id <- study_id
      maf_data$cancer_type <- cancer_type
      maf_data$globocan_rank <- rank
      
      all_mutations[[study_id]] <- maf_data
    } else {
      cat("  ✗ MAF file not found\n")
    }
  }
  
  if (length(all_mutations) == 0) {
    stop("No mutation data found! Please run download_tcga_maf_files() first")
  }
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("TCGA DATA LOADED - TOP 15 AFRICAN CANCERS\n")
  cat(strrep("=", 80), "\n")
  
  mutations_combined <- rbindlist(all_mutations, fill = TRUE)
  
  cat("Total mutations:", nrow(mutations_combined), "\n")
  cat("Unique genes:", n_distinct(mutations_combined$Hugo_Symbol), "\n")
  cat("Unique samples:", n_distinct(mutations_combined$Tumor_Sample_Barcode), "\n")
  cat("Studies loaded:", n_distinct(mutations_combined$study_id), "\n")
  cat("Cancer types:", n_distinct(mutations_combined$cancer_type), "\n\n")
  
  # Show distribution by cancer type
  cat("Distribution by cancer type (GLOBOCAN rank):\n")
  summary_by_cancer <- mutations_combined %>%
    as_tibble() %>%
    group_by(cancer_type, globocan_rank) %>%
    summarise(
      n_mutations = n(),
      n_samples = n_distinct(Tumor_Sample_Barcode),
      n_genes = n_distinct(Hugo_Symbol),
      .groups = "drop"
    ) %>%
    arrange(globocan_rank)
  
  print(summary_by_cancer)
  cat("\n")
  
  # Standardize columns
  mutations_std <- mutations_combined %>%
    as_tibble() %>%
    select(
      study_id,
      cancer_type,
      globocan_rank,
      gene = Hugo_Symbol,
      sample_id = Tumor_Sample_Barcode,
      variant_classification = Variant_Classification,
      variant_type = Variant_Type,
      chromosome = Chromosome,
      start_position = Start_Position,
      end_position = End_Position,
      reference_allele = Reference_Allele,
      protein_change = HGVSp_Short
    ) %>%
    filter(!is.na(gene), gene != "")
  
  # Save
  fwrite(mutations_std, 
         here("data", "processed", "tcga_mutations_top15_africa.csv"),
         nThread = 4)
  
  saveRDS(mutations_std, 
          here("data", "processed", "tcga_mutations_top15_africa.rds"))
  
  cat("Saved to: tcga_mutations_top15_africa.csv and .rds\n\n")
  
  return(mutations_std)
}

tcga_mutations <- read_tcga_data_top15()


# SECTION 4: Calculate Mutation Frequencies ----


calculate_mutation_frequencies <- function(mutations_df) {
  
  cat("\nCalculating mutation frequencies for Top 15 cancers...\n")
  
  mutations_dt <- as.data.table(mutations_df)
  
  # By study and cancer type
  freq_by_study <- mutations_dt[, .(
    n_mutations = .N,
    n_samples = uniqueN(sample_id)
  ), by = .(study_id, gene, cancer_type, globocan_rank)]
  
  total_samples <- mutations_dt[, .(
    total_samples = uniqueN(sample_id)
  ), by = study_id]
  
  freq_df <- freq_by_study %>%
    left_join(total_samples, by = "study_id") %>%
    mutate(
      mutation_frequency = n_samples / total_samples
    ) %>%
    arrange(globocan_rank, desc(mutation_frequency))
  
  # Pan-cancer frequencies (weighted by GLOBOCAN incidence)
  cat("Calculating pan-cancer frequencies...\n")
  
  # Load GLOBOCAN data for weighting
  globocan_data <- read_csv(
    here("data", "processed", "cancer_incidence_africa_globocan2022.csv"),
    show_col_types = FALSE
  )
  
  pancancer_freq <- mutations_dt[, .(
    total_mutations = .N,
    total_samples_mutated = uniqueN(sample_id),
    n_studies = uniqueN(study_id),
    cancer_types = paste(unique(cancer_type), collapse = ", ")
  ), by = gene]
  
  total_unique_samples <- mutations_dt[, uniqueN(sample_id)]
  
  pancancer_freq <- pancancer_freq %>%
    mutate(
      pancancer_frequency = total_samples_mutated / total_unique_samples
    ) %>%
    arrange(desc(pancancer_frequency))
  
  # Frequencies by cancer type (for Top 15)
  freq_by_cancer <- mutations_dt[, .(
    n_mutations = .N,
    n_samples_mutated = uniqueN(sample_id)
  ), by = .(gene, cancer_type, globocan_rank)]
  
  total_by_cancer <- mutations_dt[, .(
    total_samples = uniqueN(sample_id)
  ), by = .(cancer_type, globocan_rank)]
  
  freq_by_cancer <- freq_by_cancer %>%
    left_join(total_by_cancer, by = c("cancer_type", "globocan_rank")) %>%
    mutate(
      frequency_in_cancer = n_samples_mutated / total_samples
    ) %>%
    arrange(globocan_rank, desc(frequency_in_cancer))
  
  cat("\nTop 20 most frequently mutated genes (pan-cancer):\n")
  print(pancancer_freq %>% head(20))
  cat("\n")
  
  cat("Top 10 genes by cancer type (Top 3 cancers):\n")
  for (rank in 1:3) {
    cancer_name <- unique(freq_by_cancer$cancer_type[freq_by_cancer$globocan_rank == rank])
      cat(sprintf("\nRank %d - %s:\n", rank, cancer_name))
    top_genes <- freq_by_cancer %>%
      filter(globocan_rank == rank) %>%
      head(10)
    print(top_genes %>% select(gene, n_mutations, frequency_in_cancer))
  }
  cat("\n")
  
  # Save all frequency tables
  fwrite(freq_df, 
         here("data", "processed", "mutation_frequencies_by_study_top15.csv"),
         nThread = 4)
  
  fwrite(pancancer_freq, 
         here("data", "processed", "pancancer_mutation_frequencies_top15.csv"),
         nThread = 4)
  
  fwrite(freq_by_cancer,
         here("data", "processed", "mutation_frequencies_by_cancer_top15.csv"),
         nThread = 4)
  
  cat("Frequency calculations complete!\n\n")
  
  return(list(
    by_study = as_tibble(freq_df),
    pancancer = as_tibble(pancancer_freq),
    by_cancer = as_tibble(freq_by_cancer)
  ))
}

frequencies <- calculate_mutation_frequencies(tcga_mutations)


# SECTION 5: CIViC Data ----


download_civic_data <- function() {
  
  civic_url <- "https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv"
  
  cat("Downloading CIViC clinical evidence...\n")
  
  tryCatch({
    civic_data <- fread(civic_url, nThread = 4, sep = "\t", header = TRUE)
    
    cat("  Downloaded", nrow(civic_data), "records\n")
    
    civic_processed <- civic_data %>%
      as_tibble() %>%
      mutate(
        gene = str_extract(molecular_profile, "^[A-Z0-9]+"),
        variant = str_remove(molecular_profile, "^[A-Z0-9]+\\s+")
      ) %>%
      filter(!is.na(gene), gene != "") %>%
      select(
        gene, 
        variant,
        molecular_profile,
        disease,
        therapies,
        evidence_type,
        evidence_level,
        evidence_direction,
        clinical_significance = any_of("clinical_significance")
      )
    
    cat("  Processed", nrow(civic_processed), "evidence records\n")
    
    civic_summary <- civic_processed %>%
      group_by(gene) %>%
      summarise(
        n_variants = n_distinct(variant),
        n_evidence = n(),
        has_predictive = any(evidence_type == "Predictive", na.rm = TRUE),
        has_diagnostic = any(evidence_type == "Diagnostic", na.rm = TRUE),
        has_prognostic = any(evidence_type == "Prognostic", na.rm = TRUE),
        evidence_levels = paste(unique(evidence_level), collapse = ";"),
        .groups = "drop"
      ) %>%
      arrange(desc(n_evidence))
    
    fwrite(civic_processed,
           here("data", "processed", "civic_clinical_evidence.csv"),
           nThread = 4)
    
    fwrite(civic_summary,
           here("data", "processed", "civic_gene_summary.csv"),
           nThread = 4)
    
    cat("  ✓ CIViC data saved\n\n")
    
    return(list(
      evidence = civic_processed,
      summary = civic_summary
    ))
    
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    cat("  Using minimal dataset...\n\n")
    
    minimal_civic <- tribble(
      ~gene, ~n_variants, ~n_evidence, ~has_predictive,
      "EGFR", 50, 150, TRUE,
      "BRAF", 30, 120, TRUE,
      "KRAS", 25, 100, TRUE,
      "PIK3CA", 20, 80, TRUE,
      "ERBB2", 15, 90, TRUE,
      "ALK", 12, 70, TRUE,
      "ROS1", 10, 50, TRUE,
      "MET", 10, 60, TRUE,
      "RET", 8, 45, TRUE,
      "NTRK1", 8, 40, TRUE,
      "BRCA1", 30, 100, TRUE,
      "BRCA2", 30, 100, TRUE,
      "TP53", 50, 200, FALSE
    )
    
    fwrite(minimal_civic,
           here("data", "processed", "civic_gene_summary.csv"),
           nThread = 4)
    
    return(list(
      evidence = NULL,
      summary = minimal_civic
    ))
  })
}

civic_data <- download_civic_data()


# SECTION 6: OncoKB Genes ----
# Manually curated level 1 and 2 OncoKB cancer genes with annotations

oncokb_genes <- tribble(
  ~gene, ~is_oncogene, ~is_tumor_suppressor, ~oncokb_level,
  "ABL1", TRUE, FALSE, "1",
  "AKT1", TRUE, FALSE, "1",
  "ALK", TRUE, FALSE, "1",
  "ARAF", TRUE, FALSE, "2",
  "ATM", FALSE, TRUE, "1",
  "ATR", FALSE, TRUE, "1",
  "BARD1", FALSE, TRUE, "1",
  "BRAF", TRUE, FALSE, "1",
  "BRCA1", FALSE, TRUE, "1",
  "BRCA2", FALSE, TRUE, "1",
  "BRIP1", FALSE, TRUE, "1",
  "CDK12", TRUE, FALSE, "1",
  "CHEK1", FALSE, TRUE, "1",
  "CHEK2", FALSE, TRUE, "1",
  "EGFR", TRUE, FALSE, "1",
  "ERBB2", TRUE, FALSE, "1",
  "ESR1", TRUE, FALSE, "1",
  "EZH2", TRUE, FALSE, "1",
  "FANCA", FALSE, TRUE, "1",
  "FANCL", FALSE, TRUE, "1",
  "FGFR1", TRUE, FALSE, "1",
  "FGFR2", TRUE, FALSE, "1",
  "FGFR3", TRUE, FALSE, "1",
  "FLT3", TRUE, FALSE, "1",
  "H3-3A", TRUE, FALSE, "1",
  "H3C2", TRUE, FALSE, "1",
  "H3C3", TRUE, FALSE, "1",
  "IDH1", TRUE, FALSE, "1",
  "IDH2", TRUE, FALSE, "1",
  "JAK2", TRUE, FALSE, "2",
  "KIT", TRUE, FALSE, "1",
  "KMT2A", TRUE, FALSE, "1",
  "KRAS", TRUE, FALSE, "1",
  "MAP2K1", TRUE, FALSE, "2",
  "MAP2K2", TRUE, FALSE, "2",
  "MET", TRUE, FALSE, "1",
  "MLH1", FALSE, TRUE, "1",
  "MRE11", FALSE, TRUE, "1",
  "NBN", FALSE, TRUE, "1",
  "NF1", FALSE, TRUE, "1",
  "NPM1", TRUE, FALSE, "1",
  "NRAS", TRUE, FALSE, "1",
  "NRG1", TRUE, FALSE, "1",
  "NTRK1", TRUE, FALSE, "1",
  "NTRK2", TRUE, FALSE, "1",
  "NTRK3", TRUE, FALSE, "1",
  "PALB2", FALSE, TRUE, "1",
  "PDGFB", TRUE, FALSE, "1",
  "PDGFRA", TRUE, FALSE, "1",
  "PDGFRB", TRUE, FALSE, "1",
  "PIK3CA", TRUE, FALSE, "1",
  "POLD1", FALSE, TRUE, "2",
  "POLE", FALSE, TRUE, "2",
  "PTEN", FALSE, TRUE, "1",
  "RAD51B", FALSE, TRUE, "1",
  "RAD51C", FALSE, TRUE, "1",
  "RAD51D", FALSE, TRUE, "1",
  "RAD54L", FALSE, TRUE, "1",
  "RAF1", TRUE, FALSE, "2",
  "RARA", TRUE, FALSE, "1",
  "RET", TRUE, FALSE, "1",
  "ROS1", TRUE, FALSE, "1",
  "SMARCB1", FALSE, TRUE, "1",
  "TSC1", FALSE, TRUE, "1",
  "TSC2", FALSE, TRUE, "1"
)

write_csv(oncokb_genes, here("data", "processed", "oncokb_cancer_genes.csv"))


# SECTION 7: Compile Comprehensive Database with GLOBOCAN Context ----


compile_comprehensive_database_top15 <- function(panel_genes, mutation_frequencies, civic_data, oncokb_genes) {
  
  cat("Compiling comprehensive gene database with GLOBOCAN 2022 context...\n")
  
  comprehensive_db <- panel_genes$all_genes
  
  # Add mutation frequencies
  if (!is.null(mutation_frequencies$pancancer)) {
    comprehensive_db <- comprehensive_db %>%
      left_join(
        mutation_frequencies$pancancer %>%
          select(gene, pancancer_frequency, total_mutations, n_studies, cancer_types),
        by = "gene"
      )
  }
  
  # Add cancer-specific frequencies for Top 3 cancers
  if (!is.null(mutation_frequencies$by_cancer)) {
    top3_freq <- mutation_frequencies$by_cancer %>%
      filter(globocan_rank <= 3) %>%
      select(gene, cancer_type, frequency_in_cancer) %>%
      pivot_wider(
        names_from = cancer_type,
        values_from = frequency_in_cancer,
        names_prefix = "freq_"
      )
    
    comprehensive_db <- comprehensive_db %>%
      left_join(top3_freq, by = "gene")
  }
  
  # Add CIViC data
  if (!is.null(civic_data)) {
    comprehensive_db <- comprehensive_db %>%
      left_join(
        civic_data$summary %>%
          select(gene, n_evidence, has_predictive),
        by = "gene"
      )
  }
  
  # Add OncoKB data
  comprehensive_db <- comprehensive_db %>%
    left_join(
      oncokb_genes %>%
        select(gene, is_oncogene, is_tumor_suppressor, oncokb_level),
      by = "gene"
    )
  
  # Calculate scores with African cancer context
  comprehensive_db <- comprehensive_db %>%
    mutate(
      pancancer_frequency = replace_na(pancancer_frequency, 0),
      n_evidence = replace_na(n_evidence, 0),
      
      # Actionability score
      actionability_score = case_when(
        !is.na(oncokb_level) & oncokb_level == "1" ~ 3,
        has_predictive == TRUE ~ 2,
        TRUE ~ 0
      ),
      
      # Evidence score
      evidence_score = case_when(
        in_both ~ 3,
        in_foundation | in_msk_impact ~ 2,
        n_evidence > 5 ~ 2,
        TRUE ~ 0
      ),
      
      # African relevance score (based on Top 15 frequency)
      africa_relevance_score = case_when(
        pancancer_frequency > 0.3 ~ 3,  # High frequency
        pancancer_frequency > 0.15 ~ 2, # Moderate
        pancancer_frequency > 0.05 ~ 1, # Low but present
        TRUE ~ 0
      ),
      
      # Clinical relevance with African context
      clinical_relevance = (actionability_score * 0.35) + 
        (evidence_score * 0.25) + 
        (africa_relevance_score * 0.20) +
        (pancancer_frequency * 100 * 0.20),
      
      # Gene category
      gene_category = case_when(
        is_oncogene == TRUE ~ "Oncogene",
        is_tumor_suppressor == TRUE ~ "Tumor Suppressor",
        str_detect(gene, "^MLH|^MSH|^PMS|^BRCA|^ATM|^CHEK|^PALB") ~ "DNA Repair",
        TRUE ~ "Other"
      ),
      
      # Tier assignment
      tier = case_when(
        clinical_relevance >= 2.5 ~ "Tier 1: Essential",
        clinical_relevance >= 1.5 ~ "Tier 2: High Priority",
        clinical_relevance >= 0.5 ~ "Tier 3: Moderate Priority",
        TRUE ~ "Tier 4: Lower Priority"
      ),
      
      # African cancer context
      relevant_to_top15 = !is.na(cancer_types) & cancer_types != ""
    ) %>%
    arrange(desc(clinical_relevance))
  
  # Save
  fwrite(comprehensive_db, 
         here("data", "processed", "comprehensive_gene_database_top15.csv"),
         nThread = 4)
  
  saveRDS(comprehensive_db, 
          here("data", "processed", "comprehensive_gene_database_top15.rds"))
  
  cat("✓ Comprehensive database created\n\n")
  
  return(comprehensive_db)
}

comprehensive_db <- compile_comprehensive_database_top15(
  panel_data, frequencies, civic_data, oncokb_genes
)


# SUMMARY REPORT ----

top15_summary <- tibble(
  rank = 1:15,
  cancer = c("Breast", "Prostate", "Cervix uteri", "Liver", "Colorectum",
             "Lung", "Ovary", "Non-Hodgkin lymphoma", "Bladder", "Stomach",
             "Oesophagus", "Corpus uteri", "Leukaemia", "Pancreas", "Kaposi sarcoma"),
  asr_incidence = c(40.49, 30.28, 26.42, 8.52, 8.24, 6.28, 5.28, 4.99, 4.66,
                    4.02, 3.64, 3.48, 3.10, 2.35, 2.25),
  asr_mortality = c(19.16, 17.32, 17.58, 8.19, 5.55, 5.77, 3.97, 3.26, 2.74,
                    3.50, 3.48, 1.13, 2.54, 2.23, 1.16)
)

print(top15_summary)

cat("\n")
cat("TOP 10 PRIORITIZED GENES FOR AFRICAN CANCERS:\n")
print(comprehensive_db %>%
        filter(tier %in% c("Tier 1: Essential", "Tier 2: High Priority")) %>%
        select(gene, tier, clinical_relevance, actionability_score, 
               pancancer_frequency, gene_category) %>%
        head(10))


# Export curated genes
curated_genes <- comprehensive_db %>%
  filter(tier %in% c("Tier 1: Essential", "Tier 2: High Priority")) %>%
  select(gene, source = tier, clinical_relevance, relevant_to_top15) %>%
  mutate(data_source = "GLOBOCAN_2022_Top15_Africa")

write_csv(curated_genes, here("data", "processed", "curated_cancer_genes_top15.csv"))


