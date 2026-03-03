library(dplyr)
library(ggplot2)
library(tidyr)


#### Figure 1b ----

tcga_mutation_frequencies <- tcga_mutation_frequencies %>%
  arrange(desc(total_samples)) %>%
  mutate(cancer_type = factor(cancer_type, levels = unique(cancer_type)))

### Filter unique total samples
tcga_frequencies <- tcga_mutation_frequencies %>%
  distinct(cancer_type, .keep_all = TRUE)


### Make the Plot 
ggplot(tcga_frequencies) +
  aes(x = cancer_type, y = total_samples) +
  geom_col(fill = "blue") +
  geom_text(aes(label = total_samples), 
            vjust = -0.5, 
            size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Adds space for labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cancer Type", y = "Total Samples")

#### Figure 1c ----
## FIlter genes to top 20 and arrange by literature count
literature_frequency <- gene_literature_frequency %>%
  arrange(desc(literature_count)) %>%
  slice_head(n = 20) %>%
  mutate(genes = factor(genes, levels = unique(genes)))


## Make the plot
ggplot(literature_frequency) +
  aes(x = genes, y = literature_count) +
  geom_col(fill = "blue") +
  geom_text(aes(label = literature_count), 
            vjust = -0.5, 
            size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Adds space for labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Genes", y = "Literature Count")

#### Figure 2b ----
# Filter genes to top 20 and arrange by priority score
priority_score <- integrated_gene_scores %>%
  arrange(desc(priority_score)) %>%
  slice_head(n = 30) %>%
  mutate(genes = factor(gene, levels = unique(gene)))

## Make the plot
ggplot(priority_score) +
  aes(reorder(gene, -priority_score), y = priority_score) +
  geom_col(fill = "blue") +
  #geom_text(aes(label = priority_score), 
   #         vjust = -0.5, 
   #         size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Adds space for labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Genes", y = "Gene Priority Score")

#### Figure 2c ----
# Filter genes to top 20 and arrange by clinical_relevance
clinical_score <- comprehensive_gene_database_top15 %>%
  arrange(desc(clinical_relevance)) %>%
  slice_head(n = 30) %>%
  mutate(genes = factor(gene, levels = unique(gene)))

## Make the plot
ggplot(clinical_score) +
  aes(reorder(gene, -clinical_relevance), y = clinical_relevance) +
  geom_col(fill = "blue") +
  #geom_text(aes(label = priority_score), 
  #         vjust = -0.5, 
  #         size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Adds space for labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Genes", y = "Gene Clinical Relevance Score")

#### Figure 3a ----

# Find the 80% coverage inflection point
inflection_80 <- panel_size_elbow_analysis %>%
  filter(mean_coverage > 81) %>%
  arrange(panel_size) %>%
  filter(panel_size %% 10 == 0) %>%
  slice_head(n = 1)

if (nrow(inflection_80) == 0) {
  inflection_80 <- panel_size_elbow_analysis %>%
    filter(mean_coverage > 81) %>%
    arrange(panel_size) %>%
    slice_head(n = 1)
}

# Calculate elbow point from marginal benefit (as in your original code)
elbow_analysis_with_marginal <- panel_size_elbow_analysis %>%
  arrange(panel_size) %>%
  mutate(
    marginal_coverage = c(NA, diff(mean_coverage)),
    marginal_coverage_per_gene = marginal_coverage / 10  # Assuming increments of 10
  )

# Find elbow point using your threshold method
elbow_threshold <- max(elbow_analysis_with_marginal$marginal_coverage_per_gene, na.rm = TRUE) * 0.2
elbow_point_data <- elbow_analysis_with_marginal %>%
  filter(marginal_coverage_per_gene < elbow_threshold) %>%
  slice_head(n = 1)

if (nrow(elbow_point_data) == 0) {
  elbow_point_data <- panel_size_elbow_analysis %>%
    filter(panel_size == 80) %>%  # Default as in your code
    mutate(marginal_coverage = NA, marginal_coverage_per_gene = NA)
}

# Determine y-axis maximum
y_max <- max(panel_size_elbow_analysis$mean_coverage, na.rm = TRUE)
y_max <- ceiling(y_max / 20) * 20

# Create the plot
ggplot(panel_size_elbow_analysis) +
  aes(x = panel_size, y = mean_coverage) +
  geom_line(colour = "#112446") +
  geom_point(color = "#112446") +
  
  # Highlight 80% inflection point (in blue)
  geom_point(data = inflection_80,
             color = "blue", size = 4, shape = 21, fill = "white", stroke = 1.5) +
  geom_vline(data = inflection_80,
             aes(xintercept = panel_size),
             linetype = "dashed", color = "blue", alpha = 0.5) +
  
  # Highlight elbow point (in red)
  geom_point(data = elbow_point_data,
             color = "red", size = 4, shape = 24, fill = "white", stroke = 1.5) +
  geom_vline(data = elbow_point_data,
             aes(xintercept = panel_size),
             linetype = "dashed", color = "red", alpha = 0.5) +
  
  # Add horizontal line at 80%
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray", alpha = 0.7) +
  
  # Add annotations
  geom_text(data = inflection_80,
            aes(x = panel_size, y = mean_coverage,
                label = paste0("80% Coverage\n", panel_size, " genes")),
            hjust = -0.1, vjust = -0.5, size = 3, color = "blue") +
  
  geom_text(data = elbow_point_data,
            aes(x = panel_size, y = mean_coverage,
                label = paste0("Elbow Point\n", panel_size, " genes")),
            hjust = -0.1, vjust = 1.5, size = 3, color = "red") +
  
  # Y-axis with intervals of 20
  scale_y_continuous(
    limits = c(0, y_max),
    breaks = seq(0, y_max, by = 20),
    labels = function(x) paste0(x, "%")
  ) +
  
  theme_minimal() +
  labs(title = "Elbow Analysis: Coverage by Panel Size",
       subtitle = paste("Blue: First panel size >80% coverage",
                        "Red: Elbow point (marginal benefit <20% of max)"),
       x = "Panel Size",
       y = "Mean Coverage (%)",
       caption = paste("Elbow threshold:", round(elbow_threshold, 3), 
                       "% per gene increase"))

#### Figure 3b ----

# Create combined plot with coverage and marginal benefit
coverage_plot <- ggplot(panel_size_elbow_analysis) +
  aes(x = panel_size, y = mean_coverage) +
  geom_line(colour = "#112446") +
  geom_point(color = "#112446") +
  geom_point(data = inflection_80, color = "blue", size = 4) +
  geom_point(data = elbow_point_data, color = "red", size = 4) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray") +
  scale_y_continuous(
    limits = c(0, y_max),
    breaks = seq(0, y_max, by = 20),
    labels = function(x) paste0(x, "%"),
    name = "Coverage (%)"
  ) +
  labs(title = "Coverage Elbow Analysis",
       x = "Panel Size") +
  theme_minimal()

# Create marginal benefit plot
marginal_plot <- ggplot(elbow_analysis_with_marginal) +
  aes(x = panel_size, y = marginal_coverage_per_gene) +
  geom_line(color = "darkgreen") +
  geom_point(color = "darkgreen") +
  geom_hline(yintercept = elbow_threshold, 
             linetype = "dashed", color = "red") +
  geom_vline(data = elbow_point_data,
             aes(xintercept = panel_size),
             linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "Marginal Benefit per Additional 10 Genes",
       x = "Panel Size",
       y = "Coverage Increase per Gene (%)") +
  theme_minimal()

# Combine plots
library(patchwork)
coverage_plot / marginal_plot +
  plot_annotation(title = "Complete Elbow Analysis",
                  subtitle = "Top: Coverage curve | Bottom: Marginal benefit with elbow threshold")

#### Figure 3c ----
# Calculate the percentage and keep both values
panel_size_elbow_analysis <- panel_size_elbow_analysis %>%
  mutate(
    n_actionable_pct = (n_actionable / max(n_actionable, na.rm = TRUE)) * 100,
    max_actionable = max(n_actionable, na.rm = TRUE)
  )

# Find inflection point
inflection_point <- panel_size_elbow_analysis %>%
  filter(n_actionable_pct > 80) %>%
  arrange(panel_size) %>%
  filter(panel_size %% 10 == 0) %>%
  slice_head(n = 1)

if (nrow(inflection_point) == 0) {
  inflection_point <- panel_size_elbow_analysis %>%
    filter(n_actionable_pct > 80) %>%
    arrange(panel_size) %>%
    slice_head(n = 1)
}

# Get the maximum value for y-axis scaling
max_raw <- max(panel_size_elbow_analysis$n_actionable, na.rm = TRUE)
y_max_pct <- max(panel_size_elbow_analysis$n_actionable_pct, na.rm = TRUE)
y_max_pct <- ceiling(y_max_pct / 20) * 20

ggplot(panel_size_elbow_analysis) +
  aes(x = panel_size, y = n_actionable_pct) +
  geom_line(colour = "#112446") +
  geom_point(color = "#112446") +
  geom_point(data = inflection_point,
             color = "red", size = 4, shape = 21, fill = "white", stroke = 1.5) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray", alpha = 0.7) +
  geom_vline(data = inflection_point,
             aes(xintercept = panel_size),
             linetype = "dashed", color = "red", alpha = 0.7) +
  # Annotation showing both percentage and raw value
  geom_text(data = inflection_point,
            aes(x = panel_size, y = n_actionable_pct,
                label = paste0("Panel: ", panel_size, 
                               "\n", round(n_actionable_pct, 1), "% of max",
                               "\n(", n_actionable, " / ", max_actionable, ")")),
            hjust = -0.1, vjust = -0.5, size = 3) +
  scale_y_continuous(
    limits = c(0, y_max_pct),
    breaks = seq(0, y_max_pct, by = 20),
    labels = function(x) paste0(x, "%"),
    # Add secondary axis showing raw values
    sec.axis = sec_axis(
      ~ . * max_raw / 100,
      name = "Number of Actionable Variants",
      breaks = scales::pretty_breaks(n = 6)
    )
  ) +
  theme_minimal() +
  labs(title = "Elbow Analysis: Actionable Variants by Panel Size",
       subtitle = paste("Maximum actionable variants:", max_raw),
       x = "Panel Size",
       y = "Actionable Variants (% of maximum)")


#### Figure 3d ----

# Sort data and find top 2 highest efficiency ratios
panel_size_cost_effectiveness <- panel_size_cost_effectiveness %>%
  arrange(desc(efficiency_ratio))

top_2 <- panel_size_cost_effectiveness %>%
  slice_max(efficiency_ratio, n = 2)

ggplot(panel_size_cost_effectiveness) +
  aes(x = panel_size, y = efficiency_ratio) +
  geom_line(colour = "#112446") +
  geom_point(color = "#112446") +
  # Highlight top 2 points in red
  geom_point(data = top_2, 
             color = "red", 
             size = 3) +
  # Add labels for top 2 points
  geom_text(data = top_2,
            aes(label = paste0("Size: ", panel_size, "\nRatio: ", 
                               round(efficiency_ratio, 2))),
            hjust = -0.1, 
            vjust = 0.5,
            size = 3,
            color = "red") +
  # Start y-axis from 0
  scale_y_continuous(limits = c(0, NA), 
                     expand = expansion(mult = c(0, 0.1))) +
  theme_minimal() +
  labs(title = "Cost Efficiency Ratio by Panel Size",
       subtitle = "Red points highlight top 2 highest efficiency ratios",
       x = "Panel Size",
       y = "Cost Efficiency Ratio")

#### Figure 6b ----
# Boxplot by gene
ggplot(simulation_results) +
  aes(x = reorder(gene, mutation_rate, median), y = mutation_rate) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  coord_flip() +  # Flip for better gene name readability
  theme_minimal() +
  labs(title = "Mutation Rate Distribution Across Simulations",
       x = "Gene",
       y = "Mutation Rate")

# Violin plot for better distribution visualization
ggplot(simulation_results) +
  aes(x = reorder(gene, mutation_rate, median), y = mutation_rate) +
  geom_violin(fill = "steelblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Mutation Rate Distribution by Gene",
       subtitle = "Violin shows distribution, box shows quartiles")

### Top genes heatmap
# Calculate average mutation rate per gene
gene_avg <- simulation_results %>%
  group_by(gene) %>%
  summarise(
    avg_rate = mean(mutation_rate),
    avg_mutated = mean(n_mutated)
  ) %>%
  arrange(desc(avg_rate)) %>%
  slice_head(n = 20)  # Top 20 genes

# Filter for top genes
top_genes_data <- simulation_results %>%
  filter(gene %in% gene_avg$gene) %>%
  mutate(gene = factor(gene, levels = gene_avg$gene))

# Heatmap
ggplot(top_genes_data) +
  aes(x = simulation, y = gene, fill = mutation_rate) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Mutation Rate") +
  theme_minimal() +
  labs(title = "Mutation Rate Heatmap (Top 20 Genes)",
       x = "Simulation",
       y = "Gene")

## Faceted Histograms
# Top 12 genes faceted
top_12_genes <- simulation_results %>%
  group_by(gene) %>%
  summarise(avg_rate = mean(mutation_rate)) %>%
  arrange(desc(avg_rate)) %>%
  slice_head(n = 12) %>%
  pull(gene)

simulation_results %>%
  filter(gene %in% top_12_genes) %>%
  ggplot() +
  aes(x = mutation_rate) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(title = "Mutation Rate Distribution for Top 12 Genes",
       x = "Mutation Rate",
       y = "Frequency")

## Simulated Trends line plot
# Plot mutation rate trends across simulations for top genes
top_5_genes <- simulation_results %>%
  group_by(gene) %>%
  summarise(avg_rate = mean(mutation_rate)) %>%
  arrange(desc(avg_rate)) %>%
  slice_head(n = 5) %>%
  pull(gene)

simulation_results %>%
  filter(gene %in% top_5_genes) %>%
  ggplot() +
  aes(x = simulation, y = mutation_rate, color = gene, group = gene) +
  geom_line(alpha = 0.7) +
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(title = "Mutation Rate Trends Across Simulations",
       x = "Simulation",
       y = "Mutation Rate",
       color = "Gene")
