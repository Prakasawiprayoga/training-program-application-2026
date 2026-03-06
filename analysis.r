# ---------------------------------------------------------

# Melbourne Bioinformatics Training Program

# This exercise to assess your familiarity with R and git. Please follow
# the instructions on the README page and link to your repo in your application.
# If you do not link to your repo, your application will be automatically denied.

# Leave all code you used in this R script with comments as appropriate.
# Let us know if you have any questions!


# You can use the resources available on our training website for help:
# Intro to R: https://mbite.org/intro-to-r
# Version Control with Git: https://mbite.org/intro-to-git/

# ----------------------------------------------------------

# Load libraries -------------------
# You may use base R or tidyverse for this exercise

library(tidyverse)
library(here)
library(viridis)

# Load data here ----------------------
# Load each file with a meaningful variable name.
metadata <- read_csv(here("data", "GSE60450_filtered_metadata.csv"))
gene_level <- read_csv(here("data", "GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv"))


# Inspect the data -------------------------

# What are the dimensions of each data set? (How many rows/columns in each?)
# Keep the code here for each file.

## Metadata
metadata_dim <- dim(metadata)

## Expression data
gene_level_dim <- dim(gene_level)

# Prepare/combine the data for plotting ------------------------
# How can you combine this data into one data.frame?

# Pivot gene_level
pivoted_gene_level <- gene_level %>%
  pivot_longer(
    cols = starts_with("gsm"), 
    names_to = "sample", 
    values_to = "expression_level"
  )

# Rename blank column of metadata
metadata <- metadata %>%
  rename(sample = ...1)

# Combine data
combined_data <- pivoted_gene_level %>%
  left_join(metadata, by = "sample")


# Plot the data --------------------------
## Plot the expression by cell type
## Can use boxplot() or geom_boxplot() in ggplot2


# 1. Preliminary Visualization ------------------------------------------
# Initial observation: The plot appears 'flat' or 'squashed' at the bottom 
# because many genes have low or zero expression.
# Log transformation is used to manage the wide dynamic range (CPM 0 to 100,000+).

cell_type_plot_log <- combined_data %>%
  ggplot(aes(x = immunophenotype, y = log2(expression_level + 1), fill = immunophenotype)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  theme_minimal() +
  labs(y = "Log2 Expression Level")


# 2. Data Filtering -----------------------
# Pre-filtering strategy: Removing biological noise and zero-inflation.
# Criteria: Keep genes with CPM > 1 in at least 6 samples (the size of the smallest sample).
# This ensures we focus on genes consistently expressed within a cell population.

pivoted_gene_level_filtered <- pivoted_gene_level %>%
  group_by(gene_symbol) %>%
  mutate(n_samples_above_1 = sum(expression_level > 1)) %>% 
  filter(n_samples_above_1 >= 6) %>% 
  ungroup()

combined_data_filtered <- pivoted_gene_level_filtered %>%
  left_join(metadata, by = "sample")


# 3. Refined Visualization -----------------------------------------------
# Plot after filtering: The distribution shifts upward as 'noise' genes are removed.

cell_type_plot_filtered <- combined_data_filtered %>%
  ggplot(aes(x = immunophenotype, y = log2(expression_level + 1), fill = immunophenotype)) +
  geom_jitter(color = "gray60", alpha = 0.05, size = 0.3, width = 0.2) + 
  geom_boxplot(outlier.shape = NA, color = "black", size = 0.5, width = 0.6, alpha = 0.8) +
  theme_classic() +
  scale_fill_viridis(discrete = TRUE, option = "D") + 
  
  labs(
    title = "Gene Expression Distribution across Cell Types",
    subtitle = "Filtered: CPM > 1 in at least 6 samples (Log2-transformed CPM-TMM)",
    x = "Cell Type",
    y = expression(paste(Log[2], " Expression Level (CPM-TMM + 1)")),
    fill = "Immunophenotype"
  ) +
  
  theme(
    text = element_text(size = 12), 
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  )



## Save the plot
### Show code for saving the plot with ggsave() or a similar function
ggsave("results/cell_type_expression_plot.png", width = 8, height = 6, dpi = 300)
