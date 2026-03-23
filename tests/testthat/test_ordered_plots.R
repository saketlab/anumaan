# Test Updated Plot Functions with Ordering
library(dplyr)
library(ggplot2)
library(readr)

# Load functions
source("R/visualize.R")

# Load processed data
cat("Loading processed data...\n")
data <- read_csv("output_data/sheet_2024_25_processed.csv", show_col_types = FALSE)

cat("Data loaded:", nrow(data), "rows\n\n")

# Create test directory
dir.create("test_plots_ordered", showWarnings = FALSE)

# Test 1: Bar plot ordered by frequency (descending)
cat("[1] Creating ordered bar plot (Top 10 organisms)...\n")

top_orgs <- data %>%
  count(organism_normalized, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(organism_normalized)

p1 <- plot_bar(
  data = data %>% filter(organism_normalized %in% top_orgs),
  x = "organism_normalized",
  fill = "organism_normalized",
  flip_coords = TRUE,
  title = "Top 10 Organisms (Ordered by Frequency)",
  xlab = "Organism",
  ylab = "Count"
)

ggsave("test_plots_ordered/01_bar_ordered.png", p1, width = 12, height = 8, dpi = 300)
cat("  ✓ Saved: test_plots_ordered/01_bar_ordered.png\n\n")

# Test 2: Proportion plot ordered by frequency
cat("[2] Creating ordered proportion plot (Organisms vs Resistance)...\n")

p2 <- plot_proportion(
  data = data %>% filter(organism_normalized %in% top_orgs, !is.na(antibiotic_value)),
  x = "organism_normalized",
  fill = "antibiotic_value",
  palette = "resistance",
  flip_coords = TRUE,
  title = "Resistance Proportion by Organism (Ordered by Frequency)",
  xlab = "Organism",
  ylab = "Proportion"
)

ggsave("test_plots_ordered/02_proportion_ordered.png", p2, width = 12, height = 8, dpi = 300)
cat("  ✓ Saved: test_plots_ordered/02_proportion_ordered.png\n\n")

# Test 3: Stacked bar plot ordered by frequency
cat("[3] Creating ordered stacked bar plot (Antibiotic Class vs Resistance)...\n")

top_classes <- data %>%
  filter(!is.na(antibiotic_class)) %>%
  count(antibiotic_class, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(antibiotic_class)

p3 <- plot_stacked_bar(
  data = data %>% filter(antibiotic_class %in% top_classes, !is.na(antibiotic_value)),
  x = "antibiotic_class",
  fill = "antibiotic_value",
  palette = "resistance",
  flip_coords = TRUE,
  show_percentages = TRUE,
  title = "Resistance by Antibiotic Class (Ordered by Frequency)",
  xlab = "Antibiotic Class",
  ylab = "Count"
)

ggsave("test_plots_ordered/03_stacked_ordered.png", p3, width = 12, height = 10, dpi = 300)
cat("  ✓ Saved: test_plots_ordered/03_stacked_ordered.png\n\n")

cat("======================================================================\n")
cat("ALL PLOTS NOW ORDERED BY FREQUENCY (DESCENDING)\n")
cat("======================================================================\n")
cat("Check test_plots_ordered/ directory for the new ordered visualizations\n\n")
