# Test Customization Parameters for Plot Functions
library(dplyr)
library(ggplot2)
library(readr)

# Load functions
source("R/visualize.R")

# Load processed data
cat("Loading processed data...\n")
data <- read_csv("output_data/sheet_2024_25_processed.csv", show_col_types = FALSE)

# Create test directory
dir.create("test_customization", showWarnings = FALSE)

# Get top 10 organisms
top_orgs <- data %>%
  count(organism_normalized, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(organism_normalized)

cat("\n=== Test 1: Bar Plot with Descending Order (Default) ===\n")
p1 <- plot_bar(
  data = data %>% filter(organism_normalized %in% top_orgs),
  x = "organism_normalized",
  fill = "organism_normalized",
  flip_coords = TRUE,
  title = "Top 10 Organisms - Descending Order (Default)",
  xlab = "Organism",
  ylab = "Count",
  order = "desc" # Explicit descending
)
ggsave("test_customization/01_bar_descending.png", p1, width = 12, height = 8, dpi = 300)
cat("✓ Saved: test_customization/01_bar_descending.png\n")

cat("\n=== Test 2: Bar Plot with Ascending Order ===\n")
p2 <- plot_bar(
  data = data %>% filter(organism_normalized %in% top_orgs),
  x = "organism_normalized",
  fill = "organism_normalized",
  flip_coords = TRUE,
  title = "Top 10 Organisms - Ascending Order",
  xlab = "Organism",
  ylab = "Count",
  order = "asc" # Ascending order
)
ggsave("test_customization/02_bar_ascending.png", p2, width = 12, height = 8, dpi = 300)
cat("✓ Saved: test_customization/02_bar_ascending.png\n")

cat("\n=== Test 3: Bar Plot with No Ordering (Alphabetical) ===\n")
p3 <- plot_bar(
  data = data %>% filter(organism_normalized %in% top_orgs),
  x = "organism_normalized",
  fill = "organism_normalized",
  flip_coords = TRUE,
  title = "Top 10 Organisms - No Ordering",
  xlab = "Organism",
  ylab = "Count",
  order = "none" # No ordering - natural order
)
ggsave("test_customization/03_bar_none.png", p3, width = 12, height = 8, dpi = 300)
cat("✓ Saved: test_customization/03_bar_none.png\n")

cat("\n=== Test 4: Proportion Plot with Custom Text Size ===\n")
p4 <- plot_proportion(
  data = data %>% filter(organism_normalized %in% top_orgs, !is.na(antibiotic_value)),
  x = "organism_normalized",
  fill = "antibiotic_value",
  palette = "resistance",
  flip_coords = TRUE,
  title = "Resistance Proportion - Custom Text Size (2.5)",
  xlab = "Organism",
  ylab = "Proportion",
  text_size = 2.5 # Custom small text
)
ggsave("test_customization/04_proportion_custom_text.png", p4, width = 12, height = 8, dpi = 300)
cat("✓ Saved: test_customization/04_proportion_custom_text.png\n")

cat("\n=== Test 5: Stacked Bar with Custom Bar Width ===\n")
top_classes <- data %>%
  filter(!is.na(antibiotic_class)) %>%
  count(antibiotic_class, sort = TRUE) %>%
  slice_head(n = 8) %>%
  pull(antibiotic_class)

p5 <- plot_stacked_bar(
  data = data %>% filter(antibiotic_class %in% top_classes, !is.na(antibiotic_value)),
  x = "antibiotic_class",
  fill = "antibiotic_value",
  palette = "resistance",
  flip_coords = TRUE,
  show_percentages = TRUE,
  title = "Resistance by Class - Wide Bars (bar_width = 0.9)",
  xlab = "Antibiotic Class",
  ylab = "Count",
  bar_width = 0.9 # Wider bars
)
ggsave("test_customization/05_stacked_wide_bars.png", p5, width = 12, height = 10, dpi = 300)
cat("✓ Saved: test_customization/05_stacked_wide_bars.png\n")

cat("\n=== Test 6: Bar Plot with Custom Text Angle ===\n")
p6 <- plot_bar(
  data = data %>% filter(organism_normalized %in% top_orgs),
  x = "organism_normalized",
  fill = "organism_normalized",
  flip_coords = FALSE, # Vertical bars
  title = "Top 10 Organisms - 90 Degree Text Angle",
  xlab = "Organism",
  ylab = "Count",
  text_angle = 90 # Vertical text
)
ggsave("test_customization/06_bar_vertical_text.png", p6, width = 12, height = 8, dpi = 300)
cat("✓ Saved: test_customization/06_bar_vertical_text.png\n")

cat("\n=== Test 7: Auto-Adjusted Text (Many Categories) ===\n")
# Use more categories to test auto-adjustment
top_20_orgs <- data %>%
  count(organism_normalized, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(organism_normalized)

p7 <- plot_bar(
  data = data %>% filter(organism_normalized %in% top_20_orgs),
  x = "organism_normalized",
  fill = "organism_normalized",
  flip_coords = FALSE,
  title = "Top 20 Organisms - Auto Text Adjustment",
  xlab = "Organism",
  ylab = "Count"
  # text_size = NULL (default), text_angle = NULL (default) - should auto-adjust
)
ggsave("test_customization/07_bar_auto_adjust.png", p7, width = 14, height = 8, dpi = 300)
cat("✓ Saved: test_customization/07_bar_auto_adjust.png\n")

cat("\n======================================================================\n")
cat("ALL CUSTOMIZATION TESTS COMPLETED!\n")
cat("======================================================================\n")
cat("New parameters available:\n")
cat("  - order: 'desc' (default), 'asc', or 'none'\n")
cat("  - text_size: NULL (auto-adjust) or numeric value\n")
cat("  - text_angle: NULL (auto-adjust) or 0-90 degrees\n")
cat("  - bar_width: 0.7 or 0.8 (default) or any value 0-1\n")
cat("\nCheck test_customization/ directory for examples\n\n")
