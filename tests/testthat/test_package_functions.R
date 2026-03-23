# =============================================================================
# anumaan Package - Function Testing Script
# Testing all package functions on sheet 1 -2024-25.xlsx dataset
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(lubridate)

# Load the package functions directly
cat("Loading anumaan package functions...\n")
source("R/data.R")
source("R/normalize.R")
source("R/classify.R")
source("R/visualize.R")
cat("✓ Package loaded\n\n")

# =============================================================================
# STEP 1: LOAD AND PREPARE DATA
# =============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("STEP 1: LOADING DATA\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Load raw data
raw_data <- read_excel("data/sheet 1 -2024-25.xlsx", sheet = 1, skip = 1)
cat("Loaded:", nrow(raw_data), "rows\n\n")

# Clean column names
data <- raw_data %>%
  rename(
    patient_id = OrderID,
    gender = Gender,
    age = Age,
    district = `District Name`,
    facility = FacilityName,
    date_of_culture = Visitdate,
    specimen_type = `Type of Sample`,
    organism_name = `Organism Isolated`,
    colony_count = `Colony count`
  )

# Transform WIDE to LONG format
antibiotic_cols <- names(raw_data)[12:53]

long_data <- data %>%
  pivot_longer(
    cols = all_of(antibiotic_cols),
    names_to = "antibiotic_name",
    values_to = "antibiotic_value"
  ) %>%
  filter(!is.na(antibiotic_value), antibiotic_value != "", antibiotic_value != "-")

cat("Transformed to long format:", nrow(long_data), "rows\n")
cat("Sample data:\n")
print(head(long_data %>% select(patient_id, organism_name, antibiotic_name, antibiotic_value), 3))
cat("\n\n")


# =============================================================================
# STEP 2: TEST normalize.R FUNCTIONS
# =============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("STEP 2: TESTING normalize.R FUNCTIONS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# TEST 1: normalize_organism()
cat("TEST 1: normalize_organism()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    test_data <- normalize_organism(
      data = long_data,
      organism_col = "organism_name",
      add_organism_group = TRUE
    )
    cat("✓ SUCCESS\n")
    cat("  Added columns:", setdiff(names(test_data), names(long_data)), "\n")
    cat("  Unique organisms:", n_distinct(test_data$organism_normalized, na.rm = TRUE), "\n")
    cat("  Sample normalized organisms:\n")
    print(head(test_data %>% distinct(organism_name, organism_normalized, organism_group), 5))
    long_data <- test_data # Update for next test
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 2: normalize_specimen()
cat("TEST 2: normalize_specimen()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    test_data <- normalize_specimen(
      data = long_data,
      specimen_col = "specimen_type",
      add_categories = TRUE
    )
    cat("✓ SUCCESS\n")
    cat("  Added columns:", setdiff(names(test_data), names(long_data)), "\n")
    cat("  Sterile vs Non-sterile:\n")
    print(table(test_data$sterile_classification, useNA = "ifany"))
    long_data <- test_data
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 3: get_contaminant_list()
cat("TEST 3: get_contaminant_list()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    contaminants <- get_contaminant_list(
      syndrome = "Bloodstream infections",
      return_all = FALSE
    )
    cat("✓ SUCCESS\n")
    cat("  Found", length(contaminants$names), "contaminants\n")
    cat("  First 5:", paste(head(contaminants$names, 5), collapse = ", "), "\n")
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 4: is_contaminant()
cat("TEST 4: is_contaminant()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    test_organisms <- c("Staphylococcus epidermidis", "Escherichia coli", "Klebsiella pneumoniae")
    results <- is_contaminant(
      organism_name = test_organisms,
      syndrome = "Bloodstream infections"
    )
    cat("✓ SUCCESS\n")
    cat("  Test results:\n")
    for (i in 1:length(test_organisms)) {
      cat("    ", test_organisms[i], "->", results[i], "\n")
    }
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n\n")


# =============================================================================
# STEP 3: TEST classify.R FUNCTIONS
# =============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("STEP 3: TESTING classify.R FUNCTIONS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# TEST 5: derive_age()
cat("TEST 5: derive_age()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
# Add fake DOB for testing
test_data_age <- long_data %>%
  mutate(DOB = date_of_culture - lubridate::years(age))

tryCatch(
  {
    test_data_age <- derive_age(
      data = test_data_age,
      dob_col = "DOB",
      reference_date_col = "date_of_culture"
    )
    cat("✓ SUCCESS\n")
    cat("  Age summary:\n")
    print(summary(test_data_age$Age))
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 6: assign_age_bins()
cat("TEST 6: assign_age_bins()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    test_data_age <- assign_age_bins(
      data = long_data,
      age_col = "age",
      bins = "GBD_standard"
    )
    cat("✓ SUCCESS\n")
    cat("  Age bin distribution:\n")
    print(head(table(test_data_age$Age_bin, useNA = "ifany"), 10))
    long_data <- test_data_age
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 7: extract_genus()
cat("TEST 7: extract_genus()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    test_data_genus <- extract_genus(
      data = long_data,
      organism_col = "organism_normalized"
    )
    cat("✓ SUCCESS\n")
    cat("  Sample results:\n")
    print(head(test_data_genus %>% distinct(organism_normalized, org_genus), 5))
    long_data <- test_data_genus
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 8: extract_species()
cat("TEST 8: extract_species()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    test_data_species <- extract_species(
      data = long_data,
      organism_col = "organism_normalized"
    )
    cat("✓ SUCCESS\n")
    cat("  Sample results:\n")
    print(head(test_data_species %>% distinct(organism_normalized, org_species), 5))
    long_data <- test_data_species
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 9: classify_antibiotic_class()
cat("TEST 9: classify_antibiotic_class()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    who_table <- readr::read_csv("inst/extdata/WHO_aware_class.csv", show_col_types = FALSE)
    test_data_class <- classify_antibiotic_class(
      data = long_data,
      antibiotic_col = "antibiotic_name",
      who_table = who_table
    )
    cat("✓ SUCCESS\n")
    cat("  Antibiotic classes found:\n")
    print(head(table(test_data_class$antibiotic_class, useNA = "ifany"), 10))
    long_data <- test_data_class
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 10: classify_aware()
cat("TEST 10: classify_aware()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    test_data_aware <- classify_aware(
      data = long_data,
      antibiotic_col = "antibiotic_name",
      who_table = who_table
    )
    cat("✓ SUCCESS\n")
    cat("  AWaRe distribution:\n")
    print(table(test_data_aware$aware_category, useNA = "ifany"))
    long_data <- test_data_aware
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 11: calculate_los()
cat("TEST 11: calculate_los()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
# Add fake admission/discharge dates
test_data_los <- long_data %>%
  mutate(
    admission_date = date_of_culture - lubridate::days(sample(1:30, n(), replace = TRUE)),
    discharge_date = date_of_culture + lubridate::days(sample(1:14, n(), replace = TRUE))
  )

tryCatch(
  {
    test_data_los <- calculate_los(
      data = test_data_los,
      admission_col = "admission_date",
      outcome_col = "discharge_date"
    )
    cat("✓ SUCCESS\n")
    cat("  LOS summary:\n")
    print(summary(test_data_los$Length_of_stay))
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n\n")


# =============================================================================
# STEP 4: TEST visualize.R FUNCTIONS
# =============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("STEP 4: TESTING visualize.R FUNCTIONS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Create output directory
dir.create("test_plots", showWarnings = FALSE)

# TEST 12: plot_bar()
cat("TEST 12: plot_bar()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    p <- plot_bar(
      data = long_data,
      x = "organism_normalized",
      top_n = 10,
      flip_coords = TRUE,
      title = "Top 10 Organisms"
    )
    ggsave("test_plots/test_bar.png", p, width = 10, height = 6)
    cat("✓ SUCCESS - Saved to test_plots/test_bar.png\n")
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 13: plot_stacked_bar()
cat("TEST 13: plot_stacked_bar()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    p <- plot_stacked_bar(
      data = long_data,
      x = "antibiotic_name",
      fill = "antibiotic_value",
      show_percentages = TRUE,
      title = "Resistance Pattern"
    )
    ggsave("test_plots/test_stacked.png", p, width = 12, height = 6)
    cat("✓ SUCCESS - Saved to test_plots/test_stacked.png\n")
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 14: plot_proportion()
cat("TEST 14: plot_proportion()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    p <- plot_proportion(
      data = long_data %>% head(10000),
      x = "antibiotic_name",
      fill = "antibiotic_value",
      palette = "resistance",
      title = "Resistance Proportion"
    )
    ggsave("test_plots/test_proportion.png", p, width = 12, height = 6)
    cat("✓ SUCCESS - Saved to test_plots/test_proportion.png\n")
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n")

# TEST 15: plot_histogram()
cat("TEST 15: plot_histogram()\n")
cat(paste(rep("-", 50), collapse = ""), "\n")
tryCatch(
  {
    p <- plot_histogram(
      data = long_data,
      x = "age",
      bins = 30,
      title = "Age Distribution"
    )
    ggsave("test_plots/test_histogram.png", p, width = 10, height = 6)
    cat("✓ SUCCESS - Saved to test_plots/test_histogram.png\n")
  },
  error = function(e) {
    cat("✗ FAILED:", e$message, "\n")
  }
)
cat("\n\n")


# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("TESTING COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Check the test_plots/ directory for generated visualizations\n")
cat("Review any FAILED tests above and fix the issues\n\n")

# Save processed data sample
readr::write_csv(
  long_data %>% head(1000),
  "test_output_sample.csv"
)
cat("Saved sample output to test_output_sample.csv\n")
