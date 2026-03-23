# Test sp. -> spp. Normalization Fix
library(dplyr)

# Load the normalize function
source("R/normalize.R")

# Create test data with various sp. formats
test_data <- data.frame(
  patient_id = 1:10,
  organism_name = c(
    "Acinetobacter sp.", # Should normalize to "acinetobacter spp."
    "Acinetobacter sp", # Should normalize to "acinetobacter spp."
    "acinetobacter sp.", # Should normalize to "acinetobacter spp."
    "ACINETOBACTER SP", # Should normalize to "acinetobacter spp."
    "Pseudomonas sp.", # Should normalize to "pseudomonas spp."
    "Pseudomonas sp", # Should normalize to "pseudomonas spp."
    "Enterococcus sp.", # Should normalize to "enterococcus spp."
    "Proteus sp", # Should normalize to "proteus spp."
    "Acinetobacter baumannii", # Should stay as species name
    "Enterobacter sp." # Should normalize to "enterobacter spp."
  ),
  stringsAsFactors = FALSE
)

cat("=== TESTING ORGANISM NORMALIZATION: sp. -> spp. ===\n\n")

cat("Input organisms:\n")
print(test_data$organism_name)

cat("\n\nApplying normalize_organism()...\n\n")

# Apply normalization
result <- normalize_organism(test_data, organism_col = "organism_name")

# Show results
cat("Results:\n")
result_df <- result %>%
  select(organism_name, organism_normalized, organism_group) %>%
  arrange(organism_name)

print(as.data.frame(result_df))

# Check if all "sp" variants are now "spp."
cat("\n\n=== VERIFICATION ===\n")

sp_variants <- result %>%
  filter(grepl("sp\\.", organism_name, ignore.case = TRUE) |
    grepl("sp$", organism_name, ignore.case = TRUE))

cat("\nAll 'sp' variants in input:\n")
print(sp_variants$organism_name)

cat("\nNormalized to:\n")
print(unique(sp_variants$organism_normalized))

# Count unique normalizations
cat("\n\nSummary:\n")
cat(sprintf("- Total organisms tested: %d\n", nrow(test_data)))
cat(sprintf("- Unique normalized names: %d\n", length(unique(result$organism_normalized))))
cat(sprintf(
  "- Successfully grouped: %d (%.1f%%)\n",
  sum(!is.na(result$organism_group)),
  100 * sum(!is.na(result$organism_group)) / nrow(result)
))

# Check specific case: Acinetobacter
acinetobacter_cases <- result %>%
  filter(grepl("acinetobacter", organism_name, ignore.case = TRUE))

cat("\n\n=== ACINETOBACTER SPECIFIC TEST ===\n")
cat(sprintf("Input variants: %d\n", nrow(acinetobacter_cases)))
cat(sprintf("Unique normalized: %d\n", length(unique(acinetobacter_cases$organism_normalized))))
cat("\nAcinetobacter normalization:\n")
print(as.data.frame(acinetobacter_cases %>% select(organism_name, organism_normalized)))

# Success check
all_spp <- all(grepl("spp\\.", result$organism_normalized[grepl("sp", result$organism_name, ignore.case = TRUE) &
  !grepl("species|baumannii", result$organism_name, ignore.case = TRUE)]))

if (all_spp) {
  cat("\n\n✓ SUCCESS: All 'sp' variants correctly normalized to 'spp.'\n")
} else {
  cat("\n\n✗ FAILURE: Some 'sp' variants were not normalized correctly\n")
}
