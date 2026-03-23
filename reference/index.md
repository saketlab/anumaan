# Package index

## Pipeline

Main preprocessing pipeline and configuration

- [`amr_preprocess()`](https://saketlab.github.io/anumaan/reference/amr_preprocess.md)
  : AMR Data Preprocessing Pipeline
- [`amr_config()`](https://saketlab.github.io/anumaan/reference/amr_config.md)
  : Create AMR Preprocessing Configuration
- [`validate_config()`](https://saketlab.github.io/anumaan/reference/validate_config.md)
  : Validate AMR Configuration
- [`print(`*`<amr_config>`*`)`](https://saketlab.github.io/anumaan/reference/print.amr_config.md)
  : Print AMR Configuration
- [`print(`*`<amr_result>`*`)`](https://saketlab.github.io/anumaan/reference/print.amr_result.md)
  : Print Method for AMR Preprocessing Results
- [`summary(`*`<amr_result>`*`)`](https://saketlab.github.io/anumaan/reference/summary.amr_result.md)
  : Summary Method for AMR Preprocessing Results
- [`print(`*`<amr_preprocessing_report>`*`)`](https://saketlab.github.io/anumaan/reference/print.amr_preprocessing_report.md)
  : Print Preprocessing Report

## Normalization

Standardize column names, organisms, antibiotics, and specimens

- [`standardize_column_names()`](https://saketlab.github.io/anumaan/reference/standardize_column_names.md)
  : Standardize Column Names to Package Convention
- [`normalize_organism()`](https://saketlab.github.io/anumaan/reference/normalize_organism.md)
  : Normalize Organism Names
- [`normalize_antibiotic()`](https://saketlab.github.io/anumaan/reference/normalize_antibiotic.md)
  : Normalize Antibiotic Names
- [`normalize_specimen()`](https://saketlab.github.io/anumaan/reference/normalize_specimen.md)
  : Normalize Specimen/Sample Type
- [`standardize_gender()`](https://saketlab.github.io/anumaan/reference/standardize_gender.md)
  : Standardize Gender Values
- [`standardize_outcome()`](https://saketlab.github.io/anumaan/reference/standardize_outcome.md)
  : Standardize Outcome Values
- [`standardize_susceptibility()`](https://saketlab.github.io/anumaan/reference/standardize_susceptibility.md)
  : Standardize Susceptibility Values
- [`recode_intermediate()`](https://saketlab.github.io/anumaan/reference/recode_intermediate.md)
  : Recode Intermediate (I) Susceptibility Values
- [`clean_antibiotic_values()`](https://saketlab.github.io/anumaan/reference/clean_antibiotic_values.md)
  : Clean Antibiotic Susceptibility Values
- [`parse_dates()`](https://saketlab.github.io/anumaan/reference/parse_dates.md)
  : Parse Dates Safely

## Classification

MDR/XDR classification, organism grouping, and AWaRe

- [`classify_mdr()`](https://saketlab.github.io/anumaan/reference/classify_mdr.md)
  : Classify MDR (Multidrug Resistant)
- [`classify_xdr()`](https://saketlab.github.io/anumaan/reference/classify_xdr.md)
  : Classify XDR (Extensively Drug Resistant)
- [`classify_org_group()`](https://saketlab.github.io/anumaan/reference/classify_org_group.md)
  : Classify Organism Group
- [`classify_aware()`](https://saketlab.github.io/anumaan/reference/classify_aware.md)
  : Classify AWaRe Category
- [`classify_antibiotic_class()`](https://saketlab.github.io/anumaan/reference/classify_antibiotic_class.md)
  : Classify Antibiotic to WHO Class
- [`classify_mortality()`](https://saketlab.github.io/anumaan/reference/classify_mortality.md)
  : Classify Infection-Related Mortality
- [`derive_infection_type()`](https://saketlab.github.io/anumaan/reference/derive_infection_type.md)
  : Derive Infection Type (HAI / CAI) per Patient
- [`derive_infection_type_for_mortality()`](https://saketlab.github.io/anumaan/reference/derive_infection_type_for_mortality.md)
  : Derive Infection Type (HAI / CAI) for Mortality RR Model

## Enrichment

Derive missing variables from existing data

- [`enrich_age()`](https://saketlab.github.io/anumaan/reference/enrich_age.md)
  : Enrich Age
- [`enrich_los()`](https://saketlab.github.io/anumaan/reference/enrich_los.md)
  : Enrich Length of Stay
- [`enrich_infection_type()`](https://saketlab.github.io/anumaan/reference/enrich_infection_type.md)
  : Enrich Infection Type
- [`enrich_hospital_department()`](https://saketlab.github.io/anumaan/reference/enrich_hospital_department.md)
  : Enrich Hospital Department
- [`assign_age_bins()`](https://saketlab.github.io/anumaan/reference/assign_age_bins.md)
  : Assign Age Bins
- [`derive_age()`](https://saketlab.github.io/anumaan/reference/derive_age.md)
  : Derive Age from Date of Birth

## Deduplication and collapsing

Event IDs, deduplication, and data reshaping

- [`create_event_ids()`](https://saketlab.github.io/anumaan/reference/create_event_ids.md)
  : Create Event IDs from Patient-Level Data
- [`deduplicate_events()`](https://saketlab.github.io/anumaan/reference/deduplicate_events.md)
  : Deduplicate Events
- [`flag_contaminants()`](https://saketlab.github.io/anumaan/reference/flag_contaminants.md)
  : Flag Contaminant Organisms
- [`flag_polymicrobial()`](https://saketlab.github.io/anumaan/reference/flag_polymicrobial.md)
  : Flag Polymicrobial Infections (No Specimen Type, No Episode ID)
- [`collapse_to_antibiotic_level()`](https://saketlab.github.io/anumaan/reference/collapse_to_antibiotic_level.md)
  : Collapse to Antibiotic Level (OPTIONAL - Run When YOU Decide)
- [`collapse_to_class_level()`](https://saketlab.github.io/anumaan/reference/collapse_to_class_level.md)
  : Collapse to Class Level
- [`create_wide_format()`](https://saketlab.github.io/anumaan/reference/create_wide_format.md)
  : Create Wide Format Dataset
- [`pivot_wide_to_long()`](https://saketlab.github.io/anumaan/reference/pivot_wide_to_long.md)
  : Convert Wide Format to Long Format
- [`create_resistance_profile()`](https://saketlab.github.io/anumaan/reference/create_resistance_profile.md)
  : Create Resistance Profile
- [`compute_polymicrobial_weight()`](https://saketlab.github.io/anumaan/reference/compute_polymicrobial_weight.md)
  : Compute Polymicrobial Weights

## Resistance profiles

Marginal resistance, coresistance, and profile computation

- [`compute_R_kd_fatal()`](https://saketlab.github.io/anumaan/reference/compute_R_kd_fatal.md)
  : Compute Fatal Prevalence of Resistance (R_kd)
- [`compute_base_yll_from_dl()`](https://saketlab.github.io/anumaan/reference/compute_base_yll_from_dl.md)
  : Compute Base YLL Block from a Scalar D_L (Steps 1-3)
- [`compute_fatal_resistance_prevalence()`](https://saketlab.github.io/anumaan/reference/compute_fatal_resistance_prevalence.md)
  : Compute Per-Profile Fatal Prevalence of Resistance R_Kdelta\* (Eq.
  13)
- [`compute_fraction_associated()`](https://saketlab.github.io/anumaan/reference/compute_fraction_associated.md)
  : Compute Associated-Burden Fractions per Resistance Profile
- [`compute_hospital_daly()`](https://saketlab.github.io/anumaan/reference/compute_hospital_daly.md)
  : Compute Hospital-Level DALY Breakdown
- [`compute_marginal_resistance()`](https://saketlab.github.io/anumaan/reference/compute_marginal_resistance.md)
  : Compute Marginal Resistance per Pathogen and Antibiotic Class
- [`compute_p0()`](https://saketlab.github.io/anumaan/reference/compute_p0.md)
  : Compute Baseline Mortality Rate Among Fully Susceptible Patients
  (p0)
- [`compute_paf_los()`](https://saketlab.github.io/anumaan/reference/compute_paf_los.md)
  : Compute LOS Population Attributable Fraction per Resistance Profile
- [`compute_paf_rr_mortality()`](https://saketlab.github.io/anumaan/reference/compute_paf_rr_mortality.md)
  : Compute Mortality Population Attributable Fraction per Resistance
  Profile
- [`compute_pairwise_coresistance()`](https://saketlab.github.io/anumaan/reference/compute_pairwise_coresistance.md)
  : Compute Pairwise Co-resistance Matrices per Pathogen
- [`compute_patient_los()`](https://saketlab.github.io/anumaan/reference/compute_patient_los.md)
  : Compute Patient-Level Post-Infection LOS
- [`compute_polymicrobial_weight()`](https://saketlab.github.io/anumaan/reference/compute_polymicrobial_weight.md)
  : Compute Polymicrobial Weights
- [`compute_processing_statistics()`](https://saketlab.github.io/anumaan/reference/compute_processing_statistics.md)
  : Compute Processing Statistics
- [`compute_resistance_profiles()`](https://saketlab.github.io/anumaan/reference/compute_resistance_profiles.md)
  : Compute Resistance Profile Probabilities per Pathogen
- [`compute_yld_associated()`](https://saketlab.github.io/anumaan/reference/compute_yld_associated.md)
  : Compute YLDs Associated with Resistance
- [`compute_yld_attributable()`](https://saketlab.github.io/anumaan/reference/compute_yld_attributable.md)
  : Compute YLDs Attributable to Resistance
- [`compute_yll_associated()`](https://saketlab.github.io/anumaan/reference/compute_yll_associated.md)
  : Compute YLL Associated with AMR (Patient-Level, Facility-Direct)
- [`compute_yll_attributable()`](https://saketlab.github.io/anumaan/reference/compute_yll_attributable.md)
  : Compute YLL Attributable to AMR
- [`select_resistance_class()`](https://saketlab.github.io/anumaan/reference/select_resistance_class.md)
  : Select Resistance Class
- [`prioritize_resistance()`](https://saketlab.github.io/anumaan/reference/prioritize_resistance.md)
  : Prioritize Resistance

## Burden estimation (YLL/YLD)

GBD-methodology burden calculations

- [`calculate_CR_L()`](https://saketlab.github.io/anumaan/reference/calculate_CR_L.md)
  : Calculate the CFR adjustment factor (CR_L)
- [`calculate_P_Lk()`](https://saketlab.github.io/anumaan/reference/calculate_P_Lk.md)
  : Calculate Fatal Pathogen Distribution (P\_{LK})
- [`calculate_P_Lk_fatal()`](https://saketlab.github.io/anumaan/reference/calculate_P_Lk_fatal.md)
  : Calculate fatal pathogen distribution (P\_{Lk})
- [`calculate_P_Lk_prime()`](https://saketlab.github.io/anumaan/reference/calculate_P_Lk_prime.md)
  : Calculate non-fatal pathogen distribution (P'\_{Lk})
- [`calculate_Rkd_prime()`](https://saketlab.github.io/anumaan/reference/calculate_Rkd_prime.md)
  : Calculate non-fatal prevalence of resistance (R'\_{k,d})
- [`calculate_YLD()`](https://saketlab.github.io/anumaan/reference/calculate_YLD.md)
  : Calculate YLD per pathogen
- [`calculate_cfr_lk()`](https://saketlab.github.io/anumaan/reference/calculate_cfr_lk.md)
  : Calculate case fatality ratio by syndrome and pathogen (CFR\_{Lk})
- [`calculate_deaths_by_cause()`](https://saketlab.github.io/anumaan/reference/calculate_deaths_by_cause.md)
  : Calculate Deaths by Underlying Cause (D_J)
- [`calculate_distance_matrix()`](https://saketlab.github.io/anumaan/reference/calculate_distance_matrix.md)
  : Calculate Distance Matrix Between Locations
- [`calculate_incidence_L()`](https://saketlab.github.io/anumaan/reference/calculate_incidence_L.md)
  : Calculate syndrome incidence from deaths, CFR, and CR_L
  (formula-based)
- [`calculate_infection_fraction()`](https://saketlab.github.io/anumaan/reference/calculate_infection_fraction.md)
  : Calculate Infection Fraction of Deaths by Cause (S_J)
- [`calculate_los()`](https://saketlab.github.io/anumaan/reference/calculate_los.md)
  : Calculate Length of Stay
- [`calculate_spatial_autocorrelation()`](https://saketlab.github.io/anumaan/reference/calculate_spatial_autocorrelation.md)
  : Calculate Spatial Autocorrelation (Moran's I)
- [`calculate_spatial_metrics()`](https://saketlab.github.io/anumaan/reference/calculate_spatial_metrics.md)
  : Calculate AMR Metrics by Geographic Unit
- [`calculate_syndrome_deaths()`](https://saketlab.github.io/anumaan/reference/calculate_syndrome_deaths.md)
  : Calculate Deaths by Infectious Syndrome (D_L)
- [`calculate_syndrome_fraction()`](https://saketlab.github.io/anumaan/reference/calculate_syndrome_fraction.md)
  : Calculate Infectious Syndrome Fraction (M_LJ)
- [`fit_distributions()`](https://saketlab.github.io/anumaan/reference/fit_distributions.md)
  : Fit Multiple Distributions
- [`fit_los_rr_distribution()`](https://saketlab.github.io/anumaan/reference/fit_los_rr_distribution.md)
  : Estimate Per-Class LOS Relative Risk via Parametric Distribution
  Fitting
- [`fit_los_rr_nima()`](https://saketlab.github.io/anumaan/reference/fit_los_rr_nima.md)
  : Estimate Per-Profile LOS Relative Risk via Distribution Fitting
  (Nima Procedure)
- [`fit_los_rr_poisson()`](https://saketlab.github.io/anumaan/reference/fit_los_rr_poisson.md)
  : Estimate Per-Class LOS Relative Risk via Quasi-Poisson Regression
- [`fit_mortality_rr_logistic()`](https://saketlab.github.io/anumaan/reference/fit_mortality_rr_logistic.md)
  : Estimate Per-Class Mortality Odds Ratio via Mixed-Effects Logistic
  Regression
- [`assign_rr_to_profiles()`](https://saketlab.github.io/anumaan/reference/assign_rr_to_profiles.md)
  : Assign Per-Class LOS RR to Resistance Profiles (Max Rule)
- [`convert_or_to_rr()`](https://saketlab.github.io/anumaan/reference/convert_or_to_rr.md)
  : Convert Odds Ratios to Relative Risks Using Baseline Mortality (p0)
- [`count_incident_cases()`](https://saketlab.github.io/anumaan/reference/count_incident_cases.md)
  : Count incident cases by syndrome from facility data
- [`filter_profiles_to_rr_classes()`](https://saketlab.github.io/anumaan/reference/filter_profiles_to_rr_classes.md)
  : Filter Profiles to Classes with Actual RR Estimates
- [`get_top_pathogens()`](https://saketlab.github.io/anumaan/reference/get_top_pathogens.md)
  : Identify top N pathogens by occurrence
- [`load_rr_reference()`](https://saketlab.github.io/anumaan/reference/load_rr_reference.md)
  : Load RR (Relative Risk) Reference Data
- [`lookup_rr()`](https://saketlab.github.io/anumaan/reference/lookup_rr.md)
  : Lookup Relative Risk Values
- [`add_rr_mappings()`](https://saketlab.github.io/anumaan/reference/add_rr_mappings.md)
  : Add RR Pathogen and Drug Mappings
- [`map_class_to_rr()`](https://saketlab.github.io/anumaan/reference/map_class_to_rr.md)
  : Map Antibiotic Class to RR Drug Category
- [`map_to_rr_pathogen()`](https://saketlab.github.io/anumaan/reference/map_to_rr_pathogen.md)
  : Map Organism to RR Pathogen Category

## LOS modeling

Length-of-stay distribution fitting and comparison

- [`safe_fit()`](https://saketlab.github.io/anumaan/reference/safe_fit.md)
  : Safely Fit a Distribution
- [`fit_distributions()`](https://saketlab.github.io/anumaan/reference/fit_distributions.md)
  : Fit Multiple Distributions
- [`compare_distribution_aic()`](https://saketlab.github.io/anumaan/reference/compare_distribution_aic.md)
  : Compare Distribution Fits by AIC
- [`summarise_distribution()`](https://saketlab.github.io/anumaan/reference/summarise_distribution.md)
  : Summarise a Fitted Distribution
- [`plot_los_distributions()`](https://saketlab.github.io/anumaan/reference/plot_los_distributions.md)
  : Plot LOS Distribution with Fitted Overlays
- [`prepare_los_data()`](https://saketlab.github.io/anumaan/reference/prepare_los_data.md)
  : Prepare LOS Dataset
- [`get_los_by_resistance()`](https://saketlab.github.io/anumaan/reference/get_los_by_resistance.md)
  : Extract LOS Vectors by Resistance Status

## Burden plots

GBD-style and hospital-level burden visualizations

- [`compute_hospital_daly()`](https://saketlab.github.io/anumaan/reference/compute_hospital_daly.md)
  : Compute Hospital-Level DALY Breakdown
- [`plot_gbd_fig4()`](https://saketlab.github.io/anumaan/reference/plot_gbd_fig4.md)
  : GBD Figure 4 – Deaths by Organism Group (Overlapping Bars)
- [`plot_gbd_fig6()`](https://saketlab.github.io/anumaan/reference/plot_gbd_fig6.md)
  : GBD Figure 6 – Attributable Deaths Heatmap (Pathogen x Drug Class)
- [`plot_burden_comparison()`](https://saketlab.github.io/anumaan/reference/plot_burden_comparison.md)
  : Burden Comparison Bar Plot
- [`plot_burden_heatmap()`](https://saketlab.github.io/anumaan/reference/plot_burden_heatmap.md)
  : Burden Heatmap Across Hospitals
- [`plot_daly_lollipop()`](https://saketlab.github.io/anumaan/reference/plot_daly_lollipop.md)
  : DALY Lollipop Plot
- [`plot_resistance_fraction()`](https://saketlab.github.io/anumaan/reference/plot_resistance_fraction.md)
  : Resistance Fraction Stacked Bar

## Visualization

Generic AMR plotting utilities

- [`amr_theme()`](https://saketlab.github.io/anumaan/reference/amr_theme.md)
  : AMR Theme for ggplot2
- [`get_amr_palette()`](https://saketlab.github.io/anumaan/reference/get_amr_palette.md)
  : Get Color Palette
- [`plot_bar()`](https://saketlab.github.io/anumaan/reference/plot_bar.md)
  : Generic Bar Plot
- [`plot_grouped_bar()`](https://saketlab.github.io/anumaan/reference/plot_grouped_bar.md)
  : Grouped (Dodged) Bar Plot
- [`plot_stacked_bar()`](https://saketlab.github.io/anumaan/reference/plot_stacked_bar.md)
  : Stacked Bar Plot
- [`plot_heatmap()`](https://saketlab.github.io/anumaan/reference/plot_heatmap.md)
  : Heatmap
- [`plot_histogram()`](https://saketlab.github.io/anumaan/reference/plot_histogram.md)
  : Histogram
- [`plot_line()`](https://saketlab.github.io/anumaan/reference/plot_line.md)
  : Line Chart
- [`plot_proportion()`](https://saketlab.github.io/anumaan/reference/plot_proportion.md)
  : Proportion Bar Plot
- [`plot_resistance_heatmap()`](https://saketlab.github.io/anumaan/reference/plot_resistance_heatmap.md)
  : Resistance Pattern Heatmap

## Spatial analysis

Spatial metrics and mapping

- [`create_spatial_object()`](https://saketlab.github.io/anumaan/reference/create_spatial_object.md)
  : Create Spatial Object from AMR Data
- [`create_choropleth_map()`](https://saketlab.github.io/anumaan/reference/create_choropleth_map.md)
  : Create Choropleth Map
- [`create_interactive_map()`](https://saketlab.github.io/anumaan/reference/create_interactive_map.md)
  : Create Interactive Leaflet Map
- [`calculate_spatial_autocorrelation()`](https://saketlab.github.io/anumaan/reference/calculate_spatial_autocorrelation.md)
  : Calculate Spatial Autocorrelation (Moran's I)
- [`calculate_spatial_metrics()`](https://saketlab.github.io/anumaan/reference/calculate_spatial_metrics.md)
  : Calculate AMR Metrics by Geographic Unit
- [`calculate_distance_matrix()`](https://saketlab.github.io/anumaan/reference/calculate_distance_matrix.md)
  : Calculate Distance Matrix Between Locations
- [`detect_hotspots()`](https://saketlab.github.io/anumaan/reference/detect_hotspots.md)
  : Detect Spatial Hotspots (Getis-Ord Gi\*)

## Validation and reporting

Data quality checks and reports

- [`validate_data_quality()`](https://saketlab.github.io/anumaan/reference/validate_data_quality.md)
  : Validate Data Quality
- [`validate_required_fields()`](https://saketlab.github.io/anumaan/reference/validate_required_fields.md)
  : Validate Required Fields
- [`check_logical_consistency()`](https://saketlab.github.io/anumaan/reference/check_logical_consistency.md)
  : Check Logical Consistency
- [`summarize_column_mapping()`](https://saketlab.github.io/anumaan/reference/summarize_column_mapping.md)
  : Summarize Column Mapping
- [`summarize_data_quality()`](https://saketlab.github.io/anumaan/reference/summarize_data_quality.md)
  : Summarize Data Quality
- [`summarize_raw_data()`](https://saketlab.github.io/anumaan/reference/summarize_raw_data.md)
  : Summarize Raw Data
- [`summarize_transformations()`](https://saketlab.github.io/anumaan/reference/summarize_transformations.md)
  : Summarize Data Transformations
- [`generate_preprocessing_report()`](https://saketlab.github.io/anumaan/reference/generate_preprocessing_report.md)
  : Generate Preprocessing Report
- [`export_report()`](https://saketlab.github.io/anumaan/reference/export_report.md)
  : Export Report to File

## Utilities

Helper functions

- [`round_to_sum()`](https://saketlab.github.io/anumaan/reference/round_to_sum.md)
  : Largest-Remainder Rounding
- [`shorten_drug_class()`](https://saketlab.github.io/anumaan/reference/shorten_drug_class.md)
  : Shorten Antibiotic Class Names
- [`default_column_mappings`](https://saketlab.github.io/anumaan/reference/default_column_mappings.md)
  : Default Column Name Mappings
- [`get_age_bins()`](https://saketlab.github.io/anumaan/reference/get_age_bins.md)
  : Standard Age Bins (GBD Compatible)
- [`get_antimicrobial_categories()`](https://saketlab.github.io/anumaan/reference/get_antimicrobial_categories.md)
  : Get Antimicrobial Categories for MDR/XDR Classification
- [`get_beta_lactam_hierarchy()`](https://saketlab.github.io/anumaan/reference/get_beta_lactam_hierarchy.md)
  : Beta-Lactam Class Hierarchy
- [`get_contaminant_list()`](https://saketlab.github.io/anumaan/reference/get_contaminant_list.md)
  : Get Contaminant List from Reference File
- [`get_magiorakos_thresholds()`](https://saketlab.github.io/anumaan/reference/get_magiorakos_thresholds.md)
  : Get Magiorakos MDR/XDR Thresholds
- [`is_contaminant()`](https://saketlab.github.io/anumaan/reference/is_contaminant.md)
  : Check if Organism is a Contaminant
- [`groom_optional_columns()`](https://saketlab.github.io/anumaan/reference/groom_optional_columns.md)
  : Groom Optional Columns
- [`remove_duplicate_rows()`](https://saketlab.github.io/anumaan/reference/remove_duplicate_rows.md)
  : Remove Duplicate Rows
- [`extract_genus()`](https://saketlab.github.io/anumaan/reference/extract_genus.md)
  : Extract Genus from Organism Name
- [`extract_species()`](https://saketlab.github.io/anumaan/reference/extract_species.md)
  : Extract Species from Organism Name
