# prep_clean_and_standardize.R
# MDR/XDR reference lookup tables (Magiorakos 2012 criteria).
#
# All standardization functions have been moved to dedicated files:
#   find_extdata_file          -> prep_derivation.R
#   prep_standardize_column_names -> prep_schema.R
#   prep_standardize_organisms -> prep_standardize_organisms.R
#   prep_standardize_antibiotics -> prep_standardize_antibiotics.R
#   prep_standardize_specimens -> prep_standardize_specimens.R
#   prep_standardize_sex/outcome/infection_type -> prep_types.R


#' Get Beta-Lactam Antibiotic Hierarchy
#'
#' Returns the beta-lactam antibiotic class hierarchy used for MDR/XDR
#' classification (Magiorakos 2012 criteria), ordered from broadest
#' to narrowest spectrum.
#'
#' @return Character vector of antibiotic class names in hierarchy order.
#' @export
get_beta_lactam_hierarchy <- function() {
  c(
    "Carbapenems",
    "Fourth-generation-cephalosporins",
    "Third-generation-cephalosporins",
    "Beta-lactam/beta-lactamase-inhibitor_anti-pseudomonal",
    "Beta-lactam/beta-lactamase-inhibitor",
    "Aminopenicillins",
    "Penicillins"
  )
}


#' Get Magiorakos MDR/XDR Thresholds by Organism Group
#'
#' Returns the MDR and XDR category thresholds for each organism group as
#' defined by Magiorakos et al. (2012). Used by \code{prep_classify_mdr()} and
#' \code{prep_classify_xdr()}.
#'
#' @return A tibble with columns: \code{organism_group}, \code{mdr_threshold},
#'   \code{xdr_threshold}, \code{total_categories}.
#' @export
get_magiorakos_thresholds <- function() {
  tibble::tribble(
    ~organism_group,           ~mdr_threshold, ~xdr_threshold, ~total_categories,
    "Enterobacterales",         3,              "all_but_2",     9,
    "Pseudomonas aeruginosa",   3,              "all_but_2",    10,
    "Acinetobacter spp",        3,              "all_but_2",     9,
    "Staphylococcus aureus",    3,              "all_but_2",     9,
    "Enterococcus spp",         3,              "all_but_2",     7,
    "Streptococcus pneumoniae", 3,              "all_but_2",     5
  )
}


#' Get Antimicrobial Categories for MDR/XDR Classification
#'
#' Returns the antimicrobial categories used for Magiorakos MDR/XDR
#' classification for a given organism group.
#'
#' @param organism_group Character. Organism group name. Default
#'   \code{"Enterobacterales"}.
#' @return Character vector of antimicrobial category names.
#' @export
get_antimicrobial_categories <- function(organism_group = "Enterobacterales") {
  categories <- list(
    Enterobacterales = c(
      "Aminoglycosides",
      "Carbapenems",
      "Cephalosporins (3rd gen)",
      "Cephalosporins (4th gen)",
      "Fluoroquinolones",
      "Monobactams",
      "Penicillins + beta-lactamase inhibitors",
      "Polymyxins",
      "Tigecycline"
    ),
    "Pseudomonas aeruginosa" = c(
      "Aminoglycosides",
      "Antipseudomonal carbapenems",
      "Antipseudomonal cephalosporins",
      "Antipseudomonal fluoroquinolones",
      "Antipseudomonal penicillins + beta-lactamase inhibitors",
      "Monobactams",
      "Phosphonic acids",
      "Polymyxins",
      "Ceftazidime-avibactam",
      "Ceftolozane-tazobactam"
    ),
    "Acinetobacter spp" = c(
      "Aminoglycosides",
      "Carbapenems",
      "Cephalosporins (extended-spectrum)",
      "Fluoroquinolones",
      "Penicillins + beta-lactamase inhibitors",
      "Polymyxins",
      "Sulbactam",
      "Tetracyclines",
      "Tigecycline"
    ),
    "Staphylococcus aureus" = c(
      "Aminoglycosides",
      "Fluoroquinolones",
      "Folate pathway inhibitors",
      "Fusidic acid",
      "Glycopeptides",
      "Linezolid",
      "Mupirocin",
      "Rifampicin",
      "Tetracyclines"
    )
  )

  if (organism_group %in% names(categories)) {
    return(categories[[organism_group]])
  } else {
    warning(sprintf(
      "No category list for '%s', returning Enterobacterales default.",
      organism_group
    ))
    return(categories$Enterobacterales)
  }
}
