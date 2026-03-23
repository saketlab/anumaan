# spatial.R
# Spatial analysis functions for AMR data with geocomputation

#' Calculate AMR Metrics by Geographic Unit
#'
#' Aggregates AMR data by district or facility and calculates key metrics
#' including resistance rates, isolate counts, and diversity indices.
#'
#' @param data Data frame with processed AMR data
#' @param geo_level Character. "district" or "facility"
#' @param stratify_by Character vector. Additional stratification variables (e.g., "organism_group")
#' @param min_isolates Numeric. Minimum isolates required to calculate rate. Default 30.
#'
#' @return Data frame with geographic AMR metrics
#' @export
calculate_spatial_metrics <- function(data,
                                      geo_level = "district",
                                      stratify_by = NULL,
                                      min_isolates = 30) {
  if (!geo_level %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", geo_level))
  }

  # Build grouping variables
  group_vars <- c(geo_level, stratify_by)

  # Calculate metrics
  metrics <- data %>%
    dplyr::filter(!is.na(!!rlang::sym(geo_level))) %>%
    dplyr::group_by(!!!rlang::syms(group_vars)) %>%
    dplyr::summarise(
      n_isolates = dplyr::n(),
      n_tested = sum(!is.na(antibiotic_value)),
      n_resistant = sum(antibiotic_value == "R", na.rm = TRUE),
      n_susceptible = sum(antibiotic_value == "S", na.rm = TRUE),
      n_intermediate = sum(antibiotic_value == "I", na.rm = TRUE),

      # Resistance rate
      resistance_rate = ifelse(n_tested >= min_isolates,
        (n_resistant / n_tested) * 100,
        NA_real_
      ),

      # Unique organisms
      n_unique_organisms = dplyr::n_distinct(organism_normalized, na.rm = TRUE),

      # Unique patients
      n_patients = dplyr::n_distinct(patient_id, na.rm = TRUE),

      # Unique facilities (if grouping by district)
      n_facilities = if (geo_level == "district") {
        dplyr::n_distinct(facility, na.rm = TRUE)
      } else {
        NA_integer_
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Sufficient data flag
      sufficient_data = n_tested >= min_isolates
    )

  return(metrics)
}


#' Create Spatial Object from AMR Data
#'
#' Converts AMR metrics to spatial (sf) object using coordinates.
#' Requires geocoded locations or joins with spatial reference data.
#'
#' @param metrics_data Data frame from calculate_spatial_metrics()
#' @param geo_level Character. "district" or "facility"
#' @param coords_data Data frame with columns: name, longitude, latitude
#'
#' @return sf object with spatial geometries
#' @export
create_spatial_object <- function(metrics_data,
                                  geo_level = "district",
                                  coords_data = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for spatial analysis. Install with: install.packages('sf')")
  }

  # If coords_data provided, join with metrics
  if (!is.null(coords_data)) {
    # Ensure coords_data has required columns
    required_cols <- c("name", "longitude", "latitude")
    if (!all(required_cols %in% names(coords_data))) {
      stop(sprintf("coords_data must have columns: %s", paste(required_cols, collapse = ", ")))
    }

    # Join metrics with coordinates
    spatial_data <- metrics_data %>%
      dplyr::left_join(
        coords_data,
        by = setNames("name", geo_level)
      ) %>%
      dplyr::filter(!is.na(longitude), !is.na(latitude))

    # Convert to sf object (point geometry)
    sf_object <- sf::st_as_sf(
      spatial_data,
      coords = c("longitude", "latitude"),
      crs = 4326, # WGS84
      remove = FALSE
    )

    return(sf_object)
  } else {
    warning("No coords_data provided. Returning non-spatial data frame.")
    return(metrics_data)
  }
}


#' Calculate Spatial Autocorrelation (Moran's I)
#'
#' Tests for spatial clustering of AMR rates using Moran's I statistic.
#' Requires sf object with spatial geometries.
#'
#' @param spatial_obj sf object from create_spatial_object()
#' @param variable Character. Variable name to test for autocorrelation
#' @param neighbors_k Integer. Number of nearest neighbors. Default 4.
#'
#' @return List with Moran's I statistic, p-value, and interpretation
#' @export
calculate_spatial_autocorrelation <- function(spatial_obj,
                                              variable = "resistance_rate",
                                              neighbors_k = 4) {
  if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("spdep", quietly = TRUE)) {
    stop("Packages 'sf' and 'spdep' required. Install with: install.packages(c('sf', 'spdep'))")
  }

  # Remove NAs
  spatial_obj <- spatial_obj %>%
    dplyr::filter(!is.na(!!rlang::sym(variable)))

  if (nrow(spatial_obj) < 5) {
    warning("Insufficient data for spatial autocorrelation (n < 5)")
    return(NULL)
  }

  # Create spatial neighbors (k-nearest)
  coords <- sf::st_coordinates(spatial_obj)
  neighbors <- spdep::knearneigh(coords, k = min(neighbors_k, nrow(spatial_obj) - 1))
  neighbors_list <- spdep::knn2nb(neighbors)

  # Create spatial weights
  weights <- spdep::nb2listw(neighbors_list, style = "W", zero.policy = TRUE)

  # Calculate Moran's I
  moran_test <- spdep::moran.test(
    spatial_obj[[variable]],
    weights,
    zero.policy = TRUE
  )

  # Interpretation
  interpretation <- dplyr::case_when(
    moran_test$p.value > 0.05 ~ "No significant spatial autocorrelation",
    moran_test$estimate[1] > 0 ~ "Positive spatial autocorrelation (clustering)",
    moran_test$estimate[1] < 0 ~ "Negative spatial autocorrelation (dispersion)",
    TRUE ~ "Unknown pattern"
  )

  result <- list(
    statistic = moran_test$estimate[1],
    expected = moran_test$estimate[2],
    variance = moran_test$estimate[3],
    p_value = moran_test$p.value,
    interpretation = interpretation,
    significant = moran_test$p.value <= 0.05
  )

  return(result)
}


#' Detect Spatial Hotspots (Getis-Ord Gi*)
#'
#' Identifies statistically significant hotspots and coldspots of high/low AMR.
#'
#' @param spatial_obj sf object from create_spatial_object()
#' @param variable Character. Variable name to analyze
#' @param neighbors_k Integer. Number of nearest neighbors. Default 4.
#' @param alpha Numeric. Significance level. Default 0.05.
#'
#' @return sf object with hotspot classification
#' @export
detect_hotspots <- function(spatial_obj,
                            variable = "resistance_rate",
                            neighbors_k = 4,
                            alpha = 0.05) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' required. Install with: install.packages('spdep')")
  }

  # Remove NAs
  spatial_obj <- spatial_obj %>%
    dplyr::filter(!is.na(!!rlang::sym(variable)))

  if (nrow(spatial_obj) < 5) {
    warning("Insufficient data for hotspot detection (n < 5)")
    return(spatial_obj)
  }

  # Create spatial neighbors
  coords <- sf::st_coordinates(spatial_obj)
  neighbors <- spdep::knearneigh(coords, k = min(neighbors_k, nrow(spatial_obj) - 1))
  neighbors_list <- spdep::knn2nb(neighbors)
  weights <- spdep::nb2listw(neighbors_list, style = "B", zero.policy = TRUE)

  # Calculate Getis-Ord Gi* statistic
  gi_star <- spdep::localG(spatial_obj[[variable]], weights, zero.policy = TRUE)

  # Calculate p-values (two-tailed)
  p_values <- 2 * pnorm(-abs(gi_star))

  # Classify hotspots
  spatial_obj$gi_star <- as.numeric(gi_star)
  spatial_obj$gi_p_value <- p_values
  spatial_obj$hotspot_class <- dplyr::case_when(
    p_values > alpha ~ "Not Significant",
    gi_star > 0 ~ "Hotspot (High)",
    gi_star < 0 ~ "Coldspot (Low)",
    TRUE ~ "Not Significant"
  )

  spatial_obj$hotspot_class <- factor(
    spatial_obj$hotspot_class,
    levels = c("Hotspot (High)", "Not Significant", "Coldspot (Low)")
  )

  return(spatial_obj)
}


#' Create Choropleth Map
#'
#' Creates a thematic map showing AMR metrics by geographic area.
#'
#' @param spatial_obj sf object with metrics
#' @param variable Character. Variable to map
#' @param title Character. Map title
#' @param palette Character. Color palette: "sequential", "diverging", "hotcold". Default "sequential".
#' @param breaks Numeric vector. Custom breaks for classification
#' @param labels Character vector. Labels for breaks
#'
#' @return tmap object
#' @export
create_choropleth_map <- function(spatial_obj,
                                  variable = "resistance_rate",
                                  title = "AMR Resistance Rate by District",
                                  palette = "sequential",
                                  breaks = NULL,
                                  labels = NULL) {
  if (!requireNamespace("tmap", quietly = TRUE)) {
    stop("Package 'tmap' required. Install with: install.packages('tmap')")
  }

  # Select palette
  if (palette == "sequential") {
    colors <- c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")
  } else if (palette == "diverging") {
    colors <- RColorBrewer::brewer.pal(11, "RdYlGn")
    colors <- rev(colors) # Reverse so red = high
  } else if (palette == "hotcold") {
    colors <- c(
      "#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8",
      "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"
    )
  } else {
    colors <- palette # Custom palette
  }

  # Create map
  map <- tmap::tm_shape(spatial_obj) +
    tmap::tm_fill(
      col = variable,
      palette = colors,
      title = variable,
      breaks = breaks,
      labels = labels,
      textNA = "No data",
      colorNA = "grey90"
    ) +
    tmap::tm_borders(col = "white", lwd = 1.5) +
    tmap::tm_layout(
      title = title,
      title.size = 1.2,
      frame = FALSE,
      legend.outside = TRUE
    ) +
    tmap::tm_compass(position = c("right", "top")) +
    tmap::tm_scale_bar(position = c("left", "bottom"))

  return(map)
}


#' Create Interactive Leaflet Map
#'
#' Creates an interactive web map of AMR data.
#'
#' @param spatial_obj sf object with metrics
#' @param variable Character. Variable to map
#' @param popup_vars Character vector. Variables to show in popup
#' @param title Character. Map title
#' @param palette Character. Color palette function name (e.g., "YlOrRd")
#'
#' @return leaflet map object
#' @export
create_interactive_map <- function(spatial_obj,
                                   variable = "resistance_rate",
                                   popup_vars = NULL,
                                   title = "AMR Interactive Map",
                                   palette = "YlOrRd") {
  if (!requireNamespace("leaflet", quietly = TRUE)) {
    stop("Package 'leaflet' required. Install with: install.packages('leaflet')")
  }

  # Create color palette
  pal <- leaflet::colorNumeric(
    palette = palette,
    domain = spatial_obj[[variable]],
    na.color = "grey"
  )

  # Build popup content
  if (is.null(popup_vars)) {
    popup_vars <- c(variable, "n_isolates", "n_tested")
  }

  popup_content <- spatial_obj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(popup_vars)) %>%
    apply(1, function(row) {
      paste(
        paste0("<b>", names(row), ":</b> ", round(as.numeric(row), 2)),
        collapse = "<br/>"
      )
    })

  # Create map
  map <- leaflet::leaflet(spatial_obj) %>%
    leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron) %>%
    leaflet::addCircleMarkers(
      radius = ~ sqrt(n_isolates) / 2,
      color = ~ pal(get(variable)),
      fillOpacity = 0.7,
      stroke = TRUE,
      weight = 1,
      popup = popup_content
    ) %>%
    leaflet::addLegend(
      pal = pal,
      values = ~ get(variable),
      title = variable,
      position = "bottomright"
    ) %>%
    leaflet::addControl(
      html = paste0("<h4>", title, "</h4>"),
      position = "topright"
    )

  return(map)
}


#' Calculate Distance Matrix Between Locations
#'
#' Computes pairwise distances between facilities or districts.
#'
#' @param spatial_obj sf object with spatial geometries
#' @param units Character. "km" or "m". Default "km".
#'
#' @return Distance matrix
#' @export
calculate_distance_matrix <- function(spatial_obj, units = "km") {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' required")
  }

  # Calculate pairwise distances
  dist_matrix <- sf::st_distance(spatial_obj)

  # Convert units
  if (units == "km") {
    dist_matrix <- dist_matrix / 1000
  }

  # Add row/column names
  if ("district" %in% names(spatial_obj)) {
    dimnames(dist_matrix) <- list(spatial_obj$district, spatial_obj$district)
  } else if ("facility" %in% names(spatial_obj)) {
    dimnames(dist_matrix) <- list(spatial_obj$facility, spatial_obj$facility)
  }

  return(as.matrix(dist_matrix))
}
