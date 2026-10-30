# profile_fit.R
# Stan fitting and sampler diagnostics for the joint Bernoulli
# resistance-profile model (profile_build_design()'s output). This is the
# generic fitting machinery that belongs in anumaan; anumaan-analysis
# supplies the actual sampler_config (chains, warmup, seed, ...) and the
# run's output_dir -- that choice is a run-config decision, not a package
# default.
#
# profile_fit_joint_bernoulli() saves everything a later review needs: the
# fit object, a parameter summary, a config snapshot (hashed, for
# reproducibility), per-chain sampler diagnostics, and a short run log.
# Diagnostic PLOTS are deliberately a separate function,
# profile_check_sampler() -- a caller who already trusts the numeric
# diagnostics shouldn't be forced to render plots to get a fit back.
#
# Functions:
#   - profile_fit_joint_bernoulli (stage 4: compile + sample + save artifacts)
#   - profile_check_sampler       (stage 5a: pass/fail sampler-validity
#                                  verdict + trace/rank/energy diagnostic
#                                  plots, from a fit already on disk or in
#                                  memory)


#' Fit the Joint Bernoulli Resistance-Profile Model
#'
#' Compiles \code{inst/stan/joint_bernoulli_general.stan} (or a caller-
#' supplied alternative) and samples it against a Stan-ready design from
#' \code{\link{profile_build_design}}, saving every artifact a later
#' review needs: the fit object (\code{fit.rds}), a parameter summary
#' (\code{parameter_summary.csv}), per-chain sampler diagnostics
#' (\code{sampler_diagnostics.rds}), a snapshot of the exact sampler
#' config used, hashed for reproducibility (\code{config_snapshot.rds}),
#' and a short human-readable run log (\code{run_log.txt}). cmdstanr's own
#' per-chain CSV files are also written into \code{output_dir} (via
#' \code{output_dir} forwarded to \code{$sample()}).
#'
#' @param design Output of \code{\link{profile_build_design}} (reads
#'   \code{design$stan_data}; \code{design$manifest} is recorded in the
#'   config snapshot for provenance).
#' @param output_dir Character. Directory for all artifacts (created if
#'   missing).
#' @param sampler_config Named list. Any of \code{chains} (default
#'   \code{4L}), \code{iter_warmup} (\code{1000L}), \code{iter_sampling}
#'   (\code{1000L}), \code{seed} (\code{123L}), \code{parallel_chains}
#'   (default: \code{= chains}), \code{adapt_delta} (\code{0.9}),
#'   \code{max_treedepth} (\code{10L}), \code{threads_per_chain}
#'   (\code{1L} -- \code{> 1} enables \code{reduce_sum} threading in the
#'   Stan model's \code{grainsize} machinery). Unset entries use these
#'   defaults.
#' @param stan_file Character or \code{NULL}. Default: the package's
#'   \code{inst/stan/joint_bernoulli_general.stan}.
#' @param save_warmup Logical, forwarded to \code{$sample()}. Default
#'   \code{FALSE}.
#'
#' @return A list: \code{fit} (the \code{cmdstanr} \code{CmdStanMCMC}
#'   object), \code{summary} (data frame, \code{fit$summary()}),
#'   \code{diagnostics} (list from \code{fit$diagnostic_summary()}),
#'   \code{config} (the resolved sampler config, the design manifest, and
#'   \code{config_hash}), \code{output_dir}.
#' @export
profile_fit_joint_bernoulli <- function(design,
                                        output_dir,
                                        sampler_config = list(),
                                        stan_file = NULL,
                                        save_warmup = FALSE) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required for profile_fit_joint_bernoulli().", call. = FALSE)
  }
  if (!is.list(design) || !"stan_data" %in% names(design)) {
    stop("design must be the list returned by profile_build_design() (needs $stan_data).")
  }

  if (is.null(stan_file)) {
    stan_file <- system.file("stan", "joint_bernoulli_general.stan", package = "anumaan")
  }
  if (!nzchar(stan_file) || !file.exists(stan_file)) {
    stop(sprintf("Stan file not found: %s", stan_file))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  defaults <- list(
    chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L,
    seed = 123L, parallel_chains = NULL, adapt_delta = 0.9,
    max_treedepth = 10L, threads_per_chain = 1L
  )
  cfg <- utils::modifyList(defaults, sampler_config)
  if (is.null(cfg$parallel_chains)) cfg$parallel_chains <- cfg$chains

  use_threads <- cfg$threads_per_chain > 1L
  mod <- cmdstanr::cmdstan_model(
    stan_file,
    cpp_options = if (use_threads) list(stan_threads = TRUE) else list(),
    quiet = TRUE
  )

  timestamp <- as.character(Sys.time())
  config_hash <- tryCatch(
    rlang::hash(list(design$stan_data, cfg)),
    error = function(e) NA_character_
  )
  config_snapshot <- list(
    stan_file = stan_file, sampler_config = cfg,
    design_manifest = design$manifest, timestamp = timestamp,
    config_hash = config_hash
  )
  saveRDS(config_snapshot, file.path(output_dir, "config_snapshot.rds"))

  fit <- mod$sample(
    data = design$stan_data,
    chains = as.integer(cfg$chains),
    parallel_chains = as.integer(cfg$parallel_chains),
    iter_warmup = as.integer(cfg$iter_warmup),
    iter_sampling = as.integer(cfg$iter_sampling),
    adapt_delta = cfg$adapt_delta,
    max_treedepth = as.integer(cfg$max_treedepth),
    seed = as.integer(cfg$seed),
    threads_per_chain = if (use_threads) as.integer(cfg$threads_per_chain) else NULL,
    save_warmup = save_warmup,
    output_dir = output_dir
  )

  fit$save_object(file.path(output_dir, "fit.rds"))

  summary_tbl <- fit$summary()
  utils::write.csv(summary_tbl, file.path(output_dir, "parameter_summary.csv"), row.names = FALSE)

  diag <- fit$diagnostic_summary()
  saveRDS(diag, file.path(output_dir, "sampler_diagnostics.rds"))

  log_lines <- c(
    sprintf("profile_fit_joint_bernoulli run: %s", timestamp),
    sprintf("stan_file: %s", stan_file),
    sprintf("config_hash: %s", config_hash),
    sprintf(
      "N_events=%d D=%d K=%d S=%d P=%d",
      design$stan_data$N_events, design$stan_data$D, design$stan_data$K,
      design$stan_data$S, design$stan_data$P
    ),
    sprintf(
      "chains=%d warmup=%d sampling=%d seed=%d",
      cfg$chains, cfg$iter_warmup, cfg$iter_sampling, cfg$seed
    ),
    sprintf("num_divergent (total): %d", sum(diag$num_divergent)),
    sprintf("num_max_treedepth (total): %d", sum(diag$num_max_treedepth)),
    sprintf("ebfmi (per chain): %s", paste(round(diag$ebfmi, 3), collapse = ", "))
  )
  writeLines(log_lines, file.path(output_dir, "run_log.txt"))

  list(
    fit = fit, summary = summary_tbl, diagnostics = diag,
    config = config_snapshot, output_dir = output_dir
  )
}
