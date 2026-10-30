# profile_monitor_run.R
# Run-state monitor for a profile_fit_joint_bernoulli() run. See
# profile_fit.R's header for how this fits into the pipeline.
#
# profile_fit_joint_bernoulli() calls cmdstanr's $sample() SYNCHRONOUSLY --
# it blocks until every chain finishes, so there is no "live" state to
# observe from inside the SAME R call/session that is running it. What
# this function CAN do, honestly: read whatever cmdstan has ALREADY
# WRITTEN to output_dir -- cmdstan writes each chain's CSV incrementally
# DURING sampling, so calling this from a SEPARATE session/process
# pointed at the same output_dir while a fit is running elsewhere gives a
# genuine (if coarse) progress read. Called after the run completes (the
# common case in a short interactive session), it reports completed-run
# state from the same files, not fabricated live telemetry. It never
# parses a clearly-truncated trailing CSV line as a valid draw.


#' (Internal) Parse One Chain's CSV for Whatever Progress Is on Disk
#' @keywords internal
.profile_read_stan_csv_progress <- function(path) {
  comments_only <- TRUE
  lines <- readLines(path, warn = FALSE)
  is_comment <- grepl("^#", lines)
  comments <- lines[is_comment]
  data_lines <- lines[!is_comment]

  .extract_int <- function(pattern) {
    m <- grep(pattern, comments, value = TRUE)
    if (length(m) == 0L) return(NA_integer_)
    suppressWarnings(as.integer(gsub(".*=\\s*", "", m[1])))
  }
  .extract_logical <- function(pattern) {
    m <- grep(pattern, comments, value = TRUE)
    if (length(m) == 0L) return(NA)
    grepl("true", m[1], ignore.case = TRUE)
  }

  num_warmup <- .extract_int("^#\\s*num_warmup\\s*=")
  num_samples <- .extract_int("^#\\s*num_samples\\s*=")
  save_warmup <- .extract_logical("^#\\s*save_warmup\\s*=")
  completed <- any(grepl("Elapsed Time", comments))

  if (length(data_lines) == 0L) {
    return(list(
      phase = "not_started", last_iteration = 0L, n_divergent = NA_integer_,
      max_observed_treedepth = NA_integer_, num_warmup = num_warmup, num_samples = num_samples,
      completed = FALSE
    ))
  }

  header <- strsplit(data_lines[1], ",")[[1]]
  n_fields <- length(header)
  rows <- data_lines[-1]
  # A trailing row written while cmdstan is mid-flush can be truncated --
  # never trust it as a real draw. Only a row with exactly n_fields
  # comma-separated, all-numeric values counts.
  valid <- vapply(rows, function(r) {
    parts <- strsplit(r, ",", fixed = TRUE)[[1]]
    length(parts) == n_fields && !anyNA(suppressWarnings(as.numeric(parts)))
  }, logical(1))
  n_valid <- sum(valid)

  phase <- if (completed) {
    "completed"
  } else if (isTRUE(save_warmup) && !is.na(num_warmup)) {
    if (n_valid < num_warmup) "warmup" else "sampling"
  } else if (n_valid == 0L) {
    "warmup" # save_warmup = FALSE (cmdstanr's default): zero rows written means still in warmup
  } else {
    "sampling"
  }

  n_divergent <- NA_integer_
  max_observed_treedepth <- NA_integer_
  if (n_valid > 0L) {
    valid_rows <- rows[valid]
    div_col <- match("divergent__", header)
    td_col <- match("treedepth__", header)
    if (!is.na(div_col)) {
      div_vals <- vapply(valid_rows, function(r) as.numeric(strsplit(r, ",", fixed = TRUE)[[1]][div_col]), numeric(1))
      n_divergent <- sum(div_vals > 0, na.rm = TRUE)
    }
    if (!is.na(td_col)) {
      td_vals <- vapply(valid_rows, function(r) as.numeric(strsplit(r, ",", fixed = TRUE)[[1]][td_col]), numeric(1))
      max_observed_treedepth <- max(td_vals, na.rm = TRUE)
    }
  }

  list(
    phase = phase, last_iteration = n_valid, n_divergent = n_divergent,
    max_observed_treedepth = max_observed_treedepth, num_warmup = num_warmup,
    num_samples = num_samples, completed = completed
  )
}


#' Run-State Monitor for a profile_fit_joint_bernoulli() Run
#'
#' Reads whatever cmdstan has ALREADY WRITTEN to \code{output_dir} --
#' never re-runs or re-samples anything. Because
#' \code{\link{profile_fit_joint_bernoulli}} calls cmdstanr's
#' \code{$sample()} synchronously, this function cannot observe a run
#' live from inside the SAME call that started it; point it at the same
#' \code{output_dir} from a SEPARATE session/process while the fit is
#' running elsewhere for a genuine (coarse) mid-run progress read.
#'
#' \code{phase} is inferred, not read directly -- cmdstan's CSV records
#' \code{save_warmup}, so when that is \code{FALSE} (cmdstanr's default,
#' and \code{\link{profile_fit_joint_bernoulli}}'s default), a chain with
#' ZERO data rows written yet is reported \code{"warmup"} (warmup draws
#' are never written to the CSV in that mode, so no rows is the ONLY
#' signal available that warmup is still in progress -- it is
#' indistinguishable from "about to start" with this interface, and this
#' function says so rather than claiming otherwise).
#'
#' @param output_dir Character. The SAME \code{output_dir} passed to
#'   \code{\link{profile_fit_joint_bernoulli}}.
#' @param expected_chains Integer or \code{NULL}. Expected chain count --
#'   when supplied, a chain whose CSV does not exist at all is still
#'   reported as its own row (\code{phase = "not_started"}), rather than
#'   silently absent from the table.
#'
#' @return A data frame, one row per found (or, with \code{expected_chains},
#'   expected) chain: \code{chain_id}, \code{csv_file} (\code{NA} if that
#'   chain's CSV does not exist), \code{phase} (\code{"not_started"},
#'   \code{"warmup"}, \code{"sampling"}, \code{"completed"}),
#'   \code{last_iteration} (count of valid data rows written so far --
#'   NEVER counts a truncated trailing line), \code{num_warmup},
#'   \code{num_samples} (as configured, from the CSV's own comment
#'   header), \code{last_modified} (file mtime -- the closest honest
#'   proxy to "last update time" this interface exposes),
#'   \code{n_divergent}, \code{max_observed_treedepth} (both \code{NA}
#'   before any row is written).
#' @export
profile_monitor_run <- function(output_dir, expected_chains = NULL) {
  if (!dir.exists(output_dir)) {
    stop(sprintf("output_dir '%s' does not exist.", output_dir), call. = FALSE)
  }
  csv_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
  # cmdstanr's per-chain CSV naming convention: "<model>-<timestamp>-<chain_id>-<hash>.csv"
  chain_id_of <- function(f) {
    m <- regmatches(basename(f), regexpr("-([0-9]+)-[0-9a-fA-F]+\\.csv$", basename(f)))
    if (length(m) == 0L || nchar(m) == 0L) return(NA_integer_)
    as.integer(gsub("[^0-9]", "", strsplit(m, "-")[[1]][2]))
  }
  found_ids <- vapply(csv_files, chain_id_of, integer(1))
  found <- data.frame(chain_id = found_ids, csv_file = csv_files, stringsAsFactors = FALSE)
  found <- found[!is.na(found$chain_id), , drop = FALSE]

  all_ids <- if (!is.null(expected_chains)) seq_len(expected_chains) else sort(unique(found$chain_id))

  rows <- lapply(all_ids, function(cid) {
    f <- found$csv_file[found$chain_id == cid]
    if (length(f) == 0L) {
      return(data.frame(
        chain_id = cid, csv_file = NA_character_, phase = "not_started",
        last_iteration = NA_integer_, num_warmup = NA_integer_, num_samples = NA_integer_,
        last_modified = as.POSIXct(NA), n_divergent = NA_integer_, max_observed_treedepth = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }
    prog <- .profile_read_stan_csv_progress(f[1])
    data.frame(
      chain_id = cid, csv_file = f[1], phase = prog$phase,
      last_iteration = prog$last_iteration, num_warmup = prog$num_warmup, num_samples = prog$num_samples,
      last_modified = file.info(f[1])$mtime,
      n_divergent = prog$n_divergent, max_observed_treedepth = prog$max_observed_treedepth,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
