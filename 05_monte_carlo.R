
## title: "APIM meta-analysis: individual-level missing and final pipeline"
## author: "Anonymized"
## date: "15-10-2026"

# =============================================================================
# 05_monte_carlo.R
# Monte Carlo simulation of the item-level missing-data methods for the APIMeM
# meta-analysis (ADEMP; Morris, White & Crowther, 2019). 04 evaluates ONE
# realisation of the missingness per condition; this script separates method
# differences from noise.
#
#  Aims      How do the missing-data methods affect the pooled actor/partner
#            effects (rma() random effects, REML, z) of the two-study APIMeM
#            meta-analysis under MCAR, MAR and MNAR item nonresponse?
#  Data-generating mechanisms (full factorial, 9 conditions)
#            population: the toy model of simulate_toy_data.R (latent APIMeM,
#            4 Likert items per construct; Study 1 N = 130 with loadings
#            .72-.80, Study 2 N = 174 with loadings .55-.66), auxiliaries and
#            instruments as in 02; mechanism MCAR / MAR / MNAR (as in 02);
#            rate .05 / .10 / .15 of every item (exact)
#  Estimands the six partner-pooled effects (actor/partner a, b, c'), corrected
#            for unreliability (truth: latent paths) and uncorrected (truth: the
#            mean of the two studies' population scale-score paths); and, in
#            every replication, the complete-data (COM) estimate
#  Methods   COM, PRO_FIML, FIML_TS, mice_pmm (m = MC_M) and PMM_avg (the
#            same imputations averaged into one dataset: illustration of why
#            the imputations are pooled with Rubin's rules); add mice_cat or
#            heckman with MC_METHODS=...; V = ML
#  Performance  bias and RMSE (vs truth and vs COM), empirical and model SE,
#            coverage of the 95% CI, CI width relative to COM, rejection rate;
#            Monte Carlo SEs (Morris et al., 2019, Table 6). The summary also
#            makes the pre-specified mice-PMM vs FIML decision on these results
#            (choose_method() in 00_functions.R), which 08 and 09 then use.
#
#  Running
#    Rscript 05_monte_carlo.R                all conditions (resumes if interrupted)
#    Rscript 05_monte_carlo.R 7              condition 7 only (e.g. an HPC array job)
#    Rscript 05_monte_carlo.R summarise      only (re)build the summaries
#    MC_MODE=quick | pilot | full            2 / 50 / 500 replications (default pilot)
#    MC_NSIM=20                              override the number of replications
#    WORKERS=8                               replications in parallel
#  From RStudio: Sys.setenv(MC_MODE = "quick"); source("05_monte_carlo.R")
#  Every replication has its own L'Ecuyer-CMRG stream, so results do not depend
#  on WORKERS or on the order of the chunks.
#  Runtime (one core): about 1 min per replication and condition (m = 10);
#  pilot = 9 x 50 = ~8 core-hours, full = 9 x 500 = ~75 core-hours.
# =============================================================================

source("00_functions.R")

## 0. Settings

MODE <- Sys.getenv("MC_MODE", "pilot")

CFG <- switch(MODE,
  quick = list(NSIM = 2,   M = 3,  maxit = 5,  CHUNK = 1),
  pilot = list(NSIM = 50,  M = 10, maxit = 10, CHUNK = 5),
  full  = list(NSIM = 500, M = 10, maxit = 10, CHUNK = 10),
  stop("MC_MODE must be quick, pilot or full"))

CFG$M <- as.integer(Sys.getenv("MC_M", CFG$M))

CFG$NSIM <- as.integer(Sys.getenv("MC_NSIM", CFG$NSIM))

METHODS <- strsplit(Sys.getenv("MC_METHODS", "COM,PRO_FIML,FIML_TS,mice_pmm,PMM_avg"), ",")[[1]]
stopifnot(all(METHODS %in% METHODS_ALL), "COM" %in% METHODS)

MASTER  <- 20261016

WORKERS <- as.integer(Sys.getenv("WORKERS", "10"))

OUT <- file.path("output", "mc", MODE); RAW <- file.path(OUT, "raw")
dir.create(RAW, FALSE, TRUE)

# The chunks in a mode folder must come from the same settings (only the number
# of replications may grow), so results of other settings are never mixed in

MC_SETTINGS <- list(M = CFG$M, maxit = CFG$maxit, CHUNK = CFG$CHUNK, methods = METHODS, master = MASTER)

f_set <- file.path(OUT, "settings.rds")

check_mc_settings <- function() {
  have <- length(list.files(RAW, pattern = "^cond\\d+_chunk\\d+\\.rds$")) > 0
  old <- if (file.exists(f_set)) readRDS(f_set) else NULL
  if (have && !identical(old, MC_SETTINGS))
    stop(sprintf("%s holds results of other settings%s. Delete that folder (or restore its settings) and rerun.",
                 OUT, if (is.null(old)) " (or of an earlier version of this script)" else ""))
  saveRDS(MC_SETTINGS, f_set)
}

DESIGN <- expand.grid(rate = RATES, mechanism = MECHANISMS, stringsAsFactors = FALSE)[, 2:1]

DESIGN$condition <- seq_len(nrow(DESIGN))

DESIGN$label <- sprintf("%s_%02dpct", DESIGN$mechanism, round(100 * DESIGN$rate))

POP <- toy_population(file.path("data", "population_values.rds"))

SIGMA_T <- toy_implied_sigma()

## 1. Random-number streams

rep_seeds <- function(cond) {
  with_rng(c("L'Ecuyer-CMRG", "Inversion", "Rejection"), MASTER + cond, {   # restores the caller's generator
    s <- .Random.seed; out <- vector("list", CFG$NSIM)
    for (r in seq_len(CFG$NSIM)) { out[[r]] <- s; s <- parallel::nextRNGStream(s) }
    out
  })
}

## 2. One replication

one_rep <- function(cond_row, r, seed) {
  if (!exists("mice.impute.catx", envir = globalenv()))      # custom mice methods (looked up by name)
    sys.source("00_functions.R", envir = globalenv())
  old_kind <- RNGkind(); old_seed <- get0(".Random.seed", envir = globalenv())
  on.exit({ RNGkind(old_kind[1], old_kind[2], old_kind[3])      # restore the caller's generator
            if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv()) }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG"); assign(".Random.seed", seed, envir = globalenv())
  sr <- setNames(vector("list", length(METHODS)), METHODS); secs <- setNames(numeric(length(METHODS)), METHODS)
  errs <- character(0)
  for (k in seq_len(nrow(STUDIES))) {
    key <- STUDIES$key[k]
    comp <- add_aux_iv(toy_make_study(STUDIES$n[k], TOY$loadings[[key]], TOY$thresholds, SIGMA_T))
    mis  <- ampute_items(comp, cond_row$mechanism, cond_row$rate)
    imp_cache <- list()                       # PMM_avg re-uses the mice-PMM imputations
    for (mt in METHODS) {
      t0 <- proc.time()[["elapsed"]]
      res <- tryCatch(suppressWarnings({
        src <- imputations_for(mt); imps <- NULL
        if (!is.na(src)) {
          if (is.null(imp_cache[[src]]))
            imp_cache[[src]] <- impute_mice(mis, src, CFG$M, CFG$maxit, seed = NULL, donors = 5L)
          imps <- imp_cache[[src]]
        }
        run_method_study(mt, if (mt == "COM") comp else mis, STUDIES$study[k], imps, ml_only = TRUE)
      }), error = function(e) e)
      secs[mt] <- secs[mt] + proc.time()[["elapsed"]] - t0
      if (inherits(res, "error")) { errs <- c(errs, sprintf("%s/%s: %s", key, mt, conditionMessage(res))); next }
      sr[[mt]][[STUDIES$study[k]]] <- res
    }
  }
  est <- do.call(rbind, lapply(METHODS, function(mt) {
    if (length(sr[[mt]]) < nrow(STUDIES)) return(NULL)
    out <- tryCatch(meta_both(sr[[mt]], "ML", methods = c("RE_z", "CE_z"), n_mc = 2000)$direct,
                    error = function(e) NULL)
    if (is.null(out)) return(NULL)
    fmi <- if (mt %in% IMPUTATION_METHODS) mean(unlist(lapply(sr[[mt]], function(s)
      c(s$uncorrected$mi$ML$fmi, s$corrected$mi$ML$fmi)))) else NA_real_
    cbind(rep = r, missing_method = mt, out[, c("component", "method", "effect", "est", "se", "lo", "hi", "p", "tau2")],
          mean_fmi = fmi)
  }))
  list(est = est, misc = data.frame(rep = r, method = METHODS, seconds = secs,
                                    error = paste(errs, collapse = " | ")))
}

## 3. Run a condition in check-pointed chunks

run_condition <- function(cond) {
  cr <- DESIGN[cond, ]; seeds <- rep_seeds(cond)
  chunks <- split(seq_len(CFG$NSIM), ceiling(seq_len(CFG$NSIM) / CFG$CHUNK))
  todo <- Filter(function(k) !file.exists(file.path(RAW, sprintf("cond%02d_chunk%03d.rds", cond, k))),
                 seq_along(chunks))
  if (!length(todo)) { cat(sprintf("Condition %d: already complete\n", cond)); return(invisible()) }
  cat(sprintf("Condition %d (%s): %d chunk(s) to run\n", cond, cr$label, length(todo)))
  job <- function(k) {
    t0 <- proc.time()[["elapsed"]]
    out <- lapply(chunks[[k]], function(r) one_rep(cr, r, seeds[[r]]))
    res <- list(est = do.call(rbind, lapply(out, `[[`, "est")),
                misc = do.call(rbind, lapply(out, `[[`, "misc")), reps = chunks[[k]])
    saveRDS(res, file.path(RAW, sprintf("cond%02d_chunk%03d.rds", cond, k)))
    (proc.time()[["elapsed"]] - t0) / length(chunks[[k]])
  }
  secs <- unlist(par_lapply(todo, job, workers = WORKERS))
  cat(sprintf("  done: %.1f s per replication\n", mean(secs)))
}

## 4. Summaries

mcse_p <- function(p, n) sqrt(p * (1 - p) / n)

summarise_mc <- function() {
  files <- list.files(RAW, pattern = "^cond\\d+_chunk\\d+\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No raw results in ", RAW)
  raw <- lapply(files, readRDS)
  cond_of <- as.integer(sub("^cond(\\d+)_.*", "\\1", basename(files)))
  est  <- do.call(rbind, Map(function(x, c) cbind(condition = c, x$est), raw, cond_of))
  misc <- do.call(rbind, Map(function(x, c) cbind(condition = c, x$misc), raw, cond_of))
  est$truth <- ifelse(est$component == "corrected", POP$effect_truth$corrected[est$effect],
                      POP$effect_truth$uncorrected[est$effect])
  com <- est[est$missing_method == "COM", c("condition", "rep", "component", "method", "effect", "est", "lo", "hi")]
  names(com)[6:8] <- c("est_com", "lo_com", "hi_com")
  x <- merge(est, com, by = c("condition", "rep", "component", "method", "effect"))
  per <- do.call(rbind, lapply(split(x, list(x$condition, x$missing_method, x$component, x$method, x$effect),
                                     drop = TRUE), function(d) {
    d <- d[is.finite(d$est) & is.finite(d$se), ]; n <- nrow(d); if (n < 2) return(NULL)
    err <- d$est - d$truth; dc <- d$est - d$est_com; mse <- mean(err^2)
    cov <- mean(d$lo <= d$truth & d$truth <= d$hi); rej <- mean(d$p < .05)
    emp <- stats::sd(d$est); mod <- sqrt(mean(d$se^2))
    data.frame(condition = d$condition[1], missing_method = d$missing_method[1], component = d$component[1],
               method = d$method[1], effect = d$effect[1], truth = d$truth[1], n_rep = n,
               bias = mean(err), mcse_bias = stats::sd(err) / sqrt(n), rmse = sqrt(mse),
               bias_vs_com = mean(dc), rmsd_vs_com = sqrt(mean(dc^2)),
               emp_se = emp, mod_se = mod, rel_err_se = 100 * (mod / emp - 1),
               coverage = cov, mcse_coverage = mcse_p(cov, n),
               ci_width = mean(d$hi - d$lo), ci_ratio = mean((d$hi - d$lo) / (d$hi_com - d$lo_com)),
               rejection = rej, mcse_rejection = mcse_p(rej, n), mean_fmi = mean(d$mean_fmi))
  }))
  per <- merge(DESIGN[, c("condition", "mechanism", "rate", "label")], per, by = "condition")
  by_method <- aggregate(cbind(abs_bias = abs(bias), rmse, rmsd_vs_com, coverage, ci_ratio, rejection, mean_fmi) ~
                           condition + mechanism + rate + label + missing_method + component + method,
                         per, function(v) mean(v, na.rm = TRUE), na.action = na.pass)
  by_method <- by_method[order(by_method$condition, match(by_method$missing_method, METHODS_ALL),
                               match(by_method$component, COMPONENTS)), ]
  rt <- aggregate(seconds ~ condition + method, misc, mean)
  S <- list(mode = MODE, cfg = CFG, design = DESIGN, methods = METHODS, per_effect = per,
            by_method = by_method, runtime = rt,
            n_rep = tapply(misc$rep, misc$condition, function(v) length(unique(v))),
            errors = misc[nzchar(misc$error), ], truth = POP$effect_truth)
  S$decision <- choose_method(MC = S)
  saveRDS(S, file.path(OUT, "mc_summary.rds"))
  writeLines(S$decision$text, file.path(OUT, "method_decision_mc.txt"))
  cat(S$decision$text, sep = "\n")
  write_xlsx_sheets(list(design = DESIGN, by_method = by_method, per_effect = per, runtime = rt,
                         decision = data.frame(text = S$decision$text),
                         errors = if (nrow(S$errors)) S$errors else data.frame(none = TRUE)),
                    file.path(OUT, "mc_results.xlsx"))
  cat(sprintf("Summaries written to %s (replications per condition: %s)\n", OUT,
              paste(S$n_rep, collapse = ", ")))
  invisible(S)
}

## 5. Main

args <- commandArgs(trailingOnly = TRUE)

if (!length(args) || args[1] != "summarise") {
  check_mc_settings()
  conds <- if (length(args)) as.integer(args) else DESIGN$condition
  cat(sprintf("MC mode = %s: NSIM = %d, m = %d, methods = %s, %d condition(s), workers = %d\n",
              MODE, CFG$NSIM, CFG$M, paste(METHODS, collapse = ", "), length(conds), WORKERS))
  t0 <- proc.time()[["elapsed"]]
  for (cond in conds) run_condition(cond)
  cat(sprintf("Elapsed: %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
}

if (!length(args) || args[1] == "summarise") S <- summarise_mc()
