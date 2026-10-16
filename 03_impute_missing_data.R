
## title: "APIM meta-analysis: individual-level missing and final pipeline"
## author: "Anonymized"
## date: "15-10-2026"

# =============================================================================
# 03_impute_missing_data.R
# Impute the 18 incomplete datasets from 02 at the ITEM level (m = 20):
#
#   mice_pmm  mice, predictive mean matching for every item (DEFAULT METHOD)
#   mice_cat  mice, proportional-odds model for every (5-point) item, parameters
#             drawn from their posterior ("catx")
#   heckman   OPTIONAL (not run by default; METHODS=mice_pmm,mice_cat,heckman):
#             Heckman-type selection-model MI, ordered-probit outcome equation
#             with probit selection, the partner's instrument in the selection
#             equation only ("heckoprob"). It is the only MNAR model here, but
#             at the item level it estimates 24 selection models from 6-26
#             missing values each with a weak instrument (first-stage |z| ~ 2-3
#             in the MNAR data), and it takes ~10-15 min per dataset (m = 20).
#             In the m = 2 pilot it was not closer to the complete data than
#             two-stage FIML even under MNAR.
#
# Imputation model of each incomplete item (00_functions.R, section 6): the
# other items of the same scale, the scale scores of the other five scales
# (updated passively inside the chains, so own and partner scores are both
# there), and the auxiliaries aux_p1, aux_p2, aux_dyad. Heckman-MI adds the
# instrument of the item's partner. Every imputed value is one of 1-5, and the
# scale scores are recomputed from the completed items in the analysis.
#
# FIML (two-stage, item level), proration + FIML and the complete-data benchmark
# need no imputation; they are run in 04 together with the analysis of these
# imputations.
#
# Output
#   data/imputed/<dataset>__<method>.csv   stacked completed datasets (.imp = 1..m)
#   data/imputed/imputations_all.rds       list[dataset][method] of completed data
#   output/imputations/raw/<dataset>.rds   checkpoint per dataset (the script
#                                          resumes where it stopped)
#   output/imputations/imputation_log.xlsx runtime, convergence (PSRF), quality
#                                          of the imputed values, instrument relevance
#
# Runtime per dataset (one core, m = 20): mice-PMM ~0.5 min, mice-CAT ~1.5 min
# (about 40 min for the 18 datasets); Heckman-MI ~10-15 min.
# Environment variables: WORKERS (datasets in parallel), M_IMP (m),
# METHODS (comma separated, default mice_pmm,mice_cat). Fixed seeds per dataset
# and method make the imputations reproducible and independent of WORKERS; the
# checkpoints let an interrupted run resume (or add a method later).
# =============================================================================

source("00_functions.R")

## 0. Settings

IN_RDS     <- file.path("data", "amputed", "amputed_datasets.rds")

OUT_DATA   <- file.path("data", "imputed")

OUT_LOG    <- file.path("output", "imputations")

RAW        <- file.path(OUT_LOG, "raw")

M          <- as.integer(Sys.getenv("M_IMP", "20"))

WORKERS    <- as.integer(Sys.getenv("WORKERS", "15"))

METHODS    <- strsplit(Sys.getenv("METHODS", "mice_pmm,mice_cat"), ",")[[1]]

WRITE_XLSX <- FALSE                 # stacked .xlsx copies (about 1 MB each)

SEED       <- 20261013

SETTINGS   <- list(maxit_pmm = 20, donors = 5L, maxit_cat = 10, maxit_heck = 10, pred_set = PRED_SET)

stopifnot(all(METHODS %in% IMPUTATION_METHODS))

dir.create(OUT_DATA, FALSE, TRUE); dir.create(RAW, FALSE, TRUE)

if (!file.exists(IN_RDS)) source("02_ampute_data.R", local = new.env())
A <- readRDS(IN_RDS)
grid <- A$grid; datasets <- grid$dataset

# Seed of dataset i and method j (also used by 06 for delta-adjusted mice-PMM,
# so that delta = 0 reproduces the mice-PMM imputations)

imp_seed <- function(i, method) SEED + 1000 * i + match(method, IMPUTATION_METHODS)

# Cache key of dataset i and method j: the amputed data, m, the settings, the
# seed and the generator. A checkpoint is reused only for the same key, so
# imputations of other amputed data (e.g. after 01/02 were rerun) are redone.

imp_key <- function(i, method)
  cache_key("impute", method, M, SETTINGS, imp_seed(i, method), RNG_IMPUTATION,
            A$amputed[[datasets[i]]][, c(ID_VAR, ITEMS, AUX_VARS, IV_VARS)])

## 1. One job = one dataset x all methods

impute_job <- function(i) {
  if (!exists("mice.impute.heckoprob", envir = globalenv())) sys.source("00_functions.R", envir = globalenv())
  nm <- datasets[i]; f_raw <- file.path(RAW, paste0(nm, ".rds"))
  old <- if (file.exists(f_raw)) readRDS(f_raw) else list(imputations = list(), log = NULL)
  keys <- vapply(IMPUTATION_METHODS, function(mt) imp_key(i, mt), "")
  stale <- names(old$imputations)[!vapply(names(old$imputations), function(mt)
    identical(old$cache_keys[[mt]], keys[[mt]]), TRUE)]
  if (length(stale)) {                       # computed from other data or settings: discard
    old$imputations[stale] <- NULL; old$quality[stale] <- NULL; old$cache_keys[stale] <- NULL
    old$log <- old$log[!old$log$method %in% stale, , drop = FALSE]
    message(sprintf("  %s: cached %s imputations do not match the current amputed data or settings; redone",
                    nm, paste(stale, collapse = ", ")))
  }
  done <- names(old$imputations)
  dat <- A$amputed[[nm]]; comp <- A$complete[[grid$study[i]]]
  for (meth in setdiff(METHODS, done)) {
    t0 <- proc.time()[["elapsed"]]
    imps <- switch(meth,
      mice_pmm = impute_mice(dat, "mice_pmm", M, SETTINGS$maxit_pmm, imp_seed(i, meth), return_mids = TRUE,
                             pred_set = SETTINGS$pred_set, donors = SETTINGS$donors),
      mice_cat = impute_mice(dat, "mice_cat", M, SETTINGS$maxit_cat, imp_seed(i, meth), return_mids = TRUE,
                             pred_set = SETTINGS$pred_set),
      heckman  = impute_mice(dat, "heckman", M, SETTINGS$maxit_heck, imp_seed(i, meth), return_mids = TRUE,
                             pred_set = SETTINGS$pred_set))
    secs <- proc.time()[["elapsed"]] - t0
    mids <- attr(imps, "mids"); attr(imps, "mids") <- NULL
    ps <- mice_psrf(mids)
    q <- imputation_quality(comp, dat, imps)
    signed <- q$std_bias_mean * MNAR_DIRECTION[q$construct]
    old$log <- rbind(old$log[old$log$method != meth, , drop = FALSE], data.frame(
      dataset = nm, study = grid$study[i], mechanism = grid$mechanism[i], rate = grid$rate[i],
      method = meth, m = length(imps), seconds = round(secs, 1),
      max_psrf = if (length(ps)) max(ps) else NA_real_,
      nrmse = mean(q$nrmse, na.rm = TRUE), pct_exact = mean(q$pct_exact, na.rm = TRUE),
      std_bias_mean = mean(q$std_bias_mean), std_bias_mnar_direction = mean(signed),
      rel_bias_sd = mean(q$rel_bias_sd), tvd = mean(q$tvd), mabias_cor = mean(q$mabias_cor),
      logged_events = if (is.null(mids$loggedEvents)) 0L else nrow(mids$loggedEvents)))
    old$quality[[meth]] <- q
    old$imputations[[meth]] <- imps
    old$cache_keys[[meth]] <- keys[[meth]]
    old$iv <- cbind(dataset = nm, iv_first_stage(dat))
    saveRDS(old, f_raw)                                   # checkpoint
  }
  old
}

cat(sprintf("Imputing %d datasets x %d methods (m = %d, workers = %d) ...\n",
            length(datasets), length(METHODS), M, WORKERS))
t_all <- proc.time()[["elapsed"]]
out <- par_lapply(seq_along(datasets), function(i) {
  r <- impute_job(i)
  if (WORKERS == 1) cat(sprintf("  %-20s done (%s s)\n", datasets[i],
                                paste(round(r$log$seconds), collapse = "/")))
  r
}, workers = WORKERS)

names(out) <- datasets

cat(sprintf("Total imputation time: %.1f min\n", (proc.time()[["elapsed"]] - t_all) / 60))

## 2. Save
imputations <- lapply(out, function(o) o$imputations[intersect(IMPUTATION_METHODS, names(o$imputations))])
saveRDS(list(imputations = imputations, settings = SETTINGS, M = M, seed = SEED,
             methods = METHODS, grid = grid),
        file.path(OUT_DATA, "imputations_all.rds"))
for (nm in datasets) for (meth in names(imputations[[nm]])) {
  imps <- imputations[[nm]][[meth]]
  stacked <- do.call(rbind, lapply(seq_along(imps), function(l) cbind(.imp = l, imps[[l]])))
  base <- file.path(OUT_DATA, sprintf("%s__%s", nm, meth))
  utils::write.csv(stacked, paste0(base, ".csv"), row.names = FALSE)
  if (WRITE_XLSX) write_xlsx_sheets(list(imputed_stacked = stacked), paste0(base, ".xlsx"))
}

log_tab <- do.call(rbind, lapply(out, `[[`, "log")); rownames(log_tab) <- NULL

log_tab <- log_tab[order(match(log_tab$dataset, datasets), match(log_tab$method, IMPUTATION_METHODS)), ]

qual_tab <- do.call(rbind, lapply(names(out), function(nm) do.call(rbind, lapply(names(out[[nm]]$quality),
  function(meth) cbind(dataset = nm, method = meth, out[[nm]]$quality[[meth]])))))

iv_tab <- do.call(rbind, lapply(out, `[[`, "iv")); rownames(iv_tab) <- NULL

iv_wide <- reshape(iv_tab, idvar = "item", timevar = "dataset", direction = "wide")

names(iv_wide) <- sub("z_iv\\.", "", names(iv_wide)); rownames(iv_wide) <- NULL

settings_tab <- data.frame(
  setting = c("m", names(SETTINGS), "seed", "R", "mice", "lavaan"),
  value = c(M, vapply(SETTINGS, as.character, ""), SEED, R.version.string,
            as.character(utils::packageVersion("mice")), as.character(utils::packageVersion("lavaan"))))

spec_tab <- data.frame(item = ITEMS, scale = ITEM_SCALE, partner = ITEM_MEMBER,
                       predictors = vapply(ITEMS, function(it) {
                         v <- ITEM_VAR[[it]]
                         paste(c(setdiff(ITEMS[ITEM_VAR == v], it), setdiff(unname(SCALES), SCALES[[v]]),
                                 AUX_VARS), collapse = ", ") }, ""),
                       heckman_selection_only = paste0("iv_", ITEM_MEMBER), row.names = NULL)

write_xlsx_sheets(list(log = log_tab, quality_by_item = qual_tab, imputation_model = spec_tab,
                       instrument_relevance_z = iv_wide, settings = settings_tab),
                  file.path(OUT_LOG, "imputation_log.xlsx"))

## 3. Checks
for (nm in datasets) for (meth in names(imputations[[nm]])) {
  a <- A$amputed[[nm]]; obs <- !is.na(a[, ITEMS])
  for (d in imputations[[nm]][[meth]]) {
    stopifnot(!anyNA(d[, ITEMS]), all(d[, ITEMS][obs] == a[, ITEMS][obs]),
              all(as.matrix(d[, ITEMS]) %in% ITEM_VALUES),
              identical(d[, c(ID_VAR, AUX_VARS, IV_VARS)], a[, c(ID_VAR, AUX_VARS, IV_VARS)]))
  }
}

cat("Checks passed: no missing values left, every imputed value in 1-5, observed values,",
    "auxiliaries and instruments unchanged.\n\n")

print(log_tab[, c("dataset", "method", "seconds", "max_psrf", "nrmse", "std_bias_mnar_direction", "tvd")],
      digits = 3, row.names = FALSE)

cat("\nInstrument relevance (z of the partner's instrument in the response model), mean |z| by dataset:\n")

print(round(colMeans(abs(iv_wide[, -1]), na.rm = TRUE), 1))
