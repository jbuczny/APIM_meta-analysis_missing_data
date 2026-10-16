
## title: "APIM meta-analysis: individual-level missing and final pipeline"
## author: "Anonymized"
## date: "15-10-2026"

# =============================================================================
# 09_final_meta_analysis.R
# The APIMeM meta-analysis of the two studies WITH item-level missing data:
# the script to run on the real data (N = 130 and N = 174 dyads).
#
#   1. Read the two incomplete study files and describe the missingness
#      (item, scale-score and dyad level).
#   2. Handle the missing items with the method chosen on the toy data
#      (choose_method(): the Monte Carlo decision of 05 if available, else the
#      single-realisation decision of 04; default mice-PMM, m = 20, item level,
#      parcel predictors + auxiliaries; or the two-stage item-level FIML).
#   3. Fit the APIMeM in lavaan with the functions of apim_meta_analysis.R
#      (uncorrected and corrected for unreliability; V = ML and MLR) to every
#      completed dataset and aggregate the m = 20 analyses WITHIN each study
#      with Rubin's rules (Burgess et al., 2013). (The imputed data values are
#      not averaged: averaging imputations treats imputed values as known and
#      understates the uncertainty.)
#   4. Pool the studies with rma(): random effects, REML, z test (primary; the
#      specification chosen from the sensitivity report), with the common-effect
#      and Knapp-Hartung versions for reference; indirect and total effects
#      with Monte Carlo CIs and the joint-significance test.
#   5. Optional: delta-adjusted MNAR tipping-point analysis (as in 06).
#   6. Optional: rerun apim_meta_sensitivity.R (rma() models, settings of the
#      original run) on the study-level results WITH the missing data: exact
#      minimum detectable effects, critical values and indirect-effect MDEs
#      (RUN_POWER_SENSITIVITY; see 07).
#
#   Output (output/final/): stacked imputed datasets (.csv), results (.xlsx,
#   .rds), APA tables (.docx), and APIM_meta_analysis_results_MI.rds, the
#   study-level input for apim_meta_sensitivity.R.
#
# FOR THE REAL DATA: set INPUT_FILES (one wide file per study, one row per
# dyad, items named x_p1_1 ... y_p2_4 as in the toy data, values 1-5, missing
# = empty/NA) and AUX (always-observed variables that predict nonresponse or the
# scale scores; e.g. relationship duration, demographics, related scales).
# Without auxiliaries set AUX <- character(0) (MAR is then assumed given the
# items only). By default the script runs on two incomplete toy datasets.
# =============================================================================

source("00_functions.R")

## 0. Settings

INPUT_FILES <- c("Study 1" = file.path("data", "amputed", "study1_MAR_10pct.csv"),    # <- your files
                 "Study 2" = file.path("data", "amputed", "study2_MAR_10pct.csv"))

AUX     <- c("aux_p1", "aux_p2", "aux_dyad")    # <- your auxiliary variables (complete)

f_dec   <- file.path("output", "evaluation", "04_evaluation_results.rds")

f_mc    <- file.path("output", "mc", c("full", "pilot"), "mc_summary.rds")

f_mc    <- f_mc[file.exists(f_mc)][1]

DEC     <- choose_method(if (file.exists(f_dec)) readRDS(f_dec), if (!is.na(f_mc)) readRDS(f_mc))

METHOD  <- Sys.getenv("FINAL_METHOD", DEC$chosen)          # or set it: "mice_pmm" / "FIML_TS"

M       <- 20

MAXIT   <- 20

SEED    <- 20261017

RUN_DELTA <- env_flag("RUN_DELTA", TRUE)    # MNAR tipping-point analysis (~2 min per delta)

DELTAS  <- c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)

RUN_POWER_SENSITIVITY <- env_flag("RUN_POWER", TRUE)  # apim_meta_sensitivity.R rerun (~5 min)

POWER_R <- 1000                             # replications for that rerun

OUT <- file.path("output", "final"); dir.create(OUT, FALSE, TRUE)

stopifnot(METHOD %in% c("mice_pmm", "mice_cat", "FIML_TS", "PRO_FIML"), length(INPUT_FILES) >= 2)

cat("Method:", METHOD_LABELS[[METHOD]], sprintf("(decision basis: %s)\n", DEC$basis))

## 1. Data and missingness

dat <- lapply(INPUT_FILES, function(f) {
  d <- utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  if (!ID_VAR %in% names(d)) d[[ID_VAR]] <- sprintf("D%03d", seq_len(nrow(d)))   # no dyad identifier: row numbers
  miss_cols <- setdiff(c(ITEMS, AUX), names(d))
  if (length(miss_cols)) stop(f, ": missing columns ", paste(miss_cols, collapse = ", "))
  bad <- !is.na(as.matrix(d[, ITEMS])) & !(as.matrix(d[, ITEMS]) %in% ITEM_VALUES)
  if (any(bad)) stop(f, ": item values outside ", paste(range(ITEM_VALUES), collapse = "-"))
  if (length(AUX) && anyNA(d[, AUX])) stop(f, ": the auxiliaries must be complete (impute them first or drop them)")
  d
})

miss <- do.call(rbind, lapply(names(dat), function(s) {
  ms <- missing_summary(dat[[s]])
  data.frame(study = s, n_dyads = nrow(dat[[s]]), ms[ms$variable %in% c(
    "ALL item cells", "Dyads with any missing item", "Dyads complete (all 24 items)",
    "Scale scores with any missing item", "Scale scores not computable by proration (< 75% items)"), ])
}))

print(miss, row.names = FALSE)

## 2-3. Missing-data handling and study-level APIMeM

study_res <- list(); imps_all <- list(); log <- list()

for (i in seq_along(dat)) {
  s <- names(dat)[i]; t0 <- proc.time()[["elapsed"]]
  if (METHOD %in% c("mice_pmm", "mice_cat")) {
    imps <- impute_mice(dat[[s]], METHOD, m = M, maxit = if (METHOD == "mice_pmm") MAXIT else 10,
                        seed = SEED + i, return_mids = TRUE, aux = AUX, donors = 5L)
    mids <- attr(imps, "mids"); attr(imps, "mids") <- NULL
    imps_all[[s]] <- imps
    study_res[[s]] <- analyse_mi(imps, s)
    ps <- mice_psrf(mids)
    log[[s]] <- data.frame(study = s, method = METHOD, m = M, max_psrf = if (length(ps)) max(ps) else NA,
                           seconds = proc.time()[["elapsed"]] - t0)
    stacked <- do.call(rbind, lapply(seq_along(imps), function(l) cbind(.imp = l, imps[[l]][, c(ID_VAR, ITEMS, AUX)])))
    utils::write.csv(stacked, file.path(OUT, sprintf("imputed_%s__%s.csv", gsub("\\s+", "", tolower(s)), METHOD)),
                     row.names = FALSE)
  } else {
    study_res[[s]] <- if (METHOD == "FIML_TS") analyse_fiml_ts(dat[[s]], s, aux = AUX) else analyse_prorated(dat[[s]], s)
    study_res[[s]]$moments <- NULL
    log[[s]] <- data.frame(study = s, method = METHOD, m = NA, max_psrf = NA,
                           seconds = proc.time()[["elapsed"]] - t0)
  }
  cat(sprintf("  %s done (%.0f s)\n", s, log[[s]]$seconds))
}

log <- do.call(rbind, log)

## 4. Meta-analysis with rma()

MA <- meta_both(study_res, c("ML", "MLR"), n_mc = 1e5, mc_seed = SEED)

f_sens_in <- file.path(OUT, "APIM_meta_analysis_results_MI.rds")

write_sensitivity_input(study_res, f_sens_in, note = METHOD_LABELS[[METHOD]])

prim_d <- MA$direct[MA$direct$method == PRIMARY_RMA & MA$direct$V_type == PRIMARY_V, ]

cat("\nPooled direct effects (rma() RE, REML, z; V = ML):\n")

print(transform(prim_d[, c("component", "effect", "est", "se", "lo", "hi", "p", "tau2")],
                est = round(est, 3), se = round(se, 3), lo = round(lo, 3), hi = round(hi, 3),
                p = signif(p, 3), tau2 = signif(tau2, 3)), row.names = FALSE)

study_tab <- do.call(rbind, lapply(names(study_res), function(s) do.call(rbind, lapply(COMPONENTS, function(cmp) {
  r <- study_res[[s]][[cmp]]
  data.frame(study = s, component = cmp, path = PATHS$label, effect = PATHS$effect, est = r$est[PATHS$label],
             se_ML = sqrt(diag(r$V$ML))[PATHS$label], se_MLR = sqrt(diag(r$V$MLR))[PATHS$label],
             fmi = if (!is.null(r$mi)) r$mi$ML$fmi else NA_real_, row.names = NULL)
}))))

## 5. Delta-adjusted MNAR tipping points (mice-PMM)

DELTA <- NULL

if (RUN_DELTA) {
  cat("\nDelta-adjusted tipping-point analysis ...\n")
  DELTA <- do.call(rbind, lapply(DELTAS, function(dl) {
    sr <- setNames(lapply(seq_along(dat), function(i)
      analyse_mi(impute_delta(dat[[i]], dl, M, MAXIT, seed = SEED + i, aux = AUX), names(dat)[i], ml_only = TRUE)),
      names(dat))
    md <- meta_both(sr, "ML", methods = PRIMARY_RMA, n_mc = 2e4, mc_seed = SEED)
    cbind(delta = dl, md$direct)
  }))
}

## 6. Optional: apim_meta_sensitivity.R with the missing-data sampling covariances

POWER <- NULL

if (RUN_POWER_SENSITIVITY) {
  cat("\nPower sensitivity (apim_meta_sensitivity.R, rma() models) with the missing-data covariances ...\n")
  f_res <- find_sensitivity_file("APIM_sensitivity_results.rds")
  ovr <- if (!is.null(f_res)) { cf <- readRDS(f_res)$config
    cf[setdiff(names(cf), c("rds_file", "input_mode", "n_cores", "validate_mc", "out_prefix"))] } else list()
  ovr <- utils::modifyList(ovr, list(input_mode = "rds", validate_mc = FALSE, n_cores = 1L, R = POWER_R))
  POWER <- run_sensitivity_script(f_sens_in, ovr, rma_only = TRUE)
  cat(sprintf("  done (%.0f s)\n", POWER$seconds))
}

## 7. Tables

sheets <- list(missingness = miss, log = log, pooled_direct = MA$direct, pooled_indirect = MA$indirect,
               study_level = study_tab)

if (!is.null(DELTA)) sheets$delta_adjusted <- DELTA

if (!is.null(POWER)) { sheets$power_mde <- POWER$summary; sheets$power_indirect_mde <- POWER$indirect_min
                       sheets$power_critical <- POWER$observed }

sheets$method_decision <- data.frame(text = DEC$text)

write_xlsx_sheets(sheets, file.path(OUT, "09_final_results.xlsx"))

saveRDS(list(method = METHOD, decision = DEC, files = INPUT_FILES, aux = AUX, missingness = miss,
             study_results = study_res, meta = MA, delta = DELTA, power = POWER, log = log),
        file.path(OUT, "09_final_results.rds"))

pct <- "95%"

k_txt <- sprintf("Pooled across %d studies (%s dyads)", length(dat), paste(vapply(dat, nrow, 0L), collapse = " and "))

meth_txt <- switch(METHOD,
  mice_pmm = sprintf(paste0("Item-level missing values were multiply imputed with predictive mean matching (mice; m = %d; ",
                            "predictors: the other items of the scale, the five other scale scores and %s); the APIMeM was ",
                            "fitted to every completed dataset and the m analyses were pooled within each study with ",
                            "Rubin's rules."), M, if (length(AUX)) paste("the auxiliaries", paste(AUX, collapse = ", ")) else "no auxiliaries"),
  mice_cat = sprintf("Item-level multiple imputation with proportional-odds models (m = %d), pooled within each study with Rubin's rules.", M),
  FIML_TS  = "Item-level missing values were handled with two-stage FIML (saturated item-level model with the auxiliaries, then the APIMeM on the implied scale-score moments).",
  PRO_FIML = "Scale scores were prorated (\u{2265} 75% of items) and the APIMeM was fitted with FIML.")

doc <- new_landscape_doc()

## Table 1: missingness and reliability

t1 <- do.call(rbind, lapply(names(dat), function(s) {
  ms <- missing_summary(dat[[s]]); sc_any <- is.na(scale_scores(dat[[s]]))
  data.frame(v = VAR_LABELS[VARS],
             items = vapply(VARS, function(v) sprintf("%.1f", 100 * mean(is.na(dat[[s]][, ITEMS[ITEM_VAR == v]]))), ""),
             scale = sprintf("%.1f", 100 * colMeans(sc_any)),
             rel = fmt(study_res[[s]]$rel[VARS], 2, lead0 = FALSE), stringsAsFactors = FALSE)
}))

g1 <- grouped_body(t1, rep(sprintf("%s (%d dyads; %.1f%% of item values missing)", names(dat), vapply(dat, nrow, 0L),
                                   vapply(dat, function(d) 100 * mean(is.na(d[, ITEMS])), 0)), each = 6))

ft1 <- apa_table(g1$body, labels = c("Scale score", "Items missing (%)", "Scale scores incomplete (%)", "\u{03c9}"),
                 group_rows = g1$group_rows,
                 note = paste0("Scale scores incomplete = at least one of the four items missing. \u{03c9} = McDonald's omega ",
                               "total (one-factor CFA), averaged over the completed datasets for multiple imputation; used ",
                               "for the correction for unreliability. ", meth_txt),
                 first_col_width = 3.0, col_width = 1.3)

doc <- add_apa_table(doc, 1, "Item-Level Missing Data and Reliability by Study", ft1, new_page = FALSE)

## Table 2: pooled direct effects

t2 <- NULL

for (e in EFFECTS) {
  row <- data.frame(e = EFFECT_LAB[[e]], stringsAsFactors = FALSE)
  for (cmp in COMPONENTS) {
    x <- prim_d[prim_d$component == cmp & prim_d$effect == e, ]
    row[[paste0(cmp, "_b")]] <- fmt(x$est); row[[paste0(cmp, "_se")]] <- fmt(x$se)
    row[[paste0(cmp, "_ci")]] <- fmt_ci(x$lo, x$hi); row[[paste0(cmp, "_p")]] <- fmt_p(x$p)
    row[[paste0(cmp, "_t")]] <- fmt(x$tau2, 3)
  }
  t2 <- rbind(t2, row)
}

g2 <- grouped_body(t2, APIM$PATH_GROUP_LABEL[sub("^(actor|partner)_", "", EFFECTS)])

ft2 <- apa_table(g2$body, labels = c("Effect", rep(c("\u{03b2}", "SE", paste(pct, "CI"), "p", "\u{03c4}\u{00b2}"), 2)),
  spans = list(c(" " = 1L, "Uncorrected" = 5L, "Corrected for unreliability" = 5L)), group_rows = g2$group_rows,
  note = paste0(k_txt, ". \u{03b2} = pooled standardized APIMeM path: per study, the Partner 1 and Partner 2 estimates of ",
                "each effect were combined (GLS, using their sampling covariance matrix from lavaan), and the study ",
                "composites were pooled with a random-effects model (rma(); REML; Wald z). ", meth_txt,
                " With two studies \u{03c4}\u{00b2} is very imprecise; results of the common-effect and Knapp\u{2013}Hartung ",
                "versions are in the .xlsx file."),
  first_col_width = 1.2, col_width = 0.66)

ft2 <- flextable::width(ft2, j = c(4, 9), width = 1.05)

doc <- add_apa_table(doc, 2, "Pooled Direct Effects of the APIMeM With Item-Level Missing Data", ft2)

## Table 3: indirect and total effects

ind <- MA$indirect[MA$indirect$V_type == PRIMARY_V, ]

t3 <- NULL

for (lab in unique(ind$label)) {
  row <- data.frame(e = lab, stringsAsFactors = FALSE)
  for (cmp in COMPONENTS) {
    x <- ind[ind$component == cmp & ind$label == lab, ]
    row[[paste0(cmp, "_e")]] <- fmt(x$est, 3); row[[paste0(cmp, "_ci")]] <- fmt_ci(x$lo, x$hi, 3)
    row[[paste0(cmp, "_js")]] <- if (is.na(x$joint_sig)) "\u{2014}" else if (x$joint_sig) "Yes" else "No"
  }
  t3 <- rbind(t3, row)
}

g3 <- grouped_body(t3, ind$group[match(t3$e, ind$label)])

ft3 <- apa_table(g3$body, labels = c("Effect", rep(c("Estimate", paste(pct, "MC CI"), "Joint sig."), 2)),
  spans = list(c(" " = 1L, "Uncorrected" = 3L, "Corrected for unreliability" = 3L)), group_rows = g3$group_rows,
  note = paste0(k_txt, ". Products and sums of the pooled random-effects paths (Ledermann et al., 2011). MC CI = Monte ",
                "Carlo CI (100,000 draws from the pooled paths and their plug-in covariance matrix). Joint sig. = both ",
                "component paths significant (p < .05)."),
  first_col_width = 3.4, col_width = 0.9)

ft3 <- flextable::width(ft3, j = c(3, 6), width = 1.3)            # room for the CIs

doc <- add_apa_table(doc, 3, "Pooled Indirect and Total Effects", ft3)

## Table 4: study-level paths

t4 <- NULL

for (s in names(study_res)) {
  blk <- data.frame(p = sprintf("%s: %s \u{2192} %s", ifelse(PATHS$type == "actor", "Actor", "Partner"), PATHS$rhs, PATHS$lhs),
                    stringsAsFactors = FALSE)
  for (cmp in COMPONENTS) {
    x <- study_tab[study_tab$study == s & study_tab$component == cmp, ]
    blk[[paste0(cmp, "_b")]] <- fmt(x$est); blk[[paste0(cmp, "_se")]] <- fmt(x$se_ML)
    blk[[paste0(cmp, "_f")]] <- fmt(x$fmi, 2, lead0 = FALSE)
  }
  t4 <- rbind(t4, blk)
}

g4 <- grouped_body(t4, rep(sprintf("%s (%d dyads)", names(study_res), vapply(dat, nrow, 0L)), each = 12))

ft4 <- apa_table(g4$body, labels = c("Path", rep(c("\u{03b2}", "SE", "FMI"), 2)),
  spans = list(c(" " = 1L, "Uncorrected" = 3L, "Corrected for unreliability" = 3L)), group_rows = g4$group_rows,
  note = paste0("Standardized APIMeM paths of each study after the missing-data handling. SE = normal-theory (ML) ",
                "standard error including the between-imputation variance; FMI = fraction of missing information ",
                "(Barnard\u{2013}Rubin; multiple imputation only). The full 12 \u{00d7} 12 covariance matrices enter the ",
                "meta-analysis."),
  first_col_width = 2.0, col_width = 0.8, size = 9)

doc <- add_apa_table(doc, 4, "Study-Level Standardized APIMeM Paths With Item-Level Missing Data", ft4)

## Table 5: delta adjustment

if (!is.null(DELTA)) {
  t5 <- NULL
  for (cmp in COMPONENTS) for (e in EFFECTS) {
    row <- data.frame(e = EFFECT_LAB[[e]], stringsAsFactors = FALSE)
    for (dl in DELTAS) {
      x <- DELTA[DELTA$component == cmp & DELTA$effect == e & abs(DELTA$delta - dl) < 1e-9, ]
      row[[as.character(dl)]] <- paste0(fmt(x$est), ifelse(x$p < .05, "*", ""))
    }
    t5 <- rbind(t5, row)
  }
  g5 <- grouped_body(t5, rep(COMP_LAB[COMPONENTS], each = length(EFFECTS)))
  ft5 <- apa_table(g5$body, labels = c("Effect", sub("^-", MINUS, sprintf("%.2f", DELTAS))),
    spans = list(c(" " = 1L, setNames(length(DELTAS), "\u{03b4} (item-SD units; positive = hypothesized MNAR direction)"))),
    group_rows = g5$group_rows,
    note = paste0("Delta-adjusted mice-PMM: imputed item values shifted by \u{03b4} SDs of the observed item values in the ",
                  "hypothesized MNAR direction (", paste(sprintf("%s %+d", toupper(names(MNAR_DIRECTION)), MNAR_DIRECTION),
                                                          collapse = ", "), "); \u{03b4} = 0 is MAR. Pooled with rma() ",
                  "random effects; * p < .05. A conclusion is robust if it holds over the range of \u{03b4} considered plausible."),
    first_col_width = 1.3, col_width = 0.85)
  doc <- add_apa_table(doc, 5, "Sensitivity of the Pooled Direct Effects to Departures From MAR (Delta Adjustment)", ft5)
}

## Table 6: minimum detectable effects with the missing data (power sensitivity)

if (!is.null(POWER)) {
  sm <- POWER$summary[as.character(POWER$summary$method) == RMA_METHODS[[PRIMARY_RMA]] &
                        as.character(POWER$summary$k) == as.character(POWER$summary$k[1]), ]
  hl <- levels(sm$het)
  t6 <- NULL
  for (cmp in COMP_LAB) for (e in EFFECT_LAB) {
    x <- sm[as.character(sm$component) == cmp & as.character(sm$effect) == e, ]
    row <- data.frame(e = e, stringsAsFactors = FALSE)
    for (h in hl) row[[h]] <- ifelse(is.na(x$mde[as.character(x$het) == h]), sprintf("> %.2f", POWER$config$theta_max),
                                     fmt(x$mde[as.character(x$het) == h]))
    t6 <- rbind(t6, row)
  }
  g6 <- grouped_body(t6, rep(unname(COMP_LAB), each = length(EFFECT_LAB)))
  ft6 <- apa_table(g6$body, labels = c("Effect", hl), group_rows = g6$group_rows,
    note = sprintf(paste0("Minimum detectable pooled standardized effect (80%% power, two-sided \u{03b1} = .05; rma() random ",
                          "effects, REML, z test; %d replications) given the study-level sampling covariance matrices ",
                          "after the missing-data handling (%s), from apim_meta_sensitivity.R with the settings of the ",
                          "original sensitivity analysis. Compare with Table 1 of the sensitivity report."),
                   POWER$config$R, METHOD_LABELS[[METHOD]]),
    first_col_width = 2.0, col_width = 1.1)
  doc <- add_apa_table(doc, 6, "Minimum Detectable Pooled Direct Effects With the Missing Data", ft6)
}

print(doc, target = file.path(OUT, "09_final_tables.docx"))

cat("\nWritten to", OUT, ":", paste(list.files(OUT), collapse = ", "), "\n")
