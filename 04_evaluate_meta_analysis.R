
## title: "APIM meta-analysis: individual-level missing and final pipeline"
## author: "Anonymized"
## date: "15-10-2026"

# =============================================================================
# 04_evaluate_meta_analysis.R
# Evaluate the missing-data methods in the analysis they are meant for: the
# APIMeM meta-analysis of the two studies.
#
#   For every condition (MCAR/MAR/MNAR x 5/10/15%) and method:
#     1. per study, the APIMeM is fitted in lavaan with the functions of
#        apim_meta_analysis.R (12 standardized paths and their 12 x 12 sampling
#        covariance matrix; uncorrected and corrected for unreliability; V = ML
#        and MLR). Multiply imputed data: one analysis per completed dataset,
#        pooled WITHIN the study with Rubin's rules (Burgess et al., 2013);
#     2. the two studies are pooled with rma(): Partner 1/2 estimates combined
#        per study (GLS), random-effects model, REML, z test - the specification
#        chosen from the sensitivity report (APIM_sensitivity_report.docx);
#        CE and Knapp-Hartung versions are kept for reference;
#     3. indirect effects: products of the pooled paths (Monte Carlo CI and
#        joint significance).
#   The results are compared with the complete-data analysis of the same toy
#   data (COM; isolates the effect of the missing data) and with the population
#   values. The rma() code is checked against the observed-data results of the
#   sensitivity report (it must reproduce them exactly); the power analysis
#   under missing data follows in 07.
#
#   Methods: COM (benchmark), PRO_FIML (proration + scale-level FIML: what
#   apim_meta_analysis.R does by default), FIML_TS (two-stage item-level FIML
#   with the auxiliaries), mice_pmm, mice_cat, heckman (imputations from 03,
#   pooled with Rubin's rules), and PMM_avg (the 20 mice-PMM datasets averaged
#   into one and analysed once: shown only to illustrate why the imputations
#   are aggregated with Rubin's rules rather than by averaging the data).
#
#   Decision rule (pre-specified; choose_method() in 00_functions.R): mice-PMM
#   vs FIML-TS are compared on the root mean square deviation of the pooled
#   direct effects from COM (in complete-data SE units) over the MCAR and MAR
#   conditions, where both are valid. A difference within 10% is a tie, and a
#   tie goes to mice-PMM (completed datasets; needed for the delta-adjusted
#   MNAR analysis of 06). One realisation cannot separate method differences
#   from noise: when 05 has been run (pilot/full), 08 and 09 use its decision.
#
# Output (output/evaluation/)
#   raw/<dataset>__<method>.rds   study-level results (cached; delete to rerun)
#   04_evaluation_results.rds     everything below, for 06-09
#   04_evaluation.xlsx            pooled estimates, comparisons, summaries,
#                                 study level, reproduction of the report
#   method_decision.txt           the PMM vs FIML decision (single realisation)
#
# Runtime: about 2 s per completed dataset (10 lavaan fits), i.e. ~1 min per
# MI method and dataset with m = 20; ~45 min in total on one core (WORKERS=n
# runs jobs in parallel; the cache makes reruns fast).
# =============================================================================

source("00_functions.R")

## 0. Settings

OUT       <- file.path("output", "evaluation"); RAW <- file.path(OUT, "raw")

WORKERS   <- as.integer(Sys.getenv("WORKERS", "15"))

N_MC      <- 2e4             # Monte Carlo draws for the indirect-effect CIs

MC_SEED   <- 20261014

V_TYPES   <- c("ML", "MLR")

dir.create(RAW, FALSE, TRUE)

A <- readRDS(file.path("data", "amputed", "amputed_datasets.rds"))

f_imp <- file.path("data", "imputed", "imputations_all.rds")

I <- if (file.exists(f_imp)) readRDS(f_imp) else NULL

POP  <- toy_population(file.path("data", "population_values.rds"))

grid <- A$grid; datasets <- grid$dataset

conditions <- unique(grid$condition)

mi_avail <- if (is.null(I)) character(0) else
  IMPUTATION_METHODS[vapply(IMPUTATION_METHODS, function(m)
    all(vapply(datasets, function(d) !is.null(I$imputations[[d]][[m]]), TRUE)), TRUE)]

METHODS <- c("COM", "PRO_FIML", "FIML_TS", mi_avail, if ("mice_pmm" %in% mi_avail) "PMM_avg")
if (!length(mi_avail)) message("No imputations found (run 03): MI methods are skipped.")

cat("Methods evaluated:", paste(METHOD_SHORT[METHODS], collapse = ", "), "\n")

study_label <- setNames(STUDIES$study, STUDIES$key)

## 1. Study-level analyses (cached)

jobs <- c(lapply(STUDIES$key, function(k) list(dataset = paste0(k, "_complete"), study = k, method = "COM")),
          unlist(lapply(seq_len(nrow(grid)), function(i) lapply(setdiff(METHODS, "COM"), function(mt)
            list(dataset = grid$dataset[i], study = grid$study[i], method = mt))), recursive = FALSE))

raw_file <- function(j) file.path(RAW, sprintf("%s__%s.rds", j$dataset, j$method))

job_data <- function(j) if (j$method == "COM") A$complete[[j$study]] else A$amputed[[j$dataset]]

job_imps <- function(j) {
  src_imp <- imputations_for(j$method)
  if (!is.na(src_imp)) I$imputations[[j$dataset]][[src_imp]] else NULL
}

# Cache key: the method and its exact input (complete/amputed data, imputations)

jobs <- lapply(jobs, function(j) { j$key <- cache_key("analysis", j$method, job_data(j), job_imps(j)); j })

todo <- Filter(function(j) is.null(read_cached(raw_file(j), j$key)), jobs)

cat(sprintf("Study-level analyses: %d jobs (%d cached), workers = %d\n", length(jobs),
            length(jobs) - length(todo), WORKERS))

run_job <- function(j) {
  if (!exists("analyse_fiml_ts", envir = globalenv())) sys.source("00_functions.R", envir = globalenv())
  t0 <- proc.time()[["elapsed"]]
  warns <- character(0)
  res <- withCallingHandlers(run_method_study(j$method, job_data(j), study_label[[j$study]], job_imps(j)),
                             warning = function(w) { warns <<- c(warns, conditionMessage(w))
                               invokeRestart("muffleWarning") })
  res$warnings <- unique(warns)
  if (j$method == "FIML_TS") { res$delta_check <- delta_check_fiml_ts(res); res$moments <- NULL }
  res$seconds <- proc.time()[["elapsed"]] - t0
  res$cache_key <- j$key
  saveRDS(res, raw_file(j))
  res$seconds
}

t_all <- proc.time()[["elapsed"]]

if (length(todo)) invisible(par_lapply(todo, run_job, workers = WORKERS))

cat(sprintf("  done in %.1f min\n", (proc.time()[["elapsed"]] - t_all) / 60))

SR <- setNames(lapply(jobs, function(j) {
  r <- read_cached(raw_file(j), j$key)
  if (is.null(r)) stop("Study-level result missing or stale: ", raw_file(j))
  r }), vapply(jobs, function(j) paste(j$dataset, j$method, sep = "|"), ""))

get_sr <- function(dataset, method) SR[[paste(dataset, method, sep = "|")]]

## Check of the two-stage FIML (lavaan stage 2 vs delta method)

dchk <- do.call(rbind, lapply(datasets, function(d) {
  x <- get_sr(d, "FIML_TS")$delta_check
  data.frame(dataset = d, component = rownames(x), max_abs_diff_est = x[, "est"], max_abs_diff_V = x[, "V"])
}))

cat(sprintf("Two-stage FIML: lavaan stage 2 vs delta method, max |diff| est = %.1e, V = %.1e\n",
            max(dchk$max_abs_diff_est), max(dchk$max_abs_diff_V)))

## 2. Meta-analysis with rma() per condition and method

study_set <- function(cond, method) {
  if (method == "COM") return(setNames(lapply(STUDIES$key, function(k) get_sr(paste0(k, "_complete"), "COM")),
                                       STUDIES$study))
  setNames(lapply(STUDIES$key, function(k) get_sr(paste(k, cond, sep = "_"), method)), STUDIES$study)
}

MA <- list()

MA[["COM"]] <- meta_both(study_set(NULL, "COM"), V_TYPES, n_mc = N_MC, mc_seed = MC_SEED)
for (cond in conditions) for (mt in setdiff(METHODS, "COM"))
  MA[[paste(cond, mt)]] <- meta_both(study_set(cond, mt), V_TYPES, n_mc = N_MC, mc_seed = MC_SEED)

cond_info <- unique(grid[, c("condition", "mechanism", "rate")])

lab_rows <- function(tab, cond, mt) {
  ci <- if (is.null(cond)) data.frame(condition = "complete", mechanism = "none", rate = 0) else
    cond_info[cond_info$condition == cond, ]
  cbind(ci[rep(1, nrow(tab)), ], missing_method = mt, tab, row.names = NULL)
}

DIRECT <- do.call(rbind, c(list(lab_rows(MA$COM$direct, NULL, "COM")),
  lapply(names(MA)[-1], function(nm) { p <- strsplit(nm, " ")[[1]]; lab_rows(MA[[nm]]$direct, p[1], p[2]) })))

INDIRECT <- do.call(rbind, c(list(lab_rows(MA$COM$indirect, NULL, "COM")),
  lapply(names(MA)[-1], function(nm) { p <- strsplit(nm, " ")[[1]]; lab_rows(MA[[nm]]$indirect, p[1], p[2]) })))

names(INDIRECT)[names(INDIRECT) == "short"] <- "effect"

## 3. Comparison with the complete data and with the truth

com_d <- DIRECT[DIRECT$missing_method == "COM", ]

CMP <- compare_to_com(DIRECT[DIRECT$missing_method != "COM", ], com_d, c("component", "V_type", "method", "effect"))
# Yardstick for deviations: the complete-data SE of the COMMON-EFFECT estimate
# (pure sampling error). With k = 2 the RE SE depends on an unstable tau^2
# estimate, which would make the deviations of effects with tau^2 > 0 look small.

se_ce <- com_d[com_d$method == "CE_z", c("component", "V_type", "effect", "se")]

names(se_ce)[4] <- "se_com_ce"

CMP <- merge(CMP, se_ce, by = c("component", "V_type", "effect"), sort = FALSE)

CMP$dev_std_own <- CMP$dev_std                     # deviation / SE(COM) of the same rma() model

CMP$dev_std <- CMP$dev / CMP$se_com_ce             # deviation / SE(COM, common effect)

# Information loss: SE inflation of the common-effect estimate (sampling error only)

ce_ratio <- CMP[CMP$method == "CE_z", c("condition", "missing_method", "component", "V_type", "effect", "se_ratio")]

names(ce_ratio)[6] <- "se_ratio_ce"

CMP <- merge(CMP, ce_ratio, by = c("condition", "missing_method", "component", "V_type", "effect"), sort = FALSE)

CMP$truth <- ifelse(CMP$component == "corrected", POP$effect_truth$corrected[CMP$effect],
                    POP$effect_truth$uncorrected[CMP$effect])

CMP$bias_truth <- CMP$est - CMP$truth

CMP$covered <- CMP$lo <= CMP$truth & CMP$truth <= CMP$hi

COM_TRUTH <- transform(com_d, truth = ifelse(component == "corrected", POP$effect_truth$corrected[effect],
                                             POP$effect_truth$uncorrected[effect]))

COM_TRUTH$bias_truth <- COM_TRUTH$est - COM_TRUTH$truth

com_i <- INDIRECT[INDIRECT$missing_method == "COM", ]

CMP_I <- merge(INDIRECT[INDIRECT$missing_method != "COM", ],
               setNames(com_i[, c("component", "V_type", "effect", "est", "lo", "hi", "mc_sig", "joint_sig")],
                        c("component", "V_type", "effect", "est_com", "lo_com", "hi_com", "mc_sig_com", "joint_sig_com")),
               by = c("component", "V_type", "effect"), sort = FALSE)

CMP_I$dev <- CMP_I$est - CMP_I$est_com

CMP_I$ci_ratio <- (CMP_I$hi - CMP_I$lo) / (CMP_I$hi_com - CMP_I$lo_com)

CMP_I$mc_sig_agree <- CMP_I$mc_sig == CMP_I$mc_sig_com

CMP_I$joint_sig_agree <- CMP_I$joint_sig == CMP_I$joint_sig_com

## Summary over the six direct effects (primary: rma() RE z; SE inflation also
## from the common-effect model, which reflects the information loss alone)

prim <- CMP[CMP$method == PRIMARY_RMA, ]

KEY <- do.call(rbind, lapply(split(prim, list(prim$condition, prim$missing_method, prim$component, prim$V_type),
                                   drop = TRUE), function(d) {
  data.frame(
    condition = d$condition[1], mechanism = d$mechanism[1], rate = d$rate[1],
    missing_method = d$missing_method[1], component = d$component[1], V_type = d$V_type[1],
    mean_abs_dev_std = mean(abs(d$dev_std)), rmsd_std = sqrt(mean(d$dev_std^2)),
    max_abs_dev_std = max(abs(d$dev_std)), mean_ci_ratio = mean(d$ci_ratio),
    mean_se_ratio_ce = mean(d$se_ratio_ce), n_sig_agree = sum(d$sig_agree),
    mean_abs_bias_truth = mean(abs(d$bias_truth)), n_covered = sum(d$covered),
    n_tau2_pos = sum(d$tau2 > 1e-8), n_tau2_pos_com = NA_integer_)
}))

KEY$n_tau2_pos_com <- vapply(seq_len(nrow(KEY)), function(i) sum(
  com_d$method == PRIMARY_RMA & com_d$component == KEY$component[i] & com_d$V_type == KEY$V_type[i] &
    com_d$tau2 > 1e-8), 0L)

KEY <- KEY[order(match(KEY$condition, conditions), match(KEY$missing_method, METHODS),
                 match(KEY$component, COMPONENTS), KEY$V_type), ]
rownames(KEY) <- NULL

## 4. Study-level summaries: deviation from COM, FMI, reliability ------------------------

STUDY <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  com <- get_sr(paste0(grid$study[i], "_complete"), "COM")
  do.call(rbind, lapply(setdiff(METHODS, "COM"), function(mt) {
    r <- get_sr(grid$dataset[i], mt)
    do.call(rbind, lapply(COMPONENTS, function(cmp) {
      se_c <- sqrt(diag(com[[cmp]]$V$ML))[PATHS$label]
      fmi <- if (!is.null(r[[cmp]]$mi)) r[[cmp]]$mi$ML$fmi else rep(NA_real_, 12)
      data.frame(dataset = grid$dataset[i], study = grid$study[i], condition = grid$condition[i],
                 mechanism = grid$mechanism[i], rate = grid$rate[i], missing_method = mt, component = cmp,
                 path = PATHS$label, effect = PATHS$effect, est = r[[cmp]]$est[PATHS$label],
                 se = sqrt(diag(r[[cmp]]$V$ML))[PATHS$label], est_com = com[[cmp]]$est[PATHS$label],
                 se_com = se_c, fmi = fmi, admissible = r[[cmp]]$admissible, row.names = NULL)
    }))
  }))
}))

STUDY$dev_std <- (STUDY$est - STUDY$est_com) / STUDY$se_com

STUDY$se_ratio <- STUDY$se / STUDY$se_com

REL <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  com <- get_sr(paste0(grid$study[i], "_complete"), "COM")
  do.call(rbind, lapply(setdiff(METHODS, "COM"), function(mt) data.frame(
    dataset = grid$dataset[i], study = grid$study[i], condition = grid$condition[i],
    missing_method = mt, variable = VARS, omega = get_sr(grid$dataset[i], mt)$rel[VARS],
    omega_com = com$rel[VARS], row.names = NULL)))
}))

REL$diff <- REL$omega - REL$omega_com

SECS <- do.call(rbind, lapply(jobs, function(j) {
  r <- get_sr(j$dataset, j$method)
  data.frame(dataset = j$dataset, method = j$method, seconds = r$seconds,
             warnings = paste(r$warnings %||% character(0), collapse = " | "))
}))

if (any(nzchar(SECS$warnings))) {
  cat("\nWarnings in the study-level analyses (see sheet 'runtime'):\n")
  print(table(SECS$method[nzchar(SECS$warnings)]))
}

## 5. Decision: mice-PMM vs two-stage FIML (single realisation)

DECISION <- NULL

if (all(c("mice_pmm", "FIML_TS") %in% METHODS)) {
  DECISION <- choose_method(list(key = KEY, methods = METHODS))
  DECISION$text <- c(DECISION$text, "",
                     "One realisation of the missing data per condition: 05_monte_carlo.R (pilot/full) gives the",
                     "decision that 08 and 09 use when its results are available.")
  writeLines(DECISION$text, file.path(OUT, "method_decision.txt"))
  cat("\n"); cat(DECISION$text, sep = "\n")
}

## 6. Reproduction of the sensitivity report's observed-data results
# The rma() code of this pipeline, applied to the input of the sensitivity
# analysis (APIM_meta_analysis_results.rds), must reproduce the observed pooled
# estimates and SEs of APIM_sensitivity_observed_critical.csv exactly. The power
# analysis under missing data (critical values, MDE) is in 07.

SENS <- list()

f_obs <- find_sensitivity_file("APIM_sensitivity_observed_critical.csv")

f_in  <- find_sensitivity_file("APIM_meta_analysis_results.rds")

if (!is.null(f_obs) && !is.null(f_in)) {
  obs <- utils::read.csv(f_obs, stringsAsFactors = FALSE, encoding = "UTF-8")
  rr <- readRDS(f_in)
  rep_tab <- meta_both(lapply(rr$study_results, function(s) s), "ML", n_mc = 1000)$direct
  rep_tab$method_lab <- RMA_METHODS[rep_tab$method]; rep_tab$effect_lab <- EFFECT_LAB[rep_tab$effect]
  rep_tab$component_lab <- COMP_LAB[rep_tab$component]
  mm <- merge(obs, rep_tab, by.x = c("component", "method", "effect"),
              by.y = c("component_lab", "method_lab", "effect_lab"), suffixes = c("_report", "_here"))
  SENS$reproduction <- data.frame(n_compared = nrow(mm),
                                  max_abs_diff_est = max(abs(mm$est_report - mm$est_here)),
                                  max_abs_diff_se = max(abs(mm$se_report - mm$se_here)))
  cat(sprintf("\nrma() pooling reproduces the sensitivity report's observed estimates: %d comparisons, max |diff| est = %.1e, SE = %.1e\n",
              nrow(mm), SENS$reproduction$max_abs_diff_est, SENS$reproduction$max_abs_diff_se))
  if (SENS$reproduction$max_abs_diff_est > 1e-6 || SENS$reproduction$max_abs_diff_se > 1e-6)
    warning("The rma() pooling does not reproduce the sensitivity report: check the inputs.")
} else message("Sensitivity-analysis outputs not found in ", paste(SENS_DIRS, collapse = ", "),
               ": the reproduction check is skipped.")

## 7. Export

E <- list(methods = METHODS, conditions = conditions, cond_info = cond_info, direct = DIRECT,
          indirect = INDIRECT, compare = CMP, compare_indirect = CMP_I, key = KEY, com_truth = COM_TRUTH,
          study_level = STUDY, reliability = REL, decision = DECISION, sensitivity = SENS,
          fiml_ts_check = dchk, seconds = SECS, population = POP,
          settings = list(rel_margin = REL_MARGIN, n_mc = N_MC, primary_rma = PRIMARY_RMA, primary_v = PRIMARY_V))

saveRDS(E, file.path(OUT, "04_evaluation_results.rds"))

sheets <- list(summary_direct = KEY, pooled_direct = DIRECT, compare_direct = CMP,
               pooled_indirect = INDIRECT, compare_indirect = CMP_I, complete_vs_truth = COM_TRUTH,
               study_level = STUDY, reliability = REL, fiml_ts_check = dchk, runtime = SECS)

if (!is.null(SENS$reproduction)) sheets$sensitivity_reproduction <- SENS$reproduction

if (!is.null(DECISION)) sheets$decision <- data.frame(text = DECISION$text)

write_xlsx_sheets(sheets, file.path(OUT, "04_evaluation.xlsx"))

cat("\nWritten:", file.path(OUT, c("04_evaluation_results.rds", "04_evaluation.xlsx", "method_decision.txt")),
    sep = "\n")

k <- KEY[KEY$V_type == PRIMARY_V, c("condition", "missing_method", "component", "rmsd_std", "mean_ci_ratio", "n_sig_agree")]

k$missing_method <- METHOD_SHORT[k$missing_method]

cat("\nPooled direct effects vs complete data (rma() RE z, V = ML): RMSD in complete-data SE units\n")

print(reshape(k[, c("condition", "missing_method", "component", "rmsd_std")],
              idvar = c("missing_method", "component"), timevar = "condition", direction = "wide"),
      digits = 2, row.names = FALSE)
