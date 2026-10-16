
## title: "APIM meta-analysis: individual-level missing and final pipeline"
## author: "Anonymized"
## date: "15-10-2026"

# =============================================================================
# 06_delta_sensitivity.R
# Delta-adjustment (tipping-point) sensitivity analysis for MNAR, with the
# meta-analysis as the target (van Buuren, 2018, sec. 3.8; Leacy, Floyd, Yates
# & White, 2017; Ratitch, O'Kelly & Tosiello, 2013).
#
# Imputations: delta-adjusted mice-PMM. The imputed values of every item are
# shifted by delta SDs (of the observed item values) in the hypothesised MNAR
# direction of its construct (MNAR_DIRECTION in 00_functions.R: nonrespondents
# have HIGHER X and M and LOWER Y than MAR predicts), stochastically rounded and
# kept in 1-5. The shift is applied inside the chains. delta = 0 is MAR and is
# identical to mice-PMM of 03/04 (same seeds); delta > 0 moves the imputations in
# the MNAR direction, delta < 0 in the opposite direction.
#
# For every condition and delta: both studies are imputed (m = 20), the APIMeM
# is fitted in lavaan to every completed dataset (functions of
# apim_meta_analysis.R; V = ML), pooled within each study with Rubin's rules,
# and the studies are pooled with rma() (random effects, REML, z test: the
# specification chosen from the sensitivity report).
#
#   Tipping points (what you would do with the real data): for every pooled
#   effect, the smallest |delta| at which (i) its significance at alpha = .05
#   differs from MAR (direct effects: RE z test; specific indirect effects:
#   joint significance), and (ii) its 95% CI no longer contains the MAR
#   estimate. Interpolated on the z statistic between grid points.
#   Calibration against the complete data (toy data only): the delta that brings
#   the pooled direct effects closest to the complete-data results. Under MNAR it
#   shows how far from MAR the imputations must move; under MCAR and MAR it
#   should be near 0 (negative controls).
#
# Output (output/sensitivity/)
#   raw/<dataset>__delta<d>.rds    study-level results per delta (cached)
#   06_delta_results.rds           pooled results, tipping points, calibration
#   06_delta_sensitivity.xlsx
# Runtime: ~1 min per study and delta (m = 20), ~2 h on one core for the 9
# conditions x 6 non-zero deltas (WORKERS=n runs jobs in parallel).
# SENS_MODE=quick: 2 conditions, deltas 0 and 0.5, m = 5 (smoke test).
# =============================================================================

source("00_functions.R")

## ---- 0. Settings --------------------------------------------------------------------

MODE <- Sys.getenv("SENS_MODE", "full")

CFG <- switch(MODE,
  quick = list(DELTAS = c(0, 0.5), M = 5, maxit = 5, conditions = c("MAR_15pct", "MNAR_15pct")),
  full  = list(DELTAS = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), M = 20, maxit = 20, conditions = NULL),
  stop("SENS_MODE must be quick or full"))

WORKERS <- as.integer(Sys.getenv("WORKERS", "15"))

SEED_03 <- 20261013                    # seeds of 03: delta = 0 reproduces mice-PMM

OUT <- file.path("output", "sensitivity", MODE); RAW <- file.path(OUT, "raw")

N_MC <- 2e4; MC_SEED <- 20261014

dir.create(RAW, FALSE, TRUE)

A <- readRDS(file.path("data", "amputed", "amputed_datasets.rds"))

grid <- A$grid

conditions <- CFG$conditions %||% unique(grid$condition)

study_label <- setNames(STUDIES$study, STUDIES$key)

imp_seed <- function(i) SEED_03 + 1000 * i + match("mice_pmm", IMPUTATION_METHODS)   # as in 03

E04 <- if (file.exists(f <- file.path("output", "evaluation", "04_evaluation_results.rds"))) readRDS(f) else NULL

pmm04 <- function(ds) file.path("output", "evaluation", "raw", sprintf("%s__mice_pmm.rds", ds))

## 1. Study-level analyses per delta (cached)

jobs <- expand.grid(delta = CFG$DELTAS, dataset = grid$dataset[grid$condition %in% conditions],
                    stringsAsFactors = FALSE)

raw_file <- function(ds, dl) file.path(RAW, sprintf("%s__delta%+.2f.rds", ds, dl))

# Cache key: amputed data, delta and the imputation settings, seed and generator

jobs$key <- vapply(seq_len(nrow(jobs)), function(r) {
  ds <- jobs$dataset[r]; i <- match(ds, grid$dataset)
  cache_key("delta", jobs$delta[r], CFG$M, CFG$maxit, 5L, PRED_SET, MNAR_DIRECTION, imp_seed(i),
            RNG_IMPUTATION, A$amputed[[ds]][, c(ID_VAR, ITEMS, AUX_VARS, IV_VARS)])
}, "")

# delta = 0 with the settings of 03 is mice-PMM of 03 (same seed and generator),
# so 04's results are reused when they were computed from these imputations

I03 <- if (file.exists(f <- file.path("data", "imputed", "imputations_all.rds"))) readRDS(f) else NULL

same_as_03 <- MODE == "full" && !is.null(I03) && identical(as.integer(I03$M), as.integer(CFG$M)) &&
  isTRUE(I03$settings$maxit_pmm == CFG$maxit) && isTRUE(I03$seed == SEED_03) &&
  isTRUE(I03$settings$donors == 5L) && identical(I03$settings$pred_set, PRED_SET)

res04 <- function(ds) {
  imps <- I03$imputations[[ds]]$mice_pmm
  if (!same_as_03 || is.null(imps)) return(NULL)
  r <- read_cached(pmm04(ds), cache_key("analysis", "mice_pmm", A$amputed[[ds]], imps))
  if (is.null(r)) return(NULL)
  for (cmp in COMPONENTS) r[[cmp]]$V <- r[[cmp]]$V["ML"]
  r$completed_means <- colMeans(do.call(rbind, lapply(imps, function(d) colMeans(d[, ITEMS]))))
  r
}

reuse <- jobs$delta == 0 & vapply(jobs$dataset, function(ds) !is.null(res04(ds)), TRUE)

todo <- which(!reuse & !vapply(seq_len(nrow(jobs)), function(r)
  !is.null(read_cached(raw_file(jobs$dataset[r], jobs$delta[r]), jobs$key[r])), TRUE))

cat(sprintf("Delta analysis (%s): %d conditions x %d deltas x 2 studies = %d jobs (%d to run), m = %d, workers = %d\n",
            MODE, length(conditions), length(CFG$DELTAS), nrow(jobs), length(todo), CFG$M, WORKERS))

d_job <- function(r) {
  # parallel workers receive the functions used here, but not the custom mice
  # method "pmmx", which mice looks up by name: load the library there
  if (!exists("mice.impute.pmmx", envir = globalenv())) sys.source("00_functions.R", envir = globalenv())
  ds <- jobs$dataset[r]; dl <- jobs$delta[r]; i <- match(ds, grid$dataset)
  t0 <- proc.time()[["elapsed"]]
  imps <- impute_delta(A$amputed[[ds]], dl, CFG$M, CFG$maxit, seed = imp_seed(i))
  res <- analyse_mi(imps, study_label[[grid$study[i]]], ml_only = TRUE)
  res$completed_means <- colMeans(do.call(rbind, lapply(imps, function(d) colMeans(d[, ITEMS]))))
  res$seconds <- proc.time()[["elapsed"]] - t0
  res$cache_key <- jobs$key[r]
  saveRDS(res, raw_file(ds, dl))
  res$seconds
}

t0 <- proc.time()[["elapsed"]]

if (length(todo)) invisible(par_lapply(todo, d_job, workers = WORKERS))

cat(sprintf("  done in %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))

get_res <- function(ds, dl) {
  r <- read_cached(raw_file(ds, dl), jobs$key[jobs$dataset == ds & jobs$delta == dl])
  if (is.null(r) && dl == 0) r <- res04(ds)                 # delta = 0: mice-PMM from 03/04
  if (is.null(r)) stop("Result missing or stale for ", ds, ", delta = ", dl)
  r
}

## 2. Pool the studies with rma() for every condition and delta

pooled <- list(); means <- list()

for (cond in conditions) for (dl in CFG$DELTAS) {
  dss <- sprintf("%s_%s", STUDIES$key, cond)
  sr <- setNames(lapply(dss, get_res, dl = dl), STUDIES$study)
  ma <- meta_both(sr, "ML", methods = c("RE_z", "CE_z"), n_mc = N_MC, mc_seed = MC_SEED)
  pooled[[paste(cond, dl)]] <- list(direct = cbind(condition = cond, delta = dl, ma$direct),
                                    indirect = cbind(condition = cond, delta = dl, ma$indirect))
  amp <- A$amputed[dss]; comp <- A$complete[STUDIES$key]
  means[[paste(cond, dl)]] <- do.call(rbind, lapply(seq_along(dss), function(s) {
    cm <- sr[[s]]$completed_means
    if (is.null(cm)) return(NULL)
    sdv <- apply(comp[[s]][, ITEMS], 2, stats::sd)
    data.frame(condition = cond, delta = dl, study = STUDIES$key[s], construct = ITEM_CONSTRUCT,
               shift_sd = (cm - colMeans(comp[[s]][, ITEMS])) / sdv * MNAR_DIRECTION[ITEM_CONSTRUCT])
  }))
}

D <- do.call(rbind, lapply(pooled, `[[`, "direct")); rownames(D) <- NULL

IND <- do.call(rbind, lapply(pooled, `[[`, "indirect")); rownames(IND) <- NULL

MS <- do.call(rbind, means)

MS <- if (!is.null(MS)) aggregate(shift_sd ~ condition + delta + construct, MS, mean) else NULL

## 3. Tipping points
# Interpolated tipping point on the z statistic between neighbouring grid values

tip_one <- function(dl, z, sig, est, lo, hi, side) {
  i0 <- which(dl == 0); ord <- if (side > 0) which(dl > 0) else rev(which(dl < 0))
  res <- c(sig_grid = NA, sig_interp = NA, ci_grid = NA)
  prev <- i0
  for (i in ord) {
    if (is.na(res["sig_grid"]) && sig[i] != sig[i0]) {
      res["sig_grid"] <- dl[i]
      zc <- if (sig[i0]) sign(z[i0]) * stats::qnorm(.975) else sign(z[i]) * stats::qnorm(.975)
      res["sig_interp"] <- if (is.finite(z[prev]) && z[i] != z[prev])
        dl[prev] + (zc - z[prev]) * (dl[i] - dl[prev]) / (z[i] - z[prev]) else dl[i]
    }
    if (is.na(res["ci_grid"]) && (est[i0] < lo[i] || est[i0] > hi[i])) res["ci_grid"] <- dl[i]
    prev <- i
  }
  res
}

dr <- D[D$method == PRIMARY_RMA, ]

TIPS <- do.call(rbind, lapply(split(dr, list(dr$condition, dr$component, dr$effect), drop = TRUE), function(d) {
  d <- d[order(d$delta), ]
  up <- tip_one(d$delta, d$est / d$se, d$p < .05, d$est, d$lo, d$hi, +1)
  dn <- tip_one(d$delta, d$est / d$se, d$p < .05, d$est, d$lo, d$hi, -1)
  i0 <- which(d$delta == 0)
  data.frame(condition = d$condition[1], component = d$component[1], effect = d$effect[1], type = "direct",
             est_MAR = d$est[i0], p_MAR = d$p[i0], significant_MAR = d$p[i0] < .05,
             sig_tip_mnar_direction = up[["sig_grid"]], sig_tip_mnar_interp = up[["sig_interp"]],
             sig_tip_opposite = dn[["sig_grid"]], ci_tip_mnar_direction = up[["ci_grid"]],
             ci_tip_opposite = dn[["ci_grid"]])
}))
ir <- IND[IND$short %in% names(SPECIFIC_INDIRECT), ]
TIPS_I <- do.call(rbind, lapply(split(ir, list(ir$condition, ir$component, ir$short), drop = TRUE), function(d) {
  d <- d[order(d$delta), ]; i0 <- which(d$delta == 0)
  first <- function(idx, cond) { x <- idx[cond[idx]]; if (length(x)) d$delta[x[1]] else NA }
  up <- which(d$delta > 0); dn <- rev(which(d$delta < 0))
  data.frame(condition = d$condition[1], component = d$component[1], effect = d$short[1], type = "indirect",
             est_MAR = d$est[i0], joint_sig_MAR = d$joint_sig[i0], mc_sig_MAR = d$mc_sig[i0],
             sig_tip_mnar_direction = first(up, d$joint_sig != d$joint_sig[i0]),
             sig_tip_opposite = first(dn, d$joint_sig != d$joint_sig[i0]),
             ci_tip_mnar_direction = first(up, d$est[i0] < d$lo | d$est[i0] > d$hi),
             ci_tip_opposite = first(dn, d$est[i0] < d$lo | d$est[i0] > d$hi))
}))

TIPS <- TIPS[order(match(TIPS$condition, conditions), match(TIPS$component, COMPONENTS),
                   match(TIPS$effect, EFFECTS)), ]

TIPS_I <- TIPS_I[order(match(TIPS_I$condition, conditions), match(TIPS_I$component, COMPONENTS),
                       match(TIPS_I$effect, names(SPECIFIC_INDIRECT))), ]

rownames(TIPS) <- rownames(TIPS_I) <- NULL

## 4. Calibration against the complete data (toy data only)

CAL <- NULL

if (!is.null(E04)) {
  com <- E04$direct[E04$direct$missing_method == "COM" & E04$direct$V_type == "ML", ]
  ce  <- com[com$method == "CE_z", c("component", "effect", "se")]; names(ce)[3] <- "se_com_ce"
  rc  <- com[com$method == PRIMARY_RMA, c("component", "effect", "est")]; names(rc)[3] <- "est_com"
  x <- merge(merge(dr, rc, by = c("component", "effect")), ce, by = c("component", "effect"))
  x$dev_std <- (x$est - x$est_com) / x$se_com_ce
  CURVE <- aggregate(cbind(rmsd_std = dev_std^2) ~ condition + component + delta, x, mean)
  CURVE$rmsd_std <- sqrt(CURVE$rmsd_std)
  CAL <- do.call(rbind, lapply(split(CURVE, list(CURVE$condition, CURVE$component), drop = TRUE), function(d)
    data.frame(condition = d$condition[1], component = d$component[1],
               best_delta = d$delta[which.min(d$rmsd_std)], rmsd_at_best = min(d$rmsd_std),
               rmsd_at_MAR = d$rmsd_std[d$delta == 0])))
  CAL <- CAL[order(match(CAL$condition, conditions), match(CAL$component, COMPONENTS)), ]
  rownames(CAL) <- NULL
} else CURVE <- NULL

# Calibration on the item means: the delta at which the completed item means
# equal the complete-data means (linear interpolation of the mean shift, in item
# SDs signed in the MNAR direction, between the delta grid points). Under MNAR it
# says how far the nonrespondents' values lie from what MAR predicts; under MCAR
# and MAR it should be close to 0. Less noisy than the RMSD-based calibration.

CAL_MEANS <- NULL

if (!is.null(MS)) {
  cross0 <- function(dl, y) {
    o <- order(dl); dl <- dl[o]; y <- y[o]
    i <- which(diff(sign(y)) != 0)[1]
    if (is.na(i)) return(NA_real_)
    dl[i] - y[i] * (dl[i + 1] - dl[i]) / (y[i + 1] - y[i])
  }
  CAL_MEANS <- do.call(rbind, lapply(split(MS, list(MS$condition, MS$construct), drop = TRUE), function(d)
    data.frame(condition = d$condition[1], construct = d$construct[1],
               shift_at_MAR = d$shift_sd[d$delta == 0], delta_matching_means = cross0(d$delta, d$shift_sd))))
  CAL_MEANS <- CAL_MEANS[order(match(CAL_MEANS$condition, conditions), CAL_MEANS$construct), ]
  rownames(CAL_MEANS) <- NULL
}

## 5. Link to the sensitivity report: pooled estimate vs critical |beta|

f_obs <- find_sensitivity_file("APIM_sensitivity_observed_critical.csv")

CRIT <- NULL

if (!is.null(f_obs)) {
  obs <- utils::read.csv(f_obs, stringsAsFactors = FALSE, encoding = "UTF-8")
  obs <- obs[obs$method == RMA_METHODS[[PRIMARY_RMA]], c("component", "effect", "crit_est")]
  names(obs) <- c("component_lab", "effect_lab", "crit_report")
  CRIT <- dr[, c("condition", "delta", "component", "effect", "est", "se", "p")]
  CRIT$component_lab <- COMP_LAB[CRIT$component]; CRIT$effect_lab <- EFFECT_LAB[CRIT$effect]
  CRIT <- merge(CRIT, obs, by = c("component_lab", "effect_lab"), sort = FALSE)
  CRIT$crit_missing <- stats::qnorm(.975) * CRIT$se
  CRIT$exceeds_crit_report <- abs(CRIT$est) > CRIT$crit_report
  CRIT$exceeds_crit_missing <- abs(CRIT$est) > CRIT$crit_missing
  CRIT <- CRIT[order(match(CRIT$condition, conditions), CRIT$delta, match(CRIT$component, COMPONENTS),
                     match(CRIT$effect, EFFECTS)), setdiff(names(CRIT), c("component_lab", "effect_lab"))]
  rownames(CRIT) <- NULL
}

## 6. Export

S <- list(mode = MODE, cfg = CFG, conditions = conditions, direct = D, indirect = IND, tips = TIPS,
          tips_indirect = TIPS_I, calibration = CAL, calibration_curve = CURVE, completed_shift = MS,
          calibration_means = CAL_MEANS,
          critical = CRIT, mnar_direction = MNAR_DIRECTION)

saveRDS(S, file.path(OUT, "06_delta_results.rds"))

sheets <- list(settings = data.frame(
  setting = c("mode", "deltas (item-SD units; + = hypothesised MNAR direction)", "MNAR direction",
              "m", "maxit", "imputation", "pooling"),
  value = c(MODE, paste(CFG$DELTAS, collapse = ", "),
            paste(sprintf("%s %+d", names(MNAR_DIRECTION), MNAR_DIRECTION), collapse = "; "), CFG$M, CFG$maxit,
            "delta-adjusted mice-PMM (pmmx), item level, parcel predictors + auxiliaries",
            "lavaan APIMeM per completed dataset; Rubin within study; rma() RE (REML, z) across studies")),
  tipping_direct = TIPS, tipping_indirect = TIPS_I, pooled_direct = D, pooled_indirect = IND)

if (!is.null(CAL)) { sheets$calibration <- CAL; sheets$calibration_curve <- CURVE }

if (!is.null(MS)) sheets$completed_mean_shift <- MS

if (!is.null(CAL_MEANS)) sheets$calibration_means <- CAL_MEANS

if (!is.null(CRIT)) sheets$vs_critical_values <- CRIT

write_xlsx_sheets(sheets, file.path(OUT, "06_delta_sensitivity.xlsx"))

cat("Written:", file.path(OUT, c("06_delta_results.rds", "06_delta_sensitivity.xlsx")), sep = "\n")

cat("\nTipping points of the effects significant under MAR (delta, item-SD units, MNAR direction):\n")

print(TIPS[TIPS$significant_MAR, c("condition", "component", "effect", "est_MAR",
                                   "sig_tip_mnar_direction", "sig_tip_mnar_interp", "ci_tip_mnar_direction")],
      digits = 3, row.names = FALSE)

if (!is.null(CAL)) { cat("\nCalibration: delta closest to the complete data\n"); print(CAL, digits = 3, row.names = FALSE) }

if (!is.null(CAL_MEANS)) { cat("\nCalibration on the item means: delta at which the completed means equal the complete data\n")
  print(CAL_MEANS, digits = 3, row.names = FALSE) }
