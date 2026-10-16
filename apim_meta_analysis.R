###############################################################################
##  apim_meta_analysis.R
##
##  Meta-analysis of Actor-Partner Interdependence Mediation Model (APIMeM)
##  paths across dyadic studies, in two components:
##
##    Component 1  UNCORRECTED  - APIMeM fitted to observed scale scores
##    Component 2  CORRECTED    - APIMeM fitted to single-indicator latent
##                                variables whose error variances are fixed
##                                at (1 - reliability) * variance, i.e. the
##                                paths are corrected for unreliability in
##                                X, M and Y of BOTH partners simultaneously
##
##  Pipeline
##    1. Per study: scale scores + reliability (omega_total or alpha)
##    2. Per study: APIMeM in lavaan -> 12 standardized direct paths and their
##       FULL 12 x 12 sampling covariance matrix (delta method), under
##         - ML  : conventional (normal-theory) sampling covariances
##         - MLR : robust (Huber-White sandwich) sampling covariances
##         - BOOT: optional nonparametric dyad bootstrap; for the corrected
##                 component this also propagates the sampling error of the
##                 reliability coefficients themselves
##    3. Pooling with metafor::rma.mv using the block-diagonal V matrix:
##         - multilevel model (primary): the two partner-specific estimates of
##           each effect type are nested within studies
##         - common-effect model (sensitivity)
##       + cluster-robust CR2 inference with Satterthwaite df (clubSandwich);
##         reported in the tables only when df >= CFG$cr2_min_df
##    4. Indirect and total effects from the pooled paths with Monte Carlo CIs
##       (using the pooled covariance matrix of the paths)
##    5. APA-style tables (flextable + officer) -> one .docx
##
##  Why the reliability correction is NOT done on the pooled coefficients
##    APIM paths are PARTIAL regression coefficients. The classical
##    disattenuation r / sqrt(rxx * ryy) (and its variance v / (rxx * ryy))
##    is only valid for a bivariate correlation. Unreliability in one
##    predictor also biases the coefficients of the other predictors, so the
##    correction must be made on the whole model. The single-indicator
##    latent APIMeM does exactly this (it is equivalent to disattenuating the
##    6 x 6 correlation matrix and refitting the saturated APIMeM), and lavaan
##    returns the correct sampling (co)variances of the corrected paths.
##
##  Why metafor/clubSandwich rather than meta
##    meta has no arguments for artifact (reliability) correction, and its
##    three-level option (metagen(cluster = , rho = )) calls metafor::rma.mv
##    internally with ONE assumed within-cluster correlation. Here the exact
##    within-study covariance matrix of the 12 paths is known from lavaan, so
##    it is passed to rma.mv directly.
##
##  Data format (one file per study; wide, one row per dyad)
##    Items named <construct>_<member>_<item>, e.g. x_p1_1 ... y_p2_4
##    (patterns are configurable below). If a study only has scale scores,
##    name them <construct>_<member> (e.g. m_p2) and supply reliabilities
##    via CFG$reliability_supplied.
##
##  Path notation (1 = Partner 1, 2 = Partner 2; "12" = from P1 to P2)
##    a1  X1->M1  actor a    a2  X2->M2  actor a
##    a12 X1->M2  partner a  a21 X2->M1  partner a
##    b1  M1->Y1  actor b    b2  M2->Y2  actor b
##    b12 M1->Y2  partner b  b21 M2->Y1  partner b
##    c1  X1->Y1  actor c'   c2  X2->Y2  actor c'
##    c12 X1->Y2  partner c' c21 X2->Y1  partner c'
##
##  Author: <you>          R >= 4.2
###############################################################################


## =============================================================================
## 0. PACKAGES
## =============================================================================
pkgs <- c("lavaan", "metafor", "clubSandwich", "MASS", "flextable", "officer")
miss_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss_pkgs)) {
  stop("Please install: ", paste(miss_pkgs, collapse = ", "),
       "\n  install.packages(c(\"", paste(miss_pkgs, collapse = "\", \""), "\"))")
}
suppressPackageStartupMessages({
  library(lavaan)
  library(metafor)
  library(clubSandwich)
  library(flextable)
  library(officer)
})


## =============================================================================
## 1. CONFIGURATION  (the only section you should need to edit)
## =============================================================================
CFG <- list(
  ## Study files (add rows to add studies; the script is general in k)
  studies = data.frame(
    study = c("Study 1", "Study 2"),
    file  = c("toy_study1.csv", "toy_study2.csv"),
    stringsAsFactors = FALSE
  ),

  ## Construct prefixes and member codes used in the column names
  constructs    = c(X = "x", M = "m", Y = "y"),
  members       = c("p1", "p2"),
  member_labels = c("Partner 1", "Partner 2"),
  construct_labels = c(X = "Predictor (X)", M = "Mediator (M)", Y = "Outcome (Y)"),
  item_regex    = "^%s_%s_[0-9]+$",     # sprintf(item_regex, construct, member)

  ## Scale scores: mean of available items if >= this proportion is observed
  min_prop_items = 0.75,

  ## Reliability coefficient computed from the items: "omega" (omega_total
  ## from a one-factor congeneric CFA; recommended) or "alpha"
  reliability = "omega",
  ## Optional override, e.g. published coefficients or scale-score-only data:
  ## data.frame(study = "Study 1", variable = "X1", rel = .81)
  ## variable must be one of X1, X2, M1, M2, Y1, Y2
  reliability_supplied = NULL,

  ## Pooling
  pool_roles    = TRUE,          # TRUE: pool Partner 1 and Partner 2 estimates
                                 #   of each effect type (multilevel: roles
                                 #   nested in studies). FALSE: role-specific
                                 #   multivariate model with 12 pooled paths.
  re_struct     = "ID",          # between-study structure across effect types
                                 #   ("ID" = one tau^2; "DIAG" = one per type)
  primary_V     = "ML",          # "ML", "MLR" or "BOOT" (if boot_V = TRUE)
  primary_model = "multilevel",  # "multilevel" or "common"

  ## Optional bootstrap sampling covariance matrices (slow: ~8 lavaan fits
  ## per resample per study)
  boot_V = FALSE,
  n_boot = 1000,

  ## Cluster-robust (CR2) results are always computed and saved, but only
  ## shown in Table 2 when every Satterthwaite df >= this value (Tipton, 2015)
  cr2_min_df = 4,

  ## Monte Carlo CIs for indirect/total effects
  n_mc  = 1e5,
  level = 0.95,
  seed  = 20260911,

  ## Output
  out_docx = "APIM_meta_analysis_tables.docx",
  out_csv  = "APIM_meta_analysis_results.csv",
  out_rds  = "APIM_meta_analysis_results.rds",
  table_font = "Times New Roman",
  table_size = 10,
  max_table_width = 8.9            # inches (landscape, 1-inch margins)
)

## Optional: define CFG_OVERRIDE <- list(...) before source()-ing this script
## to change settings without editing the file (used for batch/sensitivity runs)
if (exists("CFG_OVERRIDE", inherits = TRUE)) {
  for (nm in names(CFG_OVERRIDE)) CFG[[nm]] <- CFG_OVERRIDE[[nm]]
}
stopifnot(CFG$reliability %in% c("omega", "alpha"),
          CFG$primary_V %in% c("ML", "MLR", "BOOT"),
          CFG$primary_model %in% c("multilevel", "common"))

set.seed(CFG$seed)


## =============================================================================
## 2. MODEL DICTIONARY
## =============================================================================
VARS <- c("X1", "X2", "M1", "M2", "Y1", "Y2")

PATHS <- data.frame(
  label = c("a1", "a2", "a12", "a21", "b1", "b2", "b12", "b21",
            "c1", "c2", "c12", "c21"),
  lhs   = c("M1", "M2", "M2", "M1", "Y1", "Y2", "Y2", "Y1",
            "Y1", "Y2", "Y2", "Y1"),
  rhs   = c("X1", "X2", "X1", "X2", "M1", "M2", "M1", "M2",
            "X1", "X2", "X1", "X2"),
  path  = rep(c("a", "b", "c"), each = 4),
  type  = rep(c("actor", "actor", "partner", "partner"), 3),
  stringsAsFactors = FALSE
)
PATHS$effect <- paste(PATHS$type, PATHS$path, sep = "_")      # e.g. actor_a
PATHS$from   <- substr(PATHS$rhs, 2, 2)
PATHS$to     <- substr(PATHS$lhs, 2, 2)
PATHS$role   <- sprintf("P%s\u2192P%s", PATHS$from, PATHS$to)  # e.g. P1->P2

EFFECTS <- c("actor_a", "partner_a", "actor_b", "partner_b",
             "actor_c", "partner_c")

PATH_GROUP_LABEL <- c(
  a = "a path (X \u2192 M)",
  b = "b path (M \u2192 Y)",
  c = "c\u2032 path (X \u2192 Y, direct)"
)

## APIMeM syntax. prefix = "" (observed) or "L" (latent, corrected model)
apim_syntax <- function(prefix = "") {
  v <- function(x) paste0(prefix, x)
  regs <- vapply(unique(PATHS$lhs), function(y) {
    ix <- which(PATHS$lhs == y)
    paste0(v(y), " ~ ", paste0(PATHS$label[ix], "*", v(PATHS$rhs[ix]),
                               collapse = " + "))
  }, character(1))
  c(regs,
    paste(v("X1"), "~~", v("X2")),     # partners' predictors
    paste(v("M1"), "~~", v("M2")),     # partners' mediator residuals
    paste(v("Y1"), "~~", v("Y2")))     # partners' outcome residuals
}

## Indirect and total effects, as expressions of the pooled coefficients.
## Decomposition follows Ledermann, Macho & Kenny (2011).
indirect_defs <- function(pool_roles) {
  if (pool_roles) {
    data.frame(
      group = c(rep("Specific indirect effects", 4),
                rep("Total indirect effects", 2),
                rep("Total effects", 2)),
      label = c("Actor\u2013actor (own X \u2192 own M \u2192 own Y)",
                "Partner\u2013partner (own X \u2192 partner's M \u2192 own Y)",
                "Actor\u2013partner (own X \u2192 own M \u2192 partner's Y)",
                "Partner\u2013actor (own X \u2192 partner's M \u2192 partner's Y)",
                "Actor (own X \u2192 own Y)",
                "Partner (own X \u2192 partner's Y)",
                "Actor (direct + indirect)",
                "Partner (direct + indirect)"),
      short = c("Actor\u2013actor", "Partner\u2013partner", "Actor\u2013partner",
                "Partner\u2013actor", "Total indirect, actor",
                "Total indirect, partner", "Total effect, actor",
                "Total effect, partner"),
      expr  = c("actor_a * actor_b",
                "partner_a * partner_b",
                "actor_a * partner_b",
                "partner_a * actor_b",
                "actor_a * actor_b + partner_a * partner_b",
                "actor_a * partner_b + partner_a * actor_b",
                "actor_c + actor_a * actor_b + partner_a * partner_b",
                "partner_c + actor_a * partner_b + partner_a * actor_b"),
      stringsAsFactors = FALSE)
  } else {
    data.frame(
      group = rep(c("X1 \u2192 Y1 (actor, Partner 1)",
                    "X2 \u2192 Y2 (actor, Partner 2)",
                    "X1 \u2192 Y2 (partner, Partner 1 \u2192 2)",
                    "X2 \u2192 Y1 (partner, Partner 2 \u2192 1)"), each = 4),
      label = rep(c("via own M", "via partner's M", "Total indirect",
                    "Total effect"), 4),
      short = paste(rep(c("X1\u2192Y1", "X2\u2192Y2", "X1\u2192Y2", "X2\u2192Y1"), each = 4),
                    rep(c("via own M", "via partner's M", "total indirect",
                          "total effect"), 4), sep = ", "),
      expr  = c("a1 * b1",   "a12 * b21", "a1 * b1 + a12 * b21",
                "c1 + a1 * b1 + a12 * b21",
                "a2 * b2",   "a21 * b12", "a2 * b2 + a21 * b12",
                "c2 + a2 * b2 + a21 * b12",
                "a1 * b12",  "a12 * b2",  "a1 * b12 + a12 * b2",
                "c12 + a1 * b12 + a12 * b2",
                "a2 * b21",  "a21 * b1",  "a2 * b21 + a21 * b1",
                "c21 + a2 * b21 + a21 * b1"),
      stringsAsFactors = FALSE)
  }
}


## =============================================================================
## 3. STUDY-LEVEL FUNCTIONS
## =============================================================================

## ---- 3a. Reliability ---------------------------------------------------------
## omega_total for a unit-weighted composite: (sum lambda)^2 /
## ((sum lambda)^2 + sum theta), from a one-factor congeneric CFA (FIML).
## alpha is a lower bound under congeneric items and therefore OVER-corrects.
reliability_coef <- function(items, type = c("omega", "alpha")) {
  type  <- match.arg(type)
  items <- as.data.frame(items)
  k     <- ncol(items)
  alpha <- function() {
    S <- stats::cov(items, use = "pairwise.complete.obs")
    (k / (k - 1)) * (1 - sum(diag(S)) / sum(S))
  }
  if (type == "alpha" || k < 3) {
    if (type == "omega") message("  < 3 items: omega not identified, using alpha")
    return(alpha())
  }
  names(items) <- paste0("i", seq_len(k))
  mod <- paste("F =~", paste(names(items), collapse = " + "))
  fit <- tryCatch(
    suppressWarnings(lavaan::cfa(mod, data = items, std.lv = TRUE,
                                 missing = "ml")),
    error = function(e) NULL)
  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
    message("  omega CFA failed; using alpha")
    return(alpha())
  }
  pe  <- lavaan::parameterEstimates(fit)
  lam <- pe$est[pe$op == "=~"]
  th  <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% names(items)]
  sum(lam)^2 / (sum(lam)^2 + sum(th))
}

## ---- 3b. Scale scores + reliabilities for one study --------------------------
prep_study <- function(raw, study, cfg, quiet = FALSE) {
  comp <- as.data.frame(matrix(NA_real_, nrow(raw), 6,
                               dimnames = list(NULL, VARS)))
  rel  <- setNames(rep(NA_real_, 6), VARS)
  n_items <- setNames(rep(NA_integer_, 6), VARS)

  for (v in VARS) {
    cons <- cfg$constructs[[substr(v, 1, 1)]]
    mem  <- cfg$members[as.integer(substr(v, 2, 2))]
    item_cols <- grep(sprintf(cfg$item_regex, cons, mem), names(raw), value = TRUE)
    score_col <- paste0(cons, "_", mem)

    if (length(item_cols) >= 2) {
      it <- as.data.frame(lapply(raw[item_cols], as.numeric))
      ok <- rowMeans(!is.na(it)) >= cfg$min_prop_items
      sc <- rowMeans(it, na.rm = TRUE)
      sc[!ok] <- NA
      comp[[v]]  <- sc
      rel[v]     <- reliability_coef(it, cfg$reliability)
      n_items[v] <- length(item_cols)
    } else if (score_col %in% names(raw)) {
      comp[[v]] <- as.numeric(raw[[score_col]])
    } else {
      stop(study, ": no items matching '", sprintf(cfg$item_regex, cons, mem),
           "' and no scale-score column '", score_col, "'")
    }
  }

  ## Supplied reliabilities override computed ones
  if (!is.null(cfg$reliability_supplied)) {
    rs <- cfg$reliability_supplied[cfg$reliability_supplied$study == study, ]
    if (nrow(rs)) rel[rs$variable] <- rs$rel
  }
  if (anyNA(rel)) stop(study, ": reliability missing for ",
                       paste(names(rel)[is.na(rel)], collapse = ", "),
                       " (supply via CFG$reliability_supplied)")
  if (any(rel <= 0 | rel > 1)) stop(study, ": reliabilities must be in (0, 1]")

  list(study = study, comp = comp, rel = rel, n_items = n_items,
       n = nrow(comp))
}

## ---- 3c. Fit the APIMeM (uncorrected or corrected) ---------------------------
## For the corrected model the error variance of each scale score is fixed at
## (1 - rel) * s^2, with s^2 the (FIML) estimate from the saturated
## uncorrected model, so both components use the same variance estimate.
fit_apim <- function(comp, rel = NULL, corrected = FALSE, estimator = "ML",
                     s2 = NULL) {
  miss <- if (anyNA(comp)) "ml" else "listwise"
  if (!corrected) {
    syntax <- apim_syntax("")
  } else {
    if (is.null(s2)) stop("s2 (observed variances) required for corrected model")
    theta  <- (1 - rel[VARS]) * s2[VARS]
    syntax <- c(paste0("L", VARS, " =~ 1*", VARS),
                paste0(VARS, " ~~ ", sprintf("%.12f", theta), "*", VARS),
                apim_syntax("L"))
  }
  fit <- lavaan::sem(paste(syntax, collapse = "\n"), data = comp,
                     estimator = estimator, missing = miss,
                     fixed.x = FALSE, meanstructure = TRUE)
  if (!lavaan::lavInspect(fit, "converged")) stop("APIMeM did not converge")
  fit
}

## Standardized paths + full sampling covariance matrix of the 12 paths
extract_paths <- function(fit) {
  ss  <- lavaan::standardizedSolution(fit, type = "std.all")
  est <- setNames(ss$est.std[match(PATHS$label, ss$label)], PATHS$label)
  V   <- lavaan::lavInspect(fit, "vcov.std.all")[PATHS$label, PATHS$label]
  admissible <- isTRUE(suppressWarnings(lavaan::lavInspect(fit, "post.check")))
  list(est = est, V = V, admissible = admissible)
}

## ---- 3d. One study: both components, ML and MLR sampling covariances --------
analyse_study <- function(ps) {
  f_u_ml  <- fit_apim(ps$comp, estimator = "ML")
  f_u_mlr <- fit_apim(ps$comp, estimator = "MLR")
  s2      <- diag(lavaan::lavInspect(f_u_ml, "cov.ov"))[VARS]
  f_c_ml  <- fit_apim(ps$comp, ps$rel, corrected = TRUE, estimator = "ML",  s2 = s2)
  f_c_mlr <- fit_apim(ps$comp, ps$rel, corrected = TRUE, estimator = "MLR", s2 = s2)

  pack <- function(f_ml, f_mlr) {
    e_ml  <- extract_paths(f_ml)
    e_mlr <- extract_paths(f_mlr)
    stopifnot(max(abs(e_ml$est - e_mlr$est)) < 1e-6)  # same point estimates
    list(est = e_ml$est, V = list(ML = e_ml$V, MLR = e_mlr$V),
         admissible = e_ml$admissible && e_mlr$admissible)
  }
  out <- list(uncorrected = pack(f_u_ml, f_u_mlr),
              corrected   = pack(f_c_ml, f_c_mlr),
              s2 = s2)
  if (!out$corrected$admissible)
    warning(ps$study, ": corrected APIMeM has an inadmissible solution ",
            "(e.g. negative residual variance). Check reliabilities.")
  out
}

## ---- 3e. Optional nonparametric dyad bootstrap of the 12 paths --------------
bootstrap_V <- function(raw, study, cfg, B) {
  draws <- list(uncorrected = matrix(NA_real_, B, 12, dimnames = list(NULL, PATHS$label)),
                corrected   = matrix(NA_real_, B, 12, dimnames = list(NULL, PATHS$label)))
  for (b in seq_len(B)) {
    idx <- sample.int(nrow(raw), replace = TRUE)
    res <- tryCatch({
      ps <- suppressMessages(prep_study(raw[idx, , drop = FALSE], study, cfg))
      fu <- fit_apim(ps$comp, estimator = "ML")
      s2 <- diag(lavaan::lavInspect(fu, "cov.ov"))[VARS]
      fc <- fit_apim(ps$comp, ps$rel, corrected = TRUE, estimator = "ML", s2 = s2)
      ec <- extract_paths(fc)
      list(u = extract_paths(fu)$est, c = if (ec$admissible) ec$est else NULL)
    }, error = function(e) NULL)
    if (!is.null(res)) {
      draws$uncorrected[b, ] <- res$u
      if (!is.null(res$c)) draws$corrected[b, ] <- res$c
    }
    if (b %% 100 == 0) message("    ", study, ": bootstrap ", b, "/", B)
  }
  lapply(draws, function(d) {
    ok <- stats::complete.cases(d)
    if (mean(ok) < .90) warning(study, ": ", round(100 * (1 - mean(ok)), 1),
                                "% of bootstrap resamples failed/inadmissible")
    stats::cov(d[ok, , drop = FALSE])
  })
}


## =============================================================================
## 4. META-ANALYTIC FUNCTIONS
## =============================================================================

## Long effect-size data + block-diagonal V for one component and V type
build_ma_data <- function(study_res, component, vtype) {
  dat <- do.call(rbind, lapply(names(study_res), function(s) {
    r <- study_res[[s]][[component]]
    data.frame(study = s, label = PATHS$label, effect = PATHS$effect,
               path = PATHS$path, type = PATHS$type, role = PATHS$role,
               yi = unname(r$est[PATHS$label]), stringsAsFactors = FALSE)
  }))
  dat$es_id <- seq_len(nrow(dat))
  V <- metafor::bldiag(lapply(names(study_res),
                              function(s) study_res[[s]][[component]]$V[[vtype]]))
  dimnames(V) <- list(dat$es_id, dat$es_id)
  list(dat = dat, V = V)
}

## Fit the pooling model
fit_ma <- function(md, model, cfg) {
  dat <- md$dat
  lev <- if (cfg$pool_roles) EFFECTS else PATHS$label
  dat$term <- factor(if (cfg$pool_roles) dat$effect else dat$label, levels = lev)
  if (model == "common") {
    fit <- metafor::rma.mv(yi, md$V, mods = ~ 0 + term, data = dat)
  } else if (cfg$pool_roles) {
    ## Level 3: study-by-effect-type deviations; level 2: partner-specific
    ## (role) deviations within study; level 1: sampling error with full V
    fit <- metafor::rma.mv(yi, md$V, mods = ~ 0 + term,
                           random = list(~ term | study, ~ 1 | es_id),
                           struct = cfg$re_struct, data = dat, method = "REML")
  } else {
    fit <- metafor::rma.mv(yi, md$V, mods = ~ 0 + term,
                           random = ~ term | study, struct = cfg$re_struct,
                           data = dat, method = "REML")
  }
  fit$dat_used <- dat
  fit
}

strip_term <- function(x) sub("^term", "", x)

## Model-based coefficient table
coef_table <- function(fit) {
  data.frame(term = strip_term(rownames(fit$beta)), est = fit$beta[, 1],
             se = fit$se, lo = fit$ci.lb, hi = fit$ci.ub, p = fit$pval,
             row.names = NULL)
}

## CR2 cluster-robust inference (clustered by study)
robust_table <- function(fit, level) {
  cl <- fit$dat_used$study
  tryCatch({
    vc <- clubSandwich::vcovCR(fit, cluster = cl, type = "CR2")
    ct <- clubSandwich::coef_test(fit, vcov = vc, test = "Satterthwaite")
    ci <- clubSandwich::conf_int(fit, vcov = vc, level = level,
                                 test = "Satterthwaite")
    list(tab = data.frame(term = strip_term(ct$Coef), se_cr2 = ct$SE,
                          df_cr2 = ct$df_Satt, p_cr2 = ct$p_Satt,
                          lo_cr2 = ci$CI_L, hi_cr2 = ci$CI_U, row.names = NULL),
         vcov = as.matrix(vc))
  }, error = function(e) {
    warning("CR2 failed: ", conditionMessage(e)); NULL
  })
}

## Heterogeneity: Q_E and multilevel I^2 (Cheung, 2014) using the typical
## sampling variance generalized to a full V (Viechtbauer's formula)
heterogeneity <- function(fit) {
  X <- fit$X
  W <- solve(fit$V)
  P <- W - W %*% X %*% solve(crossprod(X, W %*% X)) %*% crossprod(X, W)
  v_typ <- (fit$k - fit$p) / sum(diag(P))
  comps <- c(tau2_study = if (length(fit$tau2) && !all(is.na(fit$tau2)))
                            mean(fit$tau2) else NA_real_,
             sigma2_role = if (length(fit$sigma2)) fit$sigma2[1] else NA_real_)
  tot <- sum(comps, na.rm = TRUE)
  data.frame(k = fit$k, QE = fit$QE, QE_df = fit$k - fit$p, QE_p = fit$QEp,
             tau2_study = comps[["tau2_study"]],
             sigma2_role = comps[["sigma2_role"]],
             I2_study = 100 * comps[["tau2_study"]] / (tot + v_typ),
             I2_role  = 100 * comps[["sigma2_role"]] / (tot + v_typ),
             v_typical = v_typ)
}

## Monte Carlo CIs for functions (products/sums) of the pooled coefficients
mc_functions <- function(b, Vb, defs, n_mc, level) {
  names(b) <- strip_term(names(b))
  dimnames(Vb) <- list(names(b), names(b))
  draws <- tryCatch(MASS::mvrnorm(n_mc, mu = b, Sigma = Vb),
                    error = function(e) NULL)
  a <- (1 - level) / 2
  out <- lapply(seq_len(nrow(defs)), function(i) {
    ex  <- parse(text = defs$expr[i])
    est <- eval(ex, as.list(b))
    if (is.null(draws)) return(c(est = est, lo = NA, hi = NA))
    sim <- eval(ex, as.data.frame(draws))
    c(est = est, lo = unname(stats::quantile(sim, a)),
      hi = unname(stats::quantile(sim, 1 - a)))
  })
  cbind(defs, as.data.frame(do.call(rbind, out)))
}

## Run everything for one component x V type x model
run_ma <- function(study_res, component, vtype, model, cfg, do_mc = TRUE) {
  md  <- build_ma_data(study_res, component, vtype)
  fit <- fit_ma(md, model, cfg)
  ct  <- coef_table(fit)
  rb  <- robust_table(fit, cfg$level)
  if (!is.null(rb)) ct <- merge(ct, rb$tab, by = "term", sort = FALSE)
  ct  <- ct[match(strip_term(rownames(fit$beta)), ct$term), ]
  out <- list(fit = fit, coef = ct, het = heterogeneity(fit),
              robust_vcov = if (!is.null(rb)) rb$vcov else NULL)
  if (do_mc) {
    defs <- indirect_defs(cfg$pool_roles)
    b    <- fit$beta[, 1]
    out$indirect <- mc_functions(b, fit$vb, defs, cfg$n_mc, cfg$level)
  }
  out
}


## =============================================================================
## 5. RUN
## =============================================================================
message("== Reading and preparing studies")
raw_list <- setNames(lapply(CFG$studies$file, function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
  utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}), CFG$studies$study)

prepped <- setNames(lapply(names(raw_list), function(s)
  prep_study(raw_list[[s]], s, CFG)), names(raw_list))

message("== Fitting study-level APIMeM (uncorrected + corrected; ML + MLR)")
study_res <- lapply(prepped, analyse_study)

if (isTRUE(CFG$boot_V)) {
  message("== Bootstrapping sampling covariance matrices (B = ", CFG$n_boot, ")")
  for (s in names(study_res)) {
    bv <- bootstrap_V(raw_list[[s]], s, CFG, CFG$n_boot)
    study_res[[s]]$uncorrected$V$BOOT <- bv$uncorrected
    study_res[[s]]$corrected$V$BOOT   <- bv$corrected
  }
}

k_studies <- length(study_res)
if (k_studies < 5) {
  message("\n!! Only ", k_studies, " studies. Between-study variance components ",
          "are essentially\n   unidentified and CR2 Satterthwaite df will be < 4 ",
          "(robust tests unreliable;\n   Tipton, 2015). Treat the common-effect ",
          "model as at least as informative as\n   the multilevel model. CR2 results ",
          "are saved (CSV/RDS) but not tabled\n   while any df < CFG$cr2_min_df.\n")
}

message("== Meta-analysis")
V_TYPES <- intersect(c("ML", "MLR", "BOOT"), names(study_res[[1]]$uncorrected$V))
MODELS  <- c("multilevel", "common")
if (!CFG$primary_V %in% V_TYPES) stop("primary_V '", CFG$primary_V, "' not available")

MA <- list()
for (cmp in c("uncorrected", "corrected")) {
  for (vt in V_TYPES) {
    for (mdl in MODELS) {
      MA[[cmp]][[vt]][[mdl]] <- run_ma(study_res, cmp, vt, mdl, CFG)
    }
  }
}
PRIM <- lapply(MA, function(x) x[[CFG$primary_V]][[CFG$primary_model]])

## ---- Console summary ----
for (cmp in names(PRIM)) {
  cat("\n=====", toupper(cmp), "(", CFG$primary_model, ", V =", CFG$primary_V, ") =====\n")
  print(transform(PRIM[[cmp]]$coef[, c("term", "est", "se", "lo", "hi", "p")],
                  est = round(est, 3), se = round(se, 3), lo = round(lo, 3),
                  hi = round(hi, 3), p = signif(p, 3)), row.names = FALSE)
  cat("\nIndirect / total effects (MC CI)\n")
  ind <- PRIM[[cmp]]$indirect
  print(data.frame(effect = ind$expr, est = round(ind$est, 3),
                   lo = round(ind$lo, 3), hi = round(ind$hi, 3)), row.names = FALSE)
}


## =============================================================================
## 6. EXPORT RAW RESULTS (CSV + RDS)
## =============================================================================
long <- do.call(rbind, lapply(names(MA), function(cmp)
  do.call(rbind, lapply(names(MA[[cmp]]), function(vt)
    do.call(rbind, lapply(names(MA[[cmp]][[vt]]), function(mdl) {
      ct <- MA[[cmp]][[vt]][[mdl]]$coef
      ct$component <- cmp; ct$V_type <- vt; ct$model <- mdl
      ct
    }))))))
utils::write.csv(long, CFG$out_csv, row.names = FALSE)
saveRDS(list(config = CFG, prepped = lapply(prepped, `[`, c("study", "rel", "n")),
             study_results = study_res, meta = MA), CFG$out_rds)


## =============================================================================
## 7. APA TABLES (flextable + officer)
## =============================================================================

## ---- 7a. Number formatting (APA 7) ------------------------------------------
MINUS <- "\u2212"
fmt <- function(x, d = 2, lead0 = TRUE) {
  s <- formatC(round(x, d), format = "f", digits = d)
  s <- sub(paste0("^-(0\\.0+)$"), "\\1", s)             # no "-0.00"
  if (!lead0) s <- sub("^(-?)0\\.", "\\1.", s)          # values bounded by 1
  s <- sub("^-", MINUS, s)
  s[is.na(x)] <- "\u2014"
  s
}
fmt_p  <- function(p) ifelse(is.na(p), "\u2014",
                             ifelse(p < .001, "< .001", fmt(p, 3, lead0 = FALSE)))
fmt_ci <- function(lo, hi, d = 2) ifelse(is.na(lo) | is.na(hi), "\u2014",
                                         sprintf("[%s, %s]", fmt(lo, d), fmt(hi, d)))
fmt_df <- function(x) ifelse(is.na(x), "\u2014", formatC(x, format = "f", digits = 2))
pct <- paste0(round(100 * CFG$level), "%")

## ---- 7b. APA flextable builder ---------------------------------------------
##  body        data.frame (all character); column names are internal keys
##  labels      display labels for the bottom header row (same length)
##  spans       list of spanning header rows, top to bottom; each is a named
##              integer vector of column widths (names = labels, "" = blank)
##  group_rows  body row indices that are group headings (merged + italic)
##  note        character note (printed after an italic "Note.")
apa_table <- function(body, labels, spans = list(), group_rows = integer(0),
                      note = NULL, first_col_width = 2.4, col_width = 0.85,
                      font = CFG$table_font, size = CFG$table_size) {
  bd <- officer::fp_border(color = "black", width = 1)
  nc <- ncol(body)
  ft <- flextable::flextable(body)
  ft <- flextable::set_header_labels(ft, values = setNames(as.list(labels), names(body)))
  for (sp in rev(spans)) {
    ft <- flextable::add_header_row(ft, values = names(sp), colwidths = unname(sp),
                                    top = TRUE)
  }
  ft <- flextable::border_remove(ft)
  ft <- flextable::font(ft, fontname = font, part = "all")
  ft <- flextable::fontsize(ft, size = size, part = "all")
  ft <- flextable::align(ft, align = "center", part = "header")
  ft <- flextable::align(ft, j = 1, align = "left", part = "all")
  if (nc > 1) ft <- flextable::align(ft, j = 2:nc, align = "center", part = "body")
  ft <- flextable::valign(ft, valign = "bottom", part = "header")
  ft <- flextable::padding(ft, padding.top = 1.5, padding.bottom = 1.5, part = "all")
  ft <- flextable::line_spacing(ft, space = 1, part = "all")

  ## APA rules: top of header, bottom of header, bottom of body; partial
  ## rules under each non-empty spanning label
  ft <- flextable::hline_top(ft, border = bd, part = "header")
  ft <- flextable::hline_bottom(ft, border = bd, part = "header")
  ft <- flextable::hline_bottom(ft, border = bd, part = "body")
  for (r in seq_along(spans)) {
    sp <- spans[[r]]; start <- 1
    for (s in seq_along(sp)) {
      cols <- start:(start + sp[s] - 1)
      if (nzchar(names(sp)[s]))
        ft <- flextable::hline(ft, i = r, j = cols,
                               border = officer::fp_border(color = "black", width = .75),
                               part = "header")
      start <- start + sp[s]
    }
  }
  ## Group heading rows spanning the full width; indent the rows beneath
  if (length(group_rows)) {
    ind <- setdiff(seq_len(nrow(body)), group_rows)
    ft  <- flextable::padding(ft, i = ind, j = 1, padding.left = 14, part = "body")
  }
  for (g in group_rows) {
    ft <- flextable::merge_at(ft, i = g, j = 1:nc, part = "body")
    ft <- flextable::italic(ft, i = g, j = 1, part = "body")
    ft <- flextable::align(ft, i = g, j = 1, align = "left", part = "body")
  }
  if (!is.null(note)) {
    ft <- flextable::add_footer_lines(
      ft, values = flextable::as_paragraph(flextable::as_i("Note. "), note))
    ft <- flextable::font(ft, fontname = font, part = "footer")
    ft <- flextable::fontsize(ft, size = size, part = "footer")
    ft <- flextable::align(ft, align = "left", part = "footer")
  }
  ft <- flextable::width(ft, j = 1, width = first_col_width)
  if (nc > 1) ft <- flextable::width(ft, j = 2:nc, width = col_width)
  flextable::set_table_properties(ft, layout = "fixed")
}

## Raw WordprocessingML paragraph (used for table number/title so that the
## title can carry "page break before" + "keep with next": each table starts
## on a new page without ever producing blank pages, and the title never
## separates from its table)
title_par_xml <- function(text, bold = FALSE, italic = FALSE, page_break = FALSE,
                          font = CFG$table_font, size_pt = 12) {
  txt <- gsub("&", "&amp;", text, fixed = TRUE)
  txt <- gsub("<", "&lt;", txt, fixed = TRUE)
  txt <- gsub(">", "&gt;", txt, fixed = TRUE)
  rpr <- sprintf(paste0('<w:rPr><w:rFonts w:ascii="%1$s" w:hAnsi="%1$s" w:cs="%1$s"/>',
                        '%2$s%3$s<w:sz w:val="%4$d"/><w:szCs w:val="%4$d"/></w:rPr>'),
                 font, if (bold) "<w:b/>" else "", if (italic) "<w:i/>" else "",
                 as.integer(size_pt * 2))
  ppr <- sprintf('<w:pPr><w:keepNext/>%s<w:spacing w:before="0" w:after="120"/></w:pPr>',
                 if (page_break) "<w:pageBreakBefore/>" else "")
  sprintf(paste0('<w:p xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">',
                 '%s<w:r>%s<w:t xml:space="preserve">%s</w:t></w:r></w:p>'), ppr, rpr, txt)
}

## Add "Table n" (bold) + italic title + table to the document
add_apa_table <- function(doc, number, title, ft, new_page = TRUE) {
  doc <- officer::body_add_xml(doc, title_par_xml(paste("Table", number), bold = TRUE,
                                                  page_break = new_page))
  doc <- officer::body_add_xml(doc, title_par_xml(title, italic = TRUE))
  ## Never exceed the text width of a landscape page (scale columns if needed)
  w <- dim(ft)$widths
  if (sum(w) > CFG$max_table_width)
    ft <- flextable::width(ft, j = seq_along(w), width = w * CFG$max_table_width / sum(w))
  flextable::body_add_flextable(doc, ft, align = "left")
}

## Row labels for effects
effect_row_label <- function(term) {
  if (CFG$pool_roles) {
    ifelse(grepl("^actor", term), "Actor", "Partner")
  } else {
    p <- PATHS[match(term, PATHS$label), ]
    sprintf("%s (%s \u2192 %s)", ifelse(p$type == "actor", "Actor", "Partner"),
            p$rhs, p$lhs)
  }
}
term_path <- function(term) {
  if (CFG$pool_roles) sub("^(actor|partner)_", "", term)
  else PATHS$path[match(term, PATHS$label)]
}

## Insert group heading rows (by path) into a body data.frame
with_groups <- function(body, group_key, group_labels) {
  out <- NULL; grp <- integer(0)
  for (g in unique(group_key)) {
    head_row <- body[1, , drop = FALSE]; head_row[] <- ""
    head_row[[1]] <- group_labels[[g]]
    out <- rbind(out, head_row, body[group_key == g, , drop = FALSE])
    grp <- c(grp, nrow(out) - sum(group_key == g))
  }
  rownames(out) <- NULL
  list(body = out, group_rows = grp)
}

study_names <- names(study_res)
n_by_study  <- vapply(prepped, `[[`, numeric(1), "n")
rel_name    <- if (CFG$reliability == "omega") "\u03c9" else "\u03b1"
var_label   <- function(v) sprintf("%s, %s", CFG$construct_labels[[substr(v, 1, 1)]],
                                   CFG$member_labels[as.integer(substr(v, 2, 2))])
df_range <- function(obj) {
  d <- obj$coef$df_cr2
  if (is.null(d) || all(is.na(d))) return("not available")
  sprintf("%s\u2013%s", fmt_df(min(d, na.rm = TRUE)), fmt_df(max(d, na.rm = TRUE)))
}
model_txt <- if (CFG$primary_model == "multilevel") {
  if (CFG$pool_roles)
    "a multilevel random-effects model (REML) in which the Partner 1 and Partner 2 estimates of each effect are nested within studies"
  else "a multivariate random-effects model (REML)"
} else "a common-effect (fixed-effect) multivariate model"
V_txt <- c(ML = "normal-theory (ML)", MLR = "robust (MLR sandwich)",
           BOOT = "bootstrap")[[CFG$primary_V]]
num_word <- function(n) if (n >= 1 && n <= 9) c("one", "two", "three", "four",
  "five", "six", "seven", "eight", "nine")[n] else as.character(n)
k_txt <- sprintf("Pooled across %s studies (%s dyads)", num_word(k_studies),
                 if (k_studies == 2) paste(n_by_study, collapse = " and ")
                 else paste(n_by_study, collapse = ", "))

doc <- officer::read_docx()
doc <- officer::body_set_default_section(
  doc, officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1)))

## ---- Table 1: descriptives & reliabilities ---------------------------------
## Wide (studies side by side) for up to 3 studies, long (stacked) otherwise
rel_source_txt <- paste0(
  if (CFG$reliability == "omega")
    "McDonald's omega total from a one-factor congeneric CFA (FIML)"
  else "Cronbach's alpha",
  if (!is.null(CFG$reliability_supplied)) ", or as supplied by the study" else "")
note1 <- sprintf(paste0("Scale scores are item means. %s = %s. Reliabilities were ",
                        "used to fix the error variance of each scale score in the ",
                        "corrected (latent) APIMeM."), rel_name, rel_source_txt)
desc_block <- function(s) {
  cm <- prepped[[s]]$comp
  data.frame(M   = fmt(colMeans(cm, na.rm = TRUE)),
             SD  = fmt(apply(cm, 2, stats::sd, na.rm = TRUE)),
             rel = fmt(prepped[[s]]$rel, 2, lead0 = FALSE), stringsAsFactors = FALSE)
}
if (k_studies <= 3) {
  t1 <- data.frame(var = vapply(VARS, var_label, character(1)), stringsAsFactors = FALSE)
  for (s in study_names) t1 <- cbind(t1, setNames(desc_block(s), paste0(make.names(s), c("_M", "_SD", "_rel"))))
  ft1 <- apa_table(
    t1, labels = c("Variable", rep(c("M", "SD", rel_name), k_studies)),
    spans = list(c(" " = 1L, setNames(rep(3L, k_studies),
                                      sprintf("%s (%d dyads)", study_names, n_by_study)))),
    note = note1, first_col_width = 2.6, col_width = 0.7)
  hdr_row <- 2
} else {
  t1 <- NULL; g1 <- integer(0)
  for (s in study_names) {
    blk <- cbind(data.frame(var = vapply(VARS, var_label, character(1)),
                            stringsAsFactors = FALSE), desc_block(s))
    hr <- blk[1, ]; hr[] <- ""; hr$var <- sprintf("%s (%d dyads)", s, prepped[[s]]$n)
    t1 <- rbind(t1, hr, blk); g1 <- c(g1, nrow(t1) - 6)
  }
  rownames(t1) <- NULL
  ft1 <- apa_table(t1, labels = c("Variable", "M", "SD", rel_name), group_rows = g1,
                   note = note1, first_col_width = 2.8, col_width = 0.8)
  hdr_row <- 1
}
ft1 <- flextable::italic(ft1, i = hdr_row, j = which(grepl("(^|_)(M|SD)$", names(t1))),
                         part = "header")
doc <- add_apa_table(doc, 1, "Descriptive Statistics and Reliability Coefficients by Study", ft1,
                     new_page = FALSE)

## ---- Table 2: pooled direct effects -----------------------------------------
cr2_df_all <- unlist(lapply(PRIM, function(o) o$coef$df_cr2))
show_cr2   <- length(cr2_df_all) > 0 && !anyNA(cr2_df_all) &&
              all(cr2_df_all >= CFG$cr2_min_df)
direct_cols <- c("b", "se", "ci", "p", if (show_cr2) "cr2")
make_direct_body <- function(obj) {
  ct  <- obj$coef
  out <- data.frame(b = fmt(ct$est), se = fmt(ct$se), ci = fmt_ci(ct$lo, ct$hi),
                    p = fmt_p(ct$p), stringsAsFactors = FALSE)
  if (show_cr2) out$cr2 <- fmt_ci(ct$lo_cr2, ct$hi_cr2)
  out
}
terms <- PRIM$uncorrected$coef$term
t2 <- cbind(data.frame(eff = effect_row_label(terms), stringsAsFactors = FALSE),
            setNames(make_direct_body(PRIM$uncorrected), paste0("u_", direct_cols)),
            setNames(make_direct_body(PRIM$corrected),   paste0("c_", direct_cols)))
g2 <- with_groups(t2, term_path(terms), PATH_GROUP_LABEL)
nd <- length(direct_cols)
cr2_note <- if (show_cr2) {
  sprintf(paste0("CR2 = cluster-robust (bias-reduced linearization) confidence ",
                 "interval with Satterthwaite df (uncorrected: %s; corrected: %s)."),
          df_range(PRIM$uncorrected), df_range(PRIM$corrected))
} else {
  sprintf(paste0("Cluster-robust (CR2) inference was computed but is not reported ",
                 "because its Satterthwaite df (uncorrected: %s; corrected: %s) fall ",
                 "below %d, where robust tests are unreliable (Tipton, 2015)."),
          df_range(PRIM$uncorrected), df_range(PRIM$corrected), CFG$cr2_min_df)
}
ft2 <- apa_table(
  g2$body,
  labels = c("Effect", rep(c("\u03b2", "SE", paste(pct, "CI"), "p",
                             if (show_cr2) paste(pct, "CI (CR2)")), 2)),
  spans  = list(c(" " = 1L, "Uncorrected" = nd, "Corrected for unreliability" = nd)),
  group_rows = g2$group_rows,
  note = sprintf(paste0(
    "%s. \u03b2 = pooled standardized APIMeM path estimated with %s, using the full ",
    "within-study %s sampling covariance matrix of the 12 paths. Uncorrected = ",
    "APIMeM on observed scale scores; corrected = APIMeM on single-indicator latent ",
    "variables with error variances fixed at (1 \u2212 %s)s\u00b2. CI = model-based ",
    "Wald confidence interval. %s Results under robust (MLR) sampling covariances ",
    "and the common-effect model are in Table 5."),
    k_txt, model_txt, V_txt, rel_name, cr2_note),
  first_col_width = if (show_cr2) 1.3 else 1.8, col_width = if (show_cr2) 0.52 else 0.8)
ci_cols <- c(4, if (show_cr2) 6, 4 + nd, if (show_cr2) 6 + nd)
ft2 <- flextable::width(ft2, j = ci_cols, width = if (show_cr2) 1.0 else 1.1)
ft2 <- flextable::italic(ft2, i = 2, j = c(3, 5, 3 + nd, 5 + nd), part = "header")
doc <- add_apa_table(doc, 2, "Pooled Direct Effects of the Actor\u2013Partner Interdependence Mediation Model", ft2)

## ---- Table 3: pooled indirect and total effects ----------------------------
iu <- PRIM$uncorrected$indirect; ic <- PRIM$corrected$indirect
t3 <- data.frame(eff  = iu$label,
                 u_e  = fmt(iu$est, 3), u_ci = fmt_ci(iu$lo, iu$hi, 3),
                 c_e  = fmt(ic$est, 3), c_ci = fmt_ci(ic$lo, ic$hi, 3),
                 stringsAsFactors = FALSE)
grp_lab <- setNames(as.list(unique(iu$group)), unique(iu$group))
g3 <- with_groups(t3, iu$group, grp_lab)
ft3 <- apa_table(
  g3$body,
  labels = c("Effect", rep(c("Estimate", paste(pct, "MC CI")), 2)),
  spans  = list(c(" " = 1L, "Uncorrected" = 2L, "Corrected for unreliability" = 2L)),
  group_rows = g3$group_rows,
  note = sprintf(paste0(
    "%s. Estimates are products (and sums) of the pooled standardized paths in ",
    "Table 2 (Ledermann et al., 2011). MC CI = Monte Carlo confidence interval ",
    "based on %s draws from the joint sampling distribution of the pooled paths ",
    "(including their covariances). \"Own\" and \"partner's\" are relative to ",
    "the partner whose X is the predictor; e.g., partner\u2013actor = partner a ",
    "path \u00d7 actor b path. Total effects add the pooled direct (c\u2032) path. ",
    "Sensitivity analyses are in Table 6."),
    k_txt, format(CFG$n_mc, big.mark = ",", scientific = FALSE)),
  first_col_width = 3.7, col_width = 1.3)
ft3 <- flextable::width(ft3, j = c(2, 4), width = 0.9)
doc <- add_apa_table(doc, 3, "Pooled Indirect and Total Effects", ft3)

## ---- Table 4: heterogeneity -------------------------------------------------
het_row <- function(cmp, mdl) {
  h <- MA[[cmp]][[CFG$primary_V]][[mdl]]$het
  data.frame(
    cmp  = paste0(c(uncorrected = "Uncorrected", corrected = "Corrected")[[cmp]],
                  ", ", c(multilevel = "multilevel", common = "common effect")[[mdl]]),
    k    = as.character(h$k),
    Q    = sprintf("%s (%d)", fmt(h$QE, 2), h$QE_df),
    p    = fmt_p(h$QE_p),
    t2   = fmt(h$tau2_study, 4), s2 = fmt(h$sigma2_role, 4),
    I2s  = ifelse(is.na(h$I2_study), "\u2014", sprintf("%.1f", h$I2_study)),
    I2r  = ifelse(is.na(h$I2_role),  "\u2014", sprintf("%.1f", h$I2_role)),
    stringsAsFactors = FALSE)
}
t4 <- do.call(rbind, lapply(c("uncorrected", "corrected"), function(cmp)
  rbind(het_row(cmp, "multilevel"), het_row(cmp, "common"))))
ft4 <- apa_table(
  t4,
  labels = c("Model", "k", "Q\u2091 (df)", "p", "\u03c4\u00b2 (study)",
             "\u03c3\u00b2 (partner)", "I\u00b2 study (%)", "I\u00b2 partner (%)"),
  note = sprintf(paste0(
    "k = number of path estimates (%d studies \u00d7 12 paths). Q\u2091 = test of ",
    "residual heterogeneity given the full sampling covariance matrix. \u03c4\u00b2 = ",
    "between-study variance (per effect type, structure \"%s\"); \u03c3\u00b2 = variance ",
    "of Partner 1 vs. Partner 2 estimates within studies. I\u00b2 decomposed by level ",
    "(Cheung, 2014). With %s studies, \u03c4\u00b2 is very imprecisely estimated."),
    k_studies, CFG$re_struct, num_word(k_studies)),
  first_col_width = 2.6, col_width = 0.95)
ft4 <- flextable::italic(ft4, j = c(2, 4), part = "header")
doc <- add_apa_table(doc, 4, "Heterogeneity of the Pooled APIMeM Paths", ft4)

## ---- Tables 5 and 6: sensitivity (V type x pooling model) -----------------
## Columns: pooling model (spanning) x sampling covariance type; primary first.
## Rows: grouped by component.
sens_cols <- expand.grid(vt = V_TYPES, model = MODELS, stringsAsFactors = FALSE)
sens_cols <- sens_cols[order(sens_cols$model != CFG$primary_model,
                             sens_cols$vt != CFG$primary_V), ]
nsc <- nrow(sens_cols)
model_lab <- c(multilevel = "Multilevel model", common = "Common-effect model")
sens_spans <- c(" " = 1L, setNames(as.integer(table(factor(sens_cols$model,
                                     levels = unique(sens_cols$model)))),
                                   model_lab[unique(sens_cols$model)]))
V_note <- paste0(
  "ML = normal-theory and MLR = robust (Huber\u2013White sandwich) within-study ",
  "sampling covariance matrices",
  if ("BOOT" %in% V_TYPES) sprintf(paste0("; BOOT = nonparametric dyad bootstrap ",
    "(B = %d), which for the corrected component also propagates sampling error ",
    "in the reliability coefficients"), CFG$n_boot) else "",
  ". Multilevel = random effects for studies and for partner-specific estimates ",
  "within studies; common-effect = no random effects. The first column is the ",
  "primary specification (Tables 2 and 3).")

effect_full_label <- function(term) {
  if (CFG$pool_roles) {
    sym <- c(a = "a", b = "b", c = "c\u2032")[term_path(term)]
    sprintf("%s %s", ifelse(grepl("^actor", term), "Actor", "Partner"), sym)
  } else effect_row_label(term)
}
comp_lab <- c(uncorrected = "Uncorrected", corrected = "Corrected for unreliability")

sens_table <- function(get_cells, row_labels) {
  body <- NULL; grp <- integer(0)
  for (cmp in c("uncorrected", "corrected")) {
    blk <- data.frame(eff = row_labels, stringsAsFactors = FALSE)
    for (i in seq_len(nsc))
      blk[[paste0("s", i)]] <- get_cells(cmp, sens_cols$vt[i], sens_cols$model[i])
    head_row <- blk[1, , drop = FALSE]; head_row[] <- ""; head_row$eff <- comp_lab[[cmp]]
    body <- rbind(body, head_row, blk)
    grp  <- c(grp, nrow(body) - nrow(blk))
  }
  rownames(body) <- NULL
  list(body = body, group_rows = grp)
}

g5 <- sens_table(function(cmp, vt, mdl) {
  ct <- MA[[cmp]][[vt]][[mdl]]$coef
  sprintf("%s %s", fmt(ct$est), fmt_ci(ct$lo, ct$hi))
}, effect_full_label(terms))
ft5 <- apa_table(
  g5$body, labels = c("Effect", sens_cols$vt), spans = list(sens_spans),
  group_rows = g5$group_rows,
  note = sprintf("Cells: \u03b2 [%s CI]. %s", pct, V_note),
  first_col_width = if (nsc > 4) 1.5 else 2.2, col_width = min(1.65, 7.3 / nsc),
  size = if (nsc > 4) CFG$table_size - 1 else CFG$table_size)
doc <- add_apa_table(doc, 5, "Sensitivity of the Pooled Direct Effects to the Sampling Covariance Matrix and Pooling Model", ft5)

g6 <- sens_table(function(cmp, vt, mdl) {
  ind <- MA[[cmp]][[vt]][[mdl]]$indirect
  sprintf("%s %s", fmt(ind$est, 3), fmt_ci(ind$lo, ind$hi, 3))
}, iu$short)
ft6 <- apa_table(
  g6$body, labels = c("Effect", sens_cols$vt), spans = list(sens_spans),
  group_rows = g6$group_rows,
  note = sprintf("Cells: estimate [%s Monte Carlo CI]. %s", pct, V_note),
  first_col_width = if (nsc > 4) 1.7 else 2.2, col_width = min(1.65, 7.2 / nsc),
  size = if (nsc > 4) CFG$table_size - 1 else CFG$table_size)
doc <- add_apa_table(doc, 6, "Sensitivity of the Pooled Indirect and Total Effects", ft6)

## ---- Table S1: study-level estimates ----------------------------------------
## Wide for up to 2 studies, long (stacked by study) otherwise
path_lab <- sprintf("%s: %s \u2192 %s", ifelse(PATHS$type == "actor", "Actor", "Partner"),
                    PATHS$rhs, PATHS$lhs)
est_block <- function(s) {
  out <- list()
  for (cmp in c("uncorrected", "corrected")) {
    r <- study_res[[s]][[cmp]]
    out[[paste0(substr(cmp, 1, 1), "_b")]]  <- fmt(r$est)
    out[[paste0(substr(cmp, 1, 1), "_se")]] <- fmt(sqrt(diag(r$V$ML)))
    out[[paste0(substr(cmp, 1, 1), "_sr")]] <- fmt(sqrt(diag(r$V$MLR)))
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}
noteS1 <- paste0("Standardized estimates from each study's APIMeM. SE = normal-theory ",
                 "(delta-method) standard error; SE\u1d63 = robust (MLR) standard error. ",
                 "The full covariance matrices of these estimates (not only the SEs) ",
                 "enter the meta-analysis.")
if (k_studies <= 2) {
  tS1 <- data.frame(eff = path_lab, stringsAsFactors = FALSE)
  for (s in study_names) {
    eb  <- est_block(s)
    tS1 <- cbind(tS1, setNames(eb, paste0(make.names(s), "_", names(eb))))
  }
  gS1 <- with_groups(tS1, PATHS$path, PATH_GROUP_LABEL)
  ftS1 <- apa_table(
    gS1$body,
    labels = c("Path", rep(c("\u03b2", "SE", "SE\u1d63"), 2 * k_studies)),
    spans = list(
      c(" " = 1L, setNames(rep(6L, k_studies), sprintf("%s (%d dyads)", study_names, n_by_study))),
      c(" " = 1L, setNames(rep(3L, 2 * k_studies), rep(c("Uncorrected", "Corrected"), k_studies)))),
    group_rows = gS1$group_rows, note = noteS1,
    first_col_width = 1.9, col_width = 0.58)
  se_cols <- sort(c(seq(3, 1 + 6 * k_studies, by = 3), seq(4, 1 + 6 * k_studies, by = 3)))
  ftS1 <- flextable::italic(ftS1, i = 3, j = se_cols, part = "header")
} else {
  tS1 <- NULL; gS1 <- integer(0)
  for (s in study_names) {
    blk <- cbind(data.frame(eff = path_lab, stringsAsFactors = FALSE), est_block(s))
    hr <- blk[1, ]; hr[] <- ""; hr$eff <- sprintf("%s (%d dyads)", s, prepped[[s]]$n)
    tS1 <- rbind(tS1, hr, blk); gS1 <- c(gS1, nrow(tS1) - 12)
  }
  rownames(tS1) <- NULL
  ftS1 <- apa_table(
    tS1, labels = c("Path", rep(c("\u03b2", "SE", "SE\u1d63"), 2)),
    spans = list(c(" " = 1L, "Uncorrected" = 3L, "Corrected" = 3L)),
    group_rows = gS1, note = noteS1, first_col_width = 2.2, col_width = 0.8)
  ftS1 <- flextable::italic(ftS1, i = 2, j = c(3, 4, 6, 7), part = "header")
}
doc <- add_apa_table(doc, "S1", "Study-Level Standardized APIMeM Paths", ftS1)

## ---- Table S2: correlation matrices -----------------------------------------
tS2 <- NULL; grpS2 <- integer(0)
for (s in study_names) {
  R  <- stats::cor(prepped[[s]]$comp, use = "pairwise.complete.obs")
  rl <- prepped[[s]]$rel
  blk <- data.frame(var = sprintf("%d. %s", 1:6, vapply(VARS, var_label, character(1))),
                    stringsAsFactors = FALSE)
  for (j in 1:6) {
    col <- ifelse(seq_len(6) > j, fmt(R[, j], 2, lead0 = FALSE), "")
    col[j] <- sprintf("(%s)", fmt(rl[j], 2, lead0 = FALSE))
    blk[[paste0("v", j)]] <- col
  }
  head_row <- blk[1, ]; head_row[] <- ""
  head_row$var <- sprintf("%s (%d dyads)", s, prepped[[s]]$n)
  tS2 <- rbind(tS2, head_row, blk)
  grpS2 <- c(grpS2, nrow(tS2) - 6)
}
rownames(tS2) <- NULL
ftS2 <- apa_table(
  tS2, labels = c("Variable", as.character(1:6)), group_rows = grpS2,
  note = sprintf(paste0("Pairwise correlations between scale scores; reliability ",
                        "coefficients (%s) on the diagonal in parentheses. These ",
                        "matrices allow re-analysis with two-stage MASEM."), rel_name),
  first_col_width = 2.8, col_width = 0.7)
doc <- add_apa_table(doc, "S2", "Correlations Among Scale Scores by Study", ftS2)

print(doc, target = CFG$out_docx)
message("\n== Done. Tables written to ", CFG$out_docx,
        "\n   Numeric results: ", CFG$out_csv, " and ", CFG$out_rds)
