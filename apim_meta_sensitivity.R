
## title: "APIM meta-analysis: sensitivity"
## author: "Anonymized"
## date: "15-10-2026"

###############################################################################
##  apim_meta_sensitivity.R
##
##  Sensitivity (power / detectability) analysis for the pooled APIMeM paths of
##  apim_meta_analysis.R: how large do the pooled direct (actor/partner a, b,
##  c') and indirect effects have to be to become statistically significant,
##  for the uncorrected and the reliability-corrected component, under
##  different amounts of between-study heterogeneity, with rma() and rma.mv(),
##  and with standard and robust inference?
##
##  No dyad-level data are generated. Everything is simulated at the study
##  level from the ASSUMED study-level inputs: each study's 12 standardized
##  path estimates and their 12 x 12 sampling covariance matrix (by default
##  taken from the results file of apim_meta_analysis.R; alternatively typed
##  in as means and standard errors, see Section 2).
##  
##  Conservative assumptions for all the paths (see Section 2).
##
##  Data-generating model (per replication, per study i, per path p of
##  effect type t):
##
##      y_ip = mu_t + u_it + w_ip + e_ip
##      u_it ~ N(0, s * tau2_t)          study-level deviation (shared by the
##                                        Partner 1 and Partner 2 estimate)
##      w_ip ~ N(0, (1 - s) * tau2_t)    partner-specific deviation
##      e_i  ~ MVN(0, V_i)               sampling error, full within-study V
##  tau2_t is set from a target I^2 (or directly as tau), separately for the
##  uncorrected and corrected components because their V differ.
##
##  Analysis models (all fitted to every simulated dataset)
##    rma()    on the study-level composite of the two partner-specific
##             estimates (GLS average using the 2 x 2 block of V_i):
##               CE (z)          common-effect, Wald z
##               RE (z)          random-effects (REML), Wald z
##               RE (KH)         random-effects, Knapp-Hartung t (df = k - 1)
##               RE (KH ad hoc)  KH with the variance factor truncated at 1
##
##    rma.mv() the multilevel model of apim_meta_analysis.R (partner-specific
##             estimates nested in studies, full V):
##               ML (z)          model-based Wald z
##               ML (CR2)        cluster-robust CR2, Satterthwaite t
##
##  Computational shortcut (exact, not an approximation)
##    All models are intercept(-per-effect-type) models with V fixed, so every
##    variance estimate, standard error, KH factor, CR2 matrix and
##    Satterthwaite df is invariant to shifting the true effects mu. Each
##    replication is therefore simulated ONCE under mu = 0; the test for any
##    true effect mu is |mu + est0| > crit * SE. This gives complete power
##    curves (and power surfaces for indirect effects) from one set of fits.
##    Section 6 verifies the invariance numerically.
##
##  Indirect effects
##    Significance = joint significance test (both constituent paths
##    significant; Fritz & MacKinnon, 2007; Yzerbyt et al., 2018), which can
##    be evaluated on a full (a, b) grid with the shortcut above. An optional
##    check compares it with the Monte Carlo CI used in apim_meta_analysis.R.
##
##  Outputs: power curves, minimum detectable effects (MDE, e.g. 80% power),
##  distribution of the critical pooled estimate, empirical Type I error,
##  critical estimates for the observed data, figures (PNG + PDF), CSV/RDS,
##  and an APA-style .docx summary.
###############################################################################

## =============================================================================
## 0. PACKAGES
## =============================================================================

pkgs <- c("metafor", "clubSandwich", "ggplot2", "flextable", "officer", "parallel")
miss_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss_pkgs)) stop("Please install: ", paste(miss_pkgs, collapse = ", "))
suppressPackageStartupMessages({
  library(metafor)
  library(clubSandwich)
  library(ggplot2)
  library(flextable)
  library(officer)
})

## =============================================================================
## 1. CONFIGURATION
## =============================================================================

CFG <- list(
  ## ---- Inputs
  input_mode = "rds",                              # "rds" or "table" (Section 2)
  rds_file   = "APIM_meta_analysis_results.rds",   # from apim_meta_analysis.R
  V_type     = "ML",                               # "ML" or "MLR" sampling covariances

  ## ---- Heterogeneity
  het_metric = "I2",                  # "I2": values are I^2 of the study-level
                                      #       (partner-averaged) effects
                                      # "tau": values are tau (SD of true effects)
  het_values = c(0, .25, .50, .75),
  het_share_study = 1,                # s: share of tau^2 at the study level
                                      #    (1 = both partners' estimates shift
                                      #    together; < 1 adds partner-level
                                      #    heterogeneity within studies)

  ## ---- Number of studies: 1 = the observed studies only. c(1, 2, 5) also
  ##      evaluates k = 4 and k = 10 hypothetical studies with the same V.
  k_replicates = 1,

  ## ---- Simulation
  R        = 1000,                    # replications per condition
  seed     = 20260912,
  n_cores  = max(1L, parallel::detectCores() - 1L),
  tau2_method = "REML",               # tau^2 estimator in rma()
  re_struct   = "ID",                 # rma.mv random-effects structure (as in
                                      # apim_meta_analysis.R)

  ## ---- Effect grids (true pooled standardized effects, own scale)
  theta_max  = 1.00, theta_step = 0.0025,      # direct paths
  ind_max    = 0.70, ind_step   = 0.01,        # (a, b) grid for indirect paths

  ## ---- Criteria
  alpha        = .05,
  power_target = .80,

  ## ---- Optional check: joint significance vs Monte Carlo CI (indirect)
  validate_mc = TRUE, R_mc = 300, n_mc = 4000,

  ## ---- Output
  out_prefix = "APIM_sensitivity",
  fig_dpi    = 300,
  table_font = "Times New Roman", table_size = 10, max_table_width = 8.9
)

## Optional overrides: define CFG_OVERRIDE <- list(...) before source()

if (exists("CFG_OVERRIDE", inherits = TRUE)) {
  for (nm in names(CFG_OVERRIDE)) CFG[[nm]] <- CFG_OVERRIDE[[nm]]
}

stopifnot(CFG$input_mode %in% c("rds", "table"), CFG$het_metric %in% c("I2", "tau"),
          CFG$het_share_study >= 0, CFG$het_share_study <= 1,
          all(CFG$het_values >= 0), CFG$het_metric != "I2" || all(CFG$het_values < 1))

CFG_OVERRIDE <- list(k_replicates = c(1, 2, 5), validate_mc = FALSE,
                     out_prefix = "APIM_sensitivity_k")

source("apim_meta_sensitivity.R")

## =============================================================================
## 2. STUDY-LEVEL INPUTS  (assumed means and sampling (co)variances)
## =============================================================================

VARS    <- c("X1", "X2", "M1", "M2", "Y1", "Y2")

PATHS   <- data.frame(
  label  = c("a1", "a2", "a12", "a21", "b1", "b2", "b12", "b21",
             "c1", "c2", "c12", "c21"),
  effect = rep(c("actor_a", "partner_a", "actor_b", "partner_b",
                 "actor_c", "partner_c"), each = 2),
  stringsAsFactors = FALSE)

EFFECTS <- unique(PATHS$effect)

EFFECT_LAB <- c(actor_a = "Actor a", partner_a = "Partner a",
                actor_b = "Actor b", partner_b = "Partner b",
                actor_c = "Actor c\u2032", partner_c = "Partner c\u2032")

COMPONENTS <- c("uncorrected", "corrected")

COMP_LAB   <- c(uncorrected = "Uncorrected", corrected = "Corrected for unreliability")

INDIRECT <- list(
  "Actor\u2013actor"     = c("actor_a",   "actor_b"),
  "Partner\u2013partner" = c("partner_a", "partner_b"),
  "Actor\u2013partner"   = c("actor_a",   "partner_b"),
  "Partner\u2013actor"   = c("partner_a", "actor_b"))

##  Mode "table": assumed study-level means and SEs (prefilled with the
##  toy-data results of apim_meta_analysis.R; replace with your own).
##  Within-study correlations among the 12 estimates are taken from
##  R_WITHIN (identity = independent estimates) unless you supply one.

ASSUMED <- list(
  "Study 1" = list(
    n_dyads = 130,
    uncorrected = list(
      est = c(0.386, 0.268, 0.164, 0.116, 0.354, 0.306, 0.107, 0.066, 0.089, 0.087, 0.107, -0.036),
      se  = c(0.075, 0.081, 0.083, 0.081, 0.086, 0.081, 0.088, 0.086, 0.090, 0.084, 0.088, 0.086)),
    corrected = list(
      est = c(0.460, 0.310, 0.180, 0.118, 0.430, 0.353, 0.111, 0.070, 0.075, 0.086, 0.116, -0.059),
      se  = c(0.090, 0.097, 0.100, 0.099, 0.111, 0.099, 0.114, 0.105, 0.119, 0.106, 0.117, 0.108))),
  "Study 2" = list(
    n_dyads = 174,
    uncorrected = list(
      est = c(0.183, 0.207, 0.171, 0.182, 0.217, 0.261, 0.099, 0.074, 0.082, 0.086, -0.060, 0.027),
      se  = c(0.072, 0.071, 0.072, 0.072, 0.077, 0.075, 0.077, 0.078, 0.074, 0.074, 0.073, 0.075)),
    corrected = list(
      est = c(0.290, 0.306, 0.260, 0.278, 0.368, 0.403, 0.105, 0.042, 0.102, 0.089, -0.138, 0.009),
      se  = c(0.117, 0.108, 0.113, 0.114, 0.154, 0.139, 0.150, 0.154, 0.135, 0.124, 0.127, 0.133)))
)

R_WITHIN <- diag(12)   # or a 12 x 12 correlation matrix in PATHS$label order

load_inputs <- function(cfg) {
  est <- V <- setNames(vector("list", 2), COMPONENTS)
  obs_ref <- NULL
  if (cfg$input_mode == "rds") {
    if (!file.exists(cfg$rds_file))
      stop("Results file '", cfg$rds_file, "' not found. Run apim_meta_analysis.R ",
           "first, or set CFG$input_mode = \"table\".")
    rr <- readRDS(cfg$rds_file)
    if (isFALSE(rr$config$pool_roles))
      message("Note: the results file was produced with pool_roles = FALSE; ",
              "this sensitivity analysis always pools partners (6 effect types).")
    studies <- names(rr$study_results)
    for (cmp in COMPONENTS) {
      est[[cmp]] <- lapply(rr$study_results, function(s) s[[cmp]]$est[PATHS$label])
      V[[cmp]]   <- lapply(rr$study_results, function(s) {
        if (is.null(s[[cmp]]$V[[cfg$V_type]])) stop("V type ", cfg$V_type, " not in results file")
        s[[cmp]]$V[[cfg$V_type]][PATHS$label, PATHS$label]
      })
    }
    n_dyads <- vapply(rr$prepped, function(p) as.numeric(p$n), numeric(1))
    ## Primary pooled estimates of apim_meta_analysis.R (for a consistency check)
    if (!is.null(rr$meta[[1]][[cfg$V_type]][["multilevel"]]))
      obs_ref <- lapply(COMPONENTS, function(cmp) {
        ct <- rr$meta[[cmp]][[cfg$V_type]][["multilevel"]]$coef
        setNames(ct$est, ct$term)
      })
  } else {
    studies <- names(ASSUMED)
    for (cmp in COMPONENTS) {
      est[[cmp]] <- lapply(ASSUMED, function(s) setNames(s[[cmp]]$est, PATHS$label))
      V[[cmp]]   <- lapply(ASSUMED, function(s) {
        D <- diag(s[[cmp]]$se); VV <- D %*% R_WITHIN %*% D
        dimnames(VV) <- list(PATHS$label, PATHS$label); VV
      })
    }
    n_dyads <- vapply(ASSUMED, `[[`, numeric(1), "n_dyads")
  }
  for (cmp in COMPONENTS) for (s in studies) {
    ev <- eigen(V[[cmp]][[s]], symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0) stop("V for ", s, " (", cmp, ") is not positive definite")
  }
  if (!is.null(obs_ref)) names(obs_ref) <- COMPONENTS
  list(studies = studies, n_dyads = n_dyads, est = est, V = V, obs_ref = obs_ref)
}

INPUT <- load_inputs(CFG)

K_BASE <- length(INPUT$studies)

## =============================================================================
## 3. DESIGNS AND HETEROGENEITY
## =============================================================================

METHODS <- c(rma_CE    = "rma(): CE (z)",
             rma_RE    = "rma(): RE (z)",
             rma_KH    = "rma(): RE (KH)",
             rma_KHadh = "rma(): RE (KH ad hoc)",
             mv_z      = "rma.mv(): multilevel (z)",
             mv_CR2    = "rma.mv(): multilevel (CR2)")

PAL <- setNames(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00"),
                METHODS)                                 # Okabe-Ito

TYPE_IDX <- match(PATHS$effect, EFFECTS)               # path -> effect type

EFF_COLS <- lapply(EFFECTS, function(e) which(PATHS$effect == e))

## Everything that does not change across replications (based on metafor)

make_design <- function(cmp, m) {
  V_list <- rep(INPUT$V[[cmp]], m)
  k <- length(V_list)
  comp_w <- lapply(EFF_COLS, function(cols) t(vapply(V_list, function(Vi) {
    w <- solve(Vi[cols, cols], c(1, 1)); w / sum(w)       # GLS composite weights
  }, numeric(2))))
  comp_v <- sapply(EFF_COLS, function(cols) vapply(V_list, function(Vi) {
    1 / sum(solve(Vi[cols, cols], c(1, 1)))                 # composite variance
  }, numeric(1)))
  comp_v <- matrix(comp_v, nrow = k)
  dat <- data.frame(study = factor(rep(seq_len(k), each = 12)),
                    label = rep(PATHS$label, k),
                    term  = factor(rep(PATHS$effect, k), levels = EFFECTS),
                    es_id = seq_len(12 * k))
  list(k = k, m = m, cmp = cmp, V = metafor::bldiag(V_list), dat = dat,
       comp_w = comp_w, comp_v = comp_v, chol = lapply(V_list, chol))
}

## tau^2 per component x heterogeneity level x effect type

v_typical <- sapply(COMPONENTS, function(cmp) colMeans(make_design(cmp, 1)$comp_v))

dimnames(v_typical) <- list(EFFECTS, COMPONENTS)

HET_LAB <- if (CFG$het_metric == "I2") {
  sprintf("I\u00b2 = %d%%", round(100 * CFG$het_values))
} else sprintf("\u03c4 = %.2f", CFG$het_values)

TAU2 <- array(NA_real_, c(2, length(CFG$het_values), 6),
              dimnames = list(COMPONENTS, HET_LAB, EFFECTS))

for (cmp in COMPONENTS) for (h in seq_along(CFG$het_values)) {
  x <- CFG$het_values[h]
  TAU2[cmp, h, ] <- if (CFG$het_metric == "I2") x / (1 - x) * v_typical[, cmp] else x^2
}

M_VALUES <- sort(unique(CFG$k_replicates))

K_LAB    <- sprintf("k = %d", K_BASE * M_VALUES)

NH <- length(HET_LAB); NM <- length(M_VALUES)

## =============================================================================
## 4. FITTING ONE (SIMULATED OR OBSERVED) DATASET
## =============================================================================
## tau^2 for the univariate model, REML with DL fallback

tau2_uni <- function(y, v, method) {
  f <- tryCatch(suppressWarnings(metafor::rma(y, v, method = method,
                  control = list(stepadj = 0.5, maxiter = 1000))),
                error = function(e) NULL)
  if (is.null(f)) f <- metafor::rma(y, v, method = "DL")
  f$tau2
}

## Multilevel model with optimizer fallback (based on metafor)

fit_mv <- function(dat, V, struct) {
  args <- list(yi = dat$yi, V = V, mods = ~ 0 + term,
               random = list(~ term | study, ~ 1 | es_id), struct = struct,
               data = dat, method = "REML")
  f <- tryCatch(suppressWarnings(do.call(metafor::rma.mv,
                  c(args, list(control = list(optimizer = "nlm"))))),
                error = function(e) NULL)
  if (is.null(f)) f <- tryCatch(suppressWarnings(do.call(metafor::rma.mv, args)),
                                error = function(e) NULL)
  f
}

## Returns est, se, crit (6 methods x 6 effects) + vcov of the multilevel model

fit_all <- function(Y, D, cfg) {
  est <- se <- crit <- matrix(NA_real_, length(METHODS), 6,
                              dimnames = list(names(METHODS), EFFECTS))
  z <- stats::qnorm(1 - cfg$alpha / 2)
  k <- nrow(Y)

  ## ---- rma(): study-level composites per effect type
  for (j in seq_along(EFFECTS)) {
    y  <- rowSums(Y[, EFF_COLS[[j]], drop = FALSE] * D$comp_w[[j]])
    v  <- D$comp_v[, j]
    w0 <- 1 / v
    est["rma_CE", j] <- sum(w0 * y) / sum(w0)
    se["rma_CE", j]  <- sqrt(1 / sum(w0))
    t2 <- tau2_uni(y, v, cfg$tau2_method)
    w  <- 1 / (v + t2)
    mu <- sum(w * y) / sum(w)
    q  <- sum(w * (y - mu)^2) / (k - 1)                  # KH variance factor
    est[c("rma_RE", "rma_KH", "rma_KHadh"), j] <- mu
    se["rma_RE", j]    <- sqrt(1 / sum(w))
    se["rma_KH", j]    <- sqrt(q / sum(w))
    se["rma_KHadh", j] <- sqrt(max(1, q) / sum(w))
  }
  crit[c("rma_CE", "rma_RE"), ] <- z
  crit[c("rma_KH", "rma_KHadh"), ] <- stats::qt(1 - cfg$alpha / 2, df = k - 1)

  ## ---- rma.mv(): multilevel model with full V
  dat <- D$dat; dat$yi <- as.vector(t(Y))
  f <- fit_mv(dat, D$V, cfg$re_struct)
  vb <- matrix(NA_real_, 6, 6, dimnames = list(EFFECTS, EFFECTS))
  df_cr2 <- setNames(rep(NA_real_, 6), EFFECTS)
  if (!is.null(f)) {
    b <- setNames(f$beta[, 1], sub("^term", "", rownames(f$beta)))[EFFECTS]
    est["mv_z", ] <- est["mv_CR2", ] <- b
    se["mv_z", ]   <- sqrt(diag(f$vb))[paste0("term", EFFECTS)]
    crit["mv_z", ] <- z
    vb[] <- f$vb[paste0("term", EFFECTS), paste0("term", EFFECTS)]
    rb <- tryCatch({
      vc <- clubSandwich::vcovCR(f, cluster = dat$study, type = "CR2")
      clubSandwich::coef_test(f, vcov = vc, test = "Satterthwaite")
    }, error = function(e) NULL)
    if (!is.null(rb)) {
      se["mv_CR2", ]   <- rb$SE
      crit["mv_CR2", ] <- stats::qt(1 - cfg$alpha / 2, df = rb$df_Satt)
      df_cr2[]         <- rb$df_Satt
    }
  }
  list(est = est, se = se, crit = crit, vb = vb, df_cr2 = df_cr2)
}

## Simulated study-level estimates for replication r under mu = 0
## (common random numbers across components, heterogeneity levels and k)

gen_Y <- function(r, D, tau2_t, s, Z) {
  t(vapply(seq_len(D$k), function(i) {
    e <- drop(crossprod(D$chol[[i]], Z$e[r, i, ]))
    u <- sqrt(s * tau2_t)[TYPE_IDX] * Z$u[r, i, TYPE_IDX]
    w <- sqrt((1 - s) * tau2_t)[TYPE_IDX] * Z$w[r, i, ]
    e + u + w
  }, numeric(12)))
}

## One replication, all conditions
sim_rep <- function(r, DESIGNS, TAU2, Z, cfg) {
  dims <- c(NM, 2, NH, length(METHODS), 6)
  out <- list(est = array(NA_real_, dims), se = array(NA_real_, dims),
              crit = array(NA_real_, dims), vb = array(NA_real_, c(NM, 2, NH, 6, 6)),
              df_cr2 = array(NA_real_, c(NM, 2, NH, 6)))
  for (mi in seq_len(NM)) for (ci in 1:2) for (h in seq_len(NH)) {
    D <- DESIGNS[[mi]][[ci]]
    Y <- gen_Y(r, D, TAU2[ci, h, ], cfg$het_share_study, Z)
    f <- fit_all(Y, D, cfg)
    out$est[mi, ci, h, , ]  <- f$est
    out$se[mi, ci, h, , ]   <- f$se
    out$crit[mi, ci, h, , ] <- f$crit
    out$vb[mi, ci, h, , ]   <- f$vb
    out$df_cr2[mi, ci, h, ] <- f$df_cr2
  }
  out
}

## =============================================================================
## 5. OBSERVED DATA: estimates and critical values given the actual studies
## =============================================================================

DESIGNS <- lapply(M_VALUES, function(m) setNames(lapply(COMPONENTS, make_design, m = m),
                                                 COMPONENTS))
OBS <- lapply(COMPONENTS, function(cmp) {
  Y <- do.call(rbind, INPUT$est[[cmp]])
  fit_all(Y, DESIGNS[[1]][[cmp]], CFG)
})
names(OBS) <- COMPONENTS

if (!is.null(INPUT$obs_ref)) {
  d <- max(abs(unlist(lapply(COMPONENTS, function(cmp)
    OBS[[cmp]]$est["mv_z", ] - INPUT$obs_ref[[cmp]][EFFECTS]))))
  message(sprintf("Consistency with apim_meta_analysis.R (multilevel, V = %s): max |diff| = %.2e",
                  CFG$V_type, d))
}

## =============================================================================
## 6. VALIDATION CHECKS
## =============================================================================

message("== Validation checks")

## (a) Hand-computed CE / KH / KH ad hoc equal metafor's own implementations

chk_uni <- local({
  D <- DESIGNS[[1]]$uncorrected; Y <- do.call(rbind, INPUT$est$uncorrected)
  j <- 1; y <- rowSums(Y[, EFF_COLS[[j]]] * D$comp_w[[j]]); v <- D$comp_v[, j]
  f  <- fit_all(Y, D, CFG)
  mk <- function(...) suppressWarnings(metafor::rma(y, v, method = CFG$tau2_method,
                                                    control = list(stepadj = .5, maxiter = 1000), ...))
  c(CE   = abs(metafor::rma(y, v, method = "EE")$se - f$se["rma_CE", j]),
    RE   = abs(mk()$se - f$se["rma_RE", j]),
    KH   = abs(mk(test = "knha")$se - f$se["rma_KH", j]),
    KHad = abs(mk(test = "adhoc")$se - f$se["rma_KHadh", j]))
})
if (max(chk_uni) > 1e-6) warning("Hand-computed rma() SEs deviate from metafor: ",
                                 paste(signif(chk_uni, 3), collapse = ", "))
message(sprintf("   rma(): hand-computed vs metafor SEs, max |diff| = %.1e", max(chk_uni)))

## (b) Location invariance: shifting all true effects by delta shifts every
##     estimate by delta and leaves SEs / critical values unchanged

set.seed(CFG$seed + 1)

chk_shift <- local({
  D  <- DESIGNS[[1]]$corrected
  Zt <- list(e = array(rnorm(12 * D$k), c(1, D$k, 12)),
             u = array(rnorm(6 * D$k),  c(1, D$k, 6)),
             w = array(rnorm(12 * D$k), c(1, D$k, 12)))
  Y0 <- gen_Y(1, D, TAU2["corrected", NH, ], CFG$het_share_study, Zt)
  delta <- c(.30, .10, .25, .05, .15, .20)[TYPE_IDX]
  f0 <- fit_all(Y0, D, CFG)
  f1 <- fit_all(sweep(Y0, 2, delta, "+"), D, CFG)
  c(est  = max(abs(f1$est - f0$est - matrix(c(.30, .10, .25, .05, .15, .20), 6, 6,
                                            byrow = TRUE)), na.rm = TRUE),
    se   = max(abs(f1$se - f0$se), na.rm = TRUE),
    crit = max(abs(f1$crit - f0$crit), na.rm = TRUE))
})

message(sprintf("   Location invariance: max |diff| est %.1e, SE %.1e, crit %.1e",
                chk_shift["est"], chk_shift["se"], chk_shift["crit"]))
if (max(chk_shift) > 1e-3) warning("Location invariance check failed; power ",
                                   "curves may be approximate (see Section 6).")

## =============================================================================
## 7. SIMULATION
## =============================================================================

set.seed(CFG$seed)

K_MAX <- K_BASE * max(M_VALUES)

Z <- list(e = array(rnorm(CFG$R * K_MAX * 12), c(CFG$R, K_MAX, 12)),
          u = array(rnorm(CFG$R * K_MAX * 6),  c(CFG$R, K_MAX, 6)),
          w = array(rnorm(CFG$R * K_MAX * 12), c(CFG$R, K_MAX, 12)))

message(sprintf("== Simulating: R = %d x %d component(s) x %d heterogeneity level(s) x %d k value(s) on %d core(s)",
                CFG$R, 2, NH, NM, CFG$n_cores))

t_start <- Sys.time()

if (CFG$n_cores > 1) {
  cl <- parallel::makeCluster(CFG$n_cores)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages({
    library(metafor); library(clubSandwich) }))
  parallel::clusterExport(cl, c("fit_all", "fit_mv", "tau2_uni", "gen_Y", "sim_rep",
                                "METHODS", "EFFECTS", "EFF_COLS", "TYPE_IDX", "NM", "NH"),
                          envir = environment())
  SIM <- tryCatch(parallel::parLapply(cl, seq_len(CFG$R), sim_rep, DESIGNS = DESIGNS,
                                      TAU2 = TAU2, Z = Z, cfg = CFG),
                  finally = parallel::stopCluster(cl))
} else {
  SIM <- vector("list", CFG$R)
  for (r in seq_len(CFG$R)) {
    SIM[[r]] <- sim_rep(r, DESIGNS, TAU2, Z, CFG)
    if (r %% 100 == 0) message(sprintf("   %d / %d (%.1f min)", r, CFG$R,
                                       as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}

message(sprintf("   done in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

stack_rep <- function(el) simplify2array(lapply(SIM, `[[`, el))   # R = last dim
EST0 <- stack_rep("est"); SE <- stack_rep("se"); CRIT <- stack_rep("crit"); VB <- stack_rep("vb")
DF_CR2 <- stack_rep("df_cr2")                         # [m, component, het, effect, r]
THR  <- CRIT * SE                                     # critical |estimate| per rep
## dims: [m, component, het, method, effect, r]

fail_rate <- apply(is.na(THR), c(1, 2, 3, 4), mean)
if (any(fail_rate > .01))
  warning("Some fits failed (> 1% of replications in at least one condition); ",
          "see RESULTS$fail_rate")

## =============================================================================
## 8. POWER CURVES, MDE, CRITICAL ESTIMATES, TYPE I ERROR (DIRECT PATHS)
## =============================================================================

THETA <- seq(0, CFG$theta_max, by = CFG$theta_step)

power_vec <- function(e0, thr, theta) {
  ok <- !is.na(e0) & !is.na(thr)
  if (!any(ok)) return(rep(NA_real_, length(theta)))
  colMeans(abs(outer(e0[ok], theta, "+")) > thr[ok])
}
mde_from_curve <- function(theta, pow, target) {
  i <- which(pow >= target)[1]
  if (is.na(i)) return(NA_real_)
  if (i == 1) return(theta[1])
  theta[i - 1] + (target - pow[i - 1]) * (theta[i] - theta[i - 1]) / (pow[i] - pow[i - 1])
}

grid_idx <- expand.grid(mi = seq_len(NM), ci = 1:2, h = seq_len(NH),
                        me = seq_along(METHODS), j = 1:6)

lab_cols <- function(g) data.frame(
  k = factor(K_LAB[g$mi], levels = K_LAB),
  component = factor(COMP_LAB[COMPONENTS[g$ci]], levels = COMP_LAB),
  het = factor(HET_LAB[g$h], levels = HET_LAB),
  method = factor(METHODS[g$me], levels = METHODS),
  effect = factor(EFFECT_LAB[EFFECTS[g$j]], levels = EFFECT_LAB),
  stringsAsFactors = FALSE)

POWER <- do.call(rbind, lapply(seq_len(nrow(grid_idx)), function(i) {
  g <- grid_idx[i, ]
  p <- power_vec(EST0[g$mi, g$ci, g$h, g$me, g$j, ], THR[g$mi, g$ci, g$h, g$me, g$j, ], THETA)
  cbind(lab_cols(g)[rep(1, length(THETA)), ], theta = THETA, power = p)
}))

SUMMARY <- do.call(rbind, lapply(seq_len(nrow(grid_idx)), function(i) {
  g   <- grid_idx[i, ]
  e0  <- EST0[g$mi, g$ci, g$h, g$me, g$j, ]
  thr <- THR[g$mi, g$ci, g$h, g$me, g$j, ]
  p   <- power_vec(e0, thr, THETA)
  cbind(lab_cols(g),
        tau    = sqrt(TAU2[g$ci, g$h, g$j]),
        mde    = mde_from_curve(THETA, p, CFG$power_target),
        mde50  = mde_from_curve(THETA, p, .50),
        type1  = p[1],
        crit_med = stats::median(thr, na.rm = TRUE),
        crit_q10 = unname(stats::quantile(thr, .10, na.rm = TRUE)),
        crit_q90 = unname(stats::quantile(thr, .90, na.rm = TRUE)),
        fail   = mean(is.na(thr)))
}))

## Satterthwaite df of CR2 in the simulation (median and 10th/90th percentile)

CR2_DF <- do.call(rbind, lapply(seq_len(NM), function(mi) do.call(rbind, lapply(1:2, function(ci)
  do.call(rbind, lapply(seq_len(NH), function(h) {
    d <- DF_CR2[mi, ci, h, , ]
    data.frame(k = K_LAB[mi], component = unname(COMP_LAB[COMPONENTS[ci]]), het = HET_LAB[h],
               df_median = stats::median(d, na.rm = TRUE),
               df_q10 = unname(stats::quantile(d, .1, na.rm = TRUE)),
               df_q90 = unname(stats::quantile(d, .9, na.rm = TRUE)))
  }))))))

## Critical pooled estimates and significance for the OBSERVED studies

OBS_TAB <- do.call(rbind, lapply(COMPONENTS, function(cmp) {
  o <- OBS[[cmp]]
  do.call(rbind, lapply(seq_along(METHODS), function(me) data.frame(
    component = factor(unname(COMP_LAB[cmp]), levels = COMP_LAB),
    method = factor(unname(METHODS[me]), levels = METHODS),
    effect = factor(unname(EFFECT_LAB[EFFECTS]), levels = EFFECT_LAB),
    est = o$est[me, ], se = o$se[me, ], crit_est = o$crit[me, ] * o$se[me, ],
    significant = abs(o$est[me, ]) > o$crit[me, ] * o$se[me, ],
    stringsAsFactors = FALSE)))
}))
OBS_REF <- subset(OBS_TAB, method == METHODS[["mv_z"]])   # reference lines in plots

## =============================================================================
## 9. INDIRECT EFFECTS: joint-significance power surfaces over (a, b)
## =============================================================================

TH_IND <- seq(0, CFG$ind_max, by = CFG$ind_step)
sig_mat <- function(e0, thr, theta) abs(outer(e0, theta, "+")) > thr   # R x G

IND <- list(); IND_MIN <- list()
for (mi in seq_len(NM)) for (ci in 1:2) for (h in seq_len(NH)) for (me in seq_along(METHODS)) {
  for (tp in names(INDIRECT)) {
    ja <- match(INDIRECT[[tp]][1], EFFECTS); jb <- match(INDIRECT[[tp]][2], EFFECTS)
    ea <- EST0[mi, ci, h, me, ja, ]; ta <- THR[mi, ci, h, me, ja, ]
    eb <- EST0[mi, ci, h, me, jb, ]; tb <- THR[mi, ci, h, me, jb, ]
    ok <- !is.na(ea + ta + eb + tb)
    if (!any(ok)) next
    SA <- sig_mat(ea[ok], ta[ok], TH_IND) * 1
    SB <- sig_mat(eb[ok], tb[ok], TH_IND) * 1
    P  <- crossprod(SA, SB) / sum(ok)                   # G x G: rows a, cols b
    key <- paste(mi, ci, h, me, tp)
    IND[[key]] <- data.frame(
      k = K_LAB[mi], component = unname(COMP_LAB[COMPONENTS[ci]]), het = HET_LAB[h],
      method = unname(METHODS[me]), type = tp,
      a = rep(TH_IND, times = length(TH_IND)), b = rep(TH_IND, each = length(TH_IND)),
      power = as.vector(P), stringsAsFactors = FALSE)
    prod_ab <- outer(TH_IND, TH_IND)
    okc <- P >= CFG$power_target
    diag_p <- diag(P)
    IND_MIN[[key]] <- data.frame(
      k = K_LAB[mi], component = unname(COMP_LAB[COMPONENTS[ci]]), het = HET_LAB[h],
      method = unname(METHODS[me]), type = tp,
      min_ab = if (any(okc)) min(prod_ab[okc]) else NA_real_,
      at_a = if (any(okc)) TH_IND[which(okc & prod_ab == min(prod_ab[okc]), arr.ind = TRUE)[1, 1]] else NA_real_,
      at_b = if (any(okc)) TH_IND[which(okc & prod_ab == min(prod_ab[okc]), arr.ind = TRUE)[1, 2]] else NA_real_,
      diag_a = mde_from_curve(TH_IND, diag_p, CFG$power_target),
      type1_ab0 = P[1, 1], stringsAsFactors = FALSE)
  }
}

IND <- do.call(rbind, IND); IND_MIN <- do.call(rbind, IND_MIN)
rownames(IND) <- rownames(IND_MIN) <- NULL
for (nm in c("k", "component", "het", "method", "type")) {
  lev <- list(k = K_LAB, component = COMP_LAB, het = HET_LAB, method = METHODS,
              type = names(INDIRECT))[[nm]]
  IND[[nm]] <- factor(IND[[nm]], levels = lev); IND_MIN[[nm]] <- factor(IND_MIN[[nm]], levels = lev)
}

IND_MIN$diag_ab <- IND_MIN$diag_a^2

## Observed (a, b) points for the plots (multilevel z estimates)

OBS_IND <- do.call(rbind, lapply(COMPONENTS, function(cmp) do.call(rbind, lapply(names(INDIRECT), function(tp)
  data.frame(component = factor(unname(COMP_LAB[cmp]), levels = COMP_LAB),
             type = factor(tp, levels = names(INDIRECT)),
             a = abs(unname(OBS[[cmp]]$est["mv_z", INDIRECT[[tp]][1]])),
             b = abs(unname(OBS[[cmp]]$est["mv_z", INDIRECT[[tp]][2]])))))))

## Optional: joint significance vs Monte Carlo CI (rma.mv z and rma RE z)

MC_CHECK <- NULL

if (isTRUE(CFG$validate_mc)) {
  message("== Checking joint significance against Monte Carlo CIs")
  set.seed(CFG$seed + 2)
  Zmc <- matrix(rnorm(2 * CFG$n_mc), ncol = 2)
  Rm  <- min(CFG$R_mc, CFG$R)
  MC_CHECK <- do.call(rbind, lapply(c("mv_z", "rma_RE"), function(mcode) {
    me <- match(mcode, names(METHODS))
    do.call(rbind, lapply(1:2, function(ci) do.call(rbind, lapply(c(1, NH), function(h) {
      do.call(rbind, lapply(names(INDIRECT)[1:2], function(tp) {
        ja <- match(INDIRECT[[tp]][1], EFFECTS); jb <- match(INDIRECT[[tp]][2], EFFECTS)
        sub <- subset(IND, k == K_LAB[1] & component == unname(COMP_LAB[COMPONENTS[ci]]) &
                        het == HET_LAB[h] & method == unname(METHODS[mcode]) & type == tp & a == b)
        pts <- vapply(c(.3, .6, .8), function(p) {
          i <- which(sub$power >= p)[1]; if (is.na(i)) NA_real_ else sub$a[i] }, numeric(1))
        pts <- unique(pts[!is.na(pts)])
        do.call(rbind, lapply(pts, function(th) {
          rej <- vapply(seq_len(Rm), function(r) {
            mu <- c(th + EST0[1, ci, h, me, ja, r], th + EST0[1, ci, h, me, jb, r])
            S  <- if (mcode == "mv_z") VB[1, ci, h, c(ja, jb), c(ja, jb), r] else
                    diag(SE[1, ci, h, me, c(ja, jb), r]^2)
            if (anyNA(mu) || anyNA(S)) return(NA)
            L  <- tryCatch(chol(S), error = function(e) NULL)
            if (is.null(L)) return(NA)
            d  <- sweep(Zmc %*% L, 2, mu, "+")
            ci_ab <- stats::quantile(d[, 1] * d[, 2], c(CFG$alpha / 2, 1 - CFG$alpha / 2))
            ci_ab[1] > 0 || ci_ab[2] < 0
          }, logical(1))
          js <- sub$power[match(th, sub$a)]
          data.frame(method = unname(METHODS[mcode]), component = unname(COMP_LAB[COMPONENTS[ci]]),
                     het = HET_LAB[h], type = tp, a = th, b = th,
                     power_joint_sig = js, power_mc_ci = mean(rej, na.rm = TRUE),
                     R_mc = sum(!is.na(rej)))
        }))
      }))
    }))))
  }))
  print(transform(MC_CHECK, power_joint_sig = round(power_joint_sig, 3),
                  power_mc_ci = round(power_mc_ci, 3)), row.names = FALSE)
}

## =============================================================================
## 10. CONSOLE SUMMARY
## =============================================================================

cat("\n==== Minimum detectable pooled effect (", round(100 * CFG$power_target),
    "% power), first k =", K_LAB[1], "====\n")
for (cmp in COMP_LAB) {
  cat("\n--", cmp, "--\n")
  s <- subset(SUMMARY, k == K_LAB[1] & component == cmp)
  w <- reshape(s[, c("method", "effect", "het", "mde")], idvar = c("method", "effect"),
               timevar = "het", direction = "wide")
  names(w) <- sub("^mde\\.", "", names(w))
  w[, -(1:2)] <- lapply(w[, -(1:2)], function(x) ifelse(is.na(x), sprintf(">%.2f", CFG$theta_max),
                                                        sprintf("%.3f", x)))
  print(w[order(w$effect, w$method), ], row.names = FALSE)
}

## =============================================================================
## 11. FIGURES
## =============================================================================
theme_apa <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = "grey70"),
          legend.position = "bottom", legend.title = element_blank(),
          legend.key.width = unit(1.4, "lines"))
}
save_fig <- function(p, name, w, h) {
  f_png <- paste0(CFG$out_prefix, "_", name, ".png")
  ggsave(f_png, p, width = w, height = h, dpi = CFG$fig_dpi)
  if (capabilities("cairo"))
    ggsave(paste0(CFG$out_prefix, "_", name, ".pdf"), p, width = w, height = h,
           device = grDevices::cairo_pdf)
  structure(f_png, asp = h / w)
}
k1 <- K_LAB[1]
FIGS <- list()

## Figure 1-2: power curves for the direct paths

for (ci in 1:2) {
  cmp <- unname(COMP_LAB[COMPONENTS[ci]])
  d   <- subset(POWER, k == k1 & component == cmp)
  ref <- subset(OBS_REF, component == cmp)
  p <- ggplot(d, aes(theta, power, colour = method)) +
    geom_hline(yintercept = CFG$power_target, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = CFG$alpha, linetype = "dotted", colour = "grey55") +
    geom_vline(data = ref, aes(xintercept = abs(est)), linetype = "dotdash",
               colour = "grey25", linewidth = .4) +
    geom_line(linewidth = .55) +
    facet_grid(effect ~ het) +
    scale_colour_manual(values = PAL, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = c(.05, .5, .8, 1),
                       labels = c(".05", ".50", ".80", "1")) +
    labs(x = "True pooled standardized effect (\u03bc)", y = "Power (p < .05)") +
    guides(colour = guide_legend(nrow = 2)) + theme_apa()
  FIGS[[length(FIGS) + 1]] <- list(
    file = save_fig(p, sprintf("fig%d_power_direct_%s", ci, COMPONENTS[ci]), 9, 10),
    title = sprintf("Power to Detect the Pooled Direct Effects: %s Paths", cmp),
    note = sprintf(paste0("%s. Dashed line = %d%% power; dotted line = \u03b1 = .05; ",
                          "dot-dashed line = observed pooled estimate (multilevel model). ",
                          "Power under the random-effects data-generating model; %d replications per condition."),
                   k1, round(100 * CFG$power_target), CFG$R))
}

## Figure 3: minimum detectable effect vs heterogeneity

md <- subset(SUMMARY, k == k1)

md$reached  <- !is.na(md$mde)

md$mde_plot <- ifelse(md$reached, md$mde, CFG$theta_max)

pd <- position_dodge(width = .5)

p3 <- ggplot(md, aes(het, mde_plot, colour = method, group = method)) +
  geom_hline(data = transform(subset(OBS_REF, TRUE), val = abs(est)),
             aes(yintercept = val), linetype = "dotdash", colour = "grey25", linewidth = .4) +
  geom_line(data = subset(md, reached), position = pd, linewidth = .45) +
  geom_point(aes(shape = reached), position = pd, size = 1.7) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 2),
                     labels = c(`TRUE` = "MDE reached",
                                `FALSE` = sprintf("Not reached (> %.2f)", CFG$theta_max))) +
  facet_grid(component ~ effect) +
  scale_colour_manual(values = PAL, drop = FALSE) +
  coord_cartesian(ylim = c(0, CFG$theta_max)) +
  labs(x = NULL, y = sprintf("Minimum detectable pooled effect (%d%% power)",
                             round(100 * CFG$power_target))) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme_apa() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

FIGS[[length(FIGS) + 1]] <- list(
  file = save_fig(p3, "fig3_mde_direct", 11, 6.5),
  title = "Minimum Detectable Pooled Direct Effects by Heterogeneity and Method",
  note = sprintf(paste0("%s. MDE = smallest true pooled effect with %d%% power. Open triangles: ",
                        "not reached within the evaluated range (\u2264 %.2f). Dot-dashed line = ",
                        "observed pooled estimate."), k1, round(100 * CFG$power_target), CFG$theta_max))

## Figure 4: critical pooled estimate (distribution across replications)

p4 <- ggplot(md, aes(het, crit_med, colour = method, group = method)) +
  geom_hline(data = transform(OBS_REF, val = abs(est)), aes(yintercept = val),
             linetype = "dotdash", colour = "grey25", linewidth = .4) +
  geom_errorbar(aes(ymin = crit_q10, ymax = pmin(crit_q90, 2 * CFG$theta_max)),
                width = .25, position = pd, linewidth = .4) +
  geom_point(position = pd, size = 1.6) +
  facet_grid(component ~ effect) +
  scale_colour_manual(values = PAL, drop = FALSE) +
  coord_cartesian(ylim = c(0, CFG$theta_max)) +
  labs(x = NULL, y = "Critical |pooled estimate| for p < .05") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_apa() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

FIGS[[length(FIGS) + 1]] <- list(
  file = save_fig(p4, "fig4_critical_estimate", 11, 6.5),
  title = "Smallest Pooled Estimate That Would Be Statistically Significant",
  note = sprintf(paste0("%s. Points = median and bars = 10th\u201390th percentile of crit \u00d7 SE ",
                        "across replications (upper bars truncated at the axis). For z tests with ",
                        "fixed V and \u03c4\u00b2 = 0 the value is constant; for random-effects, KH ",
                        "and CR2 it depends on the observed dispersion."), k1))

## Figure 5: empirical Type I error at mu = 0

t1 <- aggregate(type1 ~ k + component + het + method, data = SUMMARY,
                FUN = function(x) c(mean = mean(x), min = min(x), max = max(x)))
t1 <- cbind(t1[, 1:4], as.data.frame(t1$type1))

p5 <- ggplot(subset(t1, k == k1), aes(het, mean, colour = method, group = method)) +
  geom_hline(yintercept = CFG$alpha, linetype = "dashed", colour = "grey40") +
  geom_line(position = pd, linewidth = .45) +
  geom_pointrange(aes(ymin = min, ymax = max), position = pd, size = .25) +
  facet_wrap(~ component) +
  scale_colour_manual(values = PAL, drop = FALSE) +
  labs(x = NULL, y = "Rejection rate when \u03bc = 0") +
  guides(colour = guide_legend(nrow = 2)) + theme_apa()

FIGS[[length(FIGS) + 1]] <- list(
  file = save_fig(p5, "fig5_type1_error", 9, 4.8),
  title = "Empirical Type I Error Rate of the Test of the Pooled Mean Effect",
  note = sprintf(paste0("%s. Mean (range) across the six direct effects. Under heterogeneity, ",
                        "the common-effect test answers a conditional question (the weighted ",
                        "average of these studies), so rejection rates above .05 reflect a ",
                        "different inferential target rather than a coding error."), k1))

## Figures 6-7: indirect-effect 80% power contours

hyper <- do.call(rbind, lapply(c(.02, .05, .10, .20), function(cc) {
  a <- seq(cc / CFG$ind_max, CFG$ind_max, length.out = 200)
  data.frame(a = a, b = cc / a, lab = sprintf("ab = %.2f", cc))
}))

hyper_lab <- do.call(rbind, lapply(split(hyper, hyper$lab), function(d) d[nrow(d), ]))
for (ci in 1:2) {
  cmp <- unname(COMP_LAB[COMPONENTS[ci]])
  d <- subset(IND, k == k1 & component == cmp)
  ## contours only where the target is reached (avoids empty-contour warnings)
  reach <- aggregate(power ~ het + method + type, data = d, FUN = max)
  reach <- reach[reach$power >= CFG$power_target, c("het", "method", "type")]
  d <- merge(d, reach, by = c("het", "method", "type"))
  p <- ggplot(d, aes(a, b)) +
    geom_line(data = hyper, aes(a, b, group = lab), colour = "grey80", linewidth = .35) +
    geom_text(data = hyper_lab, aes(a, b, label = lab), colour = "grey55", size = 2,
              hjust = 1, vjust = -0.4) +
    geom_contour(aes(z = power, colour = method), breaks = CFG$power_target,
                 linewidth = .6) +
    geom_point(data = subset(OBS_IND, component == cmp), aes(a, b), shape = 4,
               size = 2, colour = "black", inherit.aes = FALSE) +
    facet_grid(type ~ het) +
    scale_colour_manual(values = PAL, drop = FALSE) +
    coord_equal(xlim = c(0, CFG$ind_max), ylim = c(0, CFG$ind_max)) +
    labs(x = "True pooled a path", y = "True pooled b path") +
    guides(colour = guide_legend(nrow = 2)) + theme_apa()
  FIGS[[length(FIGS) + 1]] <- list(
    file = save_fig(p, sprintf("fig%d_indirect_%s", 5 + ci, COMPONENTS[ci]), 9, 10),
    title = sprintf("Combinations of Pooled a and b Paths Giving %d%% Power for the Indirect Effect: %s",
                    round(100 * CFG$power_target), cmp),
    note = sprintf(paste0("%s. Lines = %d%% power contours of the joint significance test ",
                          "(combinations above and to the right are detectable). Grey curves = ",
                          "constant indirect effect ab. \u00d7 = observed pooled (a, b). Methods ",
                          "without a contour never reach %d%% power in the plotted range."),
                   k1, round(100 * CFG$power_target), round(100 * CFG$power_target)))
}

## Figure 8 (only if several k): MDE as a function of the number of studies
if (NM > 1) {
  p8 <- ggplot(subset(SUMMARY, !is.na(mde)), aes(as.integer(sub("k = ", "", as.character(k))), mde,
                            colour = method, linetype = het)) +
    geom_line(linewidth = .45) + geom_point(size = 1) +
    facet_grid(component ~ effect) +
    scale_colour_manual(values = PAL, drop = FALSE) +
    scale_linetype_manual(values = rep(c("solid", "dotted", "dashed", "twodash", "longdash"),
                                       length.out = NH)) +
    labs(x = "Number of studies (k)", y = sprintf("MDE (%d%% power)", round(100 * CFG$power_target))) +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 1)) + theme_apa()
  FIGS[[length(FIGS) + 1]] <- list(
    file = save_fig(p8, "fig8_mde_by_k", 11, 6.5),
    title = "Minimum Detectable Pooled Direct Effect as a Function of the Number of Studies",
    note = "Additional studies are hypothetical replicates with the sampling covariance matrices of the observed studies.")
}

## =============================================================================
## 12. EXPORT: CSV / RDS / APA .docx
## =============================================================================
utils::write.csv(SUMMARY, paste0(CFG$out_prefix, "_direct_summary.csv"), row.names = FALSE)
utils::write.csv(POWER,   paste0(CFG$out_prefix, "_direct_power_curves.csv"), row.names = FALSE)
utils::write.csv(IND_MIN, paste0(CFG$out_prefix, "_indirect_summary.csv"), row.names = FALSE)
utils::write.csv(OBS_TAB, paste0(CFG$out_prefix, "_observed_critical.csv"), row.names = FALSE)

if (!is.null(MC_CHECK))
  utils::write.csv(MC_CHECK, paste0(CFG$out_prefix, "_mc_vs_joint.csv"), row.names = FALSE)

RESULTS <- list(config = CFG, inputs = INPUT, tau2 = TAU2, v_typical = v_typical,
                summary = SUMMARY, power = POWER, indirect = IND, indirect_min = IND_MIN,
                observed = OBS_TAB, cr2_df = CR2_DF, mc_check = MC_CHECK,
                fail_rate = fail_rate,
                validation = list(rma_vs_metafor = chk_uni, location_invariance = chk_shift))
saveRDS(RESULTS, paste0(CFG$out_prefix, "_results.rds"))

## APA helpers (as in apim_meta_analysis.R)

MINUS <- "\u2212"

fmt <- function(x, d = 2, lead0 = TRUE) {
  s <- formatC(round(x, d), format = "f", digits = d)
  s <- sub("^-(0\\.0+)$", "\\1", s)
  if (!lead0) s <- sub("^(-?)0\\.", "\\1.", s)
  s <- sub("^-", MINUS, s); s[is.na(x)] <- "\u2014"; s
}

apa_table <- function(body, labels, spans = list(), group_rows = integer(0), note = NULL,
                      first_col_width = 2.4, col_width = 0.85,
                      font = CFG$table_font, size = CFG$table_size) {
  bd <- officer::fp_border(color = "black", width = 1)
  nc <- ncol(body)
  names(body) <- paste0("col", seq_len(nc))          # safe internal keys
  ft <- flextable::flextable(body)
  ft <- flextable::set_header_labels(ft, values = setNames(as.list(labels), names(body)))
  for (sp in rev(spans)) ft <- flextable::add_header_row(ft, values = names(sp),
                                                         colwidths = unname(sp), top = TRUE)
  ft <- flextable::border_remove(ft)
  ft <- flextable::font(ft, fontname = font, part = "all")
  ft <- flextable::fontsize(ft, size = size, part = "all")
  ft <- flextable::align(ft, align = "center", part = "header")
  ft <- flextable::align(ft, j = 1, align = "left", part = "all")
  if (nc > 1) ft <- flextable::align(ft, j = 2:nc, align = "center", part = "body")
  pad <- if (nrow(body) > 24) 0.75 else 1.5
  ft <- flextable::padding(ft, padding.top = pad, padding.bottom = pad, part = "all")
  ft <- flextable::hline_top(ft, border = bd, part = "header")
  ft <- flextable::hline_bottom(ft, border = bd, part = "header")
  ft <- flextable::hline_bottom(ft, border = bd, part = "body")
  for (r in seq_along(spans)) {
    sp <- spans[[r]]; start <- 1
    for (s in seq_along(sp)) {
      if (nzchar(trimws(names(sp)[s])))
        ft <- flextable::hline(ft, i = r, j = start:(start + sp[s] - 1),
                               border = officer::fp_border(color = "black", width = .75),
                               part = "header")
      start <- start + sp[s]
    }
  }
  if (length(group_rows)) {
    ft <- flextable::padding(ft, i = setdiff(seq_len(nrow(body)), group_rows), j = 1,
                             padding.left = 14, part = "body")
    for (g in group_rows) {
      ft <- flextable::merge_at(ft, i = g, j = 1:nc, part = "body")
      ft <- flextable::italic(ft, i = g, j = 1, part = "body")
    }
  }
  if (!is.null(note)) {
    ft <- flextable::add_footer_lines(ft, values = flextable::as_paragraph(flextable::as_i("Note. "), note))
    ft <- flextable::font(ft, fontname = font, part = "footer")
    ft <- flextable::fontsize(ft, size = size, part = "footer")
  }
  ft <- flextable::width(ft, j = 1, width = first_col_width)
  if (nc > 1) ft <- flextable::width(ft, j = 2:nc, width = col_width)
  flextable::set_table_properties(ft, layout = "fixed")
}

title_par_xml <- function(text, bold = FALSE, italic = FALSE, page_break = FALSE,
                          font = CFG$table_font, size_pt = 12) {
  txt <- gsub("&", "&amp;", text, fixed = TRUE)
  txt <- gsub("<", "&lt;", txt, fixed = TRUE); txt <- gsub(">", "&gt;", txt, fixed = TRUE)
  rpr <- sprintf(paste0('<w:rPr><w:rFonts w:ascii="%1$s" w:hAnsi="%1$s" w:cs="%1$s"/>',
                        '%2$s%3$s<w:sz w:val="%4$d"/><w:szCs w:val="%4$d"/></w:rPr>'),
                 font, if (bold) "<w:b/>" else "", if (italic) "<w:i/>" else "",
                 as.integer(size_pt * 2))
  ppr <- sprintf('<w:pPr><w:keepNext/>%s<w:spacing w:before="0" w:after="120"/></w:pPr>',
                 if (page_break) "<w:pageBreakBefore/>" else "")
  sprintf(paste0('<w:p xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">',
                 '%s<w:r>%s<w:t xml:space="preserve">%s</w:t></w:r></w:p>'), ppr, rpr, txt)
}

add_caption <- function(doc, kind, number, title, new_page = TRUE) {
  doc <- officer::body_add_xml(doc, title_par_xml(paste(kind, number), bold = TRUE,
                                                  page_break = new_page))
  officer::body_add_xml(doc, title_par_xml(title, italic = TRUE))
}

add_apa_table <- function(doc, number, title, ft, new_page = TRUE) {
  doc <- add_caption(doc, "Table", number, title, new_page)
  w <- dim(ft)$widths
  if (sum(w) > CFG$max_table_width)
    ft <- flextable::width(ft, j = seq_along(w), width = w * CFG$max_table_width / sum(w))
  flextable::body_add_flextable(doc, ft, align = "left")
}

grouped_body <- function(df, group, label_col = 1) {
  out <- NULL; g <- integer(0)
  for (lv in unique(group)) {
    hr <- df[1, , drop = FALSE]; hr[] <- ""; hr[[label_col]] <- as.character(lv)
    out <- rbind(out, hr, df[group == lv, , drop = FALSE]); g <- c(g, nrow(out) - sum(group == lv))
  }
  rownames(out) <- NULL; list(body = out, group_rows = g)
}

max_lab <- sprintf("> %.2f", CFG$theta_max)

nh <- NH

doc <- officer::read_docx()

doc <- officer::body_set_default_section(doc, officer::prop_section(
  page_size = officer::page_size(orient = "landscape"),
  page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1)))

##  Table 1: MDE for direct paths (rows = heterogeneity within method,
##  columns = effect within component; fits one landscape page)

s1 <- subset(SUMMARY, k == k1)

t1b <- expand.grid(het = HET_LAB, method = METHODS, stringsAsFactors = FALSE)

for (cmp in COMP_LAB) for (e in EFFECT_LAB) {
  x <- s1[s1$component == cmp & s1$effect == e, ]
  v <- x$mde[match(paste(t1b$method, t1b$het), paste(x$method, x$het))]
  t1b[[paste(cmp, e)]] <- ifelse(is.na(v), max_lab, fmt(v))
}

g1 <- grouped_body(t1b[, -2], t1b$method)

ft1 <- apa_table(g1$body,
                 labels = c("Method / heterogeneity", rep(unname(EFFECT_LAB), 2)),
                 spans = list(c(" " = 1L, setNames(c(6L, 6L), COMP_LAB))),
                 group_rows = g1$group_rows,
                 note = sprintf(paste0("%s (%s dyads). Entries are the smallest true pooled ",
                                       "standardized effect detected with %d%% power (two-sided ",
                                       "\u03b1 = %.2f), %d replications; \"%s\" = not reached in ",
                                       "the evaluated range. Uncorrected and corrected effects are on ",
                                       "their own scales, so their entries are not directly comparable. ",
                                       "CE = common effect; RE = random effects (REML); KH = Knapp\u2013",
                                       "Hartung; CR2 = cluster-robust with Satterthwaite df. ",
                                       "Heterogeneity refers to the study-level effects (\u03c4\u00b2 = ",
                                       "I\u00b2/(1 \u2212 I\u00b2) \u00d7 typical sampling variance)."),
                                k1, paste(INPUT$n_dyads, collapse = " and "),
                                round(100 * CFG$power_target), CFG$alpha, CFG$R, max_lab),
                 first_col_width = 1.8, col_width = 0.59, size = CFG$table_size - 1)

doc <- add_caption(doc, "Table", 1, "Minimum Detectable Pooled Direct Effects", new_page = FALSE)

w <- dim(ft1)$widths

if (sum(w) > CFG$max_table_width) ft1 <- flextable::width(ft1, j = seq_along(w), width = w * CFG$max_table_width / sum(w))

doc <- flextable::body_add_flextable(doc, ft1, align = "left")

## Table 2: observed data - critical pooled estimates per method

t2b <- expand.grid(effect = EFFECT_LAB, component = COMP_LAB, stringsAsFactors = FALSE)

t2b$beta <- vapply(seq_len(nrow(t2b)), function(i) {
  x <- OBS_TAB[OBS_TAB$component == t2b$component[i] & OBS_TAB$effect == t2b$effect[i] &
                 OBS_TAB$method == unname(METHODS["mv_z"]), ]
  fmt(x$est[1])
}, character(1))
for (me in METHODS) {
  t2b[[me]] <- vapply(seq_len(nrow(t2b)), function(i) {
    x <- OBS_TAB[OBS_TAB$component == t2b$component[i] & OBS_TAB$effect == t2b$effect[i] &
                   OBS_TAB$method == me, ]
    paste0(fmt(x$crit_est[1]), if (isTRUE(x$significant[1])) "*" else "")
  }, character(1))
}

g2 <- grouped_body(t2b[, -2], t2b$component)

ft2 <- apa_table(g2$body, labels = c("Effect", "\u03b2", unname(METHODS)),
                 spans = list(c(" " = 1L, " " = 1L, "Critical |\u03b2| for p < .05" = 6L)),
                 group_rows = g2$group_rows,
                 note = paste0("\u03b2 = observed pooled estimate (multilevel model; the point estimate ",
                               "is nearly identical across methods, see the CSV). Remaining columns give ",
                               "the smallest absolute pooled estimate that would have been significant ",
                               "under each method with these data (critical value \u00d7 SE). ",
                               "* the observed \u03b2 exceeds that value, i.e. the effect is significant ",
                               "under that method."),
                 first_col_width = 1.5, col_width = 1.05)

doc <- add_apa_table(doc, 2, "Observed Pooled Estimates and Critical Values for Significance", ft2)

## Table 3: indirect effects - minimum ab with 80% power

i1 <- subset(IND_MIN, k == k1)

t3b <- expand.grid(method = METHODS, type = names(INDIRECT), stringsAsFactors = FALSE)

for (cmp in COMP_LAB) for (h in HET_LAB) {
  x <- i1[i1$component == cmp & i1$het == h, ]
  ix <- match(paste(t3b$method, t3b$type), paste(x$method, x$type))
  t3b[[paste(cmp, h)]] <- ifelse(is.na(x$min_ab[ix]), "\u2014", fmt(x$min_ab[ix], 3))
}

g3 <- grouped_body(t3b[, -2], t3b$type)

ft3 <- apa_table(g3$body, labels = c("Indirect effect / method", rep(HET_LAB, 2)),
                 spans = list(c(" " = 1L, setNames(c(nh, nh), COMP_LAB))),
                 group_rows = g3$group_rows,
                 note = sprintf(paste0("Smallest true indirect effect (product ab of pooled paths) ",
                                       "over the evaluated grid (a, b \u2264 %.2f) that is detected with ",
                                       "%d%% power by the joint significance test. \u2014 = not reached. ",
                                       "The full (a, b) combinations are shown in Figures 6 and 7."),
                                CFG$ind_max, round(100 * CFG$power_target)),
                 first_col_width = 2.4, col_width = 0.78)

doc <- add_apa_table(doc, 3, "Minimum Detectable Pooled Indirect Effects", ft3)

## Table 4: Type I error

t4b <- expand.grid(method = METHODS, stringsAsFactors = FALSE)

for (cmp in COMP_LAB) for (h in HET_LAB) {
  x <- subset(t1, k == k1 & component == cmp & het == h)
  t4b[[paste(cmp, h)]] <- fmt(x$mean[match(t4b$method, x$method)], 3, lead0 = FALSE)
}

ft4 <- apa_table(t4b, labels = c("Method", rep(HET_LAB, 2)),
                 spans = list(c(" " = 1L, setNames(c(nh, nh), COMP_LAB))),
                 note = sprintf(paste0("Proportion of replications with p < %.2f when the true pooled ",
                                       "effect is zero, averaged over the six direct effects. CR2 ",
                                       "Satterthwaite df (median across conditions) = %s."),
                                CFG$alpha, fmt(stats::median(CR2_DF$df_median[CR2_DF$k == k1]), 2)),
                 first_col_width = 2.4, col_width = 0.78)
doc <- add_apa_table(doc, 4, "Empirical Type I Error Rates", ft4)

## Figures

for (i in seq_along(FIGS)) {
  doc <- add_caption(doc, "Figure", i, FIGS[[i]]$title)
  asp <- attr(FIGS[[i]]$file, "asp")
  wd <- min(8.9, 6.3 / asp)
  doc <- officer::body_add_img(doc, FIGS[[i]]$file, width = wd, height = wd * asp)
  doc <- officer::body_add_fpar(doc, officer::fpar(
    officer::ftext("Note. ", officer::fp_text(italic = TRUE, font.family = CFG$table_font, font.size = 10)),
    officer::ftext(FIGS[[i]]$note, officer::fp_text(font.family = CFG$table_font, font.size = 10))))
}

print(doc, target = paste0(CFG$out_prefix, "_report.docx"))
message("\n== Done. Figures: ", CFG$out_prefix, "_fig*.png/.pdf; report: ",
        CFG$out_prefix, "_report.docx; results: ", CFG$out_prefix, "_*.csv/.rds")
