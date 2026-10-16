
## title: "APIM meta-analysis: individual-level missing and final pipeline"
## author: "Anonymized"
## date: "15-10-2026"

# =============================================================================
# 00_functions.R
# Shared library for the item-level missing-data pipeline of the APIM
# meta-analysis (two dyadic studies, N = 130 and N = 174 dyads). Single-level
# regression to the Actor-Partner Interdependence Mediation Model (APIMeM) 
# and its meta-analysis.
#
#   1. Design: studies, items, scale scores; the APIMeM functions of
#      apim_meta_analysis.R (sourced, not copied)
#   2. Toy population (simulate_toy_data.R) and the complete toy data
#   3. Auxiliary variables and instruments
#   4. Item-level amputation (MCAR / MAR / MNAR, exact counts, nested rates)
#   5. Custom mice methods: catx, pmmx, heckoprob
#   6. Item-level imputation: mice-PMM, mice-CAT, Heckman-MI, delta-adjusted PMM
#   7. Study-level APIMeM in lavaan: complete/imputed data, proration + FIML,
#      two-stage item-level FIML, Rubin's rules within each study
#   8. Pooling: rma() random-effects meta-analysis of the partner-pooled
#      effects (as in apim_meta_sensitivity.R), indirect effects
#   9. Performance measures and the method decision (mice-PMM vs FIML)
#  10. Power sensitivity: apim_meta_sensitivity.R rerun with the missing-data
#      sampling covariances, and rescaling of its published results
#  11. Output helpers (xlsx, APA-style flextables)
#
# The X, M and Y variables are scale scores (means of 4 Likert items), so the
# missing values are at the ITEM level. Every method below therefore works on
# the items and only then forms the scale scores that enter the APIMeM.
#
# Requirements: R >= 4.1; lavaan, metafor, clubSandwich, mice, MASS, pbivnorm,
# numDeriv, openxlsx, flextable, officer (future + future.apply for WORKERS > 1).
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

# Logical switch read from an environment variable (TRUE/FALSE, T/F, 1/0, yes/no, on/off)

env_flag <- function(name, default = TRUE) {
  v <- tolower(trimws(Sys.getenv(name, "")))
  if (!nzchar(v)) return(default)
  if (v %in% c("true", "t", "1", "yes", "y", "on")) return(TRUE)
  if (v %in% c("false", "f", "0", "no", "n", "off")) return(FALSE)
  stop("Environment variable ", name, " = '", v, "': use TRUE or FALSE (or 1/0)")
}

# Random numbers. The scripts set their seeds for R's default generator
# (Mersenne-Twister, Inversion, Rejection); another generator left in the session
# by other code (e.g. the L'Ecuyer-CMRG streams of the parallel workers) would
# change every result, so sourcing this file restores the default. Imputations
# draw from their own seeded L'Ecuyer-CMRG stream (with_rng(), impute_mice()),
# which makes them identical with and without parallel workers.

RNG_DEFAULT    <- c("Mersenne-Twister", "Inversion", "Rejection")

RNG_IMPUTATION <- c("L'Ecuyer-CMRG", "Inversion", "Rejection")

if (!identical(RNGkind()[1:3], RNG_DEFAULT)) RNGkind(RNG_DEFAULT[1], RNG_DEFAULT[2], RNG_DEFAULT[3])

# Evaluate expr with generator 'kind' seeded with 'seed'; the caller's generator
# and its state are restored afterwards

with_rng <- function(kind, seed, expr) {
  old_kind <- RNGkind()
  has <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (has) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (has) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind(kind[1], kind[2], kind[3])
  set.seed(seed)
  expr
}

# Labels contain non-ASCII characters (en dash, prime, superscript 2, Greek
# letters), and apim_meta_sensitivity.R uses them as list names. In a non-UTF-8
# session (e.g. R started with the C locale) such names cannot be represented,
# so the character type is switched to UTF-8 when possible.

if (!isTRUE(l10n_info()[["UTF-8"]])) {
  for (.loc in c("C.UTF-8", "C.utf8", "en_US.UTF-8", "English_United States.utf8"))
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", .loc)))) break
  if (!isTRUE(l10n_info()[["UTF-8"]]))
    message("Note: R is not running in a UTF-8 locale; labels with special characters may print as <U+....>.")
}

## 1. Design

## Files. The scripts run from the 'Missing data' folder (APIM/6 Meta-analysis/
## Missing data) and use the original scripts and data in place: each input is
## looked up in this folder first and then in the sibling folders Meta-analysis/
## (apim_meta_analysis.R, toy_study*.csv) and Sensitivity/ (apim_meta_sensitivity.R
## and its outputs). Nothing needs to be copied, and the pipeline always runs the
## current version of the meta-analysis code.

INPUT_DIRS  <- c(".", "../meta-analysis", "../sensitivity", "meta-analysis", "sensitivity")

SENS_DIRS   <- c("../sensitivity", ".", "sensitivity")   # outputs of apim_meta_sensitivity.R

find_input <- function(file, dirs = INPUT_DIRS) {
  if (file.exists(file) || dirname(file) != ".") return(file)
  hit <- file.path(dirs, file)
  hit <- hit[file.exists(hit)]
  if (length(hit)) hit[1] else file        # not found: keep the name for the error message
}

not_found <- function(file, dirs = INPUT_DIRS)
  paste0("'", basename(file), "' not found in ", paste(dirs, collapse = ", "),
         " (relative to the working directory, which must be the 'Missing data' folder)")

APIM_SCRIPT <- find_input("apim_meta_analysis.R")

STUDIES <- data.frame(key = c("study1", "study2"), study = c("Study 1", "Study 2"),
                      n = c(130L, 174L), toy_file = c("toy_study1.csv", "toy_study2.csv"),
                      stringsAsFactors = FALSE)

ID_VAR      <- "dyad_id"

N_ITEMS     <- 4L

ITEM_VALUES <- 1:5            # Likert response range (imputations are kept inside it)

# The APIMeM is estimated with the functions of apim_meta_analysis.R (sections
# 0-4: model dictionary, scale scores and reliability, lavaan fit, extraction of
# the 12 standardized paths and their 12 x 12 sampling covariance matrix). They
# are evaluated in their own environment, exactly as check_recovery.R does, so
# that this pipeline always uses the same APIMeM code as the meta-analysis.

load_apim_functions <- function(file = APIM_SCRIPT) {
  if (!file.exists(file)) stop(not_found(file), "; or set APIM_SCRIPT to its path.")
  src <- readLines(file, warn = FALSE, encoding = "UTF-8")
  stop_at <- grep("^## 5\\. RUN", src)[1] - 2
  if (is.na(stop_at) || stop_at < 50) stop("Section '## 5. RUN' not found in ", file)
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  had_ovr <- exists("CFG_OVERRIDE", envir = globalenv(), inherits = FALSE)
  if (had_ovr) { ovr <- get("CFG_OVERRIDE", envir = globalenv()); rm("CFG_OVERRIDE", envir = globalenv()) }
  env <- new.env(parent = globalenv())
  suppressPackageStartupMessages(eval(parse(text = src[seq_len(stop_at)], encoding = "UTF-8"),
                                      envir = env))
  if (had_ovr) assign("CFG_OVERRIDE", ovr, envir = globalenv())
  if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())        # the script sets a seed
  else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  needed <- c("CFG", "VARS", "PATHS", "EFFECTS", "PATH_GROUP_LABEL", "apim_syntax", "indirect_defs",
              "reliability_coef", "prep_study", "fit_apim", "extract_paths", "analyse_study")
  miss <- needed[!vapply(needed, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss)) stop(file, " does not define: ", paste(miss, collapse = ", "))
  env
}

APIM <- load_apim_functions()

APIM_CFG <- APIM$CFG                    # constructs, members, item_regex, min_prop_items (.75), omega

APIM_CFG$reliability_supplied <- NULL

VARS    <- APIM$VARS                    # X1 X2 M1 M2 Y1 Y2

PATHS   <- APIM$PATHS                   # 12 paths: label, lhs, rhs, path, type, effect

EFFECTS <- APIM$EFFECTS                 # actor_a partner_a actor_b partner_b actor_c partner_c

CONSTRUCTS <- APIM_CFG$constructs        # X = "x", M = "m", Y = "y"

MEMBERS    <- APIM_CFG$members           # "p1", "p2"

SCALES <- setNames(sprintf("%s_%s", CONSTRUCTS[substr(VARS, 1, 1)],
                           MEMBERS[as.integer(substr(VARS, 2, 2))]), VARS)   # X1 = "x_p1", ...

ITEMS  <- unlist(lapply(SCALES, function(s) sprintf("%s_%d", s, seq_len(N_ITEMS))), use.names = FALSE)

ITEM_VAR       <- setNames(rep(VARS, each = N_ITEMS), ITEMS)           # x_p1_3 -> "X1"

ITEM_SCALE     <- setNames(rep(unname(SCALES), each = N_ITEMS), ITEMS)  # x_p1_3 -> "x_p1"

ITEM_MEMBER    <- setNames(rep(MEMBERS[as.integer(substr(VARS, 2, 2))], each = N_ITEMS), ITEMS)

ITEM_CONSTRUCT <- setNames(rep(tolower(substr(VARS, 1, 1)), each = N_ITEMS), ITEMS)

VAR_LABELS <- setNames(sprintf("%s, %s", APIM_CFG$construct_labels[substr(VARS, 1, 1)],
                               APIM_CFG$member_labels[as.integer(substr(VARS, 2, 2))]), VARS)

EFFECT_LAB <- c(actor_a = "Actor a", partner_a = "Partner a", actor_b = "Actor b",
                partner_b = "Partner b", actor_c = "Actor c\u{2032}", partner_c = "Partner c\u{2032}")

COMPONENTS <- c("uncorrected", "corrected")

COMP_LAB   <- c(uncorrected = "Uncorrected", corrected = "Corrected for unreliability")

EFF_COLS   <- lapply(setNames(EFFECTS, EFFECTS), function(e) PATHS$label[PATHS$effect == e])

INDIRECT_DEFS <- APIM$indirect_defs(TRUE)  # products of the partner-pooled paths

SPECIFIC_INDIRECT <- setNames(list(c("actor_a", "actor_b"), c("partner_a", "partner_b"),
                                   c("actor_a", "partner_b"), c("partner_a", "actor_b")),
                              c("Actor\u{2013}actor", "Partner\u{2013}partner",
                                "Actor\u{2013}partner", "Partner\u{2013}actor"))   # setNames: locale-safe names

zs <- function(x) { s <- stats::sd(x); if (!is.finite(s) || s == 0) x * 0 else (x - mean(x)) / s }
zs_na <- function(x) { z <- (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE); z[is.na(z)] <- 0; z }

# Evaluate expr with R's default generator seeded with 'seed' (whatever the
# current generator), restoring the caller's generator and state afterwards

with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  with_rng(RNG_DEFAULT, seed, expr)
}

# Scale scores (item means). prorate = TRUE reproduces prep_study() of
# apim_meta_analysis.R: mean of the available items if at least min_prop of
# them are observed, NA otherwise. prorate = FALSE: NA if any item is missing.
scale_scores <- function(dat, prorate = FALSE, min_prop = APIM_CFG$min_prop_items) {
  out <- sapply(VARS, function(v) {
    it <- as.matrix(dat[, ITEMS[ITEM_VAR == v], drop = FALSE])
    if (!prorate) return(rowMeans(it))
    sc <- rowMeans(it, na.rm = TRUE); sc[rowMeans(!is.na(it)) < min_prop] <- NA; sc
  })
  if (is.null(dim(out))) out <- matrix(out, nrow = 1, dimnames = list(NULL, VARS))
  out
}

# Standardized APIMeM paths in closed form from a 6 x 6 covariance matrix (the
# APIMeM is saturated, so its ML estimates are the regressions implied by the
# saturated covariance matrix). Used for population values and for checks;
# the analyses themselves use lavaan.

std_paths_from_cov <- function(S) {
  out <- setNames(numeric(nrow(PATHS)), PATHS$label)
  sds <- sqrt(diag(S))
  for (y in unique(PATHS$lhs)) {
    ix <- which(PATHS$lhs == y); x <- PATHS$rhs[ix]
    out[ix] <- solve(S[x, x, drop = FALSE], S[x, y]) * sds[x] / sds[y]
  }
  out
}
# Partner-pooled (simple average) value of each effect type
effect_means <- function(paths) vapply(EFFECTS, function(e) mean(paths[EFF_COLS[[e]]]), 0)

## 2. Toy population and complete toy data
# Exactly the generator of simulate_toy_data.R (latent APIMeM with known
# standardized paths, 4 items per construct, 5-point Likert discretisation;
# Study 2 has lower loadings). The attached toy_study*.csv contain ~1.5% MCAR
# item-level missingness (p_miss = .015 in simulate_toy_data.R). The missing
# cells are restored by re-running the generator with the same seed and
# p_miss = 0: the random-number stream is unchanged (runif() is still drawn),
# so every observed value is reproduced exactly (checked in 01).

TOY <- list(seed = 20260911, r_x = .30, r_eM = .20, r_eY = .25,
            loadings = list(study1 = c(.80, .75, .78, .72), study2 = c(.62, .58, .66, .55)),
            thresholds = c(-1.6, -0.7, 0.3, 1.2), n = c(study1 = 130L, study2 = 174L))

# Platform independence. simulate_toy_data.R draws the latent scores with
# MASS::mvrnorm(), which uses the eigenvectors of Sigma_T. LAPACK libraries can
# return eigenvectors with opposite signs (e.g., R 4.6 on macOS vs. R 4.3 on
# Linux), which changes every draw, so the toy data could not be regenerated
# everywhere. mvrnorm_ref() is MASS::mvrnorm() with each eigenvector oriented
# as on the platform that generated the toy data (TOY_EIGEN_REF; only the signs
# are used), so the draws are the same on every platform (and bit-identical to
# MASS::mvrnorm() there).

TOY_EIGEN_REF <- matrix(c(
  -0.4029, -0.3701, -0.4228, -0.4228, -0.4048, -0.4234, -0.3564, 0.4752, -0.4552, 0.4772, -0.3780, 0.2633,
   0.5110,  0.5088,  0.0929, -0.0843, -0.4292, -0.5292,  0.0609, -0.5083, 0.4205, 0.4830, -0.5717, 0.0307,
  -0.6247,  0.2400,  0.4827,  0.1907,  0.2100, -0.4885,  0.2355, -0.2498, -0.4423, 0.5628, 0.3715, -0.4814),
  6, 6)

mvrnorm_ref <- function(n, mu, Sigma, ref = TOY_EIGEN_REF, tol = 1e-6) {
  p <- length(mu)
  eS <- eigen(Sigma, symmetric = TRUE); ev <- eS$values; V <- eS$vectors
  if (!all(ev >= -tol * abs(ev[1L]))) stop("'Sigma' is not positive definite")
  flip <- colSums(V * ref) < 0
  V[, flip] <- -V[, flip]
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + V %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  t(X)
}

# Fallback for 01 if the generator still does not reproduce the attached files:
# the values of the 118 cells that simulate_toy_data.R deleted (column-major
# order of the missing cells of toy_study1.csv and toy_study2.csv) and checksums
# that identify the attached files and the complete data.

TOY_MISSING_VALUES <- c(
  study1 = "4243334354434242523443442313344232444242522534334",
  study2 = "423243344334244243433341532344343243432115444423533422322114341433332")

TOY_CHECKSUM <- list(observed = c(study1 = 49105441, study2 = 65303981),
                     complete = c(study1 = 49890255, study2 = 66321666))
toy_checksum <- function(x) {
  v <- as.matrix(x); w <- (seq_along(v) * 7919) %% 9973 + 1; ok <- !is.na(v)
  sum(v[ok] * w[ok])
}

toy_B <- function() {
  B <- matrix(0, 6, 6, dimnames = list(VARS, VARS))
  B["M1", "X1"] <- .40; B["M1", "X2"] <- .10    # actor a (P1), partner a (X2->M1)
  B["M2", "X2"] <- .35; B["M2", "X1"] <- .15    # actor a (P2), partner a (X1->M2)
  B["Y1", "M1"] <- .30; B["Y1", "M2"] <- .10    # actor b (P1), partner b (M2->Y1)
  B["Y2", "M2"] <- .35; B["Y2", "M1"] <- .12    # actor b (P2), partner b (M1->Y2)
  B["Y1", "X1"] <- .15; B["Y1", "X2"] <- .05    # actor c' (P1), partner c' (X2->Y1)
  B["Y2", "X2"] <- .10; B["Y2", "X1"] <- .08    # actor c' (P2), partner c' (X1->Y2)
  B
}

toy_implied_sigma <- function(B = toy_B(), r_x = TOY$r_x, r_eM = TOY$r_eM, r_eY = TOY$r_eY) {
  Psi <- diag(6); dimnames(Psi) <- dimnames(B)
  Psi["X1", "X2"] <- Psi["X2", "X1"] <- r_x
  IB  <- solve(diag(6) - B)
  for (it in 1:200) {
    S <- IB %*% Psi %*% t(IB)
    for (v in c("M1", "M2", "Y1", "Y2")) Psi[v, v] <- Psi[v, v] + (1 - S[v, v])
    Psi["M1", "M2"] <- Psi["M2", "M1"] <- r_eM * sqrt(Psi["M1", "M1"] * Psi["M2", "M2"])
    Psi["Y1", "Y2"] <- Psi["Y2", "Y1"] <- r_eY * sqrt(Psi["Y1", "Y1"] * Psi["Y2", "Y2"])
  }
  S <- IB %*% Psi %*% t(IB)
  stopifnot(max(abs(diag(S) - 1)) < 1e-10)
  S
}

toy_make_study <- function(n_dyads, loadings, thresholds, Sigma_T, p_miss = 0) {
  Tscores <- mvrnorm_ref(n_dyads, mu = rep(0, 6), Sigma = Sigma_T)
  colnames(Tscores) <- VARS
  out <- data.frame(dyad_id = sprintf("D%03d", seq_len(n_dyads)))
  for (v in VARS) {
    construct <- tolower(substr(v, 1, 1)); member <- paste0("p", substr(v, 2, 2))
    for (i in seq_along(loadings)) {
      lam   <- loadings[i]
      ystar <- lam * Tscores[, v] + rnorm(n_dyads, 0, sqrt(1 - lam^2))
      item  <- as.integer(cut(ystar, breaks = c(-Inf, thresholds, Inf)))
      item[runif(n_dyads) < p_miss] <- NA          # drawn even if p_miss = 0 (same stream)
      out[[sprintf("%s_%s_%d", construct, member, i)]] <- item
    }
  }
  out
}

regenerate_toy_data <- function(p_miss = 0) {
  set.seed(TOY$seed, kind = RNG_DEFAULT[1], normal.kind = RNG_DEFAULT[2], sample.kind = RNG_DEFAULT[3])
  Sigma_T <- toy_implied_sigma()
  s1 <- toy_make_study(TOY$n[["study1"]], TOY$loadings$study1, TOY$thresholds, Sigma_T, p_miss)
  s2 <- toy_make_study(TOY$n[["study2"]], TOY$loadings$study2, TOY$thresholds, Sigma_T, p_miss)
  list(study1 = s1, study2 = s2)
}

# Population values (cached): the latent paths (the estimand of the corrected
# component; equal in both studies) and, per study, the paths among the scale
# scores (the estimand of the uncorrected component; attenuated, and more so in
# Study 2) and the population omega of each scale, from N_big simulated dyads.

toy_population <- function(file = "data/population_values.rds", N_big = 2e5, seed = 99) {
  if (!is.null(file) && file.exists(file)) return(readRDS(file))
  Sigma_T <- toy_implied_sigma()
  latent <- std_paths_from_cov(Sigma_T)
  per <- lapply(names(TOY$loadings), function(k) with_seed(seed + match(k, names(TOY$loadings)), {
    d <- toy_make_study(N_big, TOY$loadings[[k]], TOY$thresholds, Sigma_T)
    S <- stats::cov(scale_scores(d))
    om <- vapply(VARS, function(v) {
      C <- stats::cov(d[, ITEMS[ITEM_VAR == v]])
      f <- stats::factanal(covmat = C, factors = 1, n.obs = N_big)
      lam <- f$loadings[, 1] * sqrt(diag(C)); th <- f$uniquenesses * diag(C)
      sum(lam)^2 / (sum(lam)^2 + sum(th))
    }, 0)
    list(uncorrected = std_paths_from_cov(S), omega = om)
  }))
  names(per) <- names(TOY$loadings)
  out <- list(latent = latent,
              uncorrected = lapply(per, `[[`, "uncorrected"),
              omega = lapply(per, `[[`, "omega"), N_big = N_big, seed = seed)
  # Partner-pooled truth: simple average of the two partner-specific paths. For
  # the uncorrected component the two studies differ by design (different
  # reliabilities); the average of the two studies is used as a reference.
  out$effect_truth <- list(
    corrected = effect_means(latent),
    uncorrected = rowMeans(sapply(out$uncorrected, effect_means)))
  if (!is.null(file)) { dir.create(dirname(file), FALSE, TRUE); saveRDS(out, file) }
  out
}

## 3. Auxiliary variables and instruments
# aux_p1, aux_p2  partner-level AUXILIARY variables (continuous): always observed,
#                 related to the partner's own X, M and Y and (weakly) to the
#                 other partner's Y. Stand-ins for, e.g., a screening measure,
#                 relationship duration or a related questionnaire.
# aux_dyad        dyad-level AUXILIARY variable (binary), related to both
#                 partners' outcomes (stand-in for, e.g., cohabitation).
# iv_p1, iv_p2    INSTRUMENTS for item nonresponse (one per partner): independent
#                 of everything else; they shift the probability of nonresponse
#                 of that partner in the MNAR datasets only (stand-ins for, e.g.,
#                 randomised questionnaire order or reminder timing). A selection
#                 model (Heckman-MI) is identified without leaning on
#                 distributional assumptions only through such a variable
#                 (exclusion restriction; Heckman, 1979; Tchetgen Tchetgen &
#                 Wirth, 2017).
# The auxiliaries drive the MAR mechanism, so MAR is ignorable only if they are
# in the imputation / FIML model (inclusive strategy; Collins, Schafer & Kam,
# 2001). The instruments are used ONLY by Heckman-MI: among respondents an
# instrument becomes spuriously related to the values (selection on a collider),
# so a MAR model should not condition on it.

AUX_VARS <- c("aux_p1", "aux_p2", "aux_dyad")

IV_VARS  <- c("iv_p1", "iv_p2")

IV_VAR   <- IV_VARS                    # default 'iv' argument of mice.impute.heckoprob

AUX_COEF <- list(member = c(x = 0.15, m = 0.25, y = 0.40, partner_y = 0.10),
                 dyad   = c(intercept = 0, y = 0.80))

add_aux_iv <- function(dat, b = AUX_COEF) {
  n <- nrow(dat); S <- scale_scores(dat)
  stopifnot(!anyNA(S))
  for (k in 1:2) {
    o  <- 3 - k
    lp <- b$member[["x"]] * zs(S[, paste0("X", k)]) + b$member[["m"]] * zs(S[, paste0("M", k)]) +
      b$member[["y"]] * zs(S[, paste0("Y", k)]) + b$member[["partner_y"]] * zs(S[, paste0("Y", o)])
    dat[[paste0("aux_", MEMBERS[k])]] <- round(lp + rnorm(n, 0, sqrt(max(1 - stats::var(lp), 0.2))), 3)
  }
  ybar <- zs(rowMeans(cbind(zs(S[, "Y1"]), zs(S[, "Y2"]))))
  dat$aux_dyad <- stats::rbinom(n, 1, stats::plogis(b$dyad[["intercept"]] + b$dyad[["y"]] * ybar))
  for (j in MEMBERS) dat[[paste0("iv_", j)]] <- round(rnorm(n), 3)
  dat[, c(ID_VAR, ITEMS, AUX_VARS, IV_VARS)]
}

## 4. Item-level amputation
# Logistic selection model per item j of partner k:  P(M_ij = 1) = logistic(a_j + eta_ij)
#   MCAR: eta = 0
#   MAR : eta = mar_direction * strength * z_aux_ik,
#         z_aux_ik = z(z(aux_pk) + z(aux_dyad))   (partner k's auxiliary composite)
#   MNAR: eta = MNAR_DIRECTION[construct] * strength * z(y_ij) + iv_strength * z(iv_pk)
#         (self-masking: the probability depends on the item value itself)
# Exact rates: M_ij = 1 iff T_ij = logit(U_ij) - eta_ij < a_j with U_ij ~ U(0, 1),
# a_j = the n_mis-th order statistic of T_.j, so EXACTLY round(rate * N) values
# are deleted from every item. One U matrix per study is shared by all nine
# datasets: within a mechanism the missing cells are nested (5% within 10%
# within 15%) and the mechanisms differ only by their selection model.
#
# MNAR direction per construct (+1: higher item values more often missing).
# Plausible but ASSUMED (e.g. X a stressor, M strain, Y relationship quality):

MNAR_DIRECTION   <- c(x = 1, m = 1, y = -1)

AMPUTE_DEFAULTS  <- list(strength = 1, mar_direction = -1, iv_strength = 1)

MECHANISMS       <- c("MCAR", "MAR", "MNAR")

RATES            <- c(0.05, 0.10, 0.15)

aux_composite <- function(dat, member) zs(zs(dat[[paste0("aux_", member)]]) + zs(dat$aux_dyad))

amputation_driver <- function(dat, mechanism, item, strength, mar_direction, iv_strength,
                              mnar_direction = MNAR_DIRECTION) {
  mem <- ITEM_MEMBER[[item]]
  switch(mechanism,
         MCAR = rep(0, nrow(dat)),
         MAR  = mar_direction * strength * aux_composite(dat, mem),
         MNAR = mnar_direction[[ITEM_CONSTRUCT[[item]]]] * strength * zs_na(dat[[item]]) +
           iv_strength * zs(dat[[paste0("iv_", mem)]]),
         stop("Unknown mechanism: ", mechanism))
}

ampute_items <- function(dat, mechanism = c("MCAR", "MAR", "MNAR"), prop = 0.10, U = NULL,
                         items = ITEMS, strength = AMPUTE_DEFAULTS$strength,
                         mar_direction = AMPUTE_DEFAULTS$mar_direction,
                         iv_strength = AMPUTE_DEFAULTS$iv_strength) {
  mechanism <- match.arg(mechanism)
  n <- nrow(dat); J <- length(items)
  if (is.null(U)) U <- matrix(runif(n * J), n, J, dimnames = list(NULL, items))
  n_mis <- round(prop * n)
  out <- dat
  for (j in seq_len(J)) {
    v <- items[j]
    eta <- amputation_driver(dat, mechanism, v, strength, mar_direction, iv_strength)
    Tj <- stats::qlogis(U[, j]) - eta
    out[[v]][rank(Tj, ties.method = "first") <= n_mis] <- NA
  }
  attr(out, "amputation") <- list(mechanism = mechanism, prop = prop, n_mis = n_mis,
                                  strength = strength, mar_direction = mar_direction,
                                  iv_strength = if (mechanism == "MNAR") iv_strength else 0)
  out
}

# Missingness at the item, scale-score and dyad level

missing_summary <- function(dat) {
  M <- is.na(dat[, ITEMS, drop = FALSE]); r <- rowSums(M)
  sc_any <- sapply(VARS, function(v) rowSums(M[, ITEMS[ITEM_VAR == v], drop = FALSE]) > 0)
  sc_pro <- is.na(scale_scores(dat, prorate = TRUE))
  rbind(
    data.frame(variable = ITEMS, n_missing = colSums(M), pct_missing = round(100 * colMeans(M), 2)),
    data.frame(variable = c("ALL item cells", "Dyads with any missing item", "Dyads complete (all 24 items)",
                            "Scale scores with any missing item",
                            "Scale scores not computable by proration (< 75% items)",
                            "Dyads complete on the 6 prorated scale scores"),
               n_missing = c(sum(M), sum(r > 0), sum(r == 0), sum(sc_any), sum(sc_pro),
                             sum(rowSums(sc_pro) == 0)),
               pct_missing = round(100 * c(mean(M), mean(r > 0), mean(r == 0), mean(sc_any),
                                           mean(sc_pro), mean(rowSums(sc_pro) == 0)), 2)),
    make.row.names = FALSE)
}

# Standardised mean differences (missing minus observed rows) in the deleted
# value, the partner's auxiliary composite, the instrument and the scale scores

mechanism_check <- function(complete, amputed) {
  smd <- function(x, m) if (any(m) && any(!m)) (mean(x[m]) - mean(x[!m])) / stats::sd(x) else NA
  S <- scale_scores(complete)
  do.call(rbind, lapply(ITEMS, function(v) {
    m <- is.na(amputed[[v]]); mem <- ITEM_MEMBER[[v]]; k <- match(mem, MEMBERS)
    data.frame(item = v, pct_missing = round(100 * mean(m), 1),
               mnar_direction = MNAR_DIRECTION[[ITEM_CONSTRUCT[[v]]]],
               smd_own_value = smd(complete[[v]], m),
               smd_aux_composite = smd(aux_composite(complete, mem), m),
               smd_iv = smd(complete[[paste0("iv_", mem)]], m),
               smd_own_scale_score = smd(S[, ITEM_VAR[[v]]], m),
               smd_own_Y = smd(S[, paste0("Y", k)], m))
  }))
}

# Realised selection model: logistic regression of the item missingness
# indicators (stacked over items, item fixed effects) on the direction-signed
# standardised own value, the partner's auxiliary composite and instrument.
# Defining property of MAR: given the auxiliaries, the value itself does not
# predict missingness.

mechanism_regression <- function(complete, amputed) {
  st <- do.call(rbind, lapply(ITEMS, function(v) {
    mem <- ITEM_MEMBER[[v]]
    data.frame(item = v, M = as.integer(is.na(amputed[[v]])),
               z_own = MNAR_DIRECTION[[ITEM_CONSTRUCT[[v]]]] * zs(complete[[v]]),
               z_aux = aux_composite(complete, mem), z_iv = zs(complete[[paste0("iv_", mem)]]))
  }))
  fit <- suppressWarnings(stats::glm(M ~ 0 + factor(item) + z_own + z_aux + z_iv,
                                     family = stats::binomial, data = st))
  cf <- summary(fit)$coefficients[c("z_own", "z_aux", "z_iv"), 1:2]
  data.frame(term = rownames(cf), log_odds = cf[, 1], se = cf[, 2], row.names = NULL)
}

## 5. Custom mice methods
# catx     binary/ordinal variables: logistic regression / proportional-odds model
#          with parameters drawn from their approximate posterior (proper MI) and
#          mice's data augmentation against perfect prediction (White, Daniel &
#          Royston, 2010). Here: the 5-point items.
# pmmx     predictive mean matching with a delta shift of the imputed values (SD
#          units of the observed values), stochastic rounding and bounds: delta-
#          adjusted MI for the MNAR sensitivity analysis. delta = 0 is mice-PMM.
# heckoprob  Heckman-type selection-model MI for ordinal variables: ordered-probit
#          outcome equation with probit selection, joint ML (BHHH), weakly
#          informative priors; the instrument enters the selection equation only
#          (Galimard, Chevret, Curis & Resche-Rigon, 2018, for the binary case).

.drop_excluded <- function(y, ry, x, wy, exclude) {
  if (is.null(wy)) wy <- !ry
  x <- as.matrix(x)
  if (length(exclude)) {
    keep <- !(ry & y %in% exclude)
    y <- y[keep]; ry <- ry[keep]; x <- x[keep, , drop = FALSE]; wy <- wy[keep]
  }
  list(y = y, ry = ry, x = x, wy = wy)
}

.std <- function(x, rows) {
  x <- as.matrix(x)
  if (!ncol(x)) return(x)
  m <- colMeans(x[rows, , drop = FALSE]); s <- apply(x[rows, , drop = FALSE], 2, stats::sd)
  keep <- is.finite(s) & s > 1e-10
  x <- sweep(sweep(x[, keep, drop = FALSE], 2, m[keep]), 2, s[keep], "/")
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  colnames(x) <- make.names(colnames(x), unique = TRUE)
  x
}

.mvn_draw <- function(mu, V) {
  V <- (V + t(V)) / 2
  if (any(!is.finite(V))) return(mu)
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (!is.null(R)) return(mu + drop(crossprod(R, rnorm(length(mu)))))
  e <- eigen(V, symmetric = TRUE)
  mu + drop(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(length(mu))))
}

.inv <- function(H) tryCatch(solve(H), error = function(e) MASS::ginv(H))
stoch_round <- function(v) { f <- floor(v); f + (runif(length(v)) < (v - f)) }

.CAT_CACHE <- new.env()

mice.impute.catx <- function(y, ry, x, wy = NULL, exclude = NULL, delta = 0, ...) {
  s <- .drop_excluded(y, ry, x, wy, exclude); y <- s$y; ry <- s$ry; wy <- s$wy
  lev <- sort(unique(y[ry])); K <- length(lev)
  if (K == 1) return(rep(lev, sum(wy)))
  x <- .std(s$x, ry)
  yf <- factor(match(y, lev), levels = seq_len(K))
  aug <- tryCatch(mice:::augment(yf, ry, x, wy),
                  error = function(e) list(y = yf, ry = ry, x = x, wy = wy, w = rep(1, length(yf))))
  xa <- as.matrix(aug$x); xm <- x[wy, , drop = FALSE]
  if (K == 2) {
    X1 <- cbind(1, xa[aug$ry, , drop = FALSE])
    fit <- suppressWarnings(stats::glm.fit(X1, as.integer(aug$y[aug$ry]) - 1L, weights = aug$w[aug$ry],
                                           family = stats::binomial()))
    b <- fit$coefficients; ok <- !is.na(b)
    V <- stats::summary.glm(fit)$cov.unscaled
    bs <- rep(0, length(b)); bs[ok] <- .mvn_draw(b[ok], V)
    p <- stats::plogis(drop(cbind(1, xm) %*% bs) + delta)
    return(lev[1L + (runif(length(p)) < p)])
  }
  df <- data.frame(y = aug$y, xa)[aug$ry, , drop = FALSE]
  key <- paste(K, sum(ry), sum(match(y[ry], lev) * seq_len(sum(ry))), paste(colnames(x), collapse = ","))
  st <- .CAT_CACHE[[key]]
  fit <- tryCatch(suppressWarnings(MASS::polr(y ~ ., data = df, weights = aug$w[aug$ry],
                                              Hess = TRUE, method = "logistic", start = st)),
                  error = function(e) NULL)
  if (is.null(fit) && !is.null(st))
    fit <- tryCatch(suppressWarnings(MASS::polr(y ~ ., data = df, weights = aug$w[aug$ry],
                                                Hess = TRUE, method = "logistic")),
                    error = function(e) NULL)
  if (!is.null(fit)) assign(key, c(stats::coef(fit), fit$zeta), envir = .CAT_CACHE)
  if (is.null(fit)) {                              # rare: mice's polyreg (improper) as fallback
    v <- mice::mice.impute.polyreg(yf, ry, x, wy)
    return(lev[as.integer(as.character(v))])
  }
  nb <- length(stats::coef(fit))
  th <- .mvn_draw(c(stats::coef(fit), fit$zeta), .inv(fit$Hessian))
  beta <- th[seq_len(nb)]; zeta <- sort(th[-seq_len(nb)])
  eta <- drop(xm[, names(stats::coef(fit)), drop = FALSE] %*% beta) + delta
  cp <- stats::plogis(outer(-eta, zeta, `+`))
  lev[1L + rowSums(runif(length(eta)) > cp)]
}

mice.impute.pmmx <- function(y, ry, x, wy = NULL, exclude = NULL, delta = 0, bounds = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  v <- mice::mice.impute.pmm(y, ry, x, wy, exclude = if (length(exclude)) exclude else -99999999, ...)
  if (delta != 0) {
    v <- stoch_round(v + delta * stats::sd(y[ry & !(y %in% exclude)]))
    if (!is.null(bounds)) v <- pmin(bounds[2], pmax(bounds[1], v))
  }
  v
}

Phi2 <- function(a, b, r) {
  out <- (a == Inf) * stats::pnorm(b)
  ok <- is.finite(a)
  if (any(ok)) out[ok] <- pbivnorm::pbivnorm(a[ok], b[ok], r)
  out
}

RHO_MAX <- 0.95; RHO_PRIOR_SD <- 2; COEF_PRIOR_SD <- 2.5

probit_ridge <- function(W, s, prec = c(0, rep(1 / COEF_PRIOR_SD^2, ncol(W) - 1)), maxit = 50) {
  b <- rep(0, ncol(W))
  for (it in seq_len(maxit)) {
    eta <- drop(W %*% b)
    p <- pmin(pmax(stats::pnorm(eta), 1e-10), 1 - 1e-10); dn <- pmax(stats::dnorm(eta), 1e-10)
    w <- dn^2 / (p * (1 - p)); z <- eta + (s - p) / dn
    bn <- drop(solve(crossprod(W, W * w) + diag(prec, ncol(W)), crossprod(W, w * z)))
    done <- max(abs(bn - b)) < 1e-8; b <- bn
    if (done) break
  }
  b
}

heckoprob_fit <- function(yk, s, X, W, K, start = NULL, maxit = 30, tol = 1e-9) {
  pX <- ncol(X); pW <- ncol(W)
  ib <- seq_len(pX); ic <- pX + seq_len(K - 1); ig <- pX + K - 1 + seq_len(pW)
  ir <- pX + K - 1 + pW + 1
  o <- s == 1; Xo <- X[o, , drop = FALSE]; Wo <- W[o, , drop = FALSE]
  Wm <- W[!o, , drop = FALSE]; ko <- yk[o]
  Ko <- outer(ko, seq_len(K), `==`) * 1
  prec <- c(rep(1 / COEF_PRIOR_SD^2, pX), rep(0, K - 1), 0,
            rep(1 / COEF_PRIOR_SD^2, pW - 1), 1 / RHO_PRIOR_SD^2)
  parts <- function(th) {
    b <- th[ib]; cc <- th[ic]; g <- th[ig]; rho <- RHO_MAX * unname(tanh(th[ir]))
    tau <- cumsum(c(cc[1], exp(cc[-1])))
    xb <- drop(Xo %*% b); wg <- drop(Wo %*% g); wm <- drop(Wm %*% g)
    u <- c(tau, Inf)[ko] - xb; l <- c(-Inf, tau)[ko] - xb
    Lo <- pmax(Phi2(u, wg, -rho) - Phi2(l, wg, -rho), 1e-300)
    Lm <- pmax(stats::pnorm(-wm), 1e-300)
    list(th = th, cc = cc, rho = rho, wg = wg, wm = wm, u = u, l = l, Lo = Lo, Lm = Lm,
         nll = -(sum(log(Lo)) + sum(log(Lm))) + 0.5 * sum(prec * th^2))
  }
  nll <- function(th) parts(th)$nll
  scores <- function(th, q = parts(th)) {
    rho <- q$rho; wg <- q$wg; sr <- sqrt(1 - rho^2)
    fu <- is.finite(q$u); fl <- is.finite(q$l)
    uu <- ifelse(fu, q$u, 0); ll <- ifelse(fl, q$l, 0)
    fa_u <- fu * stats::dnorm(uu) * stats::pnorm((wg + rho * uu) / sr)
    fa_l <- fl * stats::dnorm(ll) * stats::pnorm((wg + rho * ll) / sr)
    fb_u <- ifelse(fu, stats::dnorm(wg) * stats::pnorm((uu + rho * wg) / sr), stats::dnorm(wg))
    fb_l <- fl * stats::dnorm(wg) * stats::pnorm((ll + rho * wg) / sr)
    k2 <- 1 / (2 * pi * sr)
    f2_u <- fu * exp(-(uu^2 + 2 * rho * uu * wg + wg^2) / (2 * sr^2)) * k2
    f2_l <- fl * exp(-(ll^2 + 2 * rho * ll * wg + wg^2) / (2 * sr^2)) * k2
    qu <- fa_u / q$Lo; ql <- fa_l / q$Lo
    Sb <- -Xo * (qu - ql)
    St <- Ko[, -K, drop = FALSE] * qu - Ko[, -1, drop = FALSE] * ql
    Sc <- St
    if (K > 2) for (m in 2:(K - 1)) Sc[, m] <- exp(q$cc[m]) * rowSums(St[, m:(K - 1), drop = FALSE])
    Sc[, 1] <- rowSums(St)
    Sg <- Wo * ((fb_u - fb_l) / q$Lo)
    Sr <- -(f2_u - f2_l) / q$Lo * RHO_MAX * (1 - tanh(q$th[ir])^2)
    So <- cbind(Sb, Sc, Sg, Sr)
    Sm <- matrix(0, nrow(Wm), ncol(So))
    Sm[, ig] <- -Wm * exp(stats::dnorm(q$wm, log = TRUE) - stats::pnorm(-q$wm, log.p = TRUE))
    rbind(So, Sm)
  }
  if (is.null(start) || length(start) != ir) {
    g0 <- probit_ridge(W, s)
    pk <- pmax(tabulate(ko, K) / length(ko), 1e-3)
    t0 <- stats::qnorm(pmin(cumsum(pk / sum(pk))[-K], 1 - 1e-3))
    start <- c(rep(0, pX), t0[1], log(pmax(diff(t0), 1e-2)), g0, 0)
  }
  q <- parts(start); conv <- FALSE; it <- 0
  while (it < maxit) {
    it <- it + 1
    S <- scores(q$th, q); B <- crossprod(S) + diag(prec + 1e-8)
    g <- colSums(S) - prec * q$th
    step <- tryCatch(solve(B, g), error = function(e) drop(MASS::ginv(B) %*% g))
    lam <- 1
    repeat {
      q1 <- parts(q$th + lam * step)
      if (is.finite(q1$nll) && q1$nll <= q$nll + 1e-10) break
      lam <- lam / 2; if (lam < 1e-8) break
    }
    if (lam < 1e-8) { conv <- TRUE; break }
    done <- abs(q$nll - q1$nll) < tol * (abs(q$nll) + tol) || max(abs(lam * step)) < 1e-6
    q <- q1
    if (done) { conv <- TRUE; break }
  }
  if (!conv) {
    opt <- stats::optim(q$th, nll, function(t) -colSums(scores(t)) + prec * t, method = "BFGS",
                        control = list(reltol = 1e-10, maxit = 500))
    if (opt$value <= q$nll) q <- parts(opt$par)
    conv <- opt$convergence == 0
  }
  S <- scores(q$th, q)
  list(par = q$th, V = .inv(crossprod(S) + diag(prec)), loglik = -q$nll, converged = conv,
       iter = it, idx = list(b = ib, c = ic, g = ig, r = ir))
}

heckoprob_probs <- function(par, idx, X, W) {
  b <- par[idx$b]; cc <- par[idx$c]; g <- par[idx$g]
  rho <- RHO_MAX * tanh(par[idx$r])
  tau <- cumsum(c(cc[1], exp(cc[-1])))
  xb <- drop(X %*% b); wg <- drop(W %*% g)
  C <- sapply(c(-Inf, tau, Inf), function(t) Phi2(t - xb, -wg, rho))
  if (is.null(dim(C))) C <- matrix(C, nrow = 1)
  P <- pmax(C[, -1, drop = FALSE] - C[, -ncol(C), drop = FALSE], 0)
  P / rowSums(P)
}

.HECK_CACHE <- new.env()

mice.impute.heckoprob <- function(y, ry, x, wy = NULL, iv = IV_VAR, exclude = NULL, ...) {
  s <- .drop_excluded(y, ry, x, wy, exclude); y <- s$y; ry <- s$ry; wy <- s$wy
  lev <- sort(unique(y[ry])); K <- length(lev)
  if (K == 1) return(rep(lev, sum(wy)))
  x <- .std(s$x, rep(TRUE, length(y)))
  ivc <- which(colnames(x) %in% make.names(iv))
  W <- cbind(1, x); X <- x[, setdiff(seq_len(ncol(x)), ivc), drop = FALSE]
  yk <- match(y, lev); yk[!ry] <- NA
  key <- paste(K, sum(ry), sum(yk[ry] * seq_len(sum(ry))), paste(colnames(W), collapse = ","))
  fit <- heckoprob_fit(yk, as.integer(ry), X, W, K, start = .HECK_CACHE[[key]])
  assign(key, fit$par, envir = .HECK_CACHE)
  th <- .mvn_draw(fit$par, fit$V)
  P <- heckoprob_probs(th, fit$idx, X[wy, , drop = FALSE], W[wy, , drop = FALSE])
  lev[1L + rowSums(runif(nrow(P)) > t(apply(P, 1, cumsum))[, -K, drop = FALSE])]
}

## 6. Item-level imputation
# Imputation model of every incomplete item (PRED_SET = "parcel", default):
#   - the other items of the same scale,
#   - the scale scores of the five other scales ("parcel summaries"; Eekhout,
#     de Vet, de Boer, Twisk & Heymans, 2018), updated passively inside the
#     chains, so that the dyadic structure (own and partner's X, M, Y) is kept,
#   - the auxiliaries (inclusive strategy; makes the MAR mechanism ignorable),
#   - Heckman-MI only: the instrument of the item's partner (selection equation).
# This keeps the imputation model congenial with an analysis of scale scores
# while using ~11 predictors instead of 26 (n = 130; Hardt, Herke & Leonhart,
# 2012, on too many predictors in small samples). PRED_SET = "items" uses all
# other items and the auxiliaries instead.

IMPUTATION_METHODS <- c("mice_pmm", "mice_cat", "heckman")

ITEM_METHOD <- c(mice_pmm = "pmm", mice_cat = "catx", heckman = "heckoprob", delta = "pmmx")

PRED_SET <- "parcel"

.mice_setup <- function(dat, item_method, include_iv = FALSE, pred_set = PRED_SET, aux = AUX_VARS) {
  passive <- pred_set == "parcel"
  inc <- ITEMS[colSums(is.na(dat[, ITEMS, drop = FALSE])) > 0]
  x <- dat[, c(ITEMS, aux, if (include_iv) IV_VARS), drop = FALSE]
  if (passive) {
    sc <- scale_scores(dat)                      # NA where any item is missing
    for (v in VARS) x[[SCALES[[v]]]] <- sc[, v]
    x <- x[, c(ITEMS, unname(SCALES), aux, if (include_iv) IV_VARS), drop = FALSE]
  }
  meth <- setNames(rep("", ncol(x)), names(x)); meth[inc] <- item_method
  if (passive) for (v in VARS) {
    its <- ITEMS[ITEM_VAR == v]
    if (any(its %in% inc))
      meth[SCALES[[v]]] <- sprintf("~ I((%s) / %d)", paste(its, collapse = " + "), length(its))
  }
  pred <- matrix(0L, ncol(x), ncol(x), dimnames = list(names(x), names(x)))
  for (it in inc) {
    v <- ITEM_VAR[[it]]
    if (passive) {
      pred[it, setdiff(ITEMS[ITEM_VAR == v], it)] <- 1L
      pred[it, setdiff(unname(SCALES), SCALES[[v]])] <- 1L
    } else pred[it, setdiff(ITEMS, it)] <- 1L
    pred[it, aux] <- 1L
    if (include_iv) pred[it, paste0("iv_", ITEM_MEMBER[[it]])] <- 1L
  }
  # Visit each scale's items, then update that scale's score
  visit <- if (passive) unlist(lapply(VARS, function(v)
    c(intersect(ITEMS[ITEM_VAR == v], inc), if (meth[SCALES[[v]]] != "") SCALES[[v]])), use.names = FALSE)
  else inc
  list(x = x, method = meth, pred = pred, visit = visit, blots = mice::make.blots(x), incomplete = inc)
}

.complete_list <- function(imp, dat) {
  lapply(seq_len(imp$m), function(l) { d <- dat; d[, ITEMS] <- mice::complete(imp, l)[, ITEMS]; d })
}

impute_mice <- function(dat, method = c("mice_pmm", "mice_cat", "heckman"), m = 20, maxit = 20,
                        seed = NULL, return_mids = FALSE, pred_set = PRED_SET, aux = AUX_VARS, ...) {
  method <- match.arg(method)
  s <- .mice_setup(dat, ITEM_METHOD[[method]], include_iv = method == "heckman", pred_set, aux)
  if (!length(s$incomplete)) return(rep(list(dat), m))
  args <- list(data = s$x, m = m, maxit = maxit, method = s$method, predictorMatrix = s$pred,
               visitSequence = s$visit, blots = s$blots, seed = NA, printFlag = FALSE, ...)
  rm(list = ls(.CAT_CACHE), envir = .CAT_CACHE); rm(list = ls(.HECK_CACHE), envir = .HECK_CACHE)
  if (method == "heckman") args$iv <- IV_VARS
  # seed: own L'Ecuyer-CMRG stream (identical in the main session and in parallel
  # workers); NULL: continue the current stream (the Monte Carlo replications)
  imp <- if (is.null(seed)) suppressWarnings(do.call(mice::mice, args)) else
    with_rng(RNG_IMPUTATION, seed, suppressWarnings(do.call(mice::mice, args)))
  out <- .complete_list(imp, dat); if (return_mids) attr(out, "mids") <- imp
  out
}

# Delta-adjusted mice-PMM for the MNAR sensitivity analysis (van Buuren, 2018,
# sec. 3.8; Leacy, Floyd, Yates & White, 2017). The imputed values of every item
# are shifted by delta SDs (of the observed item values) in the hypothesised
# MNAR direction of its construct, stochastically rounded and kept in 1-5. The
# shift is applied inside the chains, so shifted values also inform the other
# items. delta = 0 is mice-PMM (same seed -> same imputations).

impute_delta <- function(dat, delta = 0, m = 20, maxit = 20, seed = NULL,
                         direction = MNAR_DIRECTION, pred_set = PRED_SET, aux = AUX_VARS,
                         donors = 5L) {
  s <- .mice_setup(dat, "pmmx", FALSE, pred_set, aux)
  for (it in s$incomplete)
    s$blots[[it]] <- list(delta = delta * direction[[ITEM_CONSTRUCT[[it]]]], bounds = range(ITEM_VALUES))
  run <- function() suppressWarnings(mice::mice(s$x, m = m, maxit = maxit, method = s$method,
                                                predictorMatrix = s$pred, visitSequence = s$visit,
                                                blots = s$blots, seed = NA, printFlag = FALSE,
                                                donors = donors))
  imp <- if (is.null(seed)) run() else with_rng(RNG_IMPUTATION, seed, run())   # as impute_mice()
  .complete_list(imp, dat)
}

# Convergence: potential scale reduction factor of the chain means of each
# imputed item over the second half of the iterations (Gelman & Rubin, 1992)

mice_psrf <- function(mids) {
  cm <- mids$chainMean
  if (is.null(cm) || dim(cm)[2] < 4 || dim(cm)[3] < 2) return(NULL)
  it <- seq(floor(dim(cm)[2] / 2) + 1, dim(cm)[2])
  v <- intersect(dimnames(cm)[[1]], ITEMS)
  out <- sapply(v, function(vv) {
    x <- cm[vv, it, , drop = TRUE]
    if (anyNA(x) || any(apply(x, 2, stats::var) == 0)) return(NA_real_)
    n <- nrow(x); W <- mean(apply(x, 2, stats::var)); B <- n * stats::var(colMeans(x))
    sqrt(((n - 1) / n * W + B / n) / W)
  })
  out[!is.na(out)]
}

# Relevance of the instruments: z of iv_pk in a probit model for "item
# observed" with the auxiliaries as covariates (Heckman-MI needs |z| clearly > 2)
iv_first_stage <- function(dat) {
  do.call(rbind, lapply(ITEMS, function(v) {
    r <- as.integer(!is.na(dat[[v]]))
    if (all(r == 1)) return(data.frame(item = v, z_iv = NA_real_))
    iv <- dat[[paste0("iv_", ITEM_MEMBER[[v]])]]
    f <- suppressWarnings(stats::glm(r ~ dat$aux_p1 + dat$aux_p2 + dat$aux_dyad + iv,
                                     family = stats::binomial(link = "probit")))
    data.frame(item = v, z_iv = summary(f)$coefficients["iv", 3])
  }))
}

## 7. Study-level APIMeM (lavaan)
# All methods return the structure of analyse_study() in apim_meta_analysis.R:
#   list(uncorrected = list(est, V = list(ML, MLR), admissible),
#        corrected   = list(est, V = list(ML, MLR), admissible), ...)
# so that the results plug into the meta-analysis and into
# apim_meta_sensitivity.R unchanged.

# (a) Complete or completed (imputed) data, and (b) PRORATION + FIML: the
#     default handling of apim_meta_analysis.R. prep_study() prorates scale
#     scores with >= 75% observed items (NA otherwise) and computes omega from a
#     FIML CFA of the items; the APIMeM is then fitted with FIML (missing = "ml")
#     whenever a scale score is missing. No auxiliaries.

analyse_complete <- function(dat, study, ml_only = FALSE) {
  ps <- APIM$prep_study(dat[, c(ID_VAR, ITEMS)], study, APIM_CFG)
  res <- if (ml_only) analyse_study_ml(ps) else APIM$analyse_study(ps)
  res$rel <- ps$rel; res$n <- ps$n
  res
}

analyse_prorated <- analyse_complete            # same code, applied to incomplete items

# ML sampling covariances only (2 instead of 4 lavaan fits; for the delta
# analysis and the Monte Carlo study, which use V = ML)
analyse_study_ml <- function(ps) {
  f_u <- APIM$fit_apim(ps$comp, estimator = "ML")
  s2  <- diag(lavaan::lavInspect(f_u, "cov.ov"))[VARS]
  f_c <- APIM$fit_apim(ps$comp, ps$rel, corrected = TRUE, estimator = "ML", s2 = s2)
  pk <- function(f) { e <- APIM$extract_paths(f); list(est = e$est, V = list(ML = e$V), admissible = e$admissible) }
  list(uncorrected = pk(f_u), corrected = pk(f_c), s2 = s2)
}

# (c) TWO-STAGE item-level FIML (Savalei & Bentler, 2009; Savalei & Rhemtulla,
#     2017, "the analytic equivalent of item-level MI").
#   Stage 1: saturated normal model of the 24 items AND the auxiliaries by
#            FIML (EM; lavaan's h1 estimates), with the observed information
#            (valid under MAR; Kenward & Molenberghs, 1998; Savalei, 2010) for
#            V = ML and the Huber-White sandwich for V = MLR.
#   Stage 2: the scale-score moments are linear functions of the item moments
#            (mean_c = W mu, Sigma_c = W Sigma W'); their asymptotic covariance
#            matrix follows exactly. The APIMeM (uncorrected and corrected) is
#            fitted in lavaan to these moments, with the stage-1 asymptotic
#            covariance matrix as NACOV (se = "robust.sem"), so that the
#            standardized paths and their 12 x 12 covariance matrix come from
#            extract_paths() exactly as for complete data. Because the APIMeM is
#            saturated, this equals the delta method (checked in 04).
#   Omega:   from a one-factor CFA of each scale's 4 x 4 block of the stage-1
#            covariance matrix (so it, too, uses the auxiliaries).
# Tight convergence of the EM algorithm for the saturated (h1) model. lavaan
# 0.7-1 moved these settings into a sublist (em.h1.args = list(tol, max_iter))
# and rejects the old top-level em.h1.tol / em.h1.iter.max, so the argument is
# built for the installed version.

lav_em_h1_args <- function(tol = 1e-9, max_iter = 10000L) {
  o <- names(lavaan::lavOptions())
  if ("em.h1.args" %in% o) return(list(em.h1.args = list(tol = tol, max_iter = max_iter)))
  out <- list()
  if ("em.h1.tol" %in% o) out$em.h1.tol <- tol
  if ("em.h1.iter.max" %in% o) out$em.h1.iter.max <- max_iter
  out
}

fiml_stage1 <- function(dat, vars) {
  args <- c(list(model = paste(sprintf("%s ~~ %s", vars, vars), collapse = "\n"),
                 data = dat[, vars], missing = "ml", meanstructure = TRUE, int.ov.free = TRUE,
                 h1.information = "unstructured"), lav_em_h1_args())
  f0 <- suppressWarnings(do.call(lavaan::lavaan, args))
  h1 <- lavaan::lavInspect(f0, "sampstat.h1")
  A  <- lavaan::lavInspect(f0, "h1.information.observed")
  B  <- lavaan::lavInspect(f0, "h1.information.firstorder")
  N  <- lavaan::lavInspect(f0, "nobs")
  Ai <- solve(A)
  list(Sigma = as.matrix(h1$cov)[vars, vars], mu = as.numeric(h1$mean[vars]), N = N,
       names = colnames(A), acov = list(ML = Ai / N, MLR = Ai %*% B %*% Ai / N))
}

fiml_composite_moments <- function(st, items = ITEMS) {
  vech <- function(S) S[lower.tri(S, diag = TRUE)]
  names(st$mu) <- rownames(st$Sigma)
  W <- matrix(0, 6, length(items), dimnames = list(VARS, items))
  for (v in VARS) { its <- items[ITEM_VAR[items] == v]; W[v, its] <- 1 / length(its) }
  Sc <- W %*% st$Sigma[items, items] %*% t(W)
  mc <- drop(W %*% st$mu[items])
  L <- matrix(0, 6 + 21, length(st$names))       # d (mean_c, vech Sigma_c) / d stage-1 parameters
  E0 <- matrix(0, length(items), length(items), dimnames = list(items, items))
  for (p in seq_along(st$names)) {
    nm <- st$names[p]
    if (grepl("~1$", nm)) { v <- sub("~1$", "", nm); if (v %in% items) L[1:6, p] <- W[, v]; next }
    ab <- strsplit(nm, "~~", fixed = TRUE)[[1]]
    if (length(ab) != 2 || !all(ab %in% items)) next
    E <- E0; E[ab[1], ab[2]] <- 1; E[ab[2], ab[1]] <- 1
    L[7:27, p] <- vech(W %*% E %*% t(W))
  }
  list(Sigma = Sc, mean = mc, N = st$N, L = L,
       Gamma = lapply(st$acov, function(A) st$N * L %*% A %*% t(L)))
}

omega_from_cov <- function(S, N, type = APIM_CFG$reliability) {
  k <- ncol(S); nm <- paste0("i", seq_len(k)); dimnames(S) <- list(nm, nm)
  alpha <- (k / (k - 1)) * (1 - sum(diag(S)) / sum(S))
  if (type == "alpha" || k < 3) return(alpha)
  fit <- tryCatch(suppressWarnings(lavaan::cfa(paste("F =~", paste(nm, collapse = " + ")),
                                               sample.cov = S, sample.nobs = N, std.lv = TRUE,
                                               sample.cov.rescale = FALSE)),
                  error = function(e) NULL)
  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) return(alpha)
  pe <- lavaan::parameterEstimates(fit)
  lam <- pe$est[pe$op == "=~"]
  th  <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% nm]
  sum(lam)^2 / (sum(lam)^2 + sum(th))
}

fit_apim_moments <- function(Sc, mc, N, Gamma, theta = NULL) {
  syn <- if (is.null(theta)) APIM$apim_syntax("") else
    c(paste0("L", VARS, " =~ 1*", VARS), paste0(VARS, " ~~ ", sprintf("%.12f", theta[VARS]), "*", VARS),
      APIM$apim_syntax("L"))
  fit <- suppressWarnings(lavaan::sem(paste(syn, collapse = "\n"), sample.cov = Sc[VARS, VARS],
                                      sample.mean = mc[VARS], sample.nobs = N, NACOV = Gamma,
                                      se = "robust.sem", estimator = "ML", sample.cov.rescale = FALSE,
                                      meanstructure = TRUE, fixed.x = FALSE))
  if (!lavaan::lavInspect(fit, "converged")) stop("Stage-2 APIMeM did not converge")
  fit
}

analyse_fiml_ts <- function(dat, study = NULL, aux = AUX_VARS, items = ITEMS, ml_only = FALSE) {
  st <- fiml_stage1(dat, c(items, aux))
  cm <- fiml_composite_moments(st, items)
  rel <- vapply(VARS, function(v) {
    its <- items[ITEM_VAR[items] == v]; omega_from_cov(st$Sigma[its, its], st$N) }, 0)
  s2 <- diag(cm$Sigma)[VARS]
  pk <- function(theta) {
    vts <- if (ml_only) "ML" else c("ML", "MLR")
    ex <- lapply(vts, function(vt) APIM$extract_paths(fit_apim_moments(cm$Sigma, cm$mean, cm$N,
                                                                        cm$Gamma[[vt]], theta)))
    list(est = ex[[1]]$est, V = setNames(lapply(ex, `[[`, "V"), vts),
         admissible = all(vapply(ex, `[[`, TRUE, "admissible")))
  }
  list(uncorrected = pk(NULL), corrected = pk((1 - rel) * s2), s2 = s2, rel = rel, n = st$N,
       moments = cm[c("Sigma", "mean", "N", "Gamma")])
}

# Delta-method version of stage 2 (no lavaan), used only to check fit_apim_moments()

delta_check_fiml_ts <- function(res) {
  cm <- res$moments; vech_idx <- 7:27
  unvech <- function(v) { S <- matrix(0, 6, 6); S[lower.tri(S, diag = TRUE)] <- v
    S <- S + t(S) - diag(diag(S)); dimnames(S) <- list(VARS, VARS); S }
  v0 <- cm$Sigma[lower.tri(cm$Sigma, diag = TRUE)]
  one <- function(theta, cmp) {
    f <- function(v) { S <- unvech(v); if (!is.null(theta)) diag(S) <- diag(S) - theta; std_paths_from_cov(S) }
    J <- numDeriv::jacobian(f, v0)
    V <- J %*% (cm$Gamma$ML[vech_idx, vech_idx] / cm$N) %*% t(J)
    c(est = max(abs(f(v0) - res[[cmp]]$est)), V = max(abs(V - res[[cmp]]$V$ML)))
  }
  rbind(uncorrected = one(NULL, "uncorrected"),
        corrected = one((1 - res$rel) * res$s2, "corrected"))
}

# (d) Multiple imputation: Rubin's rules WITHIN each study (Burgess, White,
#     Resche-Rigon & Wood, 2013: impute and pool at the study level, then
#     meta-analyse). For each component and V type: est = mean over the m
#     analyses; V = Ubar + (1 + 1/m) B, with Ubar the mean 12 x 12 sampling
#     covariance matrix and B the between-imputation covariance matrix of the 12
#     paths. Per path: Barnard-Rubin df, FMI and lambda (mice::pool.scalar).

pool_rubin_study <- function(res_list, n, k_com = 5) {
  m <- length(res_list); lab <- PATHS$label
  out <- list()
  for (cmp in COMPONENTS) {
    E <- do.call(rbind, lapply(res_list, function(r) r[[cmp]]$est[lab]))
    qbar <- colMeans(E)
    B <- if (m > 1) stats::cov(E) else matrix(0, length(lab), length(lab))
    dimnames(B) <- list(lab, lab)
    vts <- names(res_list[[1]][[cmp]]$V)
    V <- info <- list()
    for (vt in vts) {
      Ubar <- Reduce(`+`, lapply(res_list, function(r) r[[cmp]]$V[[vt]][lab, lab])) / m
      Tm <- Ubar + (1 + 1 / m) * B
      V[[vt]] <- Tm
      ps <- t(vapply(seq_along(lab), function(j) {
        if (m == 1) return(c(df = Inf, riv = 0, fmi = 0))
        u <- vapply(res_list, function(r) r[[cmp]]$V[[vt]][lab[j], lab[j]], 0)
        p <- mice::pool.scalar(E[, j], u, n = n, k = k_com)
        c(df = p$df, riv = p$r, fmi = p$fmi)
      }, numeric(3)))
      info[[vt]] <- data.frame(label = lab, ubar = diag(Ubar), b = diag(B), t = diag(Tm), ps,
                               lambda = (1 + 1 / m) * diag(B) / diag(Tm), row.names = NULL)
    }
    out[[cmp]] <- list(est = setNames(qbar, lab), V = V, mi = info, m = m,
                       admissible = all(vapply(res_list, function(r) isTRUE(r[[cmp]]$admissible), TRUE)),
                       prop_admissible = mean(vapply(res_list, function(r) isTRUE(r[[cmp]]$admissible), TRUE)))
  }
  out$rel <- colMeans(do.call(rbind, lapply(res_list, `[[`, "rel")))
  out$n <- n
  out
}

analyse_mi <- function(imp_list, study, ml_only = FALSE) {
  fits <- lapply(imp_list, analyse_complete, study = study, ml_only = ml_only)
  pool_rubin_study(fits, n = fits[[1]]$n)
}

# (e) "Aggregated" data: the item-wise MEAN of the m completed datasets, analysed
#     once as if it were complete. This is NOT a valid way to combine multiple
#     imputations: the imputed values are treated as known (no between-
#     imputation variance), and averaging removes their spread, which shrinks
#     the variances of the scale scores, inflates omega and their correlations
#     and understates the SEs (van Buuren, 2018, sec. 5.1.2). It is included only
#     to show, on these data, why the m = 20 imputations are aggregated with
#     Rubin's rules (analyse_mi) instead.

average_imputations <- function(imps) {
  d <- imps[[1]]
  d[, ITEMS] <- Reduce(`+`, lapply(imps, function(x) as.matrix(x[, ITEMS]))) / length(imps)
  d
}

# The methods compared in 04-08. Listwise deletion is not among them: with
# 5-15% missing values per item, only 1-44% of the dyads have all 24 items.

METHODS_ALL <- c("COM", "PRO_FIML", "FIML_TS", "mice_pmm", "mice_cat", "heckman", "PMM_avg")

METHOD_LABELS <- c(COM = "Complete data (benchmark)",
                   PRO_FIML = "Proration + FIML (scale level; default of apim_meta_analysis.R)",
                   FIML_TS = "Two-stage FIML (item level, with auxiliaries)",
                   mice_pmm = "mice-PMM (item level, m = 20, Rubin's rules)",
                   mice_cat = "mice-CAT (item level, m = 20, Rubin's rules)",
                   heckman = "Heckman-MI (item level, instrument, m = 20, Rubin's rules)",
                   PMM_avg = "mice-PMM, averaged data (illustration only; not valid)")
METHOD_SHORT <- c(COM = "COM", PRO_FIML = "PRO-FIML", FIML_TS = "FIML-TS", mice_pmm = "PMM",
                  mice_cat = "CAT", heckman = "Heckman", PMM_avg = "PMM-avg")
# Which imputations a method needs (PMM_avg re-uses the mice-PMM imputations)
imputations_for <- function(method) if (method == "PMM_avg") "mice_pmm" else
  if (method %in% IMPUTATION_METHODS) method else NA_character_

# One method on one study's data (dat = complete data for COM, incomplete data
# otherwise; imps = list of completed datasets for the MI methods and PMM_avg)
run_method_study <- function(method, dat, study, imps = NULL, ml_only = FALSE) {
  if (!is.na(imputations_for(method)) && is.null(imps)) stop("Imputations are needed for ", method)
  switch(method,
    COM      = analyse_complete(dat, study, ml_only),
    PRO_FIML = analyse_prorated(dat, study, ml_only),
    FIML_TS  = analyse_fiml_ts(dat, study, ml_only = ml_only),
    PMM_avg  = analyse_complete(average_imputations(imps), study, ml_only),
    analyse_mi(imps, study, ml_only))
}

## 8. Meta-analysis with rma()
# As in apim_meta_sensitivity.R (the "rma()" family that the sensitivity report
# compared with rma.mv): per study, the Partner 1 and Partner 2 estimates of each
# effect type are combined into one GLS composite using the 2 x 2 block of the
# study's sampling covariance matrix; the composites are then pooled across
# studies with metafor::rma(). PRIMARY: random effects, REML, Wald z test
# ("rma(): RE (z)"), the specification chosen from the sensitivity report.
# The common-effect model and the Knapp-Hartung tests are kept for reference.

RMA_METHODS <- c(RE_z = "rma(): RE (z)", CE_z = "rma(): CE (z)",
                 RE_KH = "rma(): RE (KH)", RE_KHadh = "rma(): RE (KH ad hoc)")

PRIMARY_RMA <- "RE_z"

PRIMARY_V   <- "ML"            # as CFG$V_type in apim_meta_sensitivity.R

partner_composites <- function(study_res, cmp, vtype = PRIMARY_V) {
  lapply(study_res, function(r) {
    est <- r[[cmp]]$est[PATHS$label]; V <- r[[cmp]]$V[[vtype]][PATHS$label, PATHS$label]
    res <- lapply(EFFECTS, function(e) {
      cols <- EFF_COLS[[e]]; u <- solve(V[cols, cols], c(1, 1))
      list(y = sum(u / sum(u) * est[cols]), v = 1 / sum(u), w = setNames(u / sum(u), cols))
    })
    names(res) <- EFFECTS
    list(y = vapply(res, `[[`, 0, "y"), v = vapply(res, `[[`, 0, "v"), w = lapply(res, `[[`, "w"), V = V)
  })
}

.rma_re <- function(y, v, test) {
  f <- tryCatch(suppressWarnings(metafor::rma(y, v, method = "REML", test = test,
                                              control = list(stepadj = 0.5, maxiter = 1000))),
                error = function(e) NULL)
  if (is.null(f)) f <- metafor::rma(y, v, method = "DL", test = test)     # fallback as in the sensitivity script
  f
}

rma_fit <- function(y, v, method) switch(method,
  CE_z = metafor::rma(y, v, method = "EE"), RE_z = .rma_re(y, v, "z"),
  RE_KH = .rma_re(y, v, "knha"), RE_KHadh = .rma_re(y, v, "adhoc"),
  stop("Unknown rma method ", method))

# Covariance matrix of the six pooled random-effects estimates: plug-in from the
# within-study covariances of the composites (random effects of different effect
# types assumed independent). Its diagonal equals the rma() variances. Needed for
# Monte Carlo CIs of products of pooled paths.

pooled_cov_re <- function(pc, tau2) {
  Wt <- sapply(EFFECTS, function(e) vapply(pc, function(s) 1 / (s$v[[e]] + tau2[[e]]), 0))
  if (is.null(dim(Wt))) Wt <- matrix(Wt, nrow = 1, dimnames = list(NULL, EFFECTS))
  a <- sweep(Wt, 2, colSums(Wt), "/")
  C <- matrix(0, 6, 6, dimnames = list(EFFECTS, EFFECTS))
  for (i in seq_along(pc)) {
    s <- pc[[i]]
    for (j in EFFECTS) for (l in EFFECTS) {
      cij <- if (j == l) s$v[[j]] + tau2[[j]] else
        drop(crossprod(s$w[[j]], s$V[EFF_COLS[[j]], EFF_COLS[[l]]] %*% s$w[[l]]))
      C[j, l] <- C[j, l] + a[i, j] * a[i, l] * cij
    }
  }
  C
}

# Indirect and total effects from the pooled RE paths: products/sums of the
# pooled estimates (Ledermann, Macho & Kenny, 2011), Monte Carlo CIs, and the
# joint-significance test of the specific indirect effects (both component paths
# significant; Fritz & MacKinnon, 2007; Yzerbyt et al., 2018), which is the test
# evaluated in the sensitivity report.
# Multivariate normal draws that are the same on every platform: the Cholesky
# factor is unique, whereas the eigenvectors used by MASS::mvrnorm() can come
# with opposite signs from different LAPACK libraries, which changes every draw.
# A singular matrix uses its symmetric square root (also sign-free).

rmvn_portable <- function(n, mu, Sigma) {
  p <- length(mu)
  R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(R)) {
    e <- eigen(Sigma, symmetric = TRUE)
    R <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
  }
  out <- matrix(stats::rnorm(n * p), n, p) %*% R + matrix(mu, n, p, byrow = TRUE)
  colnames(out) <- names(mu)
  out
}

indirect_rma <- function(b, C, sig, n_mc = 2e4, level = .95, seed = 1) {
  draws <- with_seed(seed, tryCatch(rmvn_portable(n_mc, b[EFFECTS], C[EFFECTS, EFFECTS]),
                                    error = function(e) NULL))
  a <- (1 - level) / 2
  rows <- lapply(seq_len(nrow(INDIRECT_DEFS)), function(i) {
    ex <- parse(text = INDIRECT_DEFS$expr[i]); est <- eval(ex, as.list(b))
    if (is.null(draws)) return(c(est = est, lo = NA, hi = NA))
    sim <- eval(ex, as.data.frame(draws))
    c(est = est, lo = unname(stats::quantile(sim, a)), hi = unname(stats::quantile(sim, 1 - a)))
  })
  out <- cbind(INDIRECT_DEFS[, c("group", "label", "short", "expr")], as.data.frame(do.call(rbind, rows)))
  out$mc_sig <- out$lo > 0 | out$hi < 0
  js <- vapply(out$short, function(s) {
    tp <- SPECIFIC_INDIRECT[[s]]; if (is.null(tp)) NA else all(sig[tp]) }, NA)
  out$joint_sig <- unname(js)
  out
}

# Pool one set of study-level results (list named by study) for one component

meta_rma <- function(study_res, cmp, vtype = PRIMARY_V, methods = names(RMA_METHODS),
                     level = .95, n_mc = 2e4, mc_seed = 1) {
  pc <- partner_composites(study_res, cmp, vtype)
  Y  <- do.call(rbind, lapply(pc, `[[`, "y")); Vv <- do.call(rbind, lapply(pc, `[[`, "v"))
  rows <- list(); re <- list()
  for (e in EFFECTS) for (me in methods) {
    f <- rma_fit(Y[, e], Vv[, e], me)
    rows[[length(rows) + 1]] <- data.frame(
      component = cmp, V_type = vtype, method = me, effect = e, k = f$k,
      est = as.numeric(f$b), se = f$se, lo = f$ci.lb, hi = f$ci.ub, stat = f$zval, p = f$pval,
      tau2 = f$tau2, I2 = if (is.null(f$I2)) NA_real_ else f$I2, QE = f$QE, QEp = f$QEp,
      stringsAsFactors = FALSE)
    if (me == PRIMARY_RMA) re[[e]] <- f
  }
  tab <- do.call(rbind, rows)
  out <- list(direct = tab, composites = list(y = Y, v = Vv))
  if (PRIMARY_RMA %in% methods) {
    b <- vapply(re, function(f) as.numeric(f$b), 0); tau2 <- vapply(re, `[[`, 0, "tau2")
    sig <- vapply(re, function(f) f$pval < 1 - level, TRUE)
    C <- pooled_cov_re(pc, tau2); diag(C) <- vapply(re, function(f) f$se^2, 0)
    out$cov <- C
    out$indirect <- cbind(component = cmp, V_type = vtype,
                          indirect_rma(b, C, sig, n_mc = n_mc, level = level, seed = mc_seed))
  }
  out
}

meta_both <- function(study_res, vtypes = PRIMARY_V, ...) {
  res <- list()
  for (cmp in COMPONENTS) for (vt in vtypes) res[[paste(cmp, vt)]] <- meta_rma(study_res, cmp, vt, ...)
  list(direct = do.call(rbind, lapply(res, `[[`, "direct")),
       indirect = do.call(rbind, lapply(res, `[[`, "indirect")),
       cov = lapply(res, `[[`, "cov"))
}

# Study-level results in the format of APIM_meta_analysis_results.rds, so that
# apim_meta_sensitivity.R can be run on missing-data results
# (CFG_OVERRIDE = list(rds_file = <file>)).

write_sensitivity_input <- function(study_res, file, note = "") {
  rr <- list(config = c(APIM_CFG, list(pool_roles = TRUE, missing_data_note = note)),
             prepped = lapply(study_res, function(r) list(rel = r$rel, n = r$n)),
             study_results = lapply(study_res, function(r) r[c("uncorrected", "corrected")]),
             meta = NULL)
  saveRDS(rr, file); invisible(file)
}

## 9. Performance measures
# Quality of the imputed item values against the deleted true values, per item,
# averaged over imputations. nrmse: RMSE of the imputed values / SD of the true
# item; pct_exact: share of imputed values equal to the deleted value;
# std_bias_mean: (completed mean - complete mean) / SD; rel_bias_sd: SD ratio - 1;
# tvd: total variation distance of the 1-5 distributions; mabias_cor: mean
# |bias| of the item's correlations with all other items, scale scores and
# auxiliaries.

num_frame <- function(d) cbind(d[, ITEMS], scale_scores(d), d[, intersect(AUX_VARS, names(d)), drop = FALSE])

imputation_quality <- function(complete, amputed, imp_list) {
  M <- is.na(amputed[, ITEMS, drop = FALSE])
  Xc <- num_frame(complete); Rc <- stats::cor(Xc)
  per <- lapply(imp_list, function(d) {
    Xi <- num_frame(d); Ri <- stats::cor(Xi)
    t(vapply(ITEMS, function(v) {
      mi <- M[, v]; tv <- Xc[[v]]; iv <- Xi[[v]]
      tvd <- 0.5 * sum(abs(prop.table(table(factor(iv, ITEM_VALUES))) -
                             prop.table(table(factor(tv, ITEM_VALUES)))))
      if (!any(mi)) return(c(nrmse = NA, pct_exact = NA, std_bias_mean = 0, rel_bias_sd = 0,
                             tvd = tvd, mabias_cor = 0))
      c(nrmse = sqrt(mean((iv[mi] - tv[mi])^2)) / stats::sd(tv),
        pct_exact = 100 * mean(iv[mi] == tv[mi]),
        std_bias_mean = (mean(iv) - mean(tv)) / stats::sd(tv),
        rel_bias_sd = stats::sd(iv) / stats::sd(tv) - 1, tvd = tvd,
        mabias_cor = mean(abs((Ri[v, ] - Rc[v, ])[colnames(Rc) != v])))
    }, numeric(6)))
  })
  out <- Reduce(`+`, per) / length(per)
  data.frame(item = ITEMS, construct = ITEM_CONSTRUCT, scale = ITEM_SCALE, out, row.names = NULL)
}

# Compare pooled results with the complete-data (COM) results of the same data
compare_to_com <- function(tab, com, keys) {
  ref <- com[, c(keys, "est", "se", "lo", "hi", "p")]
  names(ref)[-seq_along(keys)] <- paste0(c("est", "se", "lo", "hi", "p"), "_com")
  x <- merge(tab, ref, by = keys, sort = FALSE)
  x$dev <- x$est - x$est_com
  x$dev_std <- x$dev / x$se_com                   # in complete-data SE units
  x$se_ratio <- x$se / x$se_com
  x$ci_ratio <- (x$hi - x$lo) / (x$hi_com - x$lo_com)
  x$sig <- x$p < .05; x$sig_com <- x$p_com < .05
  x$sig_agree <- x$sig == x$sig_com
  x
}

# Method decision (pre-specified): mice-PMM versus two-stage FIML, judged by how
# close the pooled direct effects (rma() RE z, V = ML) come to the complete-data
# results under MCAR and MAR, where both methods are valid and the complete-data
# analysis of the same dyads is the natural target. Basis: the Monte Carlo
# summary of 05 (pilot or full mode; root mean square deviation from the
# complete-data estimate over replications) when available, otherwise the
# single realisation of 04 (root mean square deviation in complete-data SE
# units). A difference within REL_MARGIN of the FIML value is a tie, and a tie
# goes to mice-PMM (completed datasets, which the delta-adjusted MNAR analysis
# needs). Safeguard (Monte Carlo only): a method whose 95% CI coverage under
# MCAR/MAR falls below .925 (Bradley's liberal criterion) is not chosen if the
# other method's coverage is at least .925.

REL_MARGIN <- 0.10

choose_method <- function(E04 = NULL, MC = NULL, rel_margin = REL_MARGIN) {
  pair <- c("mice_pmm", "FIML_TS")
  crit <- NULL; note <- character(0)
  usable <- function(cr) all(is.finite(cr$rmsd_mcar_mar))
  # a method without results (e.g. failed in every replication) gives NaN
  missing_note <- function(cr, where)
    sprintf("Note: %s produced no results in %s, so these were not used for the decision.",
            paste(METHOD_SHORT[cr$method[!is.finite(cr$rmsd_mcar_mar)]], collapse = " and "), where)
  if (!is.null(MC) && all(pair %in% MC$methods)) {
    bm <- MC$by_method[MC$by_method$method == PRIMARY_RMA, ]
    cr <- do.call(rbind, lapply(pair, function(mt) {
      x <- bm[bm$missing_method == mt, ]; a <- x$mechanism != "MNAR"
      data.frame(method = mt, rmsd_mcar_mar = mean(x$rmsd_vs_com[a]), rmsd_mnar = mean(x$rmsd_vs_com[!a]),
                 coverage_mcar_mar = mean(x$coverage[a]), ci_ratio_mcar_mar = mean(x$ci_ratio[a]))
    }))
    if (usable(cr)) {
      crit <- cr
      basis <- sprintf("Monte Carlo simulation (05, %s mode, %s replications per condition)", MC$mode,
                       paste(unique(MC$n_rep), collapse = "/"))
      units <- "pooled-beta units"
    } else note <- c(note, missing_note(cr, sprintf(
      "the Monte Carlo simulation (05, %s mode; see the 'errors' sheet of mc_results.xlsx)", MC$mode)))
  }
  if (is.null(crit) && !is.null(E04$key) && all(pair %in% E04$methods)) {
    k <- E04$key[E04$key$V_type == PRIMARY_V & E04$key$missing_method %in% pair, ]
    cr <- do.call(rbind, lapply(pair, function(mt) {
      x <- k[k$missing_method == mt, ]; a <- x$mechanism != "MNAR"
      data.frame(method = mt, rmsd_mcar_mar = mean(x$rmsd_std[a]), rmsd_mnar = mean(x$rmsd_std[!a]),
                 coverage_mcar_mar = NA_real_, ci_ratio_mcar_mar = mean(x$mean_ci_ratio[a]))
    }))
    if (usable(cr)) {
      crit <- cr
      basis <- "one realisation of the missing data per condition (04)"
      units <- "complete-data SE units"
    } else note <- c(note, missing_note(cr, "the evaluation (04; see the 'runtime' sheet of its .xlsx)"))
  }
  if (is.null(crit))
    return(list(chosen = "mice_pmm", basis = "default (no usable evaluation results)", units = NA_character_,
                criteria = NULL, difference = NA_real_, margin = NA_real_, verdict = "no evaluation",
                text = c(note, "No usable evaluation results for both methods: mice-PMM is used (default).")))
  rownames(crit) <- pair
  d <- crit["mice_pmm", "rmsd_mcar_mar"] - crit["FIML_TS", "rmsd_mcar_mar"]
  margin <- rel_margin * crit["FIML_TS", "rmsd_mcar_mar"]
  chosen <- if (d <= margin) "mice_pmm" else "FIML_TS"
  verdict <- if (abs(d) <= margin) "equivalent (difference within the margin)" else
    if (d < 0) "mice-PMM closer to the complete-data results" else "two-stage FIML closer to the complete-data results"
  guard <- ""
  cv <- crit$coverage_mcar_mar; names(cv) <- pair
  if (all(!is.na(cv))) {
    other <- setdiff(pair, chosen)
    if (cv[[chosen]] < .925 && cv[[other]] >= .925) {
      guard <- sprintf(" Coverage safeguard: %s covered %.3f < .925 under MCAR/MAR, so %s is chosen instead.",
                       METHOD_SHORT[[chosen]], cv[[chosen]], METHOD_SHORT[[other]])
      chosen <- other
    }
  }
  fmt3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
  text <- c(
    sprintf("Decision rule (pre-specified): root mean square deviation (RMSD) of the six pooled direct effects"),
    sprintf("(rma() RE z, V = %s) from the complete-data results, averaged over both components and the MCAR", PRIMARY_V),
    sprintf("and MAR conditions; basis: %s; units: %s.", basis, units),
    sprintf("A difference within %.0f%% of the FIML value is a tie, and a tie goes to mice-PMM.", 100 * rel_margin),
    "",
    sprintf("  mice-PMM : RMSD %s (MNAR %s), coverage %s, CI width ratio %s", fmt3(crit[1, 2]), fmt3(crit[1, 3]),
            fmt3(crit[1, 4]), fmt3(crit[1, 5])),
    sprintf("  FIML-TS  : RMSD %s (MNAR %s), coverage %s, CI width ratio %s", fmt3(crit[2, 2]), fmt3(crit[2, 3]),
            fmt3(crit[2, 4]), fmt3(crit[2, 5])),
    sprintf("  Difference (PMM - FIML) = %+.3f, margin %.3f: %s.%s", d, margin, verdict, guard),
    sprintf("  Chosen method: %s%s", METHOD_LABELS[[chosen]],
            if (chosen == "mice_pmm") "; the m = 20 imputations are aggregated with Rubin's rules within each study before the meta-analysis." else "."))
  text <- c(note, text)
  list(chosen = chosen, basis = basis, units = units, criteria = crit, difference = d, margin = margin,
       verdict = verdict, text = text)
}

## 10. Power sensitivity: link to apim_meta_sensitivity.R

SENS_SCRIPT <- find_input("apim_meta_sensitivity.R")      # this folder, then ../Sensitivity
# Method codes of apim_meta_sensitivity.R -> codes used here

RMA_CODES <- c(rma_CE = "CE_z", rma_RE = "RE_z", rma_KH = "RE_KH", rma_KHadh = "RE_KHadh")

# Run apim_meta_sensitivity.R (its sections 0-9: inputs, heterogeneity,
# simulation, power curves, MDE, critical values, indirect effects) unchanged,
# in its own environment, on a study-level input file in the format of
# APIM_meta_analysis_results.rds (write_sensitivity_input()). Set from outside:
#  - CFG_OVERRIDE, the script's own mechanism: rds_file and, from 07, all
#    settings of the original run (stored in APIM_sensitivity_results.rds), so
#    that a rerun on the complete-data input reproduces the report;
#  - rma_only = TRUE replaces the script's fit_mv() by a function returning NULL,
#    which skips the multilevel rma.mv() model (no longer used; most of the run
#    time). fit_all() computes the four rma() methods before and independently of
#    rma.mv(), and all random numbers are drawn up front, so every rma() result is
#    unchanged (07 checks this against the report).
# Sections 10-12 (console summary, figures, .docx) are not run. The script draws
# its random numbers with set.seed(), so the default RNG kind is set for the run
# and the caller's RNG state is restored afterwards.

run_sensitivity_script <- function(rds_file, override = list(), rma_only = TRUE, file = SENS_SCRIPT,
                                   keep_power = "rma_RE") {
  if (!file.exists(file)) stop(not_found(file))
  if (!file.exists(rds_file)) stop("Input file not found: ", rds_file)
  src <- readLines(file, warn = FALSE, encoding = "UTF-8")
  i5  <- grep("^## 5\\. OBSERVED DATA", src)[1]; i10 <- grep("^## 10\\. CONSOLE SUMMARY", src)[1]
  if (is.na(i5) || is.na(i10)) stop(file, ": section headers '## 5.' and '## 10.' not found")
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  old_kind <- RNGkind()
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  env <- new.env(parent = globalenv())
  env$CFG_OVERRIDE <- utils::modifyList(override, list(rds_file = rds_file))
  t0 <- proc.time()[["elapsed"]]
  suppressMessages(suppressPackageStartupMessages(
    eval(parse(text = src[seq_len(i5 - 2)], encoding = "UTF-8"), envir = env)))
  if (rma_only) env$fit_mv <- function(...) NULL
  # With rma_only the multilevel columns are empty by design, which the script
  # would report as failed fits: that warning is muffled (and only that one)
  out_txt <- utils::capture.output(withCallingHandlers(suppressMessages(
    eval(parse(text = src[(i5 - 1):(i10 - 2)], encoding = "UTF-8"), envir = env)),
    warning = function(w) if (rma_only && grepl("Some fits failed", conditionMessage(w)))
      invokeRestart("muffleWarning")))
  keep_lab <- unname(env$METHODS[if (rma_only) names(RMA_CODES) else names(env$METHODS)])
  sm <- env$SUMMARY[env$SUMMARY$method %in% keep_lab, ]
  pw <- env$POWER[env$POWER$method %in% unname(env$METHODS[keep_power]), ]
  im <- env$IND_MIN[env$IND_MIN$method %in% keep_lab, ]
  ob <- env$OBS_TAB[env$OBS_TAB$method %in% keep_lab, ]
  list(summary = sm, power = pw, indirect_min = im, observed = ob, tau2 = env$TAU2,
       v_typical = env$v_typical, fail_rate = env$fail_rate, config = env$CFG,
       n_dyads = env$INPUT$n_dyads, rma_only = rma_only, seconds = proc.time()[["elapsed"]] - t0)
}

# Rescaling of published power results. If every within-study sampling
# (co)variance is multiplied by the same factor c^2, and the between-study
# variance is set through I^2 (as in the report, so tau^2 scales with the typical
# sampling variance), every rma() test (CE, RE z, Knapp-Hartung) is unchanged
# when all effects are divided by c: power_missing(theta) = power(theta / c),
# MDE_missing = c * MDE and critical |estimate|_missing = c * critical |estimate|.
# With inflation that differs between studies and paths it is an approximation;
# 07 compares it with the exact rerun. c = SE of the common-effect estimate under
# missing data / its complete-data value (inverse-variance weighted inflation).

rescale_power <- function(theta, power, c, at) stats::approx(theta, power, xout = at / c, rule = 2, ties = max)$y

## 11. Output helpers

write_xlsx_sheets <- function(sheets, file) {
  wb <- openxlsx::createWorkbook()
  hs <- openxlsx::createStyle(textDecoration = "bold")
  for (nm in names(sheets)) {
    x <- sheets[[nm]]; s <- substr(nm, 1, 31)
    openxlsx::addWorksheet(wb, s); openxlsx::writeData(wb, s, x, headerStyle = hs)
    openxlsx::freezePane(wb, s, firstRow = TRUE)
    openxlsx::setColWidths(wb, s, seq_len(ncol(x)), "auto")
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  invisible(file)
}

par_lapply <- function(X, FUN, workers = 1, ...) {
  if (workers > 1) {
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    future.apply::future_lapply(X, FUN, ..., future.seed = TRUE)
  } else lapply(X, FUN, ...)
}

# Cache keys. Every cached intermediate result (output/*/raw/) stores the key of
# the inputs it was computed from and is reused only if that key equals the key
# of the current inputs, so results computed from other data (e.g. after 01 or
# 02 were rerun) or with other settings are recomputed instead of being mixed
# in. The key is the MD5 hash of a text rendering of the inputs, which is the
# same on every platform for the same values (base R only).

CACHE_VERSION <- "2"                  # increase when a code change alters cached results

cache_key <- function(...) {
  f <- tempfile(fileext = ".txt"); con <- file(f, open = "wb")
  put <- function(x, tag) {
    writeLines(tag, con)
    if (is.data.frame(x) || is.matrix(x)) {
      utils::write.table(x, con, sep = ",", eol = "\n", na = "NA", quote = FALSE, row.names = FALSE)
    } else if (is.list(x)) {
      nm <- names(x) %||% rep("", length(x))
      for (k in seq_along(x)) put(x[[k]], paste0(tag, "/", nm[k], "#", k))
    } else {
      writeLines(paste(if (is.numeric(x)) format(x, digits = 15) else as.character(x), collapse = ","), con)
    }
  }
  tryCatch(put(list(CACHE_VERSION, ...), "key"), finally = close(con))
  on.exit(unlink(f))
  unname(tools::md5sum(f))
}

# Read a cached result if it exists and carries the expected key (NULL otherwise)

read_cached <- function(file, key) {
  if (!file.exists(file)) return(NULL)
  x <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.list(x) && identical(x$cache_key, key)) x else NULL
}

find_sensitivity_file <- function(file, dirs = SENS_DIRS) {
  for (d in dirs) { f <- file.path(d, file); if (file.exists(f)) return(f) }
  NULL
}

# APA 7 number formatting and tables (as in apim_meta_analysis.R)

APA_FONT <- "Times New Roman"; APA_SIZE <- 10; APA_MAX_WIDTH <- 8.9

MINUS <- "\u{2212}"

fmt <- function(x, d = 2, lead0 = TRUE) {
  s <- formatC(round(x, d), format = "f", digits = d)
  s <- sub("^-(0\\.0+)$", "\\1", s)
  if (!lead0) s <- sub("^(-?)0\\.", "\\1.", s)
  s <- sub("^-", MINUS, s)
  s[is.na(x)] <- "\u{2014}"
  s
}

fmt_p  <- function(p) ifelse(is.na(p), "\u{2014}", ifelse(p < .001, "< .001", fmt(p, 3, lead0 = FALSE)))

fmt_ci <- function(lo, hi, d = 2) ifelse(is.na(lo) | is.na(hi), "\u{2014}",
                                         sprintf("[%s, %s]", fmt(lo, d), fmt(hi, d)))

apa_table <- function(body, labels, spans = list(), group_rows = integer(0), note = NULL,
                      first_col_width = 2.4, col_width = 0.85, font = APA_FONT, size = APA_SIZE) {
  bd <- officer::fp_border(color = "black", width = 1)
  nc <- ncol(body)
  names(body) <- paste0("col", seq_len(nc))
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
  ft <- flextable::valign(ft, valign = "bottom", part = "header")
  pad <- if (nrow(body) > 16) 0.75 else 1.5
  ft <- flextable::padding(ft, padding.top = pad, padding.bottom = pad, part = "all")
  ft <- flextable::line_spacing(ft, space = 1, part = "all")
  ft <- flextable::hline_top(ft, border = bd, part = "header")
  ft <- flextable::hline_bottom(ft, border = bd, part = "header")
  ft <- flextable::hline_bottom(ft, border = bd, part = "body")
  for (r in seq_along(spans)) {
    sp <- spans[[r]]; start <- 1
    for (s in seq_along(sp)) {
      if (nzchar(trimws(names(sp)[s])))
        ft <- flextable::hline(ft, i = r, j = start:(start + sp[s] - 1),
                               border = officer::fp_border(color = "black", width = .75), part = "header")
      start <- start + sp[s]
    }
  }
  if (length(group_rows)) {
    ft <- flextable::padding(ft, i = setdiff(seq_len(nrow(body)), group_rows), j = 1,
                             padding.left = 14, part = "body")
    for (g in group_rows) {
      ft <- flextable::merge_at(ft, i = g, j = 1:nc, part = "body")
      ft <- flextable::italic(ft, i = g, j = 1, part = "body")
      ft <- flextable::align(ft, i = g, j = 1, align = "left", part = "body")
    }
  }
  if (!is.null(note)) {
    ft <- flextable::add_footer_lines(ft, values = flextable::as_paragraph(flextable::as_i("Note. "), note))
    ft <- flextable::font(ft, fontname = font, part = "footer")
    ft <- flextable::fontsize(ft, size = size, part = "footer")
    ft <- flextable::align(ft, align = "left", part = "footer")
  }
  ft <- flextable::width(ft, j = 1, width = first_col_width)
  if (nc > 1) ft <- flextable::width(ft, j = 2:nc, width = col_width)
  # Compact rows: height at least 0.16 in (the text sets the actual height)
  ft <- flextable::height_all(ft, height = 0.16, part = "all")
  ft <- flextable::hrule(ft, rule = "atleast", part = "all")
  flextable::set_table_properties(ft, layout = "fixed")
}

title_par_xml <- function(text, bold = FALSE, italic = FALSE, page_break = FALSE,
                          font = APA_FONT, size_pt = 12) {
  txt <- gsub("&", "&amp;", text, fixed = TRUE)
  txt <- gsub("<", "&lt;", txt, fixed = TRUE); txt <- gsub(">", "&gt;", txt, fixed = TRUE)
  rpr <- sprintf(paste0('<w:rPr><w:rFonts w:ascii="%1$s" w:hAnsi="%1$s" w:cs="%1$s"/>',
                        '%2$s%3$s<w:sz w:val="%4$d"/><w:szCs w:val="%4$d"/></w:rPr>'),
                 font, if (bold) "<w:b/>" else "", if (italic) "<w:i/>" else "", as.integer(size_pt * 2))
  ppr <- sprintf('<w:pPr><w:keepNext/>%s<w:spacing w:before="0" w:after="120"/></w:pPr>',
                 if (page_break) "<w:pageBreakBefore/>" else "")
  sprintf(paste0('<w:p xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">',
                 '%s<w:r>%s<w:t xml:space="preserve">%s</w:t></w:r></w:p>'), ppr, rpr, txt)
}

add_apa_table <- function(doc, number, title, ft, new_page = TRUE) {
  doc <- officer::body_add_xml(doc, title_par_xml(paste("Table", number), bold = TRUE, page_break = new_page))
  doc <- officer::body_add_xml(doc, title_par_xml(title, italic = TRUE))
  w <- dim(ft)$widths
  if (sum(w) > APA_MAX_WIDTH) ft <- flextable::width(ft, j = seq_along(w), width = w * APA_MAX_WIDTH / sum(w))
  flextable::body_add_flextable(doc, ft, align = "left")
}

grouped_body <- function(df, group, label_col = 1) {
  out <- NULL; g <- integer(0)
  for (lv in unique(group)) {
    hr <- df[1, , drop = FALSE]; hr[] <- ""; hr[[label_col]] <- as.character(lv)
    out <- rbind(out, hr, df[group == lv, , drop = FALSE]); g <- c(g, nrow(out) - sum(group == lv))
  }
  rownames(out) <- NULL
  list(body = out, group_rows = g)
}

qnew_landscape_doc <- function() {
  doc <- officer::read_docx()
  officer::body_set_default_section(doc, officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1)))
}
