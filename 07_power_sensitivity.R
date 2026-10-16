
## title: "APIM meta-analysis: individual-level missing and final pipeline"
## author: "Anonymized"
## date: "15-10-2026"

# =============================================================================
# 07_power_sensitivity.R
# Sensitivity of the power analysis (APIM_sensitivity_report.docx) to the
# item-level missing data. This part relies on the output of
# apim_meta_sensitivity.R and on the script itself, restricted to the rma()
# models chosen on the basis of the report (random effects, REML, Wald z as
# the primary test; common effect and Knapp-Hartung kept for reference).
#
#   A  Reproduction. apim_meta_sensitivity.R is rerun (rma() methods only) on
#      its original input (../sensitivity/APIM_meta_analysis_results.rds) with
#      the settings stored in ../sensitivity/APIM_sensitivity_results.rds (seed,
#      R, heterogeneity levels, effect grids). Its minimum detectable effects
#      (MDE), Type I error rates and critical values must equal the report's.
#   B  Exact power under missing data. The same rerun on the study-level results
#      of the focal missing-data methods (default: mice-PMM and two-stage FIML)
#      in the 15% conditions (POWER_CONDITIONS=all: all 9). Their sampling
#      covariance matrices include
#      the information lost to the missing items (Rubin's total variance for
#      mice-PMM, the two-stage FIML covariance matrix), so the simulation of the
#      report now gives the MDE (80% power) at each heterogeneity level, the
#      Type I error and the MDE of the indirect effects WITH the missing data.
#   C  Rescaled results for every method. The report's MDEs, power curves and
#      indirect-effect MDEs scaled by the SE inflation c of each method and
#      condition relative to the report's input (rescale_power() in
#      00_functions.R; exact when all sampling variances are inflated by the
#      same factor, and B checks it).
#   D  Critical values for the observed (incomplete) toy data: the smallest
#      pooled |beta| that is significant at alpha = .05 under each method and
#      condition (as in Table 2 of the report), and whether the significance
#      decision agrees with the complete data (COM).
#   Note: the report analysed the toy files as attached, i.e. with their ~1.5%
#   MCAR item-level missing values (proration + FIML); COM is the same data
#   with those values restored (01). The two differ slightly; the report's
#   numbers are reproduced from the report's own input (A and D).
#
# Output (output/power/)
#   inputs/<method>_<condition>.rds  study-level inputs of the reruns (format of
#                                    APIM_meta_analysis_results.rds)
#   raw/<method>_<condition>.rds     reruns of apim_meta_sensitivity.R (cached)
#   07_power_results.rds, 07_power_sensitivity.xlsx, 07_mde_by_condition.png/.pdf
#
# Runtime: about 4-6 min per rerun on one core (R = 1000, rma() only): 7 reruns
# by default, 19 with POWER_CONDITIONS=all; WORKERS=n runs n reruns in parallel
# (each rerun itself uses one core). Environment variables: POWER_METHODS (default
# "mice_pmm,FIML_TS"), POWER_CONDITIONS (default the three 15% conditions, or
# "all", or a comma-separated list), POWER_R (default: the R of the report),
# WORKERS.
# =============================================================================

source("00_functions.R")

## 0. Settings and inputs

OUT <- file.path("output", "power"); RAW <- file.path(OUT, "raw"); INP <- file.path(OUT, "inputs")

dir.create(RAW, FALSE, TRUE); dir.create(INP, FALSE, TRUE)

WORKERS <- as.integer(Sys.getenv("WORKERS", "10"))

FOCAL   <- strsplit(Sys.getenv("POWER_METHODS", "mice_pmm,FIML_TS"), ",")[[1]]

PCOND   <- Sys.getenv("POWER_CONDITIONS", "MCAR_15pct,MAR_15pct,MNAR_15pct")

f_res <- find_sensitivity_file("APIM_sensitivity_results.rds")

f_in  <- find_sensitivity_file("APIM_meta_analysis_results.rds")

if (is.null(f_res) || is.null(f_in))
  stop("APIM_sensitivity_results.rds and APIM_meta_analysis_results.rds (outputs of the sensitivity ",
       "analysis) not found in ", paste(SENS_DIRS, collapse = ", "))

REP <- readRDS(f_res)                                  # what the report was built from

E   <- readRDS(file.path("output", "evaluation", "04_evaluation_results.rds"))

FOCAL <- intersect(FOCAL, E$methods)

conditions <- E$conditions; cond_info <- E$cond_info

B_CONDS <- if (PCOND == "all") conditions else intersect(strsplit(PCOND, ",")[[1]], conditions)

# Settings of the original run, passed through CFG_OVERRIDE (rds_file is set per run)

OVR <- REP$config[setdiff(names(REP$config), c("rds_file", "input_mode", "n_cores", "validate_mc", "out_prefix"))]

OVR <- c(OVR, list(input_mode = "rds", validate_mc = FALSE, n_cores = 1L))   # parallel over reruns instead

OVR$R <- as.integer(Sys.getenv("POWER_R", REP$config$R))

same_R <- OVR$R == REP$config$R

K1 <- as.character(REP$summary$k[1])                   # first k (the observed studies)

RMA_LAB <- RMA_METHODS[unname(RMA_CODES)]              # "rma(): CE (z)", ... (labels of the report)

names(RMA_LAB) <- unname(RMA_CODES)

lab2code <- setNames(names(RMA_LAB), RMA_LAB)

comp2code <- setNames(names(COMP_LAB), COMP_LAB); eff2code <- setNames(names(EFFECT_LAB), EFFECT_LAB)

as_codes <- function(d) {                              # report labels -> codes used in this pipeline
  d$method_code <- unname(lab2code[as.character(d$method)])
  d$component_code <- unname(comp2code[as.character(d$component)])
  if ("effect" %in% names(d)) d$effect_code <- unname(eff2code[as.character(d$effect)])
  d
}

cat(sprintf("Power sensitivity: R = %d (report: %d), focal methods: %s, workers = %d\n",
            OVR$R, REP$config$R, paste(METHOD_SHORT[FOCAL], collapse = ", "), WORKERS))

## A. Reproduction of the report (rma() methods)

f_A <- file.path(RAW, sprintf("COM_complete_R%d.rds", OVR$R))

key_A <- cache_key("power", OVR, unname(tools::md5sum(f_in)))   # settings and the report's input file

RA <- read_cached(f_A, key_A)

if (is.null(RA)) {
  RA <- run_sensitivity_script(f_in, OVR); RA$cache_key <- key_A; saveRDS(RA, f_A)
}

cat(sprintf("A. Rerun on the complete-data input: %.0f s\n", RA$seconds))

keyc <- c("k", "component", "het", "method", "effect")

rep_s <- REP$summary[as.character(REP$summary$method) %in% RMA_LAB, ]

mA <- merge(transform(rep_s, k = as.character(k), component = as.character(component), het = as.character(het),
                      method = as.character(method), effect = as.character(effect)),
            transform(RA$summary, k = as.character(k), component = as.character(component), het = as.character(het),
                      method = as.character(method), effect = as.character(effect)),
            by = keyc, suffixes = c("_report", "_rerun"))

same <- function(a, b) ifelse(is.na(a) & is.na(b), 0, ifelse(is.na(a) | is.na(b), Inf, abs(a - b)))

REPRO <- data.frame(n_compared = nrow(mA),
                    max_abs_diff_mde = max(same(mA$mde_report, mA$mde_rerun)),
                    max_abs_diff_mde50 = max(same(mA$mde50_report, mA$mde50_rerun)),
                    max_abs_diff_type1 = max(same(mA$type1_report, mA$type1_rerun)),
                    max_abs_diff_crit_med = max(same(mA$crit_med_report, mA$crit_med_rerun)),
                    max_abs_diff_obs_crit = {
                      ob <- merge(REP$observed, RA$observed, by = c("component", "method", "effect"))
                      max(abs(ob$crit_est.x - ob$crit_est.y)) },
                    same_R_as_report = same_R)

cat(sprintf("   Reproduction of the report (rma(), %d rows): max |diff| MDE = %.2g, Type I = %.2g, critical |beta| = %.2g\n",
            REPRO$n_compared, REPRO$max_abs_diff_mde, REPRO$max_abs_diff_type1, REPRO$max_abs_diff_obs_crit))

if (same_R && max(unlist(REPRO[2:6])) > 1e-8)
  warning("The rerun does not reproduce the report exactly (different package versions can change REML ",
          "iterations slightly); see sheet 'A_reproduction'.")

## B. Exact reruns with the missing-data sampling covariances

study_set <- function(cond, method) setNames(lapply(STUDIES$key, function(k)
  readRDS(file.path("output", "evaluation", "raw", sprintf("%s_%s__%s.rds", k, cond, method)))), STUDIES$study)

runs <- expand.grid(condition = B_CONDS, method = FOCAL, stringsAsFactors = FALSE)

out_file <- function(i) file.path(RAW, sprintf("%s_%s_R%d.rds", runs$method[i], runs$condition[i], OVR$R))
# Cache key: the settings and the study-level estimates and covariance matrices (04)

runs$key <- vapply(seq_len(nrow(runs)), function(i) cache_key("power", OVR,
  lapply(study_set(runs$condition[i], runs$method[i]), function(r) lapply(r[COMPONENTS], function(x)
    list(est = x$est, V = x$V[[PRIMARY_V]])))), "")

b_job <- function(i) {
  if (!exists("run_sensitivity_script", envir = globalenv())) sys.source("00_functions.R", envir = globalenv())
  cd <- runs$condition[i]; mt <- runs$method[i]
  f_inp <- file.path(INP, sprintf("%s_%s.rds", mt, cd))
  write_sensitivity_input(study_set(cd, mt), f_inp, note = sprintf("%s, toy data %s", METHOD_LABELS[[mt]], cd))
  r <- run_sensitivity_script(f_inp, OVR)
  r$cache_key <- runs$key[i]
  saveRDS(r, out_file(i))
  sprintf("B. %-9s %-11s %4.0f s", METHOD_SHORT[[mt]], cd, r$seconds)
}

todo <- which(vapply(seq_len(nrow(runs)), function(i) is.null(read_cached(out_file(i), runs$key[i])), TRUE))

if (length(todo)) cat(unlist(par_lapply(todo, b_job, workers = WORKERS)), sep = "\n")

lab_run <- function(d, cd, mt) cbind(condition = cd, missing_method = mt, d, row.names = NULL)

RB <- lapply(seq_len(nrow(runs)), function(i) {
  r <- read_cached(out_file(i), runs$key[i])
  if (is.null(r)) stop("Rerun missing or stale: ", out_file(i))
  r })

EXACT <- as_codes(do.call(rbind, c(list(lab_run(RA$summary, "complete", "COM")),
  lapply(seq_along(RB), function(i) lab_run(RB[[i]]$summary, runs$condition[i], runs$method[i])))))

EXACT_IND <- as_codes(do.call(rbind, c(list(lab_run(RA$indirect_min, "complete", "COM")),
  lapply(seq_along(RB), function(i) lab_run(RB[[i]]$indirect_min, runs$condition[i], runs$method[i])))))

EXACT <- EXACT[EXACT$k == K1, ]; EXACT_IND <- EXACT_IND[EXACT_IND$k == K1, ]

## ---- C. Rescaled report results for every method -----------------------------------------
# SE inflation c = SE of the common-effect estimate under missing data / complete
# data, per condition, method, component and effect (04; V = ML)

ce <- E$direct[E$direct$method == "CE_z" & E$direct$V_type == PRIMARY_V, ]

ce_com <- ce[ce$missing_method == "COM", c("component", "effect", "se")]; names(ce_com)[3] <- "se_com"

INFL <- merge(ce[ce$missing_method != "COM", c("condition", "mechanism", "rate", "missing_method", "component",
                                               "effect", "se")], ce_com, by = c("component", "effect"))

INFL$c_com <- INFL$se / INFL$se_com
# ... and relative to the input of the report (the V that its MDEs are based on)

rr_in  <- readRDS(f_in)

rep_ma <- meta_both(lapply(rr_in$study_results, identity), PRIMARY_V, n_mc = 1000)$direct

se_rep <- rep_ma[rep_ma$method == "CE_z", c("component", "effect", "se")]; names(se_rep)[3] <- "se_report"

INFL <- merge(INFL, se_rep, by = c("component", "effect"))

INFL$c <- INFL$se / INFL$se_report                    # used for the rescaling

INFL <- INFL[order(match(INFL$condition, conditions), match(INFL$missing_method, METHODS_ALL),
                   match(INFL$component, COMPONENTS), match(INFL$effect, EFFECTS)), ]

rownames(INFL) <- NULL

rs <- as_codes(rep_s[as.character(rep_s$k) == K1, ])

RESC <- merge(INFL[, c("condition", "mechanism", "rate", "missing_method", "component", "effect", "c")],
              rs[, c("component_code", "effect_code", "het", "method_code", "mde", "mde50", "crit_med", "type1")],
              by.x = c("component", "effect"), by.y = c("component_code", "effect_code"))

RESC$mde_rescaled <- RESC$c * RESC$mde              # NA (not reached in the grid) stays NA

RESC$mde50_rescaled <- RESC$c * RESC$mde50

RESC$crit_med_rescaled <- RESC$c * RESC$crit_med

names(RESC)[names(RESC) %in% c("mde", "mde50", "crit_med", "type1")] <- c("mde_report", "mde50_report",
                                                                          "crit_med_report", "type1_report")

# Power to detect a true pooled effect equal to the complete-data estimate (RE z)

pw <- as_codes(REP$power[as.character(REP$power$k) == K1 & as.character(REP$power$method) == RMA_LAB[["RE_z"]], ])

com_est <- E$direct[E$direct$missing_method == "COM" & E$direct$method == "RE_z" & E$direct$V_type == PRIMARY_V,
                    c("component", "effect", "est")]

POW_AT_COM <- do.call(rbind, lapply(split(INFL, seq_len(nrow(INFL))), function(x) {
  b <- abs(com_est$est[com_est$component == x$component & com_est$effect == x$effect])
  do.call(rbind, lapply(unique(as.character(pw$het)), function(h) {
    p <- pw[pw$component_code == x$component & pw$effect_code == x$effect & as.character(pw$het) == h, ]
    if (nrow(p) < 2) stop("Power curve not found for ", x$component, " / ", x$effect, " / ", h)
    data.frame(x[, c("condition", "mechanism", "rate", "missing_method", "component", "effect", "c")], het = h,
               beta_com = b, power_complete = stats::approx(p$theta, p$power, xout = b, rule = 2, ties = max)$y,
               power_missing = rescale_power(p$theta, p$power, x$c, b), row.names = NULL)
  }))
}))

# Indirect effects: MDE of ab (joint significance) scaled by c_a * c_b

ri <- as_codes(REP$indirect_min[as.character(REP$indirect_min$k) == K1 &
                                  as.character(REP$indirect_min$method) %in% RMA_LAB, ])

groups <- unique(INFL[, c("condition", "missing_method", "component")])

RESC_IND <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
  g <- groups[i, ]
  sel <- INFL$condition == g$condition & INFL$missing_method == g$missing_method & INFL$component == g$component
  cc <- setNames(INFL$c[sel], INFL$effect[sel])
  x <- ri[ri$component_code == g$component, ]
  ca <- unname(cc[vapply(as.character(x$type), function(tp) SPECIFIC_INDIRECT[[tp]][1], "")])
  cb <- unname(cc[vapply(as.character(x$type), function(tp) SPECIFIC_INDIRECT[[tp]][2], "")])
  data.frame(condition = g$condition, missing_method = g$missing_method, component = g$component,
             het = x$het, method_code = x$method_code, type = x$type, c_a = ca, c_b = cb,
             min_ab_report = x$min_ab, min_ab_rescaled = x$min_ab * ca * cb, row.names = NULL)
}))

## D. Critical values for the observed incomplete data

dd <- E$direct[E$direct$V_type == PRIMARY_V, ]

dd$crit_value <- ifelse(dd$method %in% c("RE_KH", "RE_KHadh"), stats::qt(.975, pmax(dd$k - 1, 1)), stats::qnorm(.975))

dd$crit_est <- dd$crit_value * dd$se

dd$significant <- abs(dd$est) > dd$crit_est

com_sig <- dd[dd$missing_method == "COM", c("component", "method", "effect", "est", "crit_est", "significant")]

names(com_sig)[4:6] <- c("est_com", "crit_est_com", "significant_com")

CRIT <- merge(dd[dd$missing_method != "COM", c("condition", "mechanism", "rate", "missing_method", "component",
                                              "method", "effect", "est", "se", "crit_est", "significant")],
              com_sig, by = c("component", "method", "effect"))

CRIT$crit_inflation <- CRIT$crit_est / CRIT$crit_est_com

CRIT$same_decision <- CRIT$significant == CRIT$significant_com

CRIT <- CRIT[order(match(CRIT$condition, conditions), match(CRIT$missing_method, METHODS_ALL),
                   match(CRIT$component, COMPONENTS), match(CRIT$method, names(RMA_METHODS)),
                   match(CRIT$effect, EFFECTS)), ]

rownames(CRIT) <- NULL

# Validation: the same computation on the report's input reproduces Table 2 of the report

rep_ma$crit_est <- ifelse(rep_ma$method %in% c("RE_KH", "RE_KHadh"), stats::qt(.975, pmax(rep_ma$k - 1, 1)),
                          stats::qnorm(.975)) * rep_ma$se

chk_D <- merge(as_codes(REP$observed)[, c("component_code", "method_code", "effect_code", "crit_est")],
               rep_ma[, c("component", "method", "effect", "crit_est")],
               by.x = c("component_code", "method_code", "effect_code"), by.y = c("component", "method", "effect"))

REPRO$max_abs_diff_table2 <- max(abs(chk_D$crit_est.x - chk_D$crit_est.y))

cat(sprintf("D. Critical |beta| computed here vs Table 2 of the report (report input): max |diff| = %.1e\n",
            REPRO$max_abs_diff_table2))

##E. Exact vs rescaled (focal methods)

cmpE <- merge(EXACT[EXACT$missing_method != "COM", c("condition", "missing_method", "component_code", "effect_code",
                                                     "het", "method_code", "mde", "type1")],
              RESC[, c("condition", "missing_method", "component", "effect", "het", "method_code", "c",
                       "mde_report", "mde_rescaled")],
              by.x = c("condition", "missing_method", "component_code", "effect_code", "het", "method_code"),
              by.y = c("condition", "missing_method", "component", "effect", "het", "method_code"))

names(cmpE)[names(cmpE) == "mde"] <- "mde_exact"

cmpE$ratio_exact_rescaled <- cmpE$mde_exact / cmpE$mde_rescaled

cmpE$mde_inflation_exact <- cmpE$mde_exact / cmpE$mde_report

APPROX <- aggregate(cbind(abs_log_ratio = abs(log(ratio_exact_rescaled))) ~ missing_method + method_code + het,
                    cmpE[is.finite(cmpE$ratio_exact_rescaled), ], function(v) exp(mean(v)) - 1)

names(APPROX)[4] <- "mean_relative_difference"

cat("E. Exact vs rescaled MDE: mean relative difference by method (RE z)\n")

print(APPROX[APPROX$method_code == "RE_z", ], digits = 3, row.names = FALSE)

## F. Summary table: MDE (RE z) under missing data

MDE_TAB <- EXACT[EXACT$method_code == PRIMARY_RMA, c("condition", "missing_method", "component_code", "effect_code",
                                                    "het", "mde", "mde50", "type1", "crit_med", "fail")]

MDE_TAB <- merge(MDE_TAB, setNames(rs[rs$method_code == PRIMARY_RMA, c("component_code", "effect_code", "het", "mde")],
                                   c("component_code", "effect_code", "het", "mde_report")),
                 by = c("component_code", "effect_code", "het"))

MDE_TAB$mde_inflation <- MDE_TAB$mde / MDE_TAB$mde_report

MDE_TAB <- MDE_TAB[order(match(MDE_TAB$condition, c("complete", conditions)),
                         match(MDE_TAB$missing_method, METHODS_ALL), match(MDE_TAB$component_code, COMPONENTS),
                         MDE_TAB$het, match(MDE_TAB$effect_code, EFFECTS)), ]

rownames(MDE_TAB) <- NULL

## G. Figure: MDE by condition (RE z, I^2 = 0% and 50%; focal methods)
# Lines: rescaled MDE (all conditions); open points: exact reruns (B); dotted:
# the report's complete-data MDE

if (requireNamespace("ggplot2", quietly = TRUE)) {
  hets <- intersect(c("I\u{00b2} = 0%", "I\u{00b2} = 50%"), levels(REP$summary$het))
  prep <- function(d) {
    d$effect <- factor(EFFECT_LAB[d$effect], levels = EFFECT_LAB)
    d$component <- factor(COMP_LAB[d$component], levels = COMP_LAB)
    d$het <- factor(as.character(d$het), levels = hets)
    d
  }
  fl <- RESC[RESC$missing_method %in% FOCAL & RESC$method_code == PRIMARY_RMA & as.character(RESC$het) %in% hets, ]
  fl <- prep(fl); fl$Method <- factor(METHOD_SHORT[fl$missing_method], levels = METHOD_SHORT[FOCAL])
  fl$mechanism <- factor(fl$mechanism, levels = MECHANISMS)
  fe <- cmpE[cmpE$method_code == PRIMARY_RMA & as.character(cmpE$het) %in% hets, ]
  fe <- merge(fe, cond_info, by = "condition")
  names(fe)[names(fe) == "component_code"] <- "component"; names(fe)[names(fe) == "effect_code"] <- "effect"
  fe <- prep(fe); fe$Method <- factor(METHOD_SHORT[fe$missing_method], levels = METHOD_SHORT[FOCAL])
  fe$mechanism <- factor(fe$mechanism, levels = MECHANISMS)
  ref <- unique(fl[, c("component", "effect", "het", "mde_report")])
  p <- ggplot2::ggplot(fl[!is.na(fl$mde_rescaled), ],
                       ggplot2::aes(100 * rate, mde_rescaled, colour = mechanism, linetype = Method)) +
    ggplot2::geom_hline(data = ref[!is.na(ref$mde_report), ], ggplot2::aes(yintercept = mde_report),
                        colour = "grey45", linetype = "dotted") +
    ggplot2::geom_line(linewidth = .45) +
    ggplot2::geom_point(data = fe[!is.na(fe$mde_exact), ], ggplot2::aes(100 * rate, mde_exact, shape = Method),
                        size = 1.6, inherit.aes = FALSE) +
    ggplot2::facet_grid(component + het ~ effect, labeller = ggplot2::labeller(component = ggplot2::label_wrap_gen(14))) +
    ggplot2::scale_colour_manual(values = c(MCAR = "#0072B2", MAR = "#E69F00", MNAR = "#D55E00")) +
    ggplot2::scale_shape_manual(values = c(1, 2)) +
    ggplot2::scale_x_continuous(breaks = c(5, 10, 15)) +
    ggplot2::labs(x = "Missing item values (%)", y = "Minimum detectable pooled effect (80% power)",
                  colour = NULL, linetype = "Rescaled", shape = "Exact rerun") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70"))
  ggplot2::ggsave(file.path(OUT, "07_mde_by_condition.png"), p, width = 10, height = 8, dpi = 300)
  if (capabilities("cairo")) ggplot2::ggsave(file.path(OUT, "07_mde_by_condition.pdf"), p, width = 10, height = 8,
                                             device = grDevices::cairo_pdf)
}

## H. Export

P <- list(R = OVR$R, R_report = REP$config$R, focal = FOCAL, exact_conditions = B_CONDS, k1 = K1,
          het_levels = levels(REP$summary$het),
          reproduction = REPRO, exact = EXACT, exact_indirect = EXACT_IND, inflation = INFL, rescaled = RESC,
          rescaled_indirect = RESC_IND, power_at_complete_estimate = POW_AT_COM, critical = CRIT,
          exact_vs_rescaled = cmpE, approx_summary = APPROX, mde_table = MDE_TAB,
          report_summary = rs, report_observed = as_codes(REP$observed), report_config = REP$config,
          seconds = c(A = RA$seconds, B = vapply(RB, `[[`, 0, "seconds")))

saveRDS(P, file.path(OUT, "07_power_results.rds"))

write_xlsx_sheets(list(
  settings = data.frame(setting = c("replications (R)", "R of the report", "focal methods", "exact reruns for",
                                    "seed", "heterogeneity", "power target", "alpha", "V type", "rma() methods"),
                        value = c(OVR$R, REP$config$R, paste(FOCAL, collapse = ", "), paste(B_CONDS, collapse = ", "),
                                  REP$config$seed,
                                  paste(levels(REP$summary$het), collapse = ", "), REP$config$power_target,
                                  REP$config$alpha, REP$config$V_type, paste(RMA_LAB, collapse = "; "))),
  A_reproduction = REPRO, A_rerun_vs_report = mA,
  B_mde_RE_z = MDE_TAB, B_exact_all_rma = EXACT, B_exact_indirect = EXACT_IND,
  C_se_inflation = INFL, C_rescaled = RESC, C_rescaled_indirect = RESC_IND,
  C_power_at_complete_est = POW_AT_COM, D_critical_values = CRIT,
  E_exact_vs_rescaled = cmpE, E_summary = APPROX), file.path(OUT, "07_power_sensitivity.xlsx"))

cat("Written:", file.path(OUT, c("07_power_results.rds", "07_power_sensitivity.xlsx", "07_mde_by_condition.png")),
    sep = "\n")

x <- MDE_TAB[MDE_TAB$het == levels(REP$summary$het)[1], ]

x$missing_method <- METHOD_SHORT[x$missing_method]

cat(sprintf("\nMDE (rma() RE z, %s, %s), mean over the six direct effects:\n", levels(REP$summary$het)[1], K1))

print(stats::aggregate(mde ~ condition + missing_method + component_code, x, mean), digits = 3, row.names = FALSE)
