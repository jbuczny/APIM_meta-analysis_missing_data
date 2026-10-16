# =============================================================================
# 08_apa_tables.R
# APA-style tables (Word) giving an overview of the missing-data sensitivity
# results of 02-07.
#
#   Input : data/amputed/amputed_datasets.rds                     (02)
#           output/imputations/imputation_log.xlsx                 (03)
#           output/evaluation/04_evaluation_results.rds            (04)
#           output/mc/<full|pilot>/mc_summary.rds                  (05; optional)
#           output/sensitivity/<full|quick>/06_delta_results.rds   (06; optional)
#           output/power/07_power_results.rds                      (07; optional)
#   Output: output/report/08_missing_data_sensitivity_tables.docx
#
#   A  Missing data and imputation
#      Table 1   Item-level missing data in the 18 incomplete datasets
#      Table 2   Quality of the imputed item values
#   B  Missing-data methods in the meta-analysis
#      Table 3   Method decision: mice-PMM versus two-stage FIML
#      Table 4   Pooled effects relative to the complete data: uncorrected paths
#      Table 5   The same, corrected for unreliability
#      Table 6   Pooled direct effects by condition: mice-PMM
#      Table 7   Pooled direct effects by condition: two-stage FIML
#      Table 8   Pooled specific indirect effects by condition
#   C  Power sensitivity: the sensitivity report under missing data (07)
#      Table 9   Critical values for significance (cf. Table 2 of the report)
#      Table 10  Minimum detectable direct effects (cf. Table 1 of the report)
#      Table 11  Minimum detectable indirect effects (cf. Table 3 of the report)
#   D  Departures from MAR: delta adjustment (06)
#      Table 12  Tipping points of the pooled effects
#      Table 13  Pooled direct effects by delta, MNAR 15%
#   E  Monte Carlo simulation (05, pilot or full mode only)
#      Tables 14-15  Performance of the methods (uncorrected; corrected)
#   Supplementary
#      Table S1  Reliability (omega) by missing-data method
#      Table S2  Fraction of missing information
#      Tables S3-S4  All methods side by side at 15% missing items
#      Table S5  Exact versus rescaled minimum detectable effects
#
# Tables whose input is missing are skipped (with a message) and the others
# are renumbered. Formatting follows apim_meta_analysis.R (APA 7: Times New
# Roman; rules only above and below the header and below the body; italic
# "Note."; no leading zero for values bounded by 1; landscape pages).
# =============================================================================

source("00_functions.R")

## ---- 0. Inputs -----------------------------------------------------------------------------
OUT <- file.path("output", "report"); dir.create(OUT, FALSE, TRUE)
OUT_DOCX <- file.path(OUT, "08_missing_data_sensitivity_tables.docx")
A <- readRDS(file.path("data", "amputed", "amputed_datasets.rds"))
E <- readRDS(file.path("output", "evaluation", "04_evaluation_results.rds"))
f_log <- file.path("output", "imputations", "imputation_log.xlsx")
LOG <- if (file.exists(f_log)) openxlsx::read.xlsx(f_log, sheet = "log") else NULL
pick <- function(dir, file, modes) {
  for (md in modes) { f <- file.path("output", dir, md, file); if (file.exists(f)) return(f) }
  NULL
}
f_mc <- pick("mc", "mc_summary.rds", c("full", "pilot"))        # quick mode = smoke test, not tabled
MC  <- if (!is.null(f_mc)) readRDS(f_mc) else NULL
f_se <- pick("sensitivity", "06_delta_results.rds", c("full", "quick"))
S06 <- if (!is.null(f_se)) readRDS(f_se) else NULL
f_pw <- file.path("output", "power", "07_power_results.rds")
P   <- if (file.exists(f_pw)) readRDS(f_pw) else NULL
DEC <- choose_method(E, MC)
chosen <- DEC$chosen

grid <- A$grid
conds <- E$conditions
cond_info <- E$cond_info[match(conds, E$cond_info$condition), ]
METH <- intersect(METHODS_ALL, setdiff(E$methods, "COM"))
focal <- intersect(c("mice_pmm", "FIML_TS"), METH)
spec <- names(SPECIFIC_INDIRECT)
rate_lab <- function(r) sprintf("%d%%", round(100 * r))
cond_spans <- function(first = c(" " = 1L)) c(first, MCAR = 3L, MAR = 3L, MNAR = 3L)
cond_labels <- rep(c("5%", "10%", "15%"), 3)
star <- function(est, p, d = 2) paste0(fmt(est, d), ifelse(!is.na(p) & p < .05, "*", ""))
mde_fmt <- function(v, d = 2, max_lab = "> 1.00") ifelse(is.na(v), max_lab, fmt(v, d))
mech_txt <- paste0("MCAR = missing completely at random; MAR = missingness depends on the always-observed ",
                   "auxiliaries; MNAR = missingness depends on the item value itself (higher X and M, lower Y) ",
                   "and on an instrument.")
rma_txt <- paste0("Pooled across the two studies with rma(): per study, the Partner 1 and Partner 2 estimates of ",
                  "each effect were combined (GLS, using their lavaan sampling covariance matrix), and the study ",
                  "composites were pooled with a random-effects model (REML, Wald z), the specification chosen ",
                  "on the basis of the sensitivity report.")
meth_abbr <- paste0("COM = complete data (the same toy data before deletion); PRO-FIML = proration (at least 75% ",
                    "of the items) and FIML on the scale scores, the default of apim_meta_analysis.R; FIML-TS = ",
                    "two-stage item-level FIML with the auxiliaries; PMM and CAT = item-level multiple imputation ",
                    "(predictive mean matching; proportional-odds model; m = 20) pooled within each study with ",
                    "Rubin's rules; PMM-avg = the 20 PMM datasets averaged into one dataset and analysed once ",
                    "(shown only to illustrate why this is not done).")
meth_short <- paste0("COM = complete data; PRO-FIML = proration and FIML on the scale scores; FIML-TS = two-stage ",
                     "item-level FIML with the auxiliaries; PMM, CAT = item-level multiple imputation (m = 20) with ",
                     "Rubin's rules; PMM-avg = the 20 PMM datasets averaged into one (illustration only).")
mech_short <- paste0("MCAR = missing completely at random; MAR = missingness depends on the auxiliaries; MNAR = it ",
                     "depends on the item value itself.")
prim <- function(d) d[d$method == PRIMARY_RMA & d$V_type == PRIMARY_V, ]
D_all <- prim(E$direct)
I_all <- E$indirect[E$indirect$V_type == PRIMARY_V, ]
tables <- list()
add <- function(title, ft, part = NULL, supp = FALSE)
  tables[[length(tables) + 1]] <<- list(title = title, ft = ft, part = part, supp = supp)
skipped <- character(0)

## =========================== A. Missing data and imputation ==============================
ms <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  s <- missing_summary(A$amputed[[i]]); g <- grid[i, ]
  val <- function(v) s$pct_missing[s$variable == v]
  data.frame(study = g$study, condition = g$condition, cells = val("ALL item cells"),
             dyads = val("Dyads with any missing item"),
             scales = val("Scale scores not computable by proration (< 75% items)"))
}))
sel <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  cbind(condition = grid$condition[i], mechanism_regression(A$complete[[grid$study[i]]], A$amputed[[i]]))))
sel <- aggregate(log_odds ~ condition + term, sel, mean)
t1 <- data.frame(cond = rate_lab(cond_info$rate), stringsAsFactors = FALSE)
for (k in STUDIES$key) {
  x <- ms[ms$study == k, ][match(conds, ms$condition[ms$study == k]), ]
  t1[[paste0(k, "_c")]] <- sprintf("%.1f", x$cells)
  t1[[paste0(k, "_d")]] <- sprintf("%.1f", x$dyads)
  t1[[paste0(k, "_s")]] <- sprintf("%.1f", x$scales)
}
for (tm in c("z_own", "z_aux", "z_iv"))
  t1[[tm]] <- fmt(sel$log_odds[sel$term == tm][match(conds, sel$condition[sel$term == tm])], 2)
g1 <- grouped_body(t1, cond_info$mechanism)
add("Item-Level Missing Data in the 18 Incomplete Toy Datasets", apa_table(g1$body,
  labels = c("Mechanism / rate", rep(c("Item cells", "Dyads", "Scale scores"), 2), "Own value", "Auxiliary", "Instrument"),
  spans = list(c(" " = 1L, setNames(c(3L, 3L), sprintf("%s (N = %d dyads)", STUDIES$study, STUDIES$n)),
                 "Realized selection model" = 3L)),
  group_rows = g1$group_rows,
  note = paste0("Item cells = percentage of the item values that are missing (exactly round(rate \u{00d7} N) values per ",
                "item: 6, 13 and 20 of 130; 9, 17 and 26 of 174). Dyads = percentage of dyads with at least one ",
                "missing item. Scale scores = percentage of the scale scores with fewer than 75% of their items ",
                "observed, which apim_meta_analysis.R sets to missing (proration otherwise). Realized selection model = ",
                "log-odds per SD from a logistic regression of the item missingness indicators on the item's own value ",
                "(signed in the hypothesized MNAR direction), the partner's auxiliary composite and the partner's ",
                "instrument (item fixed effects), averaged over the two studies. Within a mechanism the missing cells ",
                "are nested across rates. ", mech_txt),
  first_col_width = 1.3, col_width = 0.78), part = "A. Missing Data and Imputation")

if (!is.null(LOG)) {
  lg <- LOG[LOG$method %in% IMPUTATION_METHODS, ]
  lg$condition <- sprintf("%s_%02dpct", lg$mechanism, round(100 * lg$rate))
  agg <- aggregate(cbind(nrmse, std_bias_mnar_direction, mabias_cor, tvd) ~ condition + method, lg, mean)
  psr <- aggregate(max_psrf ~ condition + method, lg, max)
  agg <- merge(agg, psr, by = c("condition", "method"))
  meths <- intersect(IMPUTATION_METHODS, unique(agg$method))
  metr <- list("Normalized RMSE of the imputed values" = list(col = "nrmse", d = 2, lead0 = TRUE),
               "Standardized bias of the item means (MNAR direction)" = list(col = "std_bias_mnar_direction", d = 3, lead0 = TRUE),
               "Mean absolute bias of the correlations" = list(col = "mabias_cor", d = 3, lead0 = FALSE),
               "Largest potential scale reduction factor" = list(col = "max_psrf", d = 2, lead0 = TRUE))
  rows <- NULL; grp <- NULL
  for (nm in names(metr)) for (mt in meths) {
    x <- agg[agg$method == mt, ][match(conds, agg$condition[agg$method == mt]), ]
    rows <- rbind(rows, data.frame(m = METHOD_SHORT[[mt]],
                                   t(setNames(fmt(x[[metr[[nm]]$col]], metr[[nm]]$d, metr[[nm]]$lead0), conds)),
                                   check.names = FALSE, stringsAsFactors = FALSE))
    grp <- c(grp, nm)
  }
  g2 <- grouped_body(rows, grp)
  add("Quality of the Imputed Item Values by Method and Missing-Data Condition", apa_table(g2$body,
    labels = c("Measure / method", cond_labels), spans = list(cond_spans()), group_rows = g2$group_rows,
    note = paste0("Item-level imputation (m = 20) of each study; values averaged over the 24 items, the imputations and ",
                  "the two studies (potential scale reduction factor: maximum). Normalized RMSE = root mean square ",
                  "difference between the imputed and the deleted values divided by the item SD (it rewards predicting ",
                  "the most likely value and is therefore not a criterion for multiple imputation). Standardized bias = ",
                  "(mean of the completed item \u{2212} mean of the complete item) / SD, signed so that negative values fall ",
                  "short of the MNAR direction. PMM = predictive mean matching; CAT = proportional-odds model with ",
                  "parameters drawn from their posterior. ", mech_txt),
    first_col_width = 2.6, col_width = 0.68))
} else skipped <- c(skipped, "imputation quality (imputation log of 03 not found)")

## ======================= B. Missing-data methods in the meta-analysis ======================
if (!is.null(DEC$criteria)) {
  cr <- DEC$criteria; is_mc <- !is.na(cr$coverage_mcar_mar[1])
  t3 <- data.frame(
    meas = c(sprintf("RMSD from the complete data, MCAR and MAR (%s; decision criterion)", DEC$units),
             sprintf("RMSD from the complete data, MNAR (%s)", DEC$units),
             if (is_mc) "Coverage of the 95% CI, MCAR and MAR",
             "95% CI width relative to the complete data, MCAR and MAR"),
    pmm = c(fmt(cr["mice_pmm", "rmsd_mcar_mar"], 3), fmt(cr["mice_pmm", "rmsd_mnar"], 3),
            if (is_mc) fmt(cr["mice_pmm", "coverage_mcar_mar"], 3, lead0 = FALSE), fmt(cr["mice_pmm", "ci_ratio_mcar_mar"], 2)),
    fiml = c(fmt(cr["FIML_TS", "rmsd_mcar_mar"], 3), fmt(cr["FIML_TS", "rmsd_mnar"], 3),
             if (is_mc) fmt(cr["FIML_TS", "coverage_mcar_mar"], 3, lead0 = FALSE), fmt(cr["FIML_TS", "ci_ratio_mcar_mar"], 2)),
    stringsAsFactors = FALSE)
  add("Choice of the Missing-Data Method: mice-PMM Versus Two-Stage FIML", apa_table(t3,
    labels = c("Criterion", "mice-PMM (m = 20)", "Two-stage FIML"),
    note = paste0(sprintf("Basis: %s. ", DEC$basis),
                  "RMSD = root mean square deviation of the six pooled direct effects (rma() random effects, z; ",
                  "both components) from the complete-data estimates, averaged over the conditions. Pre-specified ",
                  sprintf("rule: the method with the smaller RMSD under MCAR and MAR; a difference within %.0f%% of the ",
                          100 * REL_MARGIN),
                  "FIML value is a tie, and a tie goes to mice-PMM (its completed datasets are needed for the delta ",
                  sprintf("adjustment). Difference (PMM \u{2212} FIML) = %s, margin = %s: %s. Chosen: %s.",
                          fmt(DEC$difference, 3), fmt(DEC$margin, 3), DEC$verdict, METHOD_LABELS[[chosen]]),
                  if (chosen == "mice_pmm") paste0(" The m = 20 imputations of each study are aggregated with Rubin's rules (pooled ",
                    "estimates and their full sampling covariance matrix) before the meta-analysis; the imputed data are ",
                    "not averaged (PMM-avg in Tables 4 and 5).") else ""),
    first_col_width = 4.6, col_width = 1.6), part = "B. Missing-Data Methods in the Meta-Analysis")
}

key <- E$key[E$key$V_type == PRIMARY_V, ]
ci_ag <- aggregate(cbind(joint = joint_sig_agree) ~ condition + missing_method + component,
                   E$compare_indirect[E$compare_indirect$V_type == PRIMARY_V &
                                        E$compare_indirect$effect %in% spec, ], sum)
summary_table <- function(cmp) {
  panels <- list("Root mean square deviation from the complete data (SE units)" = function(x) fmt(x$rmsd_std, 2),
                 "SE inflation relative to the complete data" = function(x) fmt(x$mean_se_ratio_ce, 2),
                 "Same significance decision as the complete data: direct effects (of 6) / specific indirect effects (of 4)" =
                   function(x) sprintf("%d/%d", x$n_sig_agree, x$joint))
  rows <- NULL; grp <- NULL
  for (nm in names(panels)) for (mt in METH) {
    x <- key[key$component == cmp & key$missing_method == mt, ]
    x <- x[match(conds, x$condition), ]
    y <- ci_ag[ci_ag$component == cmp & ci_ag$missing_method == mt, ]
    x$joint <- as.integer(y$joint[match(conds, y$condition)])
    rows <- rbind(rows, data.frame(m = METHOD_SHORT[[mt]], t(setNames(panels[[nm]](x), conds)),
                                   check.names = FALSE, stringsAsFactors = FALSE))
    grp <- c(grp, nm)
  }
  g <- grouped_body(rows, grp)
  apa_table(g$body, labels = c("Measure / method", cond_labels), spans = list(cond_spans()),
    group_rows = g$group_rows,
    note = paste0("Pooled with rma() (random effects, REML, z). Root mean square deviation of the six pooled direct ",
                  "effects from the complete-data estimates, in units of the complete-data SE of the common-effect ",
                  "estimate (with two studies the random-effects SE depends on an unstable \u{03c4}\u{00b2} estimate). SE inflation = ",
                  "mean ratio of the common-effect SE to its complete-data value (values of about 1.00 or below for ",
                  "PMM-avg show that averaging the imputed datasets hides the information loss). Significance: z test ",
                  "at \u{03b1} = .05; specific indirect effects: joint significance. ", meth_short, " ", mech_short),
    first_col_width = 3.0, col_width = 0.62, size = 9)
}
add("Pooled Effects Relative to the Complete-Data Analysis: Uncorrected Paths", summary_table("uncorrected"))
add("Pooled Effects Relative to the Complete-Data Analysis: Paths Corrected for Unreliability",
    summary_table("corrected"))

direct_by_cond <- function(mt) {
  rows <- NULL; grp <- NULL
  for (cmp in COMPONENTS) {
    com <- D_all[D_all$missing_method == "COM" & D_all$component == cmp, ]
    blk <- data.frame(e = EFFECT_LAB[EFFECTS], stringsAsFactors = FALSE)
    blk$COM <- vapply(EFFECTS, function(e) star(com$est[com$effect == e], com$p[com$effect == e]), "")
    for (cd in conds) {
      x <- D_all[D_all$missing_method == mt & D_all$component == cmp & D_all$condition == cd, ]
      blk[[cd]] <- vapply(EFFECTS, function(e) star(x$est[x$effect == e], x$p[x$effect == e]), "")
    }
    rows <- rbind(rows, blk); grp <- c(grp, rep(COMP_LAB[[cmp]], length(EFFECTS)))
  }
  g <- grouped_body(rows, grp)
  apa_table(g$body, labels = c("Effect", "COM", cond_labels), spans = list(cond_spans(c(" " = 1L, " " = 1L))),
    group_rows = g$group_rows,
    note = paste0(rma_txt, " Entries are pooled standardized paths (\u{03b2}); * p < .05 (random-effects z test). COM = ",
                  "complete data (the same toy data before deletion). ",
                  if (mt == "FIML_TS") paste0("Two-stage FIML: saturated item-level FIML model of the 24 items and the ",
                    "auxiliaries (stage 1), then the APIMeM fitted in lavaan to the implied scale-score moments, with ",
                    "the stage-1 sampling covariance matrix (stage 2). ") else
                  paste0("mice-PMM: m = 20 item-level imputations, the APIMeM fitted in lavaan to every completed ",
                         "dataset, and the 20 analyses pooled within each study with Rubin's rules. "), mech_txt),
    first_col_width = 1.3, col_width = 0.72)
}
for (mt in focal) add(sprintf("Pooled Direct Effects by Missing-Data Condition: %s",
                              c(mice_pmm = "mice-PMM (m = 20, Rubin's Rules)", FIML_TS = "Two-Stage FIML")[[mt]]),
                      direct_by_cond(mt))

ind_sig <- function(x, d) if (!nrow(x)) "\u{2014}" else
  paste0(fmt(x$est, d), ifelse(isTRUE(x$joint_sig), "*", ""), ifelse(isTRUE(x$mc_sig), "\u{2020}", ""))
rows <- NULL; grp <- NULL
for (cmp in COMPONENTS) for (mt in focal) {
  com <- I_all[I_all$missing_method == "COM" & I_all$component == cmp, ]
  blk <- data.frame(e = spec, stringsAsFactors = FALSE)
  blk$COM <- vapply(spec, function(e) ind_sig(com[com$effect == e, ], 3), "")
  for (cd in conds) {
    x <- I_all[I_all$missing_method == mt & I_all$component == cmp & I_all$condition == cd, ]
    blk[[cd]] <- vapply(spec, function(e) ind_sig(x[x$effect == e, ], 3), "")
  }
  rows <- rbind(rows, blk); grp <- c(grp, rep(sprintf("%s: %s", COMP_LAB[[cmp]], METHOD_SHORT[[mt]]), length(spec)))
}
g8 <- grouped_body(rows, grp)
add("Pooled Specific Indirect Effects by Missing-Data Condition", apa_table(g8$body,
  labels = c("Indirect effect", "COM", cond_labels), spans = list(cond_spans(c(" " = 1L, " " = 1L))),
  group_rows = g8$group_rows,
  note = paste0("Products of the pooled random-effects paths (Ledermann et al., 2011): actor\u{2013}actor = actor a \u{00d7} ",
                "actor b; partner\u{2013}partner = partner a \u{00d7} partner b; actor\u{2013}partner = actor a \u{00d7} partner b; partner\u{2013}actor ",
                "= partner a \u{00d7} actor b. * joint significance (both pooled paths p < .05), the test evaluated in the ",
                "sensitivity report; \u{2020} 95% Monte Carlo CI excludes zero. COM = complete data; PMM = mice-PMM (m = 20, ",
                "Rubin's rules); FIML-TS = two-stage FIML. ", mech_short),
  first_col_width = 1.4, col_width = 0.72, size = 9))

## ================== C. Power sensitivity: the sensitivity report under missing data ==================
if (!is.null(P)) {
  pm <- if (chosen %in% P$inflation$missing_method) chosen else focal[1]
  pm_lab <- METHOD_LABELS[[pm]]
  appr <- P$approx_summary[P$approx_summary$method_code == PRIMARY_RMA, ]
  appr_txt <- if (nrow(appr)) sprintf(paste0("Exact reruns of apim_meta_sensitivity.R with the missing-data sampling ",
                                             "covariance matrices (%s; %s) differ from the rescaled values by %s%% on ",
                                             "average (Table S5)."),
                                      paste(METHOD_SHORT[unique(appr$missing_method)], collapse = " and "),
                                      paste(sub("_", " ", sub("pct", "%", P$exact_conditions)), collapse = ", "),
                                      paste(unique(sprintf("%.1f", 100 * range(appr$mean_relative_difference))),
                                            collapse = "\u{2013}")) else ""
  ## Table: critical values ----
  cr9 <- P$critical[P$critical$missing_method == pm & P$critical$method == PRIMARY_RMA, ]
  rows <- NULL; grp <- NULL
  for (cmp in COMPONENTS) {
    x0 <- cr9[cr9$component == cmp & cr9$condition == conds[1], ]; x0 <- x0[match(EFFECTS, x0$effect), ]
    blk <- data.frame(e = EFFECT_LAB[EFFECTS], b = fmt(x0$est_com),
                      crit = paste0(fmt(x0$crit_est_com), ifelse(x0$significant_com, "*", "")), stringsAsFactors = FALSE)
    for (cd in conds) {
      x <- cr9[cr9$component == cmp & cr9$condition == cd, ]; x <- x[match(EFFECTS, x$effect), ]
      blk[[cd]] <- paste0(fmt(x$crit_est), ifelse(x$significant, "*", ""), ifelse(x$same_decision, "", "\u{2020}"))
    }
    rows <- rbind(rows, blk); grp <- c(grp, rep(COMP_LAB[[cmp]], length(EFFECTS)))
  }
  g9 <- grouped_body(rows, grp)
  add("Critical Values for the Significance of the Pooled Direct Effects Under Item-Level Missing Data", flextable::width(apa_table(g9$body,
    labels = c("Effect", "\u{03b2}", "Critical |\u{03b2}|", cond_labels),
    spans = list(c(" " = 1L, "Complete data (COM)" = 2L, MCAR = 3L, MAR = 3L, MNAR = 3L)),
    group_rows = g9$group_rows,
    note = paste0("Critical |\u{03b2}| = smallest absolute pooled estimate that is significant at \u{03b1} = .05 with these data ",
                  "(rma() random effects, REML, z test: 1.96 \u{00d7} SE), as in Table 2 of the sensitivity report (which ",
                  "analysed the toy files with their 1.5% MCAR item values missing; COM restores them, so its values ",
                  sprintf("differ slightly). Missing-data columns: %s. * the pooled estimate exceeds the critical value ", pm_lab),
                  "(significant); \u{2020} the significance decision differs from the complete data. ", mech_txt),
    first_col_width = 1.2, col_width = 0.68), j = 2:3, width = 0.86),
    part = "C. Power Sensitivity: The Sensitivity Report Under Missing Data")

  ## Table: MDE of the direct effects (report values x inflation factor) ----
  rs <- P$report_summary[P$report_summary$method_code == PRIMARY_RMA, ]
  hl <- P$het_levels
  infl <- P$inflation[P$inflation$missing_method == pm, ]
  rows <- NULL; grp <- NULL
  for (cmp in COMPONENTS) {
    blk <- data.frame(e = EFFECT_LAB[EFFECTS], stringsAsFactors = FALSE)
    for (h in hl) {
      x <- rs[rs$component_code == cmp & as.character(rs$het) == h, ]
      blk[[h]] <- mde_fmt(x$mde[match(EFFECTS, x$effect_code)])
    }
    for (cd in conds) {
      x <- infl[infl$component == cmp & infl$condition == cd, ]
      blk[[cd]] <- fmt(x$c[match(EFFECTS, x$effect)], 2)
    }
    rows <- rbind(rows, blk); grp <- c(grp, rep(COMP_LAB[[cmp]], length(EFFECTS)))
  }
  g10 <- grouped_body(rows, grp)
  add("Minimum Detectable Pooled Direct Effects (80% Power) and Their Inflation Under Item-Level Missing Data",
    apa_table(g10$body, labels = c("Effect", sub("I\u{00b2} = ", "", hl), cond_labels),
      spans = list(c(" " = 1L, setNames(length(hl), "MDE, complete data (report), by I\u{00b2}"),
                     setNames(9L, "Inflation factor under missing data (MDE under missing data = MDE \u{00d7} factor)")),
                   c(" " = 1L, setNames(length(hl), " "), MCAR = 3L, MAR = 3L, MNAR = 3L)),
      group_rows = g10$group_rows,
      note = paste0("MDE = smallest true pooled standardized effect detected with 80% power (two-sided \u{03b1} = .05; rma() ",
                    sprintf("random effects, REML, z test; k = 2 studies; %d replications), from Table 1 of the sensitivity ", P$R_report),
                    "report, at heterogeneity I\u{00b2} = 0%\u{2013}75% (\u{03c4}\u{00b2} = I\u{00b2}/(1 \u{2212} I\u{00b2}) \u{00d7} typical sampling variance). Inflation factor ",
                    sprintf("= SE of the pooled estimate under %s divided by its value in the report's input; ", pm_lab),
                    "multiplying an MDE (or a critical value) by it gives the value under missing data, exactly when ",
                    "all sampling variances are inflated by the same factor, and for every I\u{00b2} level. ", appr_txt,
                    " > 1.00 = not reached in the range evaluated in the report. ", mech_txt),
      first_col_width = 1.1, col_width = 0.56, size = 9))

  ## Table: MDE of the specific indirect effects ----
  ri <- P$rescaled_indirect[P$rescaled_indirect$missing_method == pm & P$rescaled_indirect$method_code == PRIMARY_RMA, ]
  rows <- NULL; grp <- NULL
  for (cmp in COMPONENTS) {
    blk <- data.frame(e = spec, stringsAsFactors = FALSE)
    for (h in hl) {
      x <- ri[ri$component == cmp & as.character(ri$het) == h & ri$condition == conds[1], ]
      blk[[h]] <- mde_fmt(x$min_ab_report[match(spec, as.character(x$type))], 3, "\u{2014}")
    }
    for (cd in conds) {
      x <- ri[ri$component == cmp & as.character(ri$het) == hl[1] & ri$condition == cd, ]
      blk[[cd]] <- fmt((x$c_a * x$c_b)[match(spec, as.character(x$type))], 2)
    }
    rows <- rbind(rows, blk); grp <- c(grp, rep(COMP_LAB[[cmp]], length(spec)))
  }
  g11 <- grouped_body(rows, grp)
  add("Minimum Detectable Pooled Indirect Effects (80% Power) and Their Inflation Under Item-Level Missing Data",
    apa_table(g11$body, labels = c("Indirect effect", sub("I\u{00b2} = ", "", hl), cond_labels),
      spans = list(c(" " = 1L, setNames(length(hl), "MDE of ab, complete data (report), by I\u{00b2}"),
                     setNames(9L, "Inflation factor under missing data")),
                   c(" " = 1L, setNames(length(hl), " "), MCAR = 3L, MAR = 3L, MNAR = 3L)),
      group_rows = g11$group_rows,
      note = paste0("MDE = smallest true indirect effect (product ab of the pooled paths) detected with 80% power by the ",
                    "joint significance test (rma() random effects, z), from Table 3 of the sensitivity report (grid a, b ",
                    "\u{2264} 0.70; \u{2014} = not reached). Inflation factor = product of the inflation factors of the a and the b path ",
                    sprintf("under %s; multiplied by the report's MDE it gives the MDE under missing data. ", pm_lab),
                    mech_txt),
      first_col_width = 1.3, col_width = 0.56, size = 9))
} else skipped <- c(skipped, "power sensitivity (run 07_power_sensitivity.R)")

## ======================== D. Departures from MAR (delta adjustment) ========================
if (!is.null(S06)) {
  max_d <- max(S06$cfg$DELTAS)
  rows <- NULL; grp <- NULL
  for (cmp in COMPONENTS) {
    blk <- data.frame(e = EFFECT_LAB[EFFECTS], stringsAsFactors = FALSE)
    for (cd in conds) blk[[cd]] <- vapply(EFFECTS, function(e) {
      x <- S06$tips[S06$tips$condition == cd & S06$tips$component == cmp & S06$tips$effect == e, ]
      if (!nrow(x)) return("\u{2014}")
      if (!x$significant_MAR) return("ns")
      if (is.na(x$sig_tip_mnar_interp)) sprintf("> %.2f", max_d) else fmt(x$sig_tip_mnar_interp, 2)
    }, "")
    rows <- rbind(rows, blk); grp <- c(grp, rep(sprintf("%s: direct effects", COMP_LAB[[cmp]]), length(EFFECTS)))
    blk <- data.frame(e = spec, stringsAsFactors = FALSE)
    for (cd in conds) blk[[cd]] <- vapply(spec, function(e) {
      x <- S06$tips_indirect[S06$tips_indirect$condition == cd & S06$tips_indirect$component == cmp &
                               S06$tips_indirect$effect == e, ]
      if (!nrow(x)) return("\u{2014}")
      if (!isTRUE(x$joint_sig_MAR)) return("ns")
      if (is.na(x$sig_tip_mnar_direction)) sprintf("> %.2f", max_d) else fmt(x$sig_tip_mnar_direction, 2)
    }, "")
    rows <- rbind(rows, blk); grp <- c(grp, rep(sprintf("%s: specific indirect effects", COMP_LAB[[cmp]]), length(spec)))
  }
  g12 <- grouped_body(rows, grp)
  add("Delta-Adjustment Sensitivity Analysis: Tipping Points of the Pooled Effects", apa_table(g12$body,
    labels = c("Effect", cond_labels), spans = list(cond_spans()), group_rows = g12$group_rows,
    note = paste0("Delta-adjusted mice-PMM (m = 20; \u{03b4} = 0 is the MAR analysis): imputed item values were shifted by \u{03b4} ",
                  "SDs of the observed item values in the hypothesized MNAR direction (higher X and M, lower Y among ",
                  sprintf("nonrespondents), \u{03b4} = %s. Entries: the smallest \u{03b4} > 0 at which the pooled effect ", paste(fmt(S06$cfg$DELTAS, 2), collapse = ", ")),
                  "loses significance (direct effects: rma() random effects, z test at \u{03b1} = .05, i.e. |\u{03b2}| below its ",
                  "critical value, interpolated on z; indirect effects: joint significance, grid values). ns = not ",
                  sprintf("significant under MAR; > %.2f = no tipping point in the evaluated range.", max_d)),
    first_col_width = 2.0, col_width = 0.7, size = 9), part = "D. Departures From MAR: Delta Adjustment")

  if ("MNAR_15pct" %in% S06$conditions) {
    dd <- S06$direct[S06$direct$method == PRIMARY_RMA & S06$direct$condition == "MNAR_15pct", ]
    rows <- NULL; grp <- NULL
    for (cmp in COMPONENTS) {
      blk <- data.frame(e = EFFECT_LAB[EFFECTS], stringsAsFactors = FALSE)
      for (dl in S06$cfg$DELTAS) blk[[as.character(dl)]] <- vapply(EFFECTS, function(e) {
        x <- dd[dd$component == cmp & dd$effect == e & abs(dd$delta - dl) < 1e-9, ]
        if (!nrow(x)) "\u{2014}" else star(x$est, x$p) }, "")
      com <- D_all[D_all$missing_method == "COM" & D_all$component == cmp, ]
      blk$COM <- vapply(EFFECTS, function(e) star(com$est[com$effect == e], com$p[com$effect == e]), "")
      rows <- rbind(rows, blk); grp <- c(grp, rep(COMP_LAB[[cmp]], length(EFFECTS)))
    }
    g13 <- grouped_body(rows, grp)
    cm <- S06$calibration_means
    cal_txt <- if (!is.null(cm) && "MNAR_15pct" %in% cm$condition) {
      x <- cm[cm$condition == "MNAR_15pct", ]; dlt <- setNames(x$delta_matching_means, x$construct)
      sprintf(paste0(" Calibration (toy data only): at \u{03b4} = 0 the completed item means fall %s to %s item SDs short of ",
                     "the complete data in the MNAR direction; they match at \u{03b4} = %s (X), %s (M) and %s (Y), the shift ",
                     "that corresponds to the simulated MNAR mechanism."),
              fmt(min(-x$shift_at_MAR), 2), fmt(max(-x$shift_at_MAR), 2),
              fmt(dlt[["x"]], 2), fmt(dlt[["m"]], 2), fmt(dlt[["y"]], 2))
    } else ""
    nd <- length(S06$cfg$DELTAS)
    add("Delta-Adjustment Sensitivity Analysis: Pooled Direct Effects by \u{03b4}, MNAR With 15% Missing Items",
      apa_table(g13$body, labels = c("Effect", sub("^-", MINUS, sprintf("%.2f", S06$cfg$DELTAS)), "COM"),
        spans = list(c(" " = 1L, setNames(nd, "\u{03b4} (item-SD units; positive = hypothesized MNAR direction)"), " " = 1L)),
        group_rows = g13$group_rows,
        note = paste0(rma_txt, " Entries are pooled standardized paths after delta-adjusted mice-PMM (m = 20); ",
                      "* p < .05 (random-effects z test). COM = complete data; \u{03b4} = 0 is the MAR analysis.", cal_txt),
        first_col_width = 1.3, col_width = 0.78))
  }
} else skipped <- c(skipped, "delta adjustment (run 06_delta_sensitivity.R)")

## ============================== E. Monte Carlo simulation ==============================
if (!is.null(MC)) {
  bm <- MC$by_method[MC$by_method$method == PRIMARY_RMA, ]
  mcm <- intersect(METHODS_ALL[-1], unique(bm$missing_method))
  pan <- list("Mean absolute bias relative to the population value" = function(x) fmt(x$abs_bias, 3),
              "Root mean square deviation from the complete-data estimate" = function(x) fmt(x$rmsd_vs_com, 3),
              "Coverage of the 95% CI" = function(x) fmt(x$coverage, 3, lead0 = FALSE),
              "95% CI width relative to the complete data" = function(x) fmt(x$ci_ratio, 2))
  for (cmp in COMPONENTS) {
    rows <- NULL; grp <- NULL
    for (nm in names(pan)) for (mt in c("COM", mcm)) {
      if (mt == "COM" && grepl("complete-data|complete data", nm)) next
      x <- bm[bm$component == cmp & bm$missing_method == mt, ]; x <- x[match(MC$design$label, x$label), ]
      rows <- rbind(rows, data.frame(m = METHOD_SHORT[[mt]], t(setNames(pan[[nm]](x), MC$design$label)),
                                     check.names = FALSE, stringsAsFactors = FALSE))
      grp <- c(grp, nm)
    }
    g14 <- grouped_body(rows, grp)
    add(sprintf("Monte Carlo Performance of the Missing-Data Methods: %s", c(uncorrected = "Uncorrected Paths",
                                                                          corrected = "Paths Corrected for Unreliability")[[cmp]]),
      apa_table(g14$body, labels = c("Measure / method", cond_labels), spans = list(cond_spans()),
        group_rows = g14$group_rows,
        note = sprintf(paste0("Monte Carlo simulation (%s mode: %s replications per condition; m = %d imputations). Each ",
                              "replication draws two new toy studies from the population of simulate_toy_data.R, deletes ",
                              "item values as in 02 and runs the meta-analysis. Values are means over the six pooled direct ",
                              "effects (rma() random effects, REML, z). Population value: %s. Monte Carlo SEs are in ",
                              "mc_results.xlsx. %s"),
                       MC$mode, paste(unique(MC$n_rep), collapse = "/"), MC$cfg$M,
                       if (cmp == "corrected") "the latent paths" else "the mean of the two studies' scale-score paths",
                       meth_abbr),
        first_col_width = 3.0, col_width = 0.62, size = 9),
      part = if (cmp == "uncorrected") "E. Monte Carlo Simulation" else NULL)
  }
} else skipped <- c(skipped, "Monte Carlo (run 05_monte_carlo.R in pilot or full mode)")

## ================================= Supplementary tables =================================
rel <- E$reliability
s1_conds <- intersect(c("MAR_15pct", "MNAR_15pct"), conds)
s1_meth <- intersect(c("PRO_FIML", "FIML_TS", "mice_pmm", "PMM_avg"), METH)
ts1 <- NULL
for (k in STUDIES$key) {
  com <- rel[rel$study == k & rel$condition == conds[1] & rel$missing_method == METH[1], ]
  blk <- data.frame(v = VAR_LABELS[VARS], COM = fmt(com$omega_com[match(VARS, com$variable)], 2, lead0 = FALSE),
                    stringsAsFactors = FALSE)
  for (cd in s1_conds) for (mt in s1_meth) {
    x <- rel[rel$study == k & rel$condition == cd & rel$missing_method == mt, ]
    blk[[paste(cd, mt)]] <- fmt(x$omega[match(VARS, x$variable)], 2, lead0 = FALSE)
  }
  ts1 <- rbind(ts1, blk)
}
gs1 <- grouped_body(ts1, rep(sprintf("%s (N = %d dyads)", STUDIES$study, STUDIES$n), each = 6))
add("Reliability Coefficients (\u{03c9}) by Missing-Data Method", apa_table(gs1$body,
  labels = c("Scale score", "COM", rep(METHOD_SHORT[s1_meth], length(s1_conds))),
  spans = list(c(" " = 1L, " " = 1L, setNames(rep(length(s1_meth), length(s1_conds)),
                                              sprintf("%s, 15%% missing", sub("_15pct", "", s1_conds))))),
  group_rows = gs1$group_rows,
  note = paste0("McDonald's \u{03c9} total from a one-factor CFA of the four items, used for the correction for unreliability. ",
                "COM = complete data. PRO-FIML: FIML CFA of the incomplete items (apim_meta_analysis.R); FIML-TS: CFA of ",
                "the stage-1 (FIML, with auxiliaries) item covariance matrix; PMM: mean over the 20 completed datasets; ",
                "PMM-avg: CFA of the averaged dataset (averaging removes the spread of the imputed values and inflates \u{03c9})."),
  first_col_width = 2.4, col_width = 0.62, size = 9), part = "Supplementary Tables", supp = TRUE)

fm_m <- intersect(c("mice_pmm", "mice_cat"), METH)[1]
if (!is.na(fm_m)) {
  st <- E$study_level[E$study_level$missing_method == fm_m, ]
  ts2 <- NULL
  for (k in STUDIES$key) for (cmp in COMPONENTS) {
    blk <- data.frame(s = sprintf("%s, %s", STUDIES$study[STUDIES$key == k], tolower(COMP_LAB[[cmp]])),
                      stringsAsFactors = FALSE)
    for (cd in conds) {
      x <- st$fmi[st$study == k & st$component == cmp & st$condition == cd]
      blk[[cd]] <- sprintf("%s (%s\u{2013}%s)", fmt(mean(x), 2, FALSE), fmt(min(x), 2, FALSE), fmt(max(x), 2, FALSE))
    }
    ts2 <- rbind(ts2, blk)
  }
  add("Fraction of Missing Information of the Study-Level Paths", apa_table(ts2,
    labels = c("Study, component", cond_labels), spans = list(cond_spans()),
    note = sprintf(paste0("Fraction of missing information (Barnard\u{2013}Rubin) of the 12 study-level standardized ",
                          "APIMeM paths under %s: mean (range) over the paths. %s"), METHOD_LABELS[[fm_m]], mech_txt),
    first_col_width = 2.2, col_width = 0.74, size = 9), supp = TRUE)
}

c15 <- intersect(c("MCAR_15pct", "MAR_15pct", "MNAR_15pct"), conds)
for (cmp in COMPONENTS) {
  rows <- NULL; grp <- NULL
  for (cd in c15) {
    blk <- data.frame(e = EFFECT_LAB[EFFECTS], stringsAsFactors = FALSE)
    for (mt in c("COM", METH)) {
      x <- if (mt == "COM") D_all[D_all$missing_method == "COM" & D_all$component == cmp, ] else
        D_all[D_all$missing_method == mt & D_all$component == cmp & D_all$condition == cd, ]
      blk[[mt]] <- vapply(EFFECTS, function(e) {
        y <- x[x$effect == e, ]; sprintf("%s %s", star(y$est, y$p), fmt_ci(y$lo, y$hi)) }, "")
    }
    rows <- rbind(rows, blk)
    ci <- cond_info[cond_info$condition == cd, ]
    grp <- c(grp, rep(sprintf("%s, %s missing", ci$mechanism, rate_lab(ci$rate)), length(EFFECTS)))
  }
  gs3 <- grouped_body(rows, grp)
  add(sprintf("Pooled Direct Effects of All Methods at 15%% Missing Items: %s",
              c(uncorrected = "Uncorrected Paths", corrected = "Paths Corrected for Unreliability")[[cmp]]),
    apa_table(gs3$body, labels = c("Effect", METHOD_SHORT[c("COM", METH)]), group_rows = gs3$group_rows,
      note = paste0("Cells: pooled \u{03b2} [95% CI] (rma() random effects, REML, z); * p < .05. ", meth_short, " ", mech_short),
      first_col_width = 1.0, col_width = 1.12, size = 8), supp = TRUE)
}

if (!is.null(P) && nrow(P$exact_vs_rescaled)) {
  hets <- intersect("I\u{00b2} = 0%", P$het_levels)                 # one level keeps the table on one page
  xr <- P$exact_vs_rescaled[P$exact_vs_rescaled$method_code == PRIMARY_RMA & P$exact_vs_rescaled$missing_method == pm, ]
  ex_c <- intersect(conds, unique(xr$condition))
  rows <- NULL; grp <- NULL
  for (cmp in COMPONENTS) for (h in hets) {
    blk <- data.frame(e = EFFECT_LAB[EFFECTS], stringsAsFactors = FALSE)
    for (cd in ex_c) {
      x <- xr[xr$component_code == cmp & as.character(xr$het) == h & xr$condition == cd, ]
      x <- x[match(EFFECTS, x$effect_code), ]
      blk[[paste(cd, "e")]] <- mde_fmt(x$mde_exact); blk[[paste(cd, "r")]] <- mde_fmt(x$mde_rescaled)
    }
    rows <- rbind(rows, blk); grp <- c(grp, rep(sprintf("%s, %s", COMP_LAB[[cmp]], h), length(EFFECTS)))
  }
  gs5 <- grouped_body(rows, grp)
  ci5 <- cond_info[match(ex_c, cond_info$condition), ]
  ov <- P$approx_summary[P$approx_summary$method_code == PRIMARY_RMA, ]
  ov_m <- tapply(ov$mean_relative_difference, ov$missing_method, mean)
  add("Minimum Detectable Effects Under Missing Data: Exact Reruns Versus Rescaling", apa_table(gs5$body,
    labels = c("Effect", rep(c("Exact", "Rescaled"), length(ex_c))),
    spans = list(c(" " = 1L, setNames(rep(2L, length(ex_c)), sprintf("%s, %s", ci5$mechanism, rate_lab(ci5$rate))))),
    group_rows = gs5$group_rows,
    note = paste0(sprintf("MDE (80%% power; rma() random effects, z; %d replications) under %s. Exact = ", P$R, pm_lab),
                  "apim_meta_sensitivity.R rerun (rma() models only; the multilevel model it also fits is no longer ",
                  "used) on the study-level estimates and sampling covariance matrices of the missing-data method, with ",
                  "the settings and random numbers of the report. Rescaled = the report's MDE \u{00d7} the inflation factor ",
                  sprintf("(Table 10). The rerun on the report's own input reproduces the report (largest difference in an MDE: %s). ",
                          format(signif(P$reproduction$max_abs_diff_mde, 2), scientific = FALSE)),
                  if (length(ov_m)) sprintf("Mean relative difference between exact and rescaled MDEs over all I\u{00b2} levels and effects: %s. ",
                                            paste(sprintf("%s %.1f%%", METHOD_SHORT[names(ov_m)], 100 * ov_m), collapse = ", ")) else "",
                  "> 1.00 = not reached in the evaluated range."),
    first_col_width = 1.2, col_width = 0.72, size = 9), supp = TRUE)
}

## ---- Number the tables and write the document ---------------------------------------------------
heading_xml <- function(text, page_break = TRUE, font = APA_FONT) {
  txt <- gsub("&", "&amp;", text, fixed = TRUE)
  sprintf(paste0('<w:p xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:pPr><w:keepNext/>',
                 '%s<w:jc w:val="center"/><w:spacing w:before="0" w:after="240"/></w:pPr><w:r><w:rPr>',
                 '<w:rFonts w:ascii="%2$s" w:hAnsi="%2$s" w:cs="%2$s"/><w:b/><w:sz w:val="24"/><w:szCs w:val="24"/>',
                 '</w:rPr><w:t xml:space="preserve">%3$s</w:t></w:r></w:p>'),
          if (page_break) "<w:pageBreakBefore/>" else "", font, txt)
}
num <- 0; snum <- 0
for (i in seq_along(tables)) {
  if (tables[[i]]$supp) { snum <- snum + 1; tables[[i]]$number <- paste0("S", snum) }
  else { num <- num + 1; tables[[i]]$number <- as.character(num) }
}
doc <- new_landscape_doc()
intro <- c(
  "APIM Meta-Analysis With Item-Level Missing Data: Overview of the Sensitivity Results (Toy Data)",
  sprintf(paste0("Two dyadic toy studies (N = %d and %d dyads; 24 Likert items each) were made incomplete under MCAR, MAR ",
                 "and MNAR at 5%%, 10%% and 15%% per item, handled with item-level missing-data methods, and analysed with ",
                 "the planned APIMeM meta-analysis (lavaan per study; rma() random effects across studies). Generated %s."),
          STUDIES$n[1], STUDIES$n[2], format(Sys.Date())),
  sprintf("Method decision, based on %s: %s. Chosen: %s.", DEC$basis, DEC$verdict, METHOD_LABELS[[chosen]]),
  if (length(skipped)) paste0("Not included: ", paste(skipped, collapse = "; "), ".") else "")
doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(intro[1],
  officer::fp_text(bold = TRUE, font.size = 12, font.family = APA_FONT)), fp_p = officer::fp_par(text.align = "center")))
for (tx in intro[-1][nzchar(intro[-1])])
  doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(tx, officer::fp_text(font.size = 11, font.family = APA_FONT)),
                                                   fp_p = officer::fp_par(padding.top = 6)))
for (tb in tables) {
  if (!is.null(tb$part)) doc <- officer::body_add_xml(doc, heading_xml(tb$part))
  doc <- add_apa_table(doc, tb$number, tb$title, tb$ft, new_page = is.null(tb$part))
}
print(doc, target = OUT_DOCX)
cat("Written:", OUT_DOCX, sprintf("(%d tables)\n", length(tables)))
cat("Method decision:", DEC$basis, "->", METHOD_LABELS[[chosen]], "\n")
if (length(skipped)) cat("Not included:", paste(skipped, collapse = "; "), "\n")
