# =============================================================================
# 02_ampute_data.R
# Data processing: create the 18 incomplete versions of the two toy datasets
# (2 studies x MCAR/MAR/MNAR x 5%/10%/15%), with missing values at the ITEM
# level (the X, M and Y variables are scale scores of 4 items each).
#
#   Input : data/toy_study1_complete.csv, data/toy_study2_complete.csv (from 01;
#           N = 130 and 174 dyads, 24 items, no missing values)
#   Output: data/amputed/
#     study<k>_<MECH>_<RATE>pct.csv / .xlsx   18 incomplete datasets
#     study<k>_complete_with_aux.csv / .xlsx  the complete data with the added
#                                             auxiliaries and instruments (truth)
#     amputation_summary.xlsx                 settings, codebook, missingness,
#                                             mechanism checks
#     amputed_datasets.rds                    all of the above as R objects
#
# Added variables (the toy data contain only the items):
#   aux_p1, aux_p2  AUXILIARY variables of Partner 1 / Partner 2 (continuous),
#                   aux_dyad a dyad-level AUXILIARY (binary). Always observed,
#                   related to X, M and Y; they drive the MAR mechanism. MAR
#                   needs a fully observed cause of missingness; including it in
#                   the imputation / FIML model makes the missingness ignorable
#                   (inclusive strategy; Collins, Schafer & Kam, 2001). Because
#                   they are related to the outcomes, methods that ignore them
#                   (e.g. proration + FIML on the scale scores) are biased under MAR.
#   iv_p1, iv_p2    INSTRUMENTS (one per partner). Independent of everything
#                   else; they shift the probability of item nonresponse of that
#                   partner in the MNAR datasets only. MNAR data do not need them
#                   (self-masking suffices), but an MNAR imputation model
#                   (Heckman-MI in 03) is identified without leaning on
#                   distributional assumptions only if a variable affects
#                   response and not the values (exclusion restriction; Heckman,
#                   1979; Tchetgen Tchetgen & Wirth, 2017).
#
# Mechanisms (per item; logistic selection model, see 00_functions.R):
#   MCAR  P(missing) constant
#   MAR   depends on the partner's auxiliaries only (low aux -> more missing)
#   MNAR  depends on the item value itself (self-masking, 1 log-odds per SD, in
#         the direction MNAR_DIRECTION: higher X and M, lower Y -> more missing)
#         and on the partner's instrument (1 log-odds per SD)
# Rates: EXACTLY 5%, 10% and 15% of every item (round(rate * N): 6/13/20 of 130
# and 9/17/26 of 174), hence of all item cells. One set of uniform random
# numbers per study is shared by its nine datasets, so within a mechanism the
# missing cells are nested (the 5% cells are also missing at 10% and 15%) and
# the mechanisms differ only by their selection model. The auxiliaries,
# instruments and dyad identifiers stay complete.
# =============================================================================

source("00_functions.R")

## ---- 0. Settings ----------------------------------------------------------------
IN_DIR  <- "data"
OUT_DIR <- file.path("data", "amputed")
SEED    <- 20261012
AMP     <- AMPUTE_DEFAULTS          # strength 1, MAR direction -1, instrument strength 1
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Read the complete toy data ---------------------------------------------------
files <- file.path(IN_DIR, sprintf("toy_%s_complete.csv", STUDIES$key))
if (!all(file.exists(files))) source("01_complete_toy_data.R", local = new.env())
dat0 <- setNames(lapply(files, utils::read.csv, stringsAsFactors = FALSE), STUDIES$key)
for (k in STUDIES$key) {
  stopifnot(all(c(ID_VAR, ITEMS) %in% names(dat0[[k]])), !anyNA(dat0[[k]][, ITEMS]),
            all(as.matrix(dat0[[k]][, ITEMS]) %in% ITEM_VALUES))
  cat(sprintf("%s: %d dyads, %d items, no missing values\n", k, nrow(dat0[[k]]), length(ITEMS)))
}

## ---- 2. Add the auxiliaries and the instruments, draw the uniform numbers -----------
set.seed(SEED)
complete <- list(); U <- list()
for (k in STUDIES$key) {
  complete[[k]] <- add_aux_iv(dat0[[k]])
  U[[k]] <- matrix(runif(nrow(complete[[k]]) * length(ITEMS)), ncol = length(ITEMS),
                   dimnames = list(NULL, ITEMS))          # common random numbers (9 datasets)
}
aux_check <- do.call(rbind, lapply(STUDIES$key, function(k) {
  S <- scale_scores(complete[[k]]); d <- complete[[k]]
  do.call(rbind, lapply(c(AUX_VARS, IV_VARS), function(v) data.frame(
    study = k, variable = v,
    t(setNames(round(stats::cor(d[[v]], S)[1, ], 3), paste0("r_", VARS))),
    R2_on_scale_scores = round(summary(stats::lm(d[[v]] ~ S))$r.squared, 3))))
}))
cat("\nAuxiliaries and instruments: correlations with the scale scores\n")
print(aux_check, row.names = FALSE)

## ---- 3. Amputation ------------------------------------------------------------------------
grid <- expand.grid(rate = RATES, mechanism = MECHANISMS, study = STUDIES$key,
                    stringsAsFactors = FALSE)[, 3:1]
grid$condition <- sprintf("%s_%02dpct", grid$mechanism, round(100 * grid$rate))
grid$dataset   <- sprintf("%s_%s", grid$study, grid$condition)
grid$n_dyads   <- STUDIES$n[match(grid$study, STUDIES$key)]
grid$n_missing_per_item <- round(grid$rate * grid$n_dyads)

amputed <- setNames(lapply(seq_len(nrow(grid)), function(i)
  ampute_items(complete[[grid$study[i]]], grid$mechanism[i], grid$rate[i], U = U[[grid$study[i]]],
               strength = AMP$strength, mar_direction = AMP$mar_direction,
               iv_strength = AMP$iv_strength)), grid$dataset)

## Checks: exact counts, untouched columns, observed values unchanged, nesting
for (i in seq_len(nrow(grid))) {
  a <- amputed[[i]]; cm <- complete[[grid$study[i]]]
  stopifnot(all(colSums(is.na(a[, ITEMS])) == grid$n_missing_per_item[i]),
            !anyNA(a[, c(ID_VAR, AUX_VARS, IV_VARS)]),
            identical(a[, c(ID_VAR, AUX_VARS, IV_VARS)], cm[, c(ID_VAR, AUX_VARS, IV_VARS)]),
            all(is.na(a[, ITEMS]) | a[, ITEMS] == cm[, ITEMS]))
}
for (k in STUDIES$key) for (mch in MECHANISMS) {
  ds <- grid$dataset[grid$study == k & grid$mechanism == mch]
  for (j in 2:length(ds))
    stopifnot(all(is.na(amputed[[ds[j]]][, ITEMS])[is.na(amputed[[ds[j - 1]]][, ITEMS])]))
}
cat("\nAll 18 datasets: exact missing counts per item, nesting across rates and untouched",
    "auxiliaries/instruments verified.\n")

## ---- 4. Summaries -----------------------------------------------------------------------------
miss_long <- do.call(rbind, lapply(names(amputed), function(nm)
  cbind(dataset = nm, missing_summary(amputed[[nm]]))))
miss_wide <- reshape(miss_long[, c("dataset", "variable", "pct_missing")], idvar = "variable",
                     timevar = "dataset", direction = "wide")
names(miss_wide) <- sub("pct_missing\\.", "", names(miss_wide)); rownames(miss_wide) <- NULL

mech_long <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  cbind(dataset = grid$dataset[i], mechanism_check(complete[[grid$study[i]]], amputed[[i]]))))
mech_avg <- do.call(rbind, lapply(split(mech_long, mech_long$dataset), function(d) data.frame(
  dataset = d$dataset[1],
  smd_own_value_in_mnar_direction = mean(d$smd_own_value * d$mnar_direction),
  smd_aux_composite = mean(d$smd_aux_composite), smd_iv = mean(d$smd_iv),
  smd_own_scale_score_in_mnar_direction = mean(d$smd_own_scale_score * d$mnar_direction))))
mech_avg <- mech_avg[match(grid$dataset, mech_avg$dataset), ]; rownames(mech_avg) <- NULL
mech_reg <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  r <- mechanism_regression(complete[[grid$study[i]]], amputed[[i]])
  data.frame(dataset = grid$dataset[i], study = grid$study[i], mechanism = grid$mechanism[i],
             rate = grid$rate[i], term = r$term, log_odds = r$log_odds, se = r$se)
}))
mech_reg_wide <- reshape(transform(mech_reg, value = sprintf("%.2f (%.2f)", log_odds, se))[
  , c("dataset", "term", "value")], idvar = "dataset", timevar = "term", direction = "wide")
names(mech_reg_wide) <- c("dataset", "own value (MNAR direction) | aux, iv",
                          "aux composite | own, iv", "instrument | own, aux")
rownames(mech_reg_wide) <- NULL
cat("\nRealised selection model, log-odds per SD (SE), stacked over the 24 items:\n")
print(mech_reg_wide, row.names = FALSE)

## ---- 5. Write the datasets ------------------------------------------------------------------
cb <- data.frame(
  variable = c(ID_VAR, ITEMS, AUX_VARS, IV_VARS),
  description = c("Dyad identifier",
                  sprintf("%s, item %s (missing = empty cell)", VAR_LABELS[ITEM_VAR], sub(".*_", "", ITEMS)),
                  "Auxiliary, Partner 1 (continuous): related to Partner 1's X, M, Y; always observed; drives MAR",
                  "Auxiliary, Partner 2 (continuous): related to Partner 2's X, M, Y; always observed; drives MAR",
                  "Auxiliary, dyad (binary): related to both partners' Y; always observed; drives MAR",
                  "Instrument for Partner 1's item nonresponse: independent of all other variables; MNAR datasets only",
                  "Instrument for Partner 2's item nonresponse: independent of all other variables; MNAR datasets only"),
  values = c("D001-...", rep("1-5", length(ITEMS)), "numeric (approx. standard normal)",
             "numeric (approx. standard normal)", "0/1", "numeric (standard normal)", "numeric (standard normal)"),
  role = c("id", rep("incomplete item", length(ITEMS)), rep("auxiliary", 3), rep("instrument", 2)),
  scale_score = c("", sprintf("%s = mean of %s_1-%s_%d", ITEM_VAR, ITEM_SCALE, ITEM_SCALE, N_ITEMS),
                  rep("", 5)),
  row.names = NULL)

mechanism_info <- data.frame(
  item = c("Selection model", "MCAR", "MAR", "MNAR", "MNAR direction per construct", "Strength",
           "MAR direction", "Instrument strength", "Exact rates", "Nesting", "Seed", "R"),
  value = c("P(M_ij = 1) = logistic(a_j + eta_ij), per item j; a_j set so that exactly round(rate * N) values are deleted",
            "eta = 0",
            "eta = MAR direction * strength * z(z(aux_pk) + z(aux_dyad)), k = the item's partner",
            "eta = direction(construct) * strength * z(item value) + iv_strength * z(iv_pk)  (self-masking + instrument)",
            paste(sprintf("%s %+d", names(MNAR_DIRECTION), MNAR_DIRECTION), collapse = "; "),
            AMP$strength, sprintf("%d (lower auxiliaries -> more missing)", AMP$mar_direction),
            sprintf("%s (MNAR only)", AMP$iv_strength),
            paste0(100 * RATES, "% of every item", collapse = "; "),
            "Common uniform numbers per study for its 9 datasets: missing cells nested across rates within a mechanism",
            SEED, R.version.string))

for (nm in names(amputed)) {
  a <- amputed[[nm]]; attr(a, "amputation") <- NULL
  utils::write.csv(a, file.path(OUT_DIR, paste0(nm, ".csv")), row.names = FALSE, na = "")
  write_xlsx_sheets(list(data = a, codebook = cb, missingness = missing_summary(a),
                         mechanism_check = mech_long[mech_long$dataset == nm, -1]),
                    file.path(OUT_DIR, paste0(nm, ".xlsx")))
}
for (k in STUDIES$key) {
  utils::write.csv(complete[[k]], file.path(OUT_DIR, sprintf("%s_complete_with_aux.csv", k)), row.names = FALSE)
  write_xlsx_sheets(list(data = complete[[k]], codebook = cb, aux_iv_check = aux_check[aux_check$study == k, ]),
                    file.path(OUT_DIR, sprintf("%s_complete_with_aux.xlsx", k)))
}
write_xlsx_sheets(list(datasets = grid, settings = mechanism_info, codebook = cb,
                       missingness_pct = miss_wide, missingness_long = miss_long,
                       selection_model = mech_reg_wide, selection_model_long = mech_reg,
                       smd_by_dataset = mech_avg, smd_by_item = mech_long, aux_iv_check = aux_check),
                  file.path(OUT_DIR, "amputation_summary.xlsx"))
amp_clean <- lapply(amputed, function(a) { attr(a, "amputation") <- NULL; a })
saveRDS(list(complete = complete, amputed = amp_clean, grid = grid, settings = AMP,
             mnar_direction = MNAR_DIRECTION, seed = SEED, U = U),
        file.path(OUT_DIR, "amputed_datasets.rds"))

cat("\nWritten to", OUT_DIR, ":", length(list.files(OUT_DIR, "_(MCAR|MAR|MNAR)_.*\\.csv$")),
    "incomplete datasets (.csv and .xlsx), the complete data with auxiliaries, the summary and the .rds\n")
summ <- miss_long[miss_long$variable %in% c("ALL item cells", "Dyads with any missing item",
                                            "Scale scores not computable by proration (< 75% items)"), ]
print(reshape(summ[, c("dataset", "variable", "pct_missing")], idvar = "dataset", timevar = "variable",
              direction = "wide"), row.names = FALSE)
