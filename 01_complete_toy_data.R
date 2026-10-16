
## title: "APIM meta-analysis: individual-level missing and final pipeline"
## author: "Anonymized"
## date: "15-10-2026"

# =============================================================================
# 01_complete_toy_data.R
# The complete toy datasets (no missing values) that the amputation in 02
# starts from.
#
#   Input : toy_study1.csv, toy_study2.csv (from simulate_toy_data.R)
#   Output: data/
#     toy_study1_complete.csv / .xlsx   N = 130 dyads, 24 items, no missing values
#     toy_study2_complete.csv / .xlsx   N = 174 dyads
#     population_values.rds             true paths (latent and scale-score level)
#     toy_complete_check.xlsx           restored cells, population values,
#                                       APIMeM on the complete data vs the truth
#
# Why this step: simulate_toy_data.R deletes ~1.5% of the item values (MCAR,
# p_miss = .015), so the attached toy files are NOT complete (49 and 69 missing
# cells). Those cells are restored by re-running the generator with the same
# seed and p_miss = 0. The random-number stream is unchanged (runif() is drawn
# either way), so every observed value is reproduced exactly; this is checked
# below. The generator draws the latent scores with eigenvectors oriented as on
# the platform that generated the toy files (mvrnorm_ref() in 00_functions.R),
# because MASS::mvrnorm() itself gives different draws with some LAPACK
# libraries (R 4.6 on macOS). Should the check still fail, the deleted values
# stored in 00_functions.R are used for the attached toy files, and only for
# other input files a single PMM draw fills the missing cells. The complete data
# are always compared with the reference version (checksum): all later results
# reproduce the delivered ones only if they are identical.
# =============================================================================

source("00_functions.R")

## 0. Settings
OUT_DIR <- "data"

dir.create(OUT_DIR, showWarnings = FALSE)

## 1. Read the attached toy data and regenerate it without missingness

orig <- setNames(lapply(STUDIES$toy_file, function(f) {
  f <- find_input(f)                                   # this folder, then ../meta-analysis
  if (!file.exists(f)) stop(not_found(f), " (run simulate_toy_data.R first)")
  utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}), STUDIES$key)
for (k in STUDIES$key) stopifnot(all(c(ID_VAR, ITEMS) %in% names(orig[[k]])))
regen <- regenerate_toy_data(p_miss = 0)

restore <- list(); check <- list(); complete <- list()
for (k in STUDIES$key) {
  o <- orig[[k]][, c(ID_VAR, ITEMS)]; r <- regen[[k]][, c(ID_VAR, ITEMS)]
  obs <- !is.na(o[, ITEMS])
  agree <- nrow(o) == nrow(r) && identical(o[[ID_VAR]], r[[ID_VAR]]) &&
    all(as.matrix(o[, ITEMS])[obs] == as.matrix(r[, ITEMS])[obs])
  stored <- as.integer(strsplit(TOY_MISSING_VALUES[[k]], "")[[1]])
  is_toy <- toy_checksum(o[, ITEMS]) == TOY_CHECKSUM$observed[[k]] && length(stored) == sum(!obs)
  if (agree) {
    complete[[k]] <- r
    how <- "regenerated (same seed, p_miss = 0); all observed values identical"
  } else if (is_toy) {
    message(k, ": the generator does not reproduce the attached file on this platform; ",
            "the deleted values stored in 00_functions.R are used.")
    m <- as.matrix(o[, ITEMS]); m[!obs] <- stored
    complete[[k]] <- cbind(o[ID_VAR], as.data.frame(m))
    how <- "stored values of the deleted cells (the generator did not reproduce the attached file here)"
  } else {
    warning(k, ": regeneration does not reproduce the input file; the ", sum(!obs),
            " missing cells are filled with one PMM draw instead.")
    imp <- with_seed(20261011, mice::mice(o[, ITEMS], m = 1, maxit = 20, method = "pmm", printFlag = FALSE))
    complete[[k]] <- cbind(o[ID_VAR], mice::complete(imp, 1))
    how <- "single PMM imputation of the missing cells (regeneration did not match)"
  }
  stopifnot(!anyNA(complete[[k]]), all(as.matrix(complete[[k]][, ITEMS]) %in% ITEM_VALUES))
  reference <- toy_checksum(complete[[k]][, ITEMS]) == TOY_CHECKSUM$complete[[k]]
  w <- which(!obs, arr.ind = TRUE)
  restore[[k]] <- data.frame(study = k, dyad_id = o[[ID_VAR]][w[, 1]], item = ITEMS[w[, 2]],
                             restored_value = as.matrix(complete[[k]][, ITEMS])[w])
  check[[k]] <- data.frame(study = k, n_dyads = nrow(o), n_item_cells = length(obs),
                           n_missing_in_attached_file = sum(!obs),
                           pct_missing_in_attached_file = round(100 * mean(!obs), 2),
                           observed_values_reproduced = agree, method = how,
                           identical_to_reference = reference)
}

restore <- do.call(rbind, restore); check <- do.call(rbind, check); rownames(restore) <- NULL

cat("\nAttached toy files and their completion:\n")

print(check[, c(1:6, 8)], row.names = FALSE)
if (!all(check$identical_to_reference))
  warning("The complete data differ from the reference version, so the amputed datasets (02) and all ",
          "later results will differ from the delivered ones (cached results are recomputed).")

## 2. Population values

pop <- toy_population(file.path(OUT_DIR, "population_values.rds"))

pop_tab <- data.frame(path = PATHS$label, effect = PATHS$effect,
                      latent_true = pop$latent[PATHS$label],
                      scale_score_study1 = pop$uncorrected$study1[PATHS$label],
                      scale_score_study2 = pop$uncorrected$study2[PATHS$label], row.names = NULL)

omega_tab <- data.frame(variable = VARS, scale = unname(SCALES),
                        omega_study1 = pop$omega$study1[VARS], omega_study2 = pop$omega$study2[VARS],
                        row.names = NULL)

## 3. APIMeM on the complete data (benchmark check against the truth)

com_chk <- do.call(rbind, lapply(seq_len(nrow(STUDIES)), function(i) {
  k <- STUDIES$key[i]
  r <- analyse_complete(complete[[k]], STUDIES$study[i])
  data.frame(study = STUDIES$study[i], path = PATHS$label,
             uncorrected = r$uncorrected$est[PATHS$label],
             se_uncorrected = sqrt(diag(r$uncorrected$V$ML))[PATHS$label],
             true_scale_score = pop$uncorrected[[k]][PATHS$label],
             corrected = r$corrected$est[PATHS$label],
             se_corrected = sqrt(diag(r$corrected$V$ML))[PATHS$label],
             true_latent = pop$latent[PATHS$label], row.names = NULL)
}))

com_chk$z_uncorrected <- (com_chk$uncorrected - com_chk$true_scale_score) / com_chk$se_uncorrected

com_chk$z_corrected   <- (com_chk$corrected - com_chk$true_latent) / com_chk$se_corrected

cat(sprintf("\nAPIMeM on the complete toy data vs the population values: max |z| = %.2f (uncorrected), %.2f (corrected)\n",
            max(abs(com_chk$z_uncorrected)), max(abs(com_chk$z_corrected))))

## 4. Write

codebook <- data.frame(
  variable = c(ID_VAR, ITEMS),
  description = c("Dyad identifier",
                  sprintf("%s, item %s", VAR_LABELS[ITEM_VAR], sub(".*_", "", ITEMS))),
  values = c("D001-...", rep("1-5 (Likert)", length(ITEMS))),
  scale_score = c("", sprintf("%s = mean of %s_1-%s_%d", ITEM_VAR, ITEM_SCALE, ITEM_SCALE, N_ITEMS)),
  row.names = NULL)

for (i in seq_len(nrow(STUDIES))) {
  k <- STUDIES$key[i]
  f <- file.path(OUT_DIR, sprintf("toy_%s_complete", k))
  utils::write.csv(complete[[k]], paste0(f, ".csv"), row.names = FALSE)
  write_xlsx_sheets(list(data = complete[[k]], codebook = codebook,
                         restored_cells = restore[restore$study == k, ]), paste0(f, ".xlsx"))
}

write_xlsx_sheets(list(attached_files = check, restored_cells = restore, population_paths = pop_tab,
                       population_omega = omega_tab, complete_data_apimem = com_chk),
                  file.path(OUT_DIR, "toy_complete_check.xlsx"))

cat("\nWritten:", file.path(OUT_DIR, c(sprintf("toy_%s_complete.csv/.xlsx", STUDIES$key),
                                       "population_values.rds", "toy_complete_check.xlsx")), sep = "\n")
