# APIM meta-analysis with item-level missing data: amputation, imputation and sensitivity analysis

This folder prepares the APIM (APIMeM) meta-analysis of the two dyadic studies (N = 130 and N = 174 dyads) for **missing data at the item level**. X, M and Y are scale scores (means of 4 Likert items per partner), so items, not scale scores, go missing.

The pipeline follows our usual missing-data pipeline (`LCA/3 R code/2 Analysis example/3 Demographics missing preregistration`), adapted from a regression to the APIMeM meta-analysis:

1.  The toy data are made incomplete under MCAR, MAR and MNAR at 5%, 10% and 15% (2 studies × 9 = 18 datasets).
2.  The items are imputed (mice-PMM and others, m = 20), or the missing items are handled by FIML.
3.  Every method is run through the planned meta-analysis:
    - the APIMeM is estimated per study in **lavaan** with the functions of `apim_meta_analysis.R`;
    - the two studies are pooled with **`rma()`** (random effects), the specification chosen from the sensitivity report.
4.  The sensitivity part relies on the output of the earlier power analysis (`apim_meta_sensitivity.R` and `APIM_sensitivity_report.docx`).

The scripts live in `APIM/6 Meta-analysis/Missing data/` and use the original files in place, so nothing needs to be copied: `apim_meta_analysis.R` and `toy_study1.csv`/`toy_study2.csv` are looked up in this folder and then in `../Meta-analysis/`; `apim_meta_sensitivity.R` in this folder and then in `../Sensitivity/`; the outputs of the sensitivity analysis (`APIM_sensitivity_results.rds`, `APIM_meta_analysis_results.rds`, `APIM_sensitivity_*.csv`) in `../Sensitivity/`. A copy placed in this folder takes precedence. Set the working directory to this folder and run the scripts in order.

**Running from a clone of the GitHub repository.** The repository is self-contained: `apim_meta_analysis.R`, `apim_meta_sensitivity.R` and the toy files are in the root, and the three outputs of the sensitivity analysis that the scripts read (`APIM_sensitivity_results.rds`, `APIM_meta_analysis_results.rds`, `APIM_sensitivity_observed_critical.csv`) are in `Sensitivity/`, which is searched after `../Sensitivity/`. All cached results (`output/*/raw/`) are included, so `03`–`08` rerun from these caches in seconds to minutes and reproduce the stored results.

``` r
# Open "Missing data.Rproj" (or setwd() to this folder). Everything is already computed; to rerun:
Sys.setenv(WORKERS = 6)                      # parallel jobs (optional; the default is 1)
for (f in c("01_complete_toy_data.R", "02_ampute_data.R", "03_impute_missing_data.R",
            "04_evaluate_meta_analysis.R", "06_delta_sensitivity.R", "07_power_sensitivity.R"))
  source(f, local = new.env())
Sys.setenv(MC_MODE = "pilot", MC_NSIM = 10); source("05_monte_carlo.R", local = new.env())  # optional, long
source("08_apa_tables.R", local = new.env())  # -> output/report/08_missing_data_sensitivity_tables.docx
```

The scripts reuse cached intermediate results (`output/*/raw/`). Each cached result carries a key of the exact inputs it was computed from (data, imputations, settings, seed), so it is reused only for unchanged inputs; results of other data or settings are recomputed automatically. Delete a cache folder to recompute a step anyway.

**Reproducibility.** Rerunning the scripts reproduces the delivered results exactly, on any platform and with any number of `WORKERS`:

- `simulate_toy_data.R` draws the latent scores with `MASS::mvrnorm()`, which uses eigenvectors. Some LAPACK libraries (e.g. in R 4.6 on macOS) return eigenvectors with opposite signs, which changes every draw. The toy generator here orients them as on the platform that generated the toy files. `01` also checks that the complete data equal the reference version (column `identical_to_reference`).
- Every imputation draws from its own seeded L'Ecuyer-CMRG stream, so it does not depend on whether it runs in a parallel worker.
- The Monte Carlo CIs of the indirect effects use Cholesky draws, which do not depend on the platform either.
- Sourcing `00_functions.R` resets R's random number generator to the default, so a generator left in the session by other code cannot change the results.

| File | Purpose | Runtime (1 core) |
|----|----|----|
| `00_functions.R` | Shared library, sourced by all scripts. It loads the APIMeM functions of `apim_meta_analysis.R` (read, not copied) and contains the toy population, auxiliaries and instruments, item-level amputation, item-level imputation, two-stage FIML, Rubin's rules within studies, `rma()` pooling, the method decision and the link to `apim_meta_sensitivity.R`. | – |
| `01_complete_toy_data.R` | Restores the complete toy data. The attached toy files already contain \~1.5% MCAR item-level missing values (`p_miss = .015` in `simulate_toy_data.R`: 49 and 69 cells). These are restored by re-running the generator with the same seed and `p_miss = 0`, which reproduces every observed value exactly, and the result is checked against the reference version. Also computes the population values. | seconds |
| `02_ampute_data.R` | **The data-processing step.** It adds the auxiliaries and the instruments, then writes the **18 incomplete datasets** as .csv and .xlsx, together with the complete data, a codebook and mechanism checks. | seconds |
| `03_impute_missing_data.R` | **Item-level multiple imputation** of the 18 datasets with mice-PMM and mice-CAT, m = 20 (Heckman-MI optional). It writes the stacked imputed datasets, an .rds and the imputation log. | \~40 min |
| `04_evaluate_meta_analysis.R` | Runs every method through the meta-analysis (lavaan → Rubin's rules within studies → `rma()`). It compares each method with the complete data and the population values, makes the PMM-vs-FIML decision and checks that the `rma()` code reproduces the report. | \~20 min (cached) |
| `05_monte_carlo.R` | Monte Carlo simulation (ADEMP) of the methods; quick / pilot / full modes; resumable and parallel. | \~7 min per replication |
| `06_delta_sensitivity.R` | **MNAR sensitivity analysis**: delta-adjusted mice-PMM with tipping points, judged against the `rma()` critical values of the report. | \~1.5 h |
| `07_power_sensitivity.R` | **Power sensitivity**: the sensitivity report under missing data. It reruns `apim_meta_sensitivity.R` unchanged (rma() models only) with the missing-data sampling covariances and rescales the report's MDEs, critical values and indirect-effect MDEs. | \~30 min |
| `08_apa_tables.R` | **APA tables** (Word) giving an overview of the sensitivity results of 02–07. | \~30 s |
| `09_final_meta_analysis.R` | The script for the **real data**: missing-data handling with the chosen method, lavaan per study, Rubin's rules, `rma()`, delta adjustment, power sensitivity and APA tables. | \~15 min |

Requirements: R ≥ 4.1 and the packages `lavaan`, `metafor`, `clubSandwich`, `mice`, `MASS`, `pbivnorm`, `numDeriv`, `openxlsx`, `flextable`, `officer` and `ggplot2` (plus `future` and `future.apply` for `WORKERS > 1`). The results here were computed with R 4.3.3, lavaan 0.6-17, mice 3.16.0 and metafor 4.4.0. The code also runs with lavaan 0.7 (which renamed the EM settings used by the two-stage FIML; checked with lavaan 0.7-3, where the two-stage FIML results agree with lavaan 0.6-17 to within 10⁻⁸).

``` r
install.packages(c("lavaan", "metafor", "clubSandwich", "mice", "MASS", "pbivnorm", "numDeriv", "openxlsx",
                   "flextable", "officer", "ggplot2", "future", "future.apply"))
```

## Step 1: the incomplete datasets (`02_ampute_data.R`)

**Added variables.** The toy data contain only the 24 items, so three auxiliaries and two instruments are added. All five are always observed.

- `aux_p1`, `aux_p2` (continuous) are partner-level **auxiliary variables**, related to the partner's own X, M and Y and weakly to the other partner's Y. They stand in for, e.g., a screening measure, relationship duration or a related questionnaire. Their R² on the six scale scores is .31–.45.
- `aux_dyad` (binary) is a dyad-level auxiliary related to both partners' Y.
- The auxiliaries drive the MAR mechanism. MAR is ignorable only if they are in the imputation or FIML model (inclusive strategy; Collins, Schafer & Kam, 2001).
- `iv_p1`, `iv_p2` are **instruments** for item nonresponse, one per partner. They are independent of everything else and shift the probability of nonresponse in the MNAR datasets only.
  - MNAR data do not need an instrument, because self-masking suffices.
  - A selection model, however, is identified without leaning on distributional assumptions only through a variable that affects response but not the values (exclusion restriction; Heckman, 1979; Tchetgen Tchetgen & Wirth, 2017). Heckman-MI uses it.

**Mechanisms.** Each item has its own logistic selection model, P(missing) = logistic(a_j + η_ij).

|   | η_ij | Realized log-odds per SD (own value / auxiliaries / instrument), 15% |
|----|----|----|
| MCAR | 0 | 0.02 / −0.01 / −0.02 |
| MAR | −1 · z(partner's auxiliary composite) | 0.03 / −1.05 / −0.03 |
| MNAR | direction · z(item value) + 1 · z(partner's instrument) | 1.00 / 0.04 / 0.97 |

- **MNAR directions** are assumptions: higher X and M values and lower Y values are missing more often. Think of X as a stressor, M as strain and Y as relationship quality. Change `MNAR_DIRECTION` in `00_functions.R` if your constructs differ.
- **Rates.** Exactly round(rate × N) values are deleted from **every item**: 6, 13 and 20 of 130, and 9, 17 and 26 of 174.
- **Nesting.** Within a study and mechanism, the 5% cells are also missing at 10% and 15% (common random numbers).
- **What that means for the scale scores and dyads.** At 15%, 88–99% of the dyads miss at least one item, and 11–15% of the scale scores have fewer than 75% of their items. `apim_meta_analysis.R` would set those scale scores to missing.

**Output (`data/amputed/`).**

- `study<k>_<MECH>_<RATE>pct.csv/.xlsx`: the **18 incomplete datasets**. The .xlsx files add codebook, missingness and mechanism-check sheets.
- `study<k>_complete_with_aux.csv/.xlsx`: the truth.
- `amputation_summary.xlsx` and `amputed_datasets.rds`.

## Step 2: imputation (`03_impute_missing_data.R`) and FIML

All methods work on the **items** and only then form the scale scores that enter the APIMeM.

- **mice-PMM** (m = 20, 20 iterations, 5 donors) is the default. Imputed values are always observed values 1–5.
- **The imputation model of each item** contains:
  - the other three items of its scale;
  - the **five other scale scores**, updated passively inside the chains, so both partners' X, M and Y are there and the dyadic structure of the APIMeM is preserved ("parcel summaries"; Eekhout et al., 2018);
  - the three auxiliaries.
- This keeps the model congenial with an analysis of scale scores and uses 11 instead of 26 predictors, which matters with n = 130.
- **mice-CAT**: a proportional-odds model per item, with parameters drawn from their posterior.
- **Heckman-MI** is optional (`METHODS=mice_pmm,mice_cat,heckman`). It is the only MNAR model, but at the item level it fits 24 selection models from 6–26 missing values each, with a weak instrument. It takes 10–15 min per dataset and was not closer to the complete data than FIML in a pilot.
- **Two-stage FIML (FIML-TS)** is the FIML method at the item level (Savalei & Bentler, 2009; Savalei & Rhemtulla, 2017):
  - Stage 1: a saturated normal model of the 24 items and the auxiliaries by FIML.
  - Stage 2: the APIMeM, uncorrected and corrected, fitted in lavaan to the implied scale-score moments, with the stage-1 sampling covariance matrix.
  - It is the analytic counterpart of item-level MI. The stage-2 standard errors equal the delta method to 10⁻⁸ (checked in 04).
- **PRO-FIML** is what `apim_meta_analysis.R` does by default: proration (≥ 75% of the items) and FIML on the scale scores, without auxiliaries.
- **Aggregating the m = 20 imputations.** Per study, the APIMeM is fitted to every completed dataset. The 12 standardized paths and their full 12 × 12 sampling covariance matrix are then pooled with **Rubin's rules**: T = Ū + (1 + 1/m)B, with Barnard–Rubin df and FMI (Burgess, White, Resche-Rigon & Wood, 2013). The pooled study results go into `rma()`.
- **Why not average the data.** The imputed *data* are not averaged. `PMM-avg` (the 20 datasets averaged into one and analysed once) is included only to show why: it hides the information lost and inflates ω.

## Step 3: the meta-analysis of each method (`04_evaluate_meta_analysis.R`)

- For each study, the APIMeM is estimated with `prep_study()`, `fit_apim()` and `extract_paths()` of `apim_meta_analysis.R`. This gives the uncorrected and the reliability-corrected components, with V = ML and MLR.
- The studies are pooled exactly as the `rma()` family of `apim_meta_sensitivity.R`:
  - per study, the Partner 1 and Partner 2 estimates of each effect are combined (GLS with their 2 × 2 sampling covariance block);
  - the composites are pooled with `metafor::rma()`, using random effects, REML and a Wald z test (primary).
  - The common-effect and Knapp–Hartung versions are kept for reference.
- On the input of the report, this code reproduces the report's observed estimates and SEs exactly (max \|diff\| = 6 × 10⁻¹⁶).
- Indirect effects are products of the pooled paths, tested with a Monte Carlo CI and the joint-significance test (the test of the report).

## Step 4: the sensitivity part

- **Power sensitivity (`07`).** `apim_meta_sensitivity.R` is evaluated unchanged with the settings stored in `APIM_sensitivity_results.rds`.
  - Its multilevel `rma.mv()` part is skipped, because only `rma()` is used now. That part is most of the run time, and it does not affect the `rma()` results or the random numbers.
  - Rerun on the report's own input, it **reproduces the report's MDEs and Type I error rates exactly**.
  - Rerun with the missing-data sampling covariance matrices, it gives the MDE, critical values and indirect-effect MDE **with the missing data**.
  - For every method and condition, the report's numbers are also rescaled by the SE inflation. This is exact when all sampling variances are inflated by the same factor, and the exact reruns check it.
- **MNAR sensitivity (`06`).** Delta-adjusted mice-PMM shifts the imputed item values by δ item-SDs in the hypothesized MNAR direction, with δ from −0.5 to 1. All δ values of a dataset use the same seeded random-number stream as the mice-PMM imputations of `03` (common random numbers), so δ = 0 reproduces them exactly and differences across δ reflect δ rather than imputation noise.
  - Tipping point: the smallest δ at which a pooled effect loses significance under `rma()`, i.e. falls below its critical value.
  - Calibration against the complete toy data shows which δ the simulated MNAR mechanism corresponds to.
- **Monte Carlo (`05`).** Separates method differences from the noise of one realisation. The PMM-vs-FIML decision uses it when it has been run in pilot or full mode.

## What the toy run shows

These are single realisations of the missing data, plus a Monte Carlo pilot (10 replications per condition). Run `05` in full mode before drawing firm conclusions.

**Missing data and imputation**

- The amputation worked as intended (realized selection models in Table 1 of the tables document).
- The mice chains converged: the largest PSRF was 1.18.
- Under MCAR and MAR, mice-PMM and mice-CAT reproduce the item means (\|standardized bias\| ≤ .02 SD).
- Under MNAR, every MAR-based method stays short of the MNAR direction, by 0.03 SD at 5% to 0.10 SD at 15%.
- The fraction of missing information of the study-level paths is .06, .10 and .15 on average at 5%, 10% and 15% (maximum .31).

**Missing-data methods in the meta-analysis**

- *The pipeline's `rma()` code* reproduces the observed results of the sensitivity report exactly (max \|diff\| 6 × 10⁻¹⁶).
- *Item-level methods with the auxiliaries* (FIML-TS, PMM, CAT) stay close to the complete-data meta-analysis. Their RMSD is 0.35–0.43 complete-data SEs under MCAR and MAR (single realisation).
- *PRO-FIML*, the default of `apim_meta_analysis.R`, ignores the auxiliaries and is clearly worse under MAR (0.59). In the Monte Carlo pilot its RMSD from the complete data is 1.6 times that of PMM under MAR. **For the real data, use an item-level method with auxiliaries rather than the default proration.**
- *Information loss.* The SEs of the pooled effects grow by about 2–3%, 4–6% and 7–11% at 5%, 10% and 15% missing items, most under MNAR (PRO-FIML up to 16%).
- *Significance decisions* agree with the complete data for 4–6 of the 6 direct effects and 2–4 of the 4 specific indirect effects. The disagreements concern the borderline effects (partner b and actor c′; p near .05) and the indirect effects that contain partner b.
- *Averaging the imputed datasets (PMM-avg)* gives point estimates close to Rubin's rules, but its SEs are on average 1–2% (up to 8%) *smaller* than with complete data, as if no information had been lost, and ω is inflated by up to .05. In the Monte Carlo pilot its CI coverage fell to .92–.93 under MCAR and MAR. The m = 20 imputations are therefore aggregated with Rubin's rules.

**Decision: mice-PMM versus FIML**

- *Single realisation.* FIML-TS was slightly closer to the complete data (RMSD 0.37 vs 0.40 SE units).
- *Monte Carlo pilot.* The two are equivalent: RMSD from the complete data 0.016 vs 0.016, CI coverage .949 vs .940, CI width ratio 1.06 for both.
- *Rule.* By the pre-specified rule (ties go to PMM) the method is **mice-PMM, m = 20, with Rubin's rules within each study before `rma()`**. `08` and `09` use this decision. Set `FINAL_METHOD=FIML_TS` to use FIML instead; the two give practically the same meta-analytic results.

**Power sensitivity (the sensitivity report under missing data)**

- *Reproduction.* `apim_meta_sensitivity.R`, rerun on its own input with `rma()` only, reproduces the report's MDEs and Type I error rates exactly (max \|diff\| = 0).
- *Inflation of MDEs and critical values* (mice-PMM, all I² levels): ×1.02–1.03 (5%), ×1.04–1.06 (10%) and ×1.07–1.11 (15%), largest under MNAR.
  - Example: the uncorrected actor a MDE at I² = 0% goes from .12 to .13 at 15%.
- *Indirect effects.* The MDEs of the indirect effects (products of two paths) grow by ×1.03–1.09, ×1.06–1.15 and ×1.11–1.23.
- *Check.* Exact reruns with the missing-data sampling covariance matrices (15% conditions) differ from the rescaled values by \< 1%.
- *Practical meaning.* 15% item-level missingness costs about as much power as losing 13–19% of the dyads.

**Departures from MAR (delta adjustment)**

- *Strong effects.* Actor a, partner a and actor b, and the actor–actor and partner–actor indirect effects, both components, stay significant up to δ = 1 item SD in the MNAR direction, in all nine conditions. No tipping point is found.
- *Borderline effects.* The uncorrected partner b and actor c′ (p ≈ .03–.10) switch between significant and not significant across δ and conditions; actor c′ tips at δ ≈ 0.25 under MAR 5–10%. Part of this is Monte Carlo noise of the imputations, so treat these effects as fragile whatever the missing-data method.
- *Calibration.* On the toy data, the completed item means match the complete data at δ ≈ 0 under MCAR/MAR and at δ ≈ 0.6–0.9 under the simulated MNAR mechanism. A δ range up to 1 item SD therefore covers it.

**Runtime of this run** (one core each): imputation \~40 min, evaluation \~20 min, delta analysis 84 min, power sensitivity 30 min, Monte Carlo pilot 75 min.

## Using this with the real data (`09_final_meta_analysis.R`)

1.  Put each study in one wide file (one row per dyad) with the items named `x_p1_1 … y_p2_4`, values 1–5 and missing values empty or NA. Set `INPUT_FILES`.
2.  Set `AUX` to always-observed variables that predict nonresponse or the scale scores, e.g. relationship duration, demographics or related scales. Without auxiliaries, use `AUX <- character(0)`; MAR is then assumed given the items only.
3.  The method is the one chosen on the toy data (`choose_method()`). You can force it with `FINAL_METHOD=mice_pmm` (or `FIML_TS`).
4.  Set `MNAR_DIRECTION` (in `00_functions.R`) to the direction you consider plausible for each construct. Preregister the δ range.
5.  The script writes the stacked imputed datasets, results (.xlsx/.rds) and APA tables, and `APIM_meta_analysis_results_MI.rds`, which `apim_meta_sensitivity.R` can read.
6.  The delta analysis (\~15 min) and the power rerun (\~5 min) run by default; switch them off with `Sys.setenv(RUN_DELTA = FALSE, RUN_POWER = FALSE)` (or 0/1) before sourcing the script.

## Repository, licence and contributors

* Repository: <https://github.com/jbuczny/APIM_meta-analysis_missing_data>
* Licence: Creative Commons Attribution-NonCommercial 4.0 International (`LICENSE.md`).
* Not tracked (see `.gitignore`): RStudio session files (`.Rproj.user/`, `.RData`, `.Rhistory`), system and Dropbox files, and the Monte Carlo smoke test (`output/mc/quick/`).
* Contributors: Jac ([@jbuczny](https://github.com/jbuczny)); code written with the assistance of Claude (Anthropic).

## Key references

Burgess, White, Resche-Rigon & Wood (2013), *Statistics in Medicine, 32*, 4499–4514 · Collins, Schafer & Kam (2001), *Psychological Methods, 6*, 330–351 · Eekhout, de Vet, de Boer, Twisk & Heymans (2018), *Statistical Methods in Medical Research, 27*, 1128–1140 · Eekhout et al. (2014), *Journal of Clinical Epidemiology, 67*, 335–342 · Enders (2010), *Applied Missing Data Analysis* · Fritz & MacKinnon (2007), *Psychological Science, 18*, 233–239 · Graham (2003), *Structural Equation Modeling, 10*, 80–100 · Heckman (1979), *Econometrica, 47*, 153–161 · Ledermann, Macho & Kenny (2011), *Structural Equation Modeling, 18*, 595–612 · Leacy, Floyd, Yates & White (2017), *American Journal of Epidemiology, 185*, 304–315 · Morris, White & Crowther (2019), *Statistics in Medicine, 38*, 2074–2102 · Savalei & Bentler (2009), *Structural Equation Modeling, 16*, 477–497 · Savalei & Rhemtulla (2017), *Journal of Educational and Behavioral Statistics, 42*, 405–431 · Tchetgen Tchetgen & Wirth (2017), *Biometrics, 73*, 1123–1131 · van Buuren (2018), *Flexible Imputation of Missing Data* (2nd ed.) · Viechtbauer (2010), *Journal of Statistical Software, 36*(3).
