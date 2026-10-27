# Primary Cox and restricted cubic spline analyses

This module contains the selected primary association analyses for GV, SHR and
HGI in the MIMIC-IV and NWICU sepsis cohorts. The two cohorts are analyzed
separately. It does not include raw database extraction or supplementary
landmark, hospital/ICU Logistic, subgroup or Clinical5 analyses.

## Run

Install R 4.4.2 and the packages **survival, jsonlite, ggplot2, patchwork, ragg**,
then run from any directory:

```sh
Rscript 01_cox_rcs/run_cox_rcs.R "/private/MIMIC-IV.csv" "/private/NWICU.csv" "/private/cox_results"
```

Use private local paths for inputs and generated outputs. No patient data are
distributed. The module writes aggregate tables and six individual PNG panels;
it does not export individual observations, fitted model objects or predictions.

## Input columns

Each CSV must contain these 13 columns:

| Column | Definition / encoding |
| --- | --- |
| stay_id | ICU stay identifier, unique within the cohort |
| subject_id | Patient identifier, unique within the cohort |
| age | Age in years, >18 |
| gender | F/M or Female/Male |
| race | Race category, including White as the reference |
| icd2_diabetes | Diabetes indicator, 0/1 |
| sofa | SOFA score |
| charlson | First available Charlson score within the first 24 hours |
| gv_cv_icu_mixed | Whole-ICU glucose coefficient of variation, expressed as a percentage |
| shr_first24h_chem | Early stress hyperglycemia ratio |
| hgi_first24h_chem | Cohort-referenced HGI, in HbA1c percentage points |
| day28_los | Follow-up duration in **hours**, already truncated at 28 days |
| day28_outcome | Original recorded 28-day event indicator, 0/1 |

The core input must be complete. Duplicate patients/stays, non-finite numeric
values and unexpected cohort/event counts stop execution. The locked counts are
5,195 patients / 751 events for MIMIC-IV and 1,140 / 260 for NWICU. Primary
`day28_outcome` is used unchanged; it is not replaced by hospital or ICU mortality.

Optional columns `hba1c` and `glucose_first24h_chem` must be supplied together.
When present, the code checks that supplied HGI equals the residual from the
cohort-specific OLS regression of HbA1c on first-24-hour chemistry glucose
(mg/dL), within an absolute tolerance of 1e-10. It does not replace supplied HGI
or use a MIMIC-IV equation for the primary NWICU association analysis.

## Statistical specification

- Model 1: unadjusted.
- Model 2: age, sex and race.
- Model 3: age, sex, race, diabetes, SOFA and Charlson.
- Separate continuous and cohort-specific quartile models are fitted for each
  biomarker. Continuous HRs are per 10 percentage points of GV, 0.1 SHR and one
  percentage point of HGI; Q1 is the quartile reference.
- Cox models use Efron ties and robust variance clustered by the unique ICU stay.
  The shared inference functions preserve the historical finite-degree-of-freedom
  confidence intervals and coefficient tests. Twenty identical copies of each
  complete-data fit are used only to preserve this calculation. There is no new
  multiple imputation and no added information from these repeated fit copies.
- The Cox table reports individual coefficient/contrast P values. Quartile
  omnibus P values are not included.
- Model 3 restricted cubic splines use three knots at the cohort-specific 10th,
  50th and 90th percentiles and the exposure median as reference. Holm corrections
  are applied separately to the three overall spline tests and the three
  nonlinearity tests within each cohort.
- Spline plots display the 1st–99th percentile range. This is a display limit,
  not a sample exclusion; spline models use the full cohort. Each panel includes
  its own exposure histogram. The two cohort panels for a biomarker share the
  same HR-axis limits.

## Outputs

- `cox_models.csv`: continuous and individual quartile contrasts, Model 1–3.
- `rcs_tests.csv`: knot/reference values, overall and nonlinear spline tests.
- `rcs_curves.csv` / `rcs_histograms.csv`: aggregate plot data.
- `cohort_summary.csv`, `cohort_quartile_cutpoints.csv`,
  `hgi_verification.csv`: aggregate input checks.
- `figures/RCS_*.png`: six individual 600-dpi panels, Times New Roman.
- `sessionInfo.txt` / `run_summary.json`: runtime information.

The selected module assumes an already curated, eligible first-patient cohort
and precomputed GV/SHR/HGI. It is not an implementation of Sepsis-3 phenotype
curation or a claim that all manuscript supplementary analyses are included.

