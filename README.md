# Sepsis glycemic indices: primary association and machine-learning analyses

Analysis code for **Glycemic variability, stress hyperglycemia ratio and hemoglobin glycation index in relation to 28-day mortality in sepsis**.

This repository contains two selected analysis modules:

| Module | Contents |
| --- | --- |
| `01_cox_rcs/` | Separate MIMIC-IV and NWICU Cox Models 1–3, continuous and quartile contrasts, and restricted cubic splines |
| `02_machine_learning/` | Training feature selection, five algorithms, internal/external evaluation, calibration, decision curves and XGBoost TreeSHAP |

Patient data, individual predictions, fitted models, database extraction credentials and workstation configuration are not included. This is a selected code release; raw-database extraction and the remaining supplementary analyses are outside its scope.

## Study specification

- One first ICU admission per patient, adult age >18 years, eligible sepsis cohort.
- MIMIC-IV: 5,195 patients and 751 recorded 28-day events. NWICU: 1,140 patients and 260 recorded 28-day events.
- The original `day28_outcome` is retained. Core survival follow-up `day28_los` is supplied in **hours**, already truncated at 28 days, and converted to days by the Cox runner. Hospital mortality is a separate outcome and is not substituted.
- GV is the whole-ICU glucose coefficient of variation (%). SHR uses early chemistry glucose and HbA1c; HGI is cohort-referenced for association analyses. The ML HGI equation is estimated from its training data and fixed for validation.
- Cox Model 1 is unadjusted; Model 2 adjusts for age, sex and race; Model 3 additionally adjusts for diabetes, SOFA and the first Charlson score within 24 hours. No missing core-model values are imputed.
- Continuous HRs are per 10 percentage points of GV, 0.1 SHR and one percentage point of HGI. Quartiles are cohort-specific, with Q1 as reference. Quartile omnibus P values are not reported.
- RCS uses three cohort-specific knots (10th, 50th and 90th percentiles), Model 3 adjustment and the median as reference. Overall and nonlinear tests receive separate within-cohort Holm corrections.
- ML uses a patient-level, outcome-stratified MIMIC-IV 70:30 split, with 3,635 training and 1,560 internal-validation patients. Feature selection, tuning, preprocessing, Platt calibration and thresholds derive from training; NWICU is evaluated without local refitting.

Whole-ICU GV is unavailable at ICU admission. These analyses should not be interpreted as a validated admission-time clinical prediction tool or as evidence of causal treatment effects.

## Data access and preparation

Authorized researchers must obtain MIMIC-IV v3.1 and NWICU v0.1.0 through the applicable PhysioNet credentialing, training and data-use processes:

- [MIMIC-IV v3.1](https://physionet.org/content/mimiciv/3.1/)
- [NWICU v0.1.0](https://physionet.org/content/nwicu-northwestern-icu/0.1.0/)

Prepare the curated first-patient cohort CSVs in a private location **outside this repository**. The two files in `input_templates/` have headers only and zero patient rows. They document the combined input schema for the two modules; they are not example patient datasets. Field details are in each module's README and `02_machine_learning/input_fields.json`. Use the final NWICU export with the four raw-record re-extracted first-day vital signs.

The runners enforce the study cohort sizes and event counts. They assume that eligibility and biomarker derivation have already been completed; they do not reconstruct a complete Sepsis-3 phenotype from raw hospital tables.

## Environment

The analysis used R 4.4.2. `requirements-R.txt` lists the direct package requirements; `software_environment.json` records observed package versions and is not a complete dependency lockfile. Use a UTF-8-capable environment. Times New Roman is required to reproduce figure typography.

Install the direct R packages if needed:

```r
packages <- readLines("requirements-R.txt", warn = FALSE)
install.packages(packages[nzchar(packages)])
```

For the original XGBoost interface, use the recorded package version when a newer release changes its API. Runtime fitting checks stop when expected features, cohort counts or calibration invariants do not match.

## Run

From the repository root, using quoted private paths:

```sh
Rscript 01_cox_rcs/run_cox_rcs.R "/private/MIMIC-IV.csv" "/private/NWICU.csv" "/private/cox_results"
Rscript 02_machine_learning/run_machine_learning.R "/private/MIMIC-IV.csv" "/private/NWICU.csv" "/private/ml_results"
```

The output parent directory must exist. Keep inputs and all generated outputs outside the repository. The Cox module produces aggregate result tables and six individual PNG spline panels. The ML module produces local feature-selection, tuning, performance, calibration, decision-curve and SHAP outputs, including private patient-level predictions and model checkpoints. Never upload those patient-level outputs or checkpoints.

ML training can take hours. Use a fresh output directory when changing inputs or analytical settings; existing completed checkpoints may be resumed only for the unchanged analysis.

## Release validation

The reorganized Cox/RCS module was run against the original private final-cohort inputs. Its 90 Cox table rows, six RCS tests, 1,806 curve points and 180 histogram bins matched the authoritative analysis outputs numerically, with maximum absolute differences of zero in the compared numerical fields.

The ML release passed R syntax and CLI/input-boundary checks. Its statistical, preprocessing and SHAP functions were compared with the original source. Full ML feature selection and training were **not rerun during release preparation**. This repository does not claim fresh end-to-end reproducibility from raw database installation.

## Citation and archival status

Please cite the associated manuscript when published. A repository URL or GitHub commit is not a DOI. No DOI is claimed in this release; a reviewed release can subsequently be archived in a DOI-assigning repository.
