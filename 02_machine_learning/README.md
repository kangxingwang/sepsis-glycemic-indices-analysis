# Five-algorithm machine learning

Run from a terminal with R 4.4.2 and the packages in `requirements-R.txt` installed:

```text
Rscript 02_machine_learning/run_machine_learning.R PRIVATE_MIMIC.csv PRIVATE_NWICU.csv PRIVATE_OUTPUT_DIR
```

Use quoted absolute paths when paths contain spaces. Both input CSVs and the output directory must be outside this public repository. The output directory's parent must already exist; the entrypoint creates its analysis subdirectories. Input templates in the repository root contain headers only and cannot be used as patient data. `input_fields.json` lists all six mandatory identity/outcome/derivation fields and all 49 candidate predictor source fields. Numeric glucose is mg/dL and HbA1c is percent. Identifiers remain character values; each input contains one ICU admission per unique patient.

The study assertions require 5,195 MIMIC-IV patients/751 original 28-day events and 1,140 NWICU patients/260 events. MIMIC-IV is outcome-stratified into 3,635 training patients/525 events and 1,560 internal validation patients/226 events using seed 1117. NWICU uses the final raw-record re-extracted first-day vital signs. The original `day28_outcome` is used throughout; hospital mortality is a different outcome and is not substituted here.

## Methods preserved from the actual analysis

Candidate fields must be present in both datasets and have <=20% missingness in each, with the four prespecified vital-sign exceptions retained. Missingness/availability screening uses predictor availability in both cohorts; NWICU outcomes are not used for training or tuning. Forty-three study candidates remain. Numeric median imputation, categorical mode imputation, a training HbA1c-on-glucose HGI equation, first/99th percentile clipping of vital/laboratory predictors, categorical encoding and scaling are estimated using the training data only. These same steps, including HGI estimation and feature selection, are repeated within each training CV fold.

Feature selection uses the intersection of binomial LASSO (`lambda.1se`, 10-fold binomial deviance) and Boruta Confirmed variables (p=.01, multiplicity adjustment, maximum 200 iterations, 500 trees per iteration). Tentative Boruta variables are excluded, and no glucose index is forced into the model. The final intersection contains 25 original clinical variables, including GV and SHR; HGI was not selected.

The five algorithms are logistic regression, XGBoost, random forest, radial-kernel SVM and k-nearest neighbors. Hyperparameters and the preferred algorithm are selected by mean training-fold Average Precision over 10 folds. Random search uses 30 XGBoost settings and 20 each for random forest, SVM and k-NN; logistic regression has one specification. The final random forest uses 1,500 trees. All final models are trained only on the 70% MIMIC-IV training set, then applied unchanged to the internal and NWICU cohorts.

Platt calibration and an actual finite training-OOF Youden maximizer are estimated from training OOF predictions and fixed for evaluation. Among tied maximizers, the lower central returned threshold is used. Evaluation calibration intercepts/slopes are diagnostic quantities and do not refit NWICU predictions. AUROC intervals use DeLong; AP and Brier intervals use 1,000 outcome-stratified bootstrap resamples of fixed evaluation predictions, without retraining. Pairwise AUROC comparisons use paired DeLong tests with Holm adjustment within cohort. Exact TreeSHAP explains the prespecified XGBoost model, grouping dummy contributions at patient level before ranking original variables. Raw and grouped additivity are checked.

## Files and private outputs

- `feature_engineering.R`: data validation, locked split, training feature engineering and final LASSO/Boruta intersection.
- `train_validate_shap.R`: fold-local processing, five-algorithm tuning/training, fixed calibration/thresholds, evaluation and SHAP.
- `plots_and_audits.R`: independent training-isolation, prediction, threshold, metric and SHAP checks, followed by PNG selection/ROC/calibration/DCA/SHAP plots.
- `run_machine_learning.R`: input/privacy boundaries, directories and execution order.

All generated CSV, RDS, fitted models, individual predictions, individual SHAP values and plot data are private outputs. They must never be uploaded to the code repository. Plots are PNG only. Generated source filenames retain UTF-8 Chinese labels from the actual analysis; use a UTF-8-capable R environment. Presentation styling here is an analysis output and does not replace the separately assembled submission figures.

An unchanged private output directory can resume completed checkpoints. Use a fresh output directory for a complete rerun, especially after changing source code, preprocessing, selection rules or input data. Full feature selection and training can take hours. This release preparation performed syntax and CLI/input-boundary checks only; it did not rerun the full ML analysis. Package versions used in the original analysis are documented in the repository environment file. The code does not include patient data or fitted patient models.
