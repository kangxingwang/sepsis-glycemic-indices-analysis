#!/usr/bin/env Rscript
# Run from any directory. Input CSVs and output directory are supplied explicitly.
# This module contains selected primary Cox and restricted cubic spline analyses.
options(stringsAsFactors = FALSE, width = 140)
if (.Platform$OS.type == "windows") {
  try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE)
}
suppressPackageStartupMessages({
  library(survival)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript run_cox_rcs.R PRIVATE_MIMIC.csv PRIVATE_NWICU.csv PRIVATE_OUTPUT_DIR")
}
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this entry point with Rscript.")
module_dir <- dirname(normalizePath(sub("^--file=", "", script_arg),
                                    winslash = "/", mustWork = TRUE))
source(file.path(module_dir, "cox_functions.R"), encoding = "UTF-8")
repository_dir <- normalizePath(dirname(module_dir), winslash = "/", mustWork = TRUE)
inside_repository <- function(path) {
  p <- sub("/+$", "", path); r <- sub("/+$", "", repository_dir)
  if (.Platform$OS.type == "windows") { p <- tolower(p); r <- tolower(r) }
  identical(p, r) || startsWith(p, paste0(r, "/"))
}
inputs <- setNames(vapply(args[1:2], normalizePath, character(1),
                         winslash = "/", mustWork = TRUE), c("MIMIC-IV", "NWICU"))
if (identical(unname(inputs[1]), unname(inputs[2]))) stop("Cohort inputs must be different files.")
if (any(vapply(inputs, inside_repository, logical(1)))) stop("Patient inputs must be outside this repository.")
out_arg <- path.expand(args[3])
if (dir.exists(out_arg)) {
  out_dir <- normalizePath(out_arg, winslash = "/", mustWork = TRUE)
} else {
  parent <- dirname(out_arg)
  if (!dir.exists(parent)) stop("The private output parent directory must exist.")
  out_dir <- file.path(normalizePath(parent, winslash = "/", mustWork = TRUE), basename(out_arg))
}
if (inside_repository(out_dir)) stop("Generated outputs must be outside this repository.")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create output directory.")
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

core <- c("Age", "Gender", "Race", "Diabetes", "SOFA", "Charlson")
exposures <- c("GV", "SHR", "HGI")
units <- c(GV = 10, SHR = 0.1, HGI = 1)
unit_labels <- c(GV = "Per 10 percentage points", SHR = "Per 0.1",
                 HGI = "Per 1 percentage point")
models <- list(Model1 = character(), Model2 = c("Age", "Gender", "Race"),
               Model3 = core)
mapping <- c(stay_id = "stay_id", subject_id = "subject_id", age = "Age",
             gender = "Gender", race = "Race", icd2_diabetes = "Diabetes",
             sofa = "SOFA", charlson = "Charlson", gv_cv_icu_mixed = "GV",
             shr_first24h_chem = "SHR", hgi_first24h_chem = "HGI",
             day28_los = "Day28_LOS", day28_outcome = "Day28_Outcome")
write_csv <- function(x, filename) {
  write.csv(x, file.path(out_dir, filename), row.names = FALSE, na = "",
            fileEncoding = "UTF-8")
}
repeat_fit <- function(fit) rep(list(fit), 20L)
# These complete-case fits use the historical finite-df inference in
# cox_functions.R. Repeated copies preserve that calculation; no imputation occurs.
rcs_basis <- function(x, k) {
  h <- (pmax(x-k[1], 0)^3 -
          pmax(x-k[2], 0)^3*(k[3]-k[1])/(k[3]-k[2]) +
          pmax(x-k[3], 0)^3*(k[2]-k[1])/(k[3]-k[2]))/(k[3]-k[1])^2
  cbind(Exposure = x, RCS_N1 = h)
}
read_cohort <- function(path, database) {
  raw <- read.csv(path, fileEncoding = "UTF-8-BOM", check.names = FALSE,
                  colClasses = "character", na.strings = c("", "NA", "NULL"))
  absent <- setdiff(names(mapping), names(raw))
  if (length(absent)) stop(database, ": missing columns: ", paste(absent, collapse=", "))
  d <- raw[names(mapping)]
  names(d) <- unname(mapping)
  numeric_fields <- c("Age", "Diabetes", "SOFA", "Charlson", exposures,
                      "Day28_LOS", "Day28_Outcome")
  for (field in numeric_fields) {
    d[[field]] <- suppressWarnings(as.numeric(d[[field]]))
    if (any(!is.finite(d[[field]]))) stop(database, ": non-finite ", field)
  }
  if (anyNA(d) || any(trimws(d$stay_id) == "") ||
      any(trimws(d$subject_id) == "")) stop(database, ": incomplete core input.")
  if (anyDuplicated(d$subject_id) || anyDuplicated(d$stay_id)) {
    stop(database, ": each subject and stay must occur once.")
  }
  sex <- toupper(trimws(d$Gender))
  if (!all(sex %in% c("F", "FEMALE", "M", "MALE"))) {
    stop(database, ": gender must be F/FEMALE or M/MALE.")
  }
  if (!all(d$Diabetes %in% c(0,1)) || !all(d$Day28_Outcome %in% c(0,1))) {
    stop(database, ": diabetes and day28_outcome must be binary 0/1.")
  }
  d$Gender <- factor(ifelse(sex %in% c("F","FEMALE"), "Female", "Male"),
                      levels = c("Female", "Male"))
  d$Race <- factor(trimws(d$Race))
  if (!"White" %in% levels(d$Race) || any(d$Race == "")) {
    stop(database, ": race requires the White reference category.")
  }
  d$Race <- relevel(d$Race, "White")
  d$Diabetes <- factor(ifelse(d$Diabetes == 1, "Yes", "No"),
                        levels = c("No", "Yes"))
  # The analysis export stores follow-up in HOURS, already truncated at 28 days.
  # Retain the supplied event flag; do not replace it by hospital/ICU mortality.
  d$Day28_LOS <- d$Day28_LOS/24
  expected <- if (database == "MIMIC-IV") c(5195L, 751L) else c(1140L, 260L)
  if (nrow(d) != expected[1] || sum(d$Day28_Outcome) != expected[2]) {
    stop(database, ": cohort/event counts differ from the locked analysis.")
  }
  if (!all(d$Day28_LOS > 0 & d$Day28_LOS <= 28) || !all(d$Age > 18)) {
    stop(database, ": invalid adult age or 28-day follow-up; follow-up must be in hours.")
  }
  hgi_check <- data.frame(Database = database, Raw_pair_available = FALSE,
                          Intercept = NA_real_, Slope = NA_real_,
                          Maximum_absolute_HGI_difference = NA_real_)
  raw_pair <- c("hba1c", "glucose_first24h_chem")
  if (any(raw_pair %in% names(raw))) {
    if (!all(raw_pair %in% names(raw))) {
      stop(database, ": supply both hba1c and glucose_first24h_chem for HGI verification.")
    }
    hba1c <- suppressWarnings(as.numeric(raw$hba1c))
    glucose <- suppressWarnings(as.numeric(raw$glucose_first24h_chem))
    if (!all(is.finite(hba1c)) || !all(is.finite(glucose))) {
      stop(database, ": non-finite HGI calibration inputs.")
    }
    equation <- coef(lm(hba1c ~ glucose))
    calculated <- hba1c-equation[1]-equation[2]*glucose
    error <- max(abs(d$HGI-calculated))
    if (!is.finite(error) || error > 1e-10) {
      stop(database, ": HGI differs from the cohort-specific HbA1c-on-glucose OLS residual.")
    }
    hgi_check <- data.frame(Database = database, Raw_pair_available = TRUE,
                            Intercept = unname(equation[1]), Slope = unname(equation[2]),
                            Maximum_absolute_HGI_difference = error)
  }
  cuts <- list()
  for (e in exposures) {
    q <- unname(quantile(d[[e]], c(.25, .5, .75)))
    if (anyDuplicated(q)) stop(database, ": repeated quartile boundaries for ", e)
    d[[paste0(e, "_Q")]] <- cut(d[[e]], c(-Inf, q, Inf), labels = paste0("Q",1:4))
    cuts[[e]] <- data.frame(Database=database, Exposure=e, Q25=q[1], Q50=q[2], Q75=q[3])
  }
  list(Data = d, HGI_check = hgi_check, Quartiles = do.call(rbind, cuts))
}

main_results <- rcs_results <- rcs_curves <- rcs_histograms <- list()
cohort_summaries <- hgi_checks <- quartiles <- list()
for (database in names(inputs)) {
  message("Cox/RCS: ", database)
  prepared <- read_cohort(inputs[[database]], database)
  d <- prepared$Data
  hgi_checks[[database]] <- prepared$HGI_check
  quartiles[[database]] <- prepared$Quartiles
  cohort_summaries[[database]] <- data.frame(
    Database=database, N=nrow(d), Events=sum(d$Day28_Outcome),
    Unique_patients=length(unique(d$subject_id)),
    Unique_stays=length(unique(d$stay_id)), Missing_core_values=sum(is.na(d)),
    Followup_unit="days after converting input hours", New_imputation=FALSE)
  rows <- list()
  for (e in exposures) for (form in c("continuous", "quartile")) {
    z <- d
    z$Exposure <- z[[e]]/units[e]
    z$Exposure_Q <- z[[paste0(e, "_Q")]]
    term <- if (form == "continuous") "Exposure" else "Exposure_Q"
    for (model in names(models)) {
      fit <- fit_cox(z, c(term, models[[model]]))
      result <- cox_table(repeat_fit(fit), e, form, model, "Day 0-28")
      # Only individual coefficient P values are reported for quartile contrasts.
      # Quartile omnibus P values are outside the reporting scope.
      result$P_overall <- NULL
      result$Database <- database
      rows[[length(rows)+1L]] <- result
    }
  }
  main_results[[database]] <- do.call(rbind, rows)
  tests <- curves <- histograms <- list()
  for (e in exposures) {
    x <- d[[e]]
    knots <- unname(quantile(x, c(.1, .5, .9)))
    reference <- median(x)
    basis <- rcs_basis(x/units[e], knots/units[e])
    z <- d
    z$Exposure <- basis[,1]
    z$RCS_N1 <- basis[,2]
    fit <- fit_cox(z, c("Exposure", "RCS_N1", core))
    pooled <- pool_mi(repeat_fit(fit))
    tests[[e]] <- data.frame(
      Database=database, Exposure=e, N=nrow(d), Events=sum(d$Day28_Outcome),
      Knot10=knots[1], Knot50=knots[2], Knot90=knots[3], Reference=reference,
      P_overall=joint_p(pooled, c("Exposure", "RCS_N1")),
      P_nonlinear=joint_p(pooled, "RCS_N1"))
    displayed <- unname(quantile(x, c(.01, .99)))
    xx <- sort(unique(c(seq(displayed[1], displayed[2], length.out=300), reference)))
    bb <- rcs_basis(xx/units[e], knots/units[e])
    b0 <- rcs_basis(reference/units[e], knots/units[e])
    contrast <- matrix(0, length(xx), length(pooled$q))
    contrast[, match(c("Exposure","RCS_N1"), pooled$names)] <-
      sweep(bb, 2, as.numeric(b0), "-")
    curves[[e]] <- cbind(Database=database, Exposure=e, Value=xx,
                         pool_contrast(pooled, contrast))
    h <- hist(x[x >= displayed[1] & x <= displayed[2]],
               breaks=seq(displayed[1], displayed[2], length.out=31),
               plot=FALSE, include.lowest=TRUE)
    histograms[[e]] <- data.frame(
      Database=database, Exposure=e, Left=head(h$breaks,-1),
      Right=tail(h$breaks,-1), Count=h$counts,
      Excluded_left=sum(x < displayed[1]), Excluded_right=sum(x > displayed[2]))
  }
  tests <- do.call(rbind, tests)
  tests$P_overall_Holm <- p.adjust(tests$P_overall, "holm")
  tests$P_nonlinear_Holm <- p.adjust(tests$P_nonlinear, "holm")
  rcs_results[[database]] <- tests
  rcs_curves[[database]] <- do.call(rbind, curves)
  rcs_histograms[[database]] <- do.call(rbind, histograms)
}
main <- do.call(rbind, main_results)
tests <- do.call(rbind, rcs_results)
curves <- do.call(rbind, rcs_curves)
histograms <- do.call(rbind, rcs_histograms)
write_csv(main, "cox_models.csv")
write_csv(tests, "rcs_tests.csv")
write_csv(curves, "rcs_curves.csv")
write_csv(histograms, "rcs_histograms.csv")
write_csv(do.call(rbind, cohort_summaries), "cohort_summary.csv")
write_csv(do.call(rbind, hgi_checks), "hgi_verification.csv")
write_csv(do.call(rbind, quartiles), "cohort_quartile_cutpoints.csv")
capture.output(sessionInfo(), file=file.path(out_dir, "sessionInfo.txt"))
# No individual records, fitted model objects or predictions are written.
source(file.path(module_dir, "plot_rcs.R"), encoding="UTF-8")
plot_rcs_panels(curves, tests, histograms, file.path(out_dir, "figures"))
write_json(list(
  status="complete", scope="primary Cox Model1-3 and Model3 RCS",
  cohorts=c("MIMIC-IV","NWICU"), cox_rows=nrow(main), rcs_tests=nrow(tests),
  individual_quartile_contrasts=TRUE, quartile_omnibus_P_reported=FALSE,
  new_imputation_performed=FALSE, individual_records_written=FALSE,
  fitted_objects_written=FALSE, R_version=as.character(getRversion())
), file.path(out_dir, "run_summary.json"), auto_unbox=TRUE, pretty=TRUE)
message("Completed primary Cox/RCS analyses. Aggregate results: ", out_dir)

