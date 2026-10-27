# Five-algorithm sepsis glycemia analysis; code-only publication entrypoint.
# Usage: Rscript run_machine_learning.R PRIVATE_MIMIC.csv PRIVATE_NW.csv PRIVATE_OUTPUT_DIR
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) stop("Usage: Rscript run_machine_learning.R <private_MIMIC_CSV> <private_NWICU_CSV> <private_output_dir>", call.=FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
if (length(script_arg) != 1L) stop("Run this entrypoint using Rscript.", call.=FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), winslash="/", mustWork=TRUE))
repository_dir <- normalizePath(dirname(script_dir), winslash="/", mustWork=TRUE)
inside_repository <- function(path) {
  p <- sub("/+$", "", path); r <- sub("/+$", "", repository_dir)
  if (.Platform$OS.type == "windows") {p <- tolower(p); r <- tolower(r)}
  identical(p,r) || startsWith(p,paste0(r,"/"))
}
for (i in 1:2) if (!file.exists(args[i]) || dir.exists(args[i])) stop("Input CSV does not exist: ", args[i], call.=FALSE)
inputs <- vapply(args[1:2], normalizePath, character(1), winslash="/", mustWork=TRUE)
if (identical(unname(inputs[1]),unname(inputs[2]))) stop("MIMIC-IV and NWICU input CSVs must be different files.", call.=FALSE)
if (any(vapply(inputs,inside_repository,logical(1)))) stop("Patient input CSVs must be stored outside this public repository.", call.=FALSE)
# Resolve a real existing parent before creating an output directory, including junctions.
out_arg <- args[3]
if (dir.exists(out_arg)) {
  output_dir <- normalizePath(out_arg,winslash="/",mustWork=TRUE)
} else {
  if (file.exists(out_arg)) stop("Output path names a file.",call.=FALSE)
  parent <- dirname(out_arg)
  if (!dir.exists(parent)) stop("The parent of the requested output directory must already exist.",call.=FALSE)
  parent <- normalizePath(parent,winslash="/",mustWork=TRUE)
  leaf <- basename(out_arg)
  if (!nzchar(leaf) || leaf %in% c(".","..")) stop("Invalid output directory name.",call.=FALSE)
  output_dir <- file.path(parent,leaf)
}
if (inside_repository(output_dir)) stop("Generated patient outputs must be stored outside this public repository.",call.=FALSE)
required <- c("stay_id","subject_id","hadm_id","day28_outcome","hba1c","glucose_first24h_chem")
for (f in inputs) {
  header <- names(read.csv(f,nrows=0,check.names=FALSE,fileEncoding="UTF-8-BOM"))
  missing <- setdiff(required,header)
  if (length(missing)) stop("Input CSV lacks required fields: ",paste(missing,collapse=", "),call.=FALSE)
  if (length(readLines(f,n=2L,warn=FALSE,encoding="UTF-8")) < 2L) stop("Input CSV contains a header only; authorized patient data are required.",call.=FALSE)
}
packages <- scan(file.path(script_dir,"requirements-R.txt"),what=character(),quiet=TRUE)
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop("Install required R packages before running: ",paste(missing,collapse=", "),call.=FALSE)
for (d in c("analysis_notes","results","audit","private_model_inputs","figures","private_checkpoints","private_plot_data")) dir.create(file.path(output_dir,d),recursive=TRUE,showWarnings=FALSE)
Sys.setenv(SEPSIS_ML_MIMIC_CSV=inputs[1],SEPSIS_ML_NW_CSV=inputs[2],SEPSIS_ML_OUTPUT=output_dir)
if (.Platform$OS.type == "windows") try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
Sys.setenv(OMP_NUM_THREADS="4",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
options(warn=1)
fe_file <- file.path(output_dir,"private_checkpoints/特征工程结果.rds")
reuse_fe <- FALSE
if (file.exists(fe_file)) {
  cached_fe <- readRDS(fe_file)
  reuse_fe <- identical(unname(cached_fe$paths),unname(inputs)) && identical(unname(cached_fe$source_md5),unname(tools::md5sum(inputs))) && cached_fe$settings$seed==1117L && cached_fe$settings$boruta_runs==200L && cached_fe$settings$boruta_trees==500L && all(c("MAP","Heart_Rate","Resp_Rate","Temperature") %in% cached_fe$candidates)
}
if (reuse_fe) cat("Reuse verified formal feature-engineering checkpoint\n") else source(file.path(script_dir,"feature_engineering.R"),encoding="UTF-8")
source(file.path(script_dir,"train_validate_shap.R"),encoding="UTF-8")
cat("FORMAL_ML_TRAINING_COMPLETED\n")
source(file.path(script_dir,"plots_and_audits.R"),encoding="UTF-8")
cat("FORMAL_FULL_R_WORKFLOW_COMPLETED\n")
