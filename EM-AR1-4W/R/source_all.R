.cursor <- normalizePath(getwd(),winslash="/",mustWork=TRUE); .project <- NULL
repeat {
  if (file.exists(file.path(.cursor,"EM-AR1-4W","R","ar1_4w.R"))) { .project <- .cursor; break }
  .parent <- dirname(.cursor); if (identical(.parent,.cursor)) break; .cursor <- .parent
}
if (is.null(.project)) stop("Cannot locate EM-AR1-4W module")
source(file.path(.project,"EM-AR1-4W","R","ar1_4w.R"),local=FALSE)
source(file.path(.project,"EM-AR1-4W","R","fmm_ar1_4w.R"),local=FALSE)
source(file.path(.project,"EM-AR1-4W","R","fmm_common_transitions_type_error_4w.R"),local=FALSE)
covinc_file <- file.path(.project,"EM-AR1-4W","R","fmm_covariates_inconsistency_4w.R")
if (file.exists(covinc_file)) source(covinc_file,local=FALSE)
rm(list=intersect(c(".cursor",".project",".parent"),ls()))
