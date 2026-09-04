# Update only the two new Table 9 rows in Table 8. Existing estimates unchanged.
source("EM-tenure/R/source_all.R")
source("EM-tenure/R/four_wave_duration_implications.R")
path <- "EM-AR2/output/results/set4_reliability/apparent_transition_decomposition.csv"
before <- read.csv(path,check.names=FALSE)
new <- build_duration_implications_4w()
old <- before[!before$model %in% new$model,]
stopifnot(nrow(old)==27L,identical(names(old),names(new)))
updated <- rbind(old,new)
stopifnot(nrow(updated)==33L,
  !anyDuplicated(updated[c("direction","model")]),
  max(abs(rowSums(updated[c("classification_only","true_reversal","true_persistent")])-1))<1e-8)
write.csv(updated,path,row.names=FALSE)
out <- "paper/generated/four_wave_duration"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(new,file.path(out,"table8_duration_rows.csv"),row.names=FALSE)
write.csv(updated,file.path(out,"table8_all_models.csv"),row.names=FALSE)
inputs <- attr(new,"input_md5")
write.csv(data.frame(path=names(inputs),md5=unname(inputs)),
  file.path(out,"table8_input_provenance.csv"),row.names=FALSE)
# Execute just the authoritative Table 8 chunk; do not rebuild other tables.
document <- readLines("paper/branch_changes_summary_short.Rmd",warn=FALSE)
start <- grep("^```\\{r apparent-transition-table,",document)
stopifnot(length(start)==1L)
end <- start + which(document[seq.int(start+1L,length(document))]=="```")[1L]
render_env <- new.env(parent=globalenv())
render_env$project_root <- normalizePath(".",winslash="/")
tex <- capture.output(eval(parse(text=document[seq.int(start+1L,end-1L)]),
  envir=render_env))
stopifnot(sum(grepl("Table 9 (1)",tex,fixed=TRUE))==3L,
          sum(grepl("Table 9 (2)",tex,fixed=TRUE))==3L,
          !any(grepl(" NA ",tex,fixed=TRUE)))
writeLines(tex,file.path(out,"table8_all_models.tex"))
display <- updated
display[4:6] <- lapply(display[4:6],function(x) sprintf("%.2f",100*x))
markdown <- c("# Table 8: Implications for apparent transitions", "",
  "| Direction | Model | Classification only (%) | Genuine reversal (%) | Genuine persistence (%) |",
  "|:--|:--|--:|--:|--:|",
  vapply(seq_len(nrow(display)),function(i) paste0("| ",
    paste(display[i,c(1,2,4,5,6)],collapse=" | ")," |"),character(1)),
  "", "Table 7 and Table 9 rows use full-duration posteriors; other rows condition on the status pair only. Models use their own estimation samples. Classification only is not the unconditional status-error probability. See the LaTeX table for full notes.")
writeLines(markdown,file.path(out,"table8_all_models.md"))
print(transform(new,classification_only=100*classification_only,
  true_reversal=100*true_reversal,true_persistent=100*true_persistent),row.names=FALSE)
cat("Preserved all 27 existing rows; added 6 rows for the two duration models.\n")
