# Table 8 definitions, intentionally different from "not confirmed in the
# reported direction" in summarise_four_wave_ar2.R. Opposite-direction latent
# changes belong to reversal/persistence, not classification_only.
posterior_duration_implications_4w <- function(df, fit, model) {
  h <- latent_histories_eps_4w()
  gamma <- fit$gamma
  stopifnot(nrow(gamma) == nrow(df), ncol(gamma) == nrow(h),
    all(is.finite(gamma)), all(gamma >= 0),
    max(abs(rowSums(gamma)-1)) < 1e-8,
    all(is.finite(df$weight)), all(df$weight >= 0))
  directions <- c("All apparent transitions", "Apparent entries", "Apparent exits")
  numerator <- matrix(0,3,3,dimnames=list(directions,
    c("classification_only","true_reversal","true_persistent")))
  denominator <- setNames(numeric(3),directions)
  for (t in 1:2) {
    masks <- cbind(
      classification_only=h[,t]==h[,t+1],
      true_reversal=h[,t]!=h[,t+1] & h[,t+2]==h[,t],
      true_persistent=h[,t]!=h[,t+1] & h[,t+2]==h[,t+1])
    stopifnot(all(rowSums(masks)==1L))
    probability <- gamma %*% masks
    y0 <- df[[paste0("y",t)]]; y1 <- df[[paste0("y",t+1)]]
    selected <- list(y0!=y1, y0==0 & y1==1, y0==1 & y1==0)
    for (d in 1:3) {
      w <- df$weight * selected[[d]]
      numerator[d,] <- numerator[d,] + colSums(probability*w)
      denominator[d] <- denominator[d] + sum(w)
    }
  }
  stopifnot(all(denominator>0),
    max(abs(numerator[1,]-colSums(numerator[2:3,,drop=FALSE]))) < 1e-6,
    abs(denominator[1]-sum(denominator[2:3])) < 1e-6)
  shares <- numerator/denominator
  stopifnot(max(abs(rowSums(shares)-1)) < 1e-8)
  result <- data.frame(direction=directions,model=model,
    conditioning="full_duration_posterior",shares,row.names=NULL)
  attr(result,"apparent_weight") <- denominator
  result
}

build_duration_implications_4w <- function() {
  paths <- c(
    ar1="EM-tenure/output/results/job_change_monthly/four_wave_ar1/continuous_clock_empirical/comparison_latest.rds",
    ar2="EM-tenure/output/results/job_change_monthly/four_wave_ar2/empirical/comparison_latest.rds",
    data="EM-tenure/output/results/job_change_monthly/four_wave_ar1/prepared_cells_latest.rds")
  ar1 <- readRDS(paths["ar1"])$best
  saved <- readRDS(paths["ar2"])
  stopifnot(identical(saved$source_md5,four_wave_ar2_source_fingerprint()),
    isTRUE(ar1$converged),isTRUE(saved$best$converged))
  df <- readRDS(paths["data"])$df4
  rows <- rbind(
    posterior_duration_implications_4w(df,ar1,"Table 9 col 1: duration AR(1)"),
    posterior_duration_implications_4w(df,saved$best,"Table 9 col 2: duration AR(2)"))
  attr(rows,"input_md5") <- tools::md5sum(paths)
  rows
}
