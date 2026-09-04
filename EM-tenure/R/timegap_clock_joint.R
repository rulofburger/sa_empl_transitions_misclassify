# One continuous nonemployment clock D0 is drawn at the latent spell start.
# A clean category at offset t restricts D0 to its interval shifted by -t/4.
# Missing reports impose no restriction; gross reports are independent draws
# from the unchanged marginal category distribution. Sum over report types.
log_emission_timegap_clock_joint <- function(categories, lambda_d,
    beta_d=0, eps_d, offsets=NULL, observed=NULL) {
  categories <- as.matrix(categories)
  n <- nrow(categories); k <- ncol(categories)
  if (is.null(offsets)) offsets <- seq_len(k)-1L
  if (length(offsets)!=k || any(!is.finite(offsets) | offsets<0) ||
      any(diff(offsets)<=0)) stop("Invalid timegap clock offsets")
  if (is.null(observed)) observed <- !is.na(categories)
  observed <- as.matrix(observed)
  if (!identical(dim(observed),dim(categories)) || anyNA(observed) ||
      any(observed & (is.na(categories) | !categories %in% 1:7)))
    stop("Invalid observed timegap categories")
  if (length(eps_d)!=1L || !is.finite(eps_d) || eps_d<0 || eps_d>1)
    stop("eps_d must be in [0,1]")
  if (!k) return(rep(0,n))
  intervals <- vapply(1:7,.timegap_interval,numeric(2))
  marginal <- log_emission_interval_d(1:7,lambda_d,beta_d)
  group <- as.vector(observed %*% 2^(seq_len(k)-1L))
  result <- numeric(n)
  for (code in unique(group)) {
    rows <- which(group==code)
    use <- which(observed[rows[1],])
    if (!length(use)) next
    pattern <- as.matrix(expand.grid(rep(list(0:1),length(use))))
    terms <- matrix(-Inf,length(rows),nrow(pattern))
    for (s in seq_len(nrow(pattern))) {
      dirty <- use[pattern[s,]==1L]; clean <- use[pattern[s,]==0L]
      if ((length(dirty)>0 && eps_d==0) ||
          (length(clean)>0 && eps_d==1)) next
      value <- rep(0,length(rows))
      if (length(dirty)) {
        value <- value+length(dirty)*log(eps_d)
        for (t in dirty) value <- value+marginal[categories[rows,t]]
      }
      if (length(clean)) {
        value <- value+length(clean)*log1p(-eps_d)
        lower <- rep(0,length(rows)); upper <- rep(Inf,length(rows))
        for (t in clean) {
          lower <- pmax(lower,intervals[1,categories[rows,t]]-
            offsets[t]*.QUARTER_YEARS)
          upper <- pmin(upper,intervals[2,categories[rows,t]]-
            offsets[t]*.QUARTER_YEARS)
        }
        valid <- lower<upper
        mass <- rep(-Inf,length(rows))
        if (any(valid)) {
          hl <- .duration_cumhaz(lower[valid],lambda_d,beta_d)
          hu <- .duration_cumhaz(upper[valid],lambda_d,beta_d)
          mass[valid] <- -hl+log(-expm1(-(hu-hl)))
        }
        value <- value+mass
      }
      terms[,s] <- value
    }
    result[rows] <- .row_logsumexp(terms)
  }
  result
}

.four_wave_timegap_clock_tables <- function(lambda_d,components) {
  lapply(1:4,function(k) {
    categories <- as.matrix(expand.grid(rep(list(0:7),k)))
    categories[categories==0L] <- NA_integer_
    vapply(components,function(p) log_emission_timegap_clock_joint(
      categories,lambda_d,eps_d=p$eps_d),numeric(8^k))
  })
}
