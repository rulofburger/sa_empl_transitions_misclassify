# ==============================================================================
# Four-wave AR(1) employment model with misclassification
# ==============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

latent_histories_ar1_4w <- function() {
  as.matrix(expand.grid(h1 = 0:1, h2 = 0:1, h3 = 0:1, h4 = 0:1))
}

stationary_alpha_ar1 <- function(theta0, theta1) {
  theta0 / (theta0 + 1 - theta1)
}

collapse_ar1_4w_cells <- function(df) {
  required <- c("y1", "y2", "y3", "y4", "weight")
  missing <- setdiff(required, names(df))
  if (length(missing)) stop("collapse_ar1_4w_cells: missing ", paste(missing, collapse = ", "))
  y <- as.matrix(df[, paste0("y", 1:4), drop = FALSE])
  if (anyNA(y) || !all(y %in% 0:1)) stop("y1-y4 must be non-missing binary values")
  if (anyNA(df$weight) || any(!is.finite(df$weight)) || any(df$weight <= 0))
    stop("weights must be finite and strictly positive")
  index <- 1L + y[, 1] + 2L*y[, 2] + 4L*y[, 3] + 8L*y[, 4]
  histories <- as.data.frame(expand.grid(y1=0:1, y2=0:1, y3=0:1, y4=0:1))
  weight <- weight_sq <- numeric(16L)
  sw <- rowsum(df$weight, index, reorder = FALSE)
  sw2 <- rowsum(df$weight^2, index, reorder = FALSE)
  weight[as.integer(rownames(sw))] <- sw[, 1]
  weight_sq[as.integer(rownames(sw2))] <- sw2[, 1]
  out <- cbind(histories, weight = weight, weight_sq = weight_sq,
               count = tabulate(index, nbins = 16L))
  attr(out, "n_obs") <- nrow(df)
  attr(out, "weight_sum") <- sum(df$weight)
  out
}

.ar1_4w_cache <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      observed <- as.matrix(expand.grid(y1=0:1, y2=0:1, y3=0:1, y4=0:1))
      latent <- latent_histories_ar1_4w()
      mismatch <- matrix(0L, 16L, 16L)
      for (tt in 1:4) mismatch <- mismatch + outer(observed[, tt], latent[, tt], "!=")
      cache <<- list(observed = observed, latent = latent, mismatch = mismatch)
    }
    cache
  }
})

ar1_4w_cell_probabilities <- function(params, model_type = "symmetric") {
  if (!model_type %in% c("none", "symmetric"))
    stop("model_type must be 'none' or 'symmetric'")
  cache <- .ar1_4w_cache()
  h <- cache$latent
  alpha <- params$alpha
  prior <- ifelse(h[, 1] == 1, alpha, 1 - alpha)
  for (tt in 2:4) {
    p1 <- ifelse(h[, tt-1] == 1, params$theta1, params$theta0)
    prior <- prior * ifelse(h[, tt] == 1, p1, 1 - p1)
  }
  pi <- if (model_type == "symmetric") params$pi else 0
  emission <- (1 - pi)^(4L - cache$mismatch) * pi^cache$mismatch
  probs <- as.vector(emission %*% prior)
  probs / sum(probs)
}

.ar1_4w_eta_names <- function(model_type, stationary) {
  out <- c("theta0", "theta1")
  if (!stationary) out <- c(out, "alpha")
  if (model_type == "symmetric") out <- c(out, "pi")
  out
}

pack_ar1_4w_params <- function(params, model_type, stationary) {
  clamp <- function(x, eps=1e-8) pmin(pmax(x, eps), 1-eps)
  nms <- .ar1_4w_eta_names(model_type, stationary)
  out <- setNames(numeric(length(nms)), nms)
  out["theta0"] <- qlogis(clamp(params$theta0))
  out["theta1"] <- qlogis(clamp(params$theta1))
  if (!stationary) out["alpha"] <- qlogis(clamp(params$alpha))
  if (model_type == "symmetric") out["pi"] <- qlogis(clamp(params$pi / .5))
  out
}

unpack_ar1_4w_eta <- function(eta, model_type, stationary) {
  eta <- setNames(as.numeric(eta), .ar1_4w_eta_names(model_type, stationary))
  theta0 <- plogis(eta["theta0"])
  theta1 <- plogis(eta["theta1"])
  alpha <- if (stationary) stationary_alpha_ar1(theta0, theta1) else plogis(eta["alpha"])
  out <- list(theta0=unname(theta0), theta1=unname(theta1), alpha=unname(alpha))
  if (model_type == "symmetric") out$pi <- unname(.5 * plogis(eta["pi"]))
  out
}

.ar1_4w_jacobian <- function(fn, x, rel_step=1e-5) {
  f0 <- fn(x)
  J <- matrix(NA_real_, length(f0), length(x), dimnames=list(names(f0), names(x)))
  for (j in seq_along(x)) {
    h <- rel_step * max(1, abs(x[j])); xp <- xm <- x
    xp[j] <- xp[j]+h; xm[j] <- xm[j]-h
    J[, j] <- (fn(xp)-fn(xm))/(2*h)
  }
  J
}

.ar1_4w_gradient <- function(fn, x, rel_step=1e-5) {
  as.numeric(.ar1_4w_jacobian(function(z) setNames(fn(z), "value"), x, rel_step))
}

initial_ar1_4w_params <- function(df, model_type="symmetric", stationary=FALSE) {
  y <- as.matrix(df[, paste0("y",1:4)])
  w <- rep(df$weight, 3L)
  from <- as.vector(y[,1:3]); to <- as.vector(y[,2:4])
  theta0 <- weighted.mean(to[from==0], w[from==0])
  theta1 <- weighted.mean(to[from==1], w[from==1])
  alpha <- weighted.mean(y[,1], df$weight)
  if (stationary) alpha <- stationary_alpha_ar1(theta0, theta1)
  out <- list(theta0=theta0, theta1=theta1, alpha=alpha)
  if (model_type == "symmetric") out$pi <- .02
  out
}

fit_ar1_4w_mle <- function(df, model_type="symmetric", stationary=FALSE,
                            starts=NULL, random_starts=12L, seed=20260818L,
                            maxit=4000L, reltol=1e-12, verbose=1L) {
  if (!model_type %in% c("none","symmetric")) stop("invalid model_type")
  cells <- collapse_ar1_4w_cells(df)
  if (is.null(starts)) starts <- list(initial_ar1_4w_params(df, model_type, stationary))
  base_eta <- pack_ar1_4w_params(starts[[1]], model_type, stationary)
  set.seed(seed)
  if (random_starts > 0) for (i in seq_len(random_starts))
    starts[[length(starts)+1L]] <- unpack_ar1_4w_eta(base_eta+rnorm(length(base_eta),0,.8), model_type, stationary)
  fn <- function(z) {
    p <- unpack_ar1_4w_eta(z, model_type, stationary)
    probs <- pmax(ar1_4w_cell_probabilities(p, model_type), 1e-300)
    -sum(cells$weight*log(probs))/sum(cells$weight)
  }
  candidates <- lapply(starts, function(start) tryCatch(
    optim(pack_ar1_4w_params(start,model_type,stationary), fn, method="BFGS",
          control=list(maxit=maxit,reltol=reltol)), error=function(e) NULL))
  candidates <- Filter(function(x) !is.null(x) && is.finite(x$value), candidates)
  if (!length(candidates)) stop("all optimization starts failed")
  opt <- candidates[[which.min(vapply(candidates, `[[`, numeric(1), "value"))]]
  eta <- setNames(opt$par, .ar1_4w_eta_names(model_type,stationary))
  params <- unpack_ar1_4w_eta(eta,model_type,stationary)
  probs <- ar1_4w_cell_probabilities(params,model_type)
  loglik <- sum(cells$weight*log(pmax(probs,1e-300)))
  score <- .ar1_4w_gradient(fn,eta)
  H <- optimHess(eta,fn); H <- (H+t(H))/2
  eig <- eigen(H,symmetric=TRUE,only.values=TRUE)$values
  prob_jac <- .ar1_4w_jacobian(function(z) setNames(ar1_4w_cell_probabilities(
    unpack_ar1_4w_eta(z,model_type,stationary),model_type)[1:15],paste0("cell",1:15)),eta)
  rank <- qr(prob_jac,tol=1e-8)$rank
  converged <- opt$convergence==0L && max(abs(score))<1e-6 && min(eig)>1e-9
  if (verbose) cat(sprintf("AR1-4W [%s, stationary=%s]: ll=%.3f score=%.2e rank=%d/%d minEig=%.2e\n",
    model_type,stationary,loglik,max(abs(score)),rank,length(eta),min(eig)))
  list(params=params,eta=eta,loglik=loglik,converged=converged,
       identified=rank==length(eta),model_type=model_type,stationary=stationary,
       cell_probabilities=probs,diagnostics=list(optimizer_code=opt$convergence,
       max_abs_score=max(abs(score)),information_min_eigenvalue=min(eig),
       information_condition=max(eig)/min(eig),probability_jacobian_rank=rank,
       parameter_count=length(eta),start_count=length(candidates)),
       sample=list(n=attr(cells,"n_obs"),weight_sum=attr(cells,"weight_sum")),
       estimator="direct_sixteen_cell_mle")
}

.ar1_4w_quantities <- function(eta,model_type,stationary) {
  p <- unpack_ar1_4w_eta(eta,model_type,stationary)
  out <- c(initial_employment=p$alpha,entry_rate=p$theta0,
           exit_rate=1-p$theta1,
           steady_state_employment=stationary_alpha_ar1(p$theta0,p$theta1),
           initial_minus_steady_state=p$alpha-stationary_alpha_ar1(p$theta0,p$theta1))
  if (model_type=="symmetric") out <- c(out,pi=p$pi)
  out
}

analytical_se_ar1_4w <- function(df,fit,finite_sample=TRUE) {
  cells <- collapse_ar1_4w_cells(df); eta <- fit$eta
  type <- fit$model_type; stationary <- fit$stationary
  logcells <- function(z) setNames(log(pmax(ar1_4w_cell_probabilities(
    unpack_ar1_4w_eta(z,type,stationary),type),1e-300)),paste0("cell",1:16))
  scores <- .ar1_4w_jacobian(logcells,eta)
  avg_nll <- function(z) -sum(cells$weight*logcells(z))/sum(cells$weight)
  bread <- optimHess(eta,avg_nll); bread <- (bread+t(bread))/2
  meat <- matrix(0,length(eta),length(eta),dimnames=list(names(eta),names(eta)))
  for (j in 1:16) meat <- meat + cells$weight_sq[j]/sum(cells$weight)^2*tcrossprod(scores[j,])
  inv <- solve(bread); vcov_eta <- inv%*%meat%*%inv
  n <- attr(cells,"n_obs"); K <- length(eta)
  if (finite_sample && n>K) vcov_eta <- vcov_eta*n/(n-K)
  qfun <- function(z) .ar1_4w_quantities(z,type,stationary)
  estimates <- qfun(eta); delta <- .ar1_4w_jacobian(qfun,eta)
  vcov_q <- delta%*%vcov_eta%*%t(delta)
  list(summary=data.frame(quantity=names(estimates),estimate=unname(estimates),
       se=sqrt(pmax(diag(vcov_q),0))),vcov_eta=vcov_eta,vcov_quantities=vcov_q,
       bread=bread,meat=meat,delta_jacobian=delta,
       method="individual_survey_weighted_sandwich_delta")
}
