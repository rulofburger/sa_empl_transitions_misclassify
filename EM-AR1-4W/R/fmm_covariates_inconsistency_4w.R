# Two-type four-wave AR(1) FMM with common demographic, observed-duration,
# and job-characteristic slopes plus inconsistency-dependent misclassification.

.clamp_probability <- function(x, lo = 1e-9, hi = 1 - 1e-9)
  pmin(pmax(x, lo), hi)

compute_inconsistencies_4w <- function(df) {
  needed <- c(paste0("age", 1:4), paste0("educ", 1:4),
              paste0("race", 1:4), paste0("female", 1:4))
  missing <- setdiff(needed, names(df))
  if (length(missing)) stop("Missing inconsistency variables: ", paste(missing, collapse = ", "))
  edge_flags <- function(prefix, admissible_change = TRUE) {
    x <- as.matrix(df[paste0(prefix, 1:4)])
    if (admissible_change) {
      d <- x[, 2:4, drop = FALSE] - x[, 1:3, drop = FALSE]
      !is.na(d) & !(abs(d) < .01 | abs(d - 1) < .01)
    } else {
      !is.na(x[, 2:4, drop = FALSE]) & !is.na(x[, 1:3, drop = FALSE]) &
        x[, 2:4, drop = FALSE] != x[, 1:3, drop = FALSE]
    }
  }
  attribute_flags <- function(v) cbind(
    v[, 1] & !v[, 2], v[, 1] & v[, 2],
    (!v[, 1] & v[, 2]) | (v[, 2] & v[, 3]), !v[, 2] & v[, 3]) * 1L
  age <- attribute_flags(edge_flags("age"))
  educ <- attribute_flags(edge_flags("educ"))
  race <- attribute_flags(edge_flags("race", FALSE))
  gender <- attribute_flags(edge_flags("female", FALSE))
  colnames(age) <- paste0("Y_age_", 1:4)
  colnames(educ) <- paste0("Y_edu_", 1:4)
  colnames(race) <- paste0("Y_race_", 1:4)
  colnames(gender) <- paste0("Y_gender_", 1:4)
  components <- cbind(age, educ, race, gender)
  counts <- do.call(cbind, lapply(1:4, function(tt) {
    count <- rowSums(components[, paste0("Y_", c("age", "edu", "race", "gender"),
                                                "_", tt), drop = FALSE])
    out <- cbind(as.integer(count == 2L), as.integer(count == 3L),
                 as.integer(count == 4L))
    colnames(out) <- paste0("Y_exactly_", 2:4, "_", tt)
    out
  }))
  cbind(components, counts)
}

prepare_fmm_covariates_inconsistency_4w <- function(
    panel_path = "data/raw/df_qlfs_A.rds",
    sector_path = "data/raw/QLFSmerged_mapped.rds", collapse = TRUE,
    error_design = c("age_education", "intercept_only", "table3_column1")) {
  error_design <- match.arg(error_design)
  panel <- readRDS(panel_path)
  required <- c("hhnr", "pnr", paste0("period", 1:4), paste0("employed", 1:4),
    paste0("age", 1:4), paste0("educ", 1:4), paste0("race", 1:4),
    paste0("female", 1:4),
    paste0("contracttype", 1:3), paste0("weight", 1:4),
    paste0("tenure", 1:3), paste0("timegap", 1:3), paste0("neverworked", 1:3))
  missing <- setdiff(required, names(panel))
  if (length(missing)) stop("Four-wave panel is missing: ", paste(missing, collapse = ", "))
  keep <- panel$age1 > 17 & panel$age1 < 56 &
    complete.cases(panel[c(paste0("employed", 1:4), paste0("age", 1:4),
      paste0("educ", 1:4), "race1", "female1", paste0("weight", 1:4))])
  panel <- panel[keep, required, drop = FALSE]
  names(panel)[match(paste0("employed", 1:4), names(panel))] <- paste0("y", 1:4)
  for (nm in c(paste0("y", 1:4), paste0("age", 1:4), paste0("educ", 1:4),
               paste0("race", 1:4), paste0("female", 1:4),
               paste0("contracttype", 1:3),
               paste0("tenure", 1:3), paste0("timegap", 1:3),
               paste0("neverworked", 1:3)))
    panel[[nm]] <- as.numeric(unclass(panel[[nm]]))
  panel$weight <- with(panel, (weight1 * weight2 * weight3 * weight4) ^ .25)
  if (any(!is.finite(panel$weight) | panel$weight <= 0)) stop("Invalid survey weights")
  panel$weight <- nrow(panel) * panel$weight / sum(panel$weight)

  # Match the three-wave observed-duration specification. Values are retained
  # only in the reported origin state; missingness indicators keep the full
  # four-wave sample without interpreting missing durations as genuine zeros.
  timegap_months <- c(`1`=1.5, `2`=4.5, `3`=7.5, `4`=10.5,
                      `5`=24, `6`=48, `7`=90)
  for (tt in 1:3) {
    y <- panel[[paste0("y", tt)]]
    tenure <- panel[[paste0("tenure", tt)]]
    tenure[tenure < 0] <- NA_real_
    tenure[y == 0] <- 0
    gap_code <- panel[[paste0("timegap", tt)]]
    gap <- unname(timegap_months[as.character(gap_code)])
    gap[gap_code == 0] <- 0
    gap[y == 1] <- 0
    never <- panel[[paste0("neverworked", tt)]]
    gap[y == 0 & never == 1] <- pmax(
      panel[[paste0("age", tt)]][y == 0 & never == 1] -
        panel[[paste0("educ", tt)]][y == 0 & never == 1] - 6, 0)
    panel[[paste0("tenure_missing_cov", tt)]] <- as.integer(y == 1 & is.na(tenure))
    panel[[paste0("timegap_missing_cov", tt)]] <- as.integer(y == 0 & is.na(gap))
    panel[[paste0("tenure_cov", tt)]] <- ifelse(y == 1, ifelse(is.na(tenure), 0, tenure / 12), 0)
    panel[[paste0("timegap_cov", tt)]] <- ifelse(y == 0, ifelse(is.na(gap), 0, gap / 12), 0)
    panel[[paste0("neverworked_cov", tt)]] <- ifelse(
      y == 0, ifelse(is.na(never), 0, never), 0)
  }

  # Keep base-data-frame subsetting semantics even when the saved source was
  # produced as a data.table (data.table treats a variable in `j` differently).
  sector <- as.data.frame(readRDS(sector_path))
  sector_required <- c("hhnr", "pnr_methodA", "wave", "sector2")
  missing <- setdiff(sector_required, names(sector))
  if (length(missing)) stop("Sector source is missing: ", paste(missing, collapse = ", "))
  relevant <- sector$wave %in% unique(unlist(panel[paste0("period", 1:3)])) &
    !is.na(sector$pnr_methodA)
  sector <- sector[relevant, sector_required, drop = FALSE]
  skey <- paste(sector$hhnr, sector$pnr_methodA, sector$wave, sep = "|")
  if (anyDuplicated(skey)) stop("Sector merge key is not unique")
  sector2 <- as.numeric(unclass(sector$sector2))
  for (tt in 1:3) {
    pkey <- paste(panel$hhnr, panel$pnr, panel[[paste0("period", tt)]], sep = "|")
    idx <- match(pkey, skey)
    if (anyNA(idx)) stop("Sector merge failed at wave ", tt)
    employed <- panel[[paste0("y", tt)]] == 1
    if (anyNA(sector2[idx][employed])) stop("Sector missing for employed respondents at wave ", tt)
    panel[[paste0("informal", tt)]] <- as.integer(employed & sector2[idx] == 2)
    contract <- panel[[paste0("contracttype", tt)]]
    panel[[paste0("permanent", tt)]] <- as.integer(employed & !is.na(contract) & contract == 1)
  }
  rm(sector)
  std <- function(x) {
    s <- sd(x); if (!is.finite(s) || s < 1e-10) stop("Constant covariate")
    list(value = (x - mean(x)) / s, center = mean(x), scale = s)
  }
  age <- std(panel$age1); age_sq <- std(panel$age1^2); educ <- std(panel$educ1)
  races <- sort(unique(panel$race1))
  race_dummies <- sapply(races[-1], function(z) as.integer(panel$race1 == z))
  if (is.null(dim(race_dummies))) race_dummies <- matrix(race_dummies, ncol = 1)
  colnames(race_dummies) <- paste0("race_", races[-1])
  common <- cbind(age = age$value, age_sq = age_sq$value, educ = educ$value,
                  race_dummies, female = panel$female1)
  duration_std <- function(prefix, transform = identity) {
    values <- transform(unlist(panel[paste0(prefix, 1:3)], use.names = FALSE))
    z <- std(values)
    list(values = split(z$value, rep(1:3, each = nrow(panel))),
         center = z$center, scale = z$scale)
  }
  log_tenure <- duration_std("tenure_cov", log1p)
  log_timegap <- duration_std("timegap_cov", log1p)
  X <- lapply(1:3, function(tt) cbind(common,
    log_tenure = log_tenure$values[[tt]],
    log_time_since_work = log_timegap$values[[tt]],
    never_worked = panel[[paste0("neverworked_cov", tt)]],
    tenure_missing = panel[[paste0("tenure_missing_cov", tt)]],
    timegap_missing = panel[[paste0("timegap_missing_cov", tt)]],
    permanent_contract = panel[[paste0("permanent", tt)]],
    informal_sector = panel[[paste0("informal", tt)]]))
  # Drop a missingness indicator if the final sample has no such reports.
  varying <- vapply(colnames(X[[1]]), function(nm)
    any(vapply(X, function(x) any(x[, nm] != 0), logical(1))), logical(1))
  optional_missing <- colnames(X[[1]]) %in% c("tenure_missing", "timegap_missing")
  keep_x <- varying | !optional_missing
  X <- lapply(X, function(x) x[, keep_x, drop = FALSE])
  entry_active <- !colnames(X[[1]]) %in%
    c("log_tenure", "tenure_missing", "permanent_contract", "informal_sector")
  persistence_active <- !colnames(X[[1]]) %in%
    c("log_time_since_work", "never_worked", "timegap_missing")
  inc <- compute_inconsistencies_4w(panel)
  Z <- lapply(1:4, function(tt) {
    if (error_design == "intercept_only")
      return(cbind(intercept = rep(1, nrow(inc))))
    out <- cbind(intercept = rep(1, nrow(inc)),
      age_inconsistent = inc[, paste0("Y_age_", tt)],
      education_inconsistent = inc[, paste0("Y_edu_", tt)])
    if (error_design == "table3_column1") out <- cbind(out,
      race_inconsistent = inc[, paste0("Y_race_", tt)],
      gender_inconsistent = inc[, paste0("Y_gender_", tt)],
      two_inconsistencies = inc[, paste0("Y_exactly_2_", tt)],
      three_inconsistencies = inc[, paste0("Y_exactly_3_", tt)],
      four_inconsistencies = inc[, paste0("Y_exactly_4_", tt)])
    storage.mode(out) <- "double"
    out
  })
  out <- list(y = as.matrix(panel[paste0("y", 1:4)]), weight = panel$weight,
    weight_sq = panel$weight^2, X = X, Z = Z, entry_active = entry_active,
    persistence_active = persistence_active,
    covariate_names = colnames(X[[1]]), error_design = error_design,
    scaling = list(age = age[c("center", "scale")],
      age_sq = age_sq[c("center", "scale")], educ = educ[c("center", "scale")],
      log_tenure = log_tenure[c("center", "scale")],
      log_time_since_work = log_timegap[c("center", "scale")]),
    n_original = nrow(panel))
  if (!collapse) return(out)
  z_blocks <- lapply(Z, function(z) z[, -1, drop = FALSE])
  z_no_intercept <- if (sum(vapply(z_blocks, ncol, integer(1))) > 0L)
    do.call(cbind, z_blocks) else matrix(numeric(0), nrow(out$y), 0L)
  key_parts <- cbind(out$y, do.call(cbind, X), z_no_intercept)
  key <- do.call(paste, c(as.data.frame(key_parts), sep = "\r"))
  first <- !duplicated(key); group <- match(key, key[first])
  collapse_sum <- function(x) as.vector(rowsum(x, group, reorder = FALSE))
  out$y <- out$y[first, , drop = FALSE]
  out$X <- lapply(out$X, function(x) x[first, , drop = FALSE])
  out$Z <- lapply(out$Z, function(x) x[first, , drop = FALSE])
  out$weight <- collapse_sum(out$weight); out$weight_sq <- collapse_sum(out$weight_sq)
  out$n_cells <- length(out$weight)
  out
}

.initial_fmm_covinc_4w <- function(data, base_fit = NULL) {
  base <- if (is.null(base_fit)) list(theta0_A=.025,theta1_A=.975,alpha_A=.48,
    theta0_B=.42,theta1_B=.53,alpha_B=.48,phi=.92,pi=.014) else base_fit$params
  nms <- data$covariate_names
  b0 <- setNames(rep(0, 2 + sum(data$entry_active)),
    c("intercept_A", "intercept_B", nms[data$entry_active]))
  active1 <- data$persistence_active %||% rep(TRUE, length(nms))
  b1 <- setNames(rep(0, 2 + sum(active1)),
    c("intercept_A", "intercept_B", nms[active1]))
  b0[1:2] <- qnorm(c(base$theta0_A, base$theta0_B))
  b1[1:2] <- qnorm(c(base$theta1_A, base$theta1_B))
  delta <- setNames(c(qlogis(2 * base$pi), rep(0, ncol(data$Z[[1]]) - 1L)),
                    colnames(data$Z[[1]]))
  list(alpha=c(A=base$alpha_A,B=base$alpha_B),phi=base$phi,beta0=b0,beta1=b1,
    delta=delta)
}

.expand_fmm_covinc_start_4w <- function(data, old_params) {
  out <- .initial_fmm_covinc_4w(data)
  out$alpha <- old_params$alpha
  out$phi <- old_params$phi
  for (block in c("beta0", "beta1", "delta")) {
    common <- intersect(names(out[[block]]), names(old_params[[block]]))
    out[[block]][common] <- old_params[[block]][common]
  }
  out
}

.fmm_covinc_type_specific_delta <- function(delta)
  all(c("intercept_A", "intercept_B") %in% names(delta))

.fmm_covinc_pi <- function(Z, delta, type=1L) {
  if (.fmm_covinc_type_specific_delta(delta)) {
    intercept <- delta[c("intercept_A", "intercept_B")[type]]
    slope_names <- colnames(Z)[-1L]
    if (!all(slope_names %in% names(delta)))
      stop("Type-specific error vector does not match the error design")
    eta <- intercept
    if (length(slope_names))
      eta <- eta + as.vector(Z[, slope_names, drop=FALSE] %*% delta[slope_names])
  } else {
    if (!all(colnames(Z) %in% names(delta)))
      stop("Common error vector does not match the error design")
    eta <- as.vector(Z %*% delta[colnames(Z)])
  }
  .clamp_probability(.5*plogis(eta),1e-8,.49)
}

.fmm_covinc_add_delta_score <- function(gradient, Z, score, delta, type) {
  if (.fmm_covinc_type_specific_delta(delta)) {
    intercept_name <- c("intercept_A", "intercept_B")[type]
    gradient[intercept_name] <- gradient[intercept_name] + sum(score)
    slope_names <- colnames(Z)[-1L]
    if (length(slope_names)) gradient[slope_names] <- gradient[slope_names] +
      as.vector(crossprod(Z[, slope_names, drop=FALSE], score))
  } else {
    gradient[colnames(Z)] <- gradient[colnames(Z)] +
      as.vector(crossprod(Z, score))
  }
  gradient
}

.fmm_covinc_transition <- function(X, beta, type, active=rep(TRUE,ncol(X))) {
  intercept <- beta[if(type==1L) 1L else 2L]
  .clamp_probability(pnorm(intercept + as.vector(X[,active,drop=FALSE] %*% beta[-(1:2)])))
}

e_step_fmm_covinc_4w <- function(data, params, retain=TRUE) {
  n <- nrow(data$y)
  pi <- lapply(1:2,function(k)lapply(1:4,function(tt)
    .fmm_covinc_pi(data$Z[[tt]],params$delta,k)))
  emission <- lapply(1:2,function(k)lapply(1:4,function(tt)cbind(
    ifelse(data$y[,tt]==0,1-pi[[k]][[tt]],pi[[k]][[tt]]),
    ifelse(data$y[,tt]==1,1-pi[[k]][[tt]],pi[[k]][[tt]]))))
  prior_type <- c(params$phi,1-params$phi)
  forward <- backward <- vector("list",2); likelihood_type <- matrix(NA_real_,n,2)
  q0 <- q1 <- vector("list",2)
  for(k in 1:2) {
    q0[[k]] <- lapply(1:3,function(tt).fmm_covinc_transition(
      data$X[[tt]],params$beta0,k,data$entry_active))
    active1 <- data$persistence_active %||% rep(TRUE, length(data$entry_active))
    q1[[k]] <- lapply(1:3,function(tt).fmm_covinc_transition(
      data$X[[tt]],params$beta1,k,active1))
    forward[[k]] <- vector("list",4)
    initial <- matrix(c(1-params$alpha[k],params$alpha[k]),n,ncol=2,byrow=TRUE)
    forward[[k]][[1]] <- initial*emission[[k]][[1]]
    for(tt in 2:4) {
      old <- forward[[k]][[tt-1]]
      pred0 <- old[,1]*(1-q0[[k]][[tt-1]])+old[,2]*(1-q1[[k]][[tt-1]])
      pred1 <- old[,1]*q0[[k]][[tt-1]]+old[,2]*q1[[k]][[tt-1]]
      forward[[k]][[tt]] <- cbind(pred0,pred1)*emission[[k]][[tt]]
    }
    likelihood_type[,k] <- rowSums(forward[[k]][[4]])
    backward[[k]] <- vector("list",4); backward[[k]][[4]] <- matrix(1,n,2)
    for(tt in 3:1) {
      nxt <- emission[[k]][[tt+1]]*backward[[k]][[tt+1]]
      backward[[k]][[tt]] <- cbind(
        (1-q0[[k]][[tt]])*nxt[,1]+q0[[k]][[tt]]*nxt[,2],
        (1-q1[[k]][[tt]])*nxt[,1]+q1[[k]][[tt]]*nxt[,2])
    }
  }
  likelihood <- pmax(as.vector(likelihood_type%*%prior_type),1e-300)
  loglik <- sum(data$weight*log(likelihood))
  if(!retain)return(loglik)
  state <- lapply(1:2,function(k)lapply(1:4,function(tt)
    prior_type[k]*forward[[k]][[tt]]*backward[[k]][[tt]]/likelihood))
  xi <- lapply(1:2,function(k)lapply(1:3,function(tt){
    old <- forward[[k]][[tt]]; nxt <- emission[[k]][[tt+1]]*backward[[k]][[tt+1]]
    cbind(prior_type[k]*old[,1]*(1-q0[[k]][[tt]])*nxt[,1]/likelihood,
      prior_type[k]*old[,1]*q0[[k]][[tt]]*nxt[,2]/likelihood,
      prior_type[k]*old[,2]*(1-q1[[k]][[tt]])*nxt[,1]/likelihood,
      prior_type[k]*old[,2]*q1[[k]][[tt]]*nxt[,2]/likelihood)}))
  list(loglik=loglik,state=state,xi=xi,q0=q0,q1=q1,
    posterior_type=cbind(rowSums(state[[1]][[1]]),rowSums(state[[2]][[1]])))
}

.fit_fractional_probit <- function(data,suff,old_beta,entry) {
  active <- if(entry)data$entry_active else
    (data$persistence_active %||% rep(TRUE,length(data$entry_active)))
  objective <- function(beta) {
    q <- 0
    for(k in 1:2)for(tt in 1:3) {
      x <- data$X[[tt]][,active,drop=FALSE]
      eta <- beta[k]+as.vector(x%*%beta[-(1:2)]); p <- .clamp_probability(pnorm(eta))
      w <- data$weight*suff[[k]][[tt]]$risk
      s <- data$weight*suff[[k]][[tt]]$success
      q <- q+sum(s*log(p)+(w-s)*log1p(-p))
    }
    q/sum(data$weight)
  }
  # Fisher scoring for grouped fractional probit observations. Building the
  # small cross-products directly is much faster than repeatedly evaluating a
  # generic optimizer over the six type-transition blocks.
  beta <- old_beta; pcoef <- length(beta); common <- seq.int(3L,pcoef)
  for(iter in 1:6) {
    A <- matrix(0,pcoef,pcoef); b <- numeric(pcoef)
    for(k in 1:2)for(tt in 1:3) {
      x <- data$X[[tt]][,active,drop=FALSE]
      eta <- beta[k]+as.vector(x%*%beta[common]); mu <- .clamp_probability(pnorm(eta))
      risk <- data$weight*suff[[k]][[tt]]$risk
      success <- data$weight*suff[[k]][[tt]]$success
      ybar <- ifelse(risk>1e-14,success/risk,mu)
      deriv <- pmax(dnorm(eta),1e-10)
      ww <- risk*deriv^2/(mu*(1-mu)); wz <- ww*(eta+(ybar-mu)/deriv)
      sx <- as.vector(crossprod(x,ww))
      A[k,k] <- A[k,k]+sum(ww); A[k,common] <- A[k,common]+sx
      A[common,k] <- A[common,k]+sx
      A[common,common] <- A[common,common]+crossprod(x,x*ww)
      b[k] <- b[k]+sum(wz); b[common] <- b[common]+as.vector(crossprod(x,wz))
    }
    proposal <- tryCatch(solve(A,b),error=function(e)beta)
    if(max(abs(proposal-beta))<1e-7){beta<-proposal;break}
    beta <- proposal
  }
  q_old <- objective(old_beta); q_new <- objective(beta); fraction <- 1
  while((!is.finite(q_new)||q_new<q_old-1e-12)&&fraction>2^-15) {
    fraction <- fraction/2; beta <- old_beta+fraction*(beta-old_beta)
    q_new <- objective(beta)
  }
  if(!is.finite(q_new)||q_new<q_old-1e-12)old_beta else beta
}

m_step_fmm_covinc_4w <- function(data,e,params) {
  out <- params
  type_mass <- sapply(1:2,function(k)sum(data$weight*e$posterior_type[,k]))
  out$phi <- .clamp_probability(type_mass[1]/sum(type_mass))
  out$alpha <- sapply(1:2,function(k).clamp_probability(
    sum(data$weight*e$state[[k]][[1]][,2])/type_mass[k]))
  names(out$alpha) <- c("A","B")
  entry <- persist <- lapply(1:2,function(k)lapply(1:3,function(tt)list()))
  for(k in 1:2)for(tt in 1:3) {
    z <- e$xi[[k]][[tt]]
    entry[[k]][[tt]] <- list(risk=z[,1]+z[,2],success=z[,2])
    persist[[k]][[tt]] <- list(risk=z[,3]+z[,4],success=z[,4])
  }
  out$beta0 <- .fit_fractional_probit(data,entry,params$beta0,TRUE)
  out$beta1 <- .fit_fractional_probit(data,persist,params$beta1,FALSE)
  names(out$beta0) <- names(params$beta0); names(out$beta1) <- names(params$beta1)
  mismatch <- lapply(1:2,function(k)lapply(1:4,function(tt){
    h1 <- e$state[[k]][[tt]][,2]; pt <- e$posterior_type[,k]
    ifelse(data$y[,tt]==1,pt-h1,h1)}))
  delta_q <- function(delta,gradient=FALSE) {
    q <- 0; grad <- setNames(numeric(length(delta)), names(delta))
    for(k in 1:2)for(tt in 1:4) {
      z <- data$Z[[tt]]; pi <- .fmm_covinc_pi(z,delta,k)
      m <- mismatch[[k]][[tt]]; pt <- e$posterior_type[,k]; w <- data$weight
      q <- q+sum(w*(m*log(pi)+(pt-m)*log1p(-pi)))
      if(gradient) {
        dpi <- pi*(1-2*pi)
        score <- w*(m/pi-(pt-m)/(1-pi))*dpi
        grad <- .fmm_covinc_add_delta_score(grad, z, score, delta, k)
      }
    }
    scale <- sum(data$weight)
    if(gradient)grad/scale else q/scale
  }
  out$delta <- optim(params$delta,function(d)-delta_q(d),function(d)-delta_q(d,TRUE),
    method="BFGS",control=list(maxit=20L,reltol=1e-10))$par
  names(out$delta) <- names(params$delta)
  out
}

.pack_fmm_covinc_4w <- function(p)c(
  alpha_A=qlogis(unname(p$alpha[1])),alpha_B=qlogis(unname(p$alpha[2])),
  phi=qlogis(p$phi),setNames(p$beta0,paste0("entry_",names(p$beta0))),
  setNames(p$beta1,paste0("persistence_",names(p$beta1))),
  setNames(p$delta,paste0("misclassification_",names(p$delta))))

.unpack_fmm_covinc_4w <- function(z,template) {
  z <- as.numeric(z); n0 <- length(template$beta0); n1 <- length(template$beta1)
  nd <- length(template$delta)
  out <- template
  out$alpha <- setNames(plogis(z[1:2]),c("A","B")); out$phi <- plogis(z[3])
  out$beta0 <- setNames(z[3+seq_len(n0)],names(template$beta0))
  out$beta1 <- setNames(z[3+n0+seq_len(n1)],names(template$beta1))
  out$delta <- setNames(z[3+n0+n1+seq_len(nd)],names(template$delta)); out
}

.direct_fmm_covinc_objective <- function(data,template) {
  cache_z <- cache <- NULL; W <- sum(data$weight)
  evaluate <- function(z) {
    if(!is.null(cache_z)&&identical(as.numeric(z),cache_z))return(cache)
    p <- .unpack_fmm_covinc_4w(z,template); e <- e_step_fmm_covinc_4w(data,p)
    score <- numeric(length(z)); score[1] <- sum(data$weight*(
      e$state[[1]][[1]][,2]-e$posterior_type[,1]*p$alpha[1]))/W
    score[2] <- sum(data$weight*(
      e$state[[2]][[1]][,2]-e$posterior_type[,2]*p$alpha[2]))/W
    score[3] <- sum(data$weight*(e$posterior_type[,1]-p$phi))/W
    n0 <- length(p$beta0); n1 <- length(p$beta1)
    s0 <- numeric(n0); s1 <- numeric(n1)
    for(k in 1:2)for(tt in 1:3) {
      xi <- e$xi[[k]][[tt]]
      r0 <- xi[,1]+xi[,2]; q0 <- e$q0[[k]][[tt]]
      sc0 <- data$weight*dnorm(qnorm(q0))*(xi[,2]/q0-(r0-xi[,2])/(1-q0))
      s0[k] <- s0[k]+sum(sc0)
      s0[-(1:2)] <- s0[-(1:2)]+as.vector(crossprod(
        data$X[[tt]][,data$entry_active,drop=FALSE],sc0))
      r1 <- xi[,3]+xi[,4]; q1 <- e$q1[[k]][[tt]]
      sc1 <- data$weight*dnorm(qnorm(q1))*(xi[,4]/q1-(r1-xi[,4])/(1-q1))
      s1[k] <- s1[k]+sum(sc1)
      active1 <- data$persistence_active %||% rep(TRUE, length(data$entry_active))
      s1[-(1:2)] <- s1[-(1:2)]+as.vector(crossprod(
        data$X[[tt]][,active1,drop=FALSE],sc1))
    }
    score[3+seq_len(n0)] <- s0/W; score[3+n0+seq_len(n1)] <- s1/W
    nd <- length(p$delta); sd <- setNames(numeric(nd), names(p$delta))
    for(k in 1:2)for(tt in 1:4) {
      h1 <- e$state[[k]][[tt]][,2]; pt <- e$posterior_type[,k]
      m <- ifelse(data$y[,tt]==1,pt-h1,h1)
      pi <- .fmm_covinc_pi(data$Z[[tt]],p$delta,k); dpi <- pi*(1-2*pi)
      sc <- data$weight*(m/pi-(pt-m)/(1-pi))*dpi
      sd <- .fmm_covinc_add_delta_score(sd, data$Z[[tt]], sc, p$delta, k)
    }
    score[3+n0+n1+seq_len(nd)] <- sd/W
    cache_z <<- as.numeric(z); cache <<- list(value=-e$loglik/W,gradient=-score,e=e,params=p)
    cache
  }
  list(fn=function(z)evaluate(z)$value,gr=function(z)evaluate(z)$gradient,
    details=function(z)evaluate(z))
}

.normalize_fmm_covinc_labels <- function(p,data) {
  active1 <- data$persistence_active %||% rep(TRUE, length(data$entry_active))
  avg <- sapply(1:2,function(k)mean(sapply(1:3,function(tt)
    weighted.mean(.fmm_covinc_transition(
      data$X[[tt]],p$beta1,k,active1),data$weight))))
  if(avg[1]>=avg[2])return(p)
  p$alpha <- rev(p$alpha); p$phi <- 1-p$phi
  p$beta0[1:2] <- rev(p$beta0[1:2]); p$beta1[1:2] <- rev(p$beta1[1:2])
  if (.fmm_covinc_type_specific_delta(p$delta))
    p$delta[c("intercept_A", "intercept_B")] <-
      rev(p$delta[c("intercept_A", "intercept_B")])
  p
}

fit_fmm_covariates_inconsistency_4w <- function(data,starts=NULL,random_starts=6L,
    seed=20260818L,max_iter=600L,ll_tol=1e-10,param_tol=2e-6,verbose=1L,
    em_warm_start=TRUE) {
  if(is.null(starts))starts <- list(.initial_fmm_covinc_4w(data))
  set.seed(seed); base <- starts[[1]]
  for(s in seq_len(random_starts)) {
    p <- base
    p$alpha <- plogis(qlogis(p$alpha)+rnorm(2,0,.5))
    p$phi <- plogis(qlogis(p$phi)+rnorm(1,0,.8))
    p$beta0[1:2] <- p$beta0[1:2]+rnorm(2,0,.8)
    p$beta1[1:2] <- p$beta1[1:2]+rnorm(2,0,.8)
    delta_sd <- rep(.3, length(p$delta)); names(delta_sd) <- names(p$delta)
    intercepts <- intersect(names(delta_sd), c("intercept", "intercept_A", "intercept_B"))
    delta_sd[intercepts] <- .6
    p$delta <- p$delta + rnorm(length(p$delta), 0, delta_sd)
    starts[[length(starts)+1L]] <- p
  }
  run_one <- function(p,loose=FALSE) {
    # One EM update supplies good common-slope starting values; analytical
    # observed-likelihood scores then converge much faster than nested EM.
    if (em_warm_start) {
      e0 <- e_step_fmm_covinc_4w(data,p)
      p <- m_step_fmm_covinc_4w(data,e0,p)
    }
    objective <- .direct_fmm_covinc_objective(data,p)
    opt <- optim(.pack_fmm_covinc_4w(p),objective$fn,objective$gr,method="BFGS",
      control=list(maxit=if(loose)15L else max_iter,reltol=if(loose)1e-6 else ll_tol))
    details <- objective$details(opt$par)
    p <- .normalize_fmm_covinc_labels(details$params,data)
    e <- e_step_fmm_covinc_4w(data,p)
    list(params=p,e=e,loglik=e$loglik,iterations=opt$counts[[1]],
      converged=opt$convergence==0L,fixed_point_error=max(abs(objective$gr(opt$par))),
      optimizer_code=opt$convergence)
  }
  screened <- lapply(seq_along(starts),function(i){
    if(verbose)cat(sprintf("Screening start %d/%d\n",i,length(starts)))
    tryCatch(run_one(starts[[i]],TRUE),error=function(e)NULL)})
  screened <- Filter(Negate(is.null),screened)
  if(!length(screened))stop("All starts failed")
  screened <- screened[order(vapply(screened,`[[`,numeric(1),"loglik"),decreasing=TRUE)]
  finalists <- lapply(screened[seq_len(1L)],function(x) {
    if (verbose) cat("Refining best screened start\n")
    run_one(x$params,FALSE)
  })
  best <- finalists[[which.max(vapply(finalists,`[[`,numeric(1),"loglik"))]]
  best$candidates <- data.frame(loglik=vapply(finalists,`[[`,numeric(1),"loglik"),
    converged=vapply(finalists,`[[`,logical(1),"converged"))
  best$sample <- list(n=data$n_original,cells=nrow(data$y),weight_sum=sum(data$weight))
  best$model <- "two_type_4w_ar1_common_covariate_slopes_inconsistency_pi"
  if(verbose)cat(sprintf("Final log likelihood %.3f; converged=%s; fixed-point error=%.2e\n",
    best$loglik,best$converged,best$fixed_point_error))
  best
}

summarize_fmm_covariates_inconsistency_4w <- function(data,fit) {
  p <- fit$params; e <- fit$e
  type_rows <- lapply(1:2,function(k) {
    risk0 <- risk1 <- entry_flow <- exit_flow <- 0
    for(tt in 1:3) {
      z <- e$xi[[k]][[tt]]; r0 <- z[,1]+z[,2]; r1 <- z[,3]+z[,4]
      risk0 <- risk0+sum(data$weight*r0); risk1 <- risk1+sum(data$weight*r1)
      entry_flow <- entry_flow+sum(data$weight*r0*e$q0[[k]][[tt]])
      exit_flow <- exit_flow+sum(data$weight*r1*(1-e$q1[[k]][[tt]]))
    }
    data.frame(type=c("A (stable)","B (mobile)")[k],share=c(p$phi,1-p$phi)[k],
      initial_employment=p$alpha[k],risk_weighted_entry=entry_flow/risk0,
      risk_weighted_exit=exit_flow/risk1)
  })
  type_summary <- do.call(rbind,type_rows)
  pattern_names <- c("consistent", "age_inconsistent",
                     "education_inconsistent", "both_inconsistent")
  patterns <- matrix(0, 4L, ncol(data$Z[[1]]),
    dimnames=list(pattern_names, colnames(data$Z[[1]])))
  patterns[, 1L] <- 1
  if ("age_inconsistent" %in% colnames(patterns))
    patterns[c("age_inconsistent", "both_inconsistent"), "age_inconsistent"] <- 1
  if ("education_inconsistent" %in% colnames(patterns))
    patterns[c("education_inconsistent", "both_inconsistent"),
             "education_inconsistent"] <- 1
  pi_table <- if (.fmm_covinc_type_specific_delta(p$delta)) {
    do.call(rbind,lapply(1:2,function(k)data.frame(
      type=c("A (stable)","B (mobile)")[k],pattern=rownames(patterns),
      probability=vapply(seq_len(nrow(patterns)), function(jj)
        .fmm_covinc_pi(patterns[jj,,drop=FALSE], p$delta, k), numeric(1)))))
  } else {
    data.frame(type="common",pattern=rownames(patterns),
      probability=vapply(seq_len(nrow(patterns)), function(jj)
        .fmm_covinc_pi(patterns[jj,,drop=FALSE], p$delta, 1L), numeric(1)))
  }
  list(type_summary=type_summary,entry_coefficients=p$beta0,
    persistence_coefficients=p$beta1,misclassification_coefficients=p$delta,
    misclassification_probabilities=pi_table)
}

.score_rows_fmm_covinc_4w <- function(data,p,e) {
  n <- nrow(data$y); n0 <- length(p$beta0); n1 <- length(p$beta1)
  nd <- length(p$delta); S <- matrix(0,n,3+n0+n1+nd)
  colnames(S) <- names(.pack_fmm_covinc_4w(p))
  S[,1] <- e$state[[1]][[1]][,2]-e$posterior_type[,1]*p$alpha[1]
  S[,2] <- e$state[[2]][[1]][,2]-e$posterior_type[,2]*p$alpha[2]
  S[,3] <- e$posterior_type[,1]-p$phi
  for(k in 1:2)for(tt in 1:3) {
    xi <- e$xi[[k]][[tt]]
    r0 <- xi[,1]+xi[,2]; q0 <- e$q0[[k]][[tt]]
    sc0 <- dnorm(qnorm(q0))*(xi[,2]/q0-(r0-xi[,2])/(1-q0))
    S[,3+k] <- S[,3+k]+sc0
    idx0 <- 3+seq_len(n0)
    S[,idx0[-(1:2)]] <- S[,idx0[-(1:2)],drop=FALSE]+
      data$X[[tt]][,data$entry_active,drop=FALSE]*sc0
    r1 <- xi[,3]+xi[,4]; q1 <- e$q1[[k]][[tt]]
    sc1 <- dnorm(qnorm(q1))*(xi[,4]/q1-(r1-xi[,4])/(1-q1))
    idx1 <- 3+n0+seq_len(n1)
    S[,idx1[k]] <- S[,idx1[k]]+sc1
    active1 <- data$persistence_active %||% rep(TRUE, length(data$entry_active))
    S[,idx1[-(1:2)]] <- S[,idx1[-(1:2)],drop=FALSE]+
      data$X[[tt]][,active1,drop=FALSE]*sc1
  }
  idxd <- 3+n0+n1+seq_len(nd)
  for(k in 1:2)for(tt in 1:4) {
    h1 <- e$state[[k]][[tt]][,2]; pt <- e$posterior_type[,k]
    m <- ifelse(data$y[,tt]==1,pt-h1,h1)
    pi <- .fmm_covinc_pi(data$Z[[tt]],p$delta,k); dpi <- pi*(1-2*pi)
    sc <- (m/pi-(pt-m)/(1-pi))*dpi
    if (.fmm_covinc_type_specific_delta(p$delta)) {
      intercept_name <- c("intercept_A", "intercept_B")[k]
      S[,idxd[match(intercept_name, names(p$delta))]] <-
        S[,idxd[match(intercept_name, names(p$delta))]] + sc
      slope_names <- colnames(data$Z[[tt]])[-1L]
      if (length(slope_names)) {
        slope_idx <- idxd[match(slope_names, names(p$delta))]
        S[,slope_idx] <- S[,slope_idx,drop=FALSE] +
          data$Z[[tt]][,slope_names,drop=FALSE] * sc
      }
    } else {
      design_idx <- idxd[match(colnames(data$Z[[tt]]), names(p$delta))]
      S[,design_idx] <- S[,design_idx,drop=FALSE] + data$Z[[tt]] * sc
    }
  }
  S
}

analytical_se_fmm_covinc_4w <- function(data,fit,finite_sample=TRUE) {
  p <- fit$params; z <- .pack_fmm_covinc_4w(p)
  objective <- .direct_fmm_covinc_objective(data,p)
  bread <- optimHess(z,objective$fn,objective$gr); bread <- (bread+t(bread))/2
  eig <- eigen(bread,symmetric=TRUE,only.values=TRUE)$values
  scores <- .score_rows_fmm_covinc_4w(data,p,fit$e); W <- sum(data$weight)
  meat <- crossprod(scores,scores*data$weight_sq)/W^2
  inv <- solve(bread); vcov <- inv%*%meat%*%inv
  K <- length(z); n <- data$n_original
  if(finite_sample&&n>K)vcov <- vcov*n/(n-K)
  list(coefficient_table=data.frame(parameter=names(z),estimate=z,
    se=sqrt(pmax(diag(vcov),0)),row.names=NULL),vcov=vcov,bread=bread,meat=meat,
    diagnostics=list(rank=qr(bread,tol=1e-8)$rank,parameters=K,
      minimum_eigenvalue=min(eig),condition=max(eig)/min(eig),
      max_abs_score=max(abs(colSums(scores*data$weight)/W))),
    method="survey_weighted_sandwich_analytical_scores")
}

derived_se_fmm_covinc_4w <- function(data,fit,inference,rel_step=1e-4) {
  template <- fit$params; z <- .pack_fmm_covinc_4w(template)
  quantities <- function(zz) {
    p <- .unpack_fmm_covinc_4w(zz,template)
    e <- e_step_fmm_covinc_4w(data,p)
    s <- summarize_fmm_covariates_inconsistency_4w(data,list(params=p,e=e))
    posterior_mean_pi <- sapply(1:2, function(k) {
      pt <- e$posterior_type[, k]
      numerator <- sum(vapply(1:4, function(tt)
        sum(data$weight * pt * .fmm_covinc_pi(data$Z[[tt]], p$delta, k)),
        numeric(1)))
      numerator / (4 * sum(data$weight * pt))
    })
    c(type_A_share=p$phi,type_B_share=1-p$phi,
      A_initial_employment=unname(p$alpha[1]),B_initial_employment=unname(p$alpha[2]),
      A_risk_weighted_entry=s$type_summary$risk_weighted_entry[1],
      A_risk_weighted_exit=s$type_summary$risk_weighted_exit[1],
      B_risk_weighted_entry=s$type_summary$risk_weighted_entry[2],
      B_risk_weighted_exit=s$type_summary$risk_weighted_exit[2],
      A_posterior_mean_misclassification=posterior_mean_pi[1],
      B_posterior_mean_misclassification=posterior_mean_pi[2],
      setNames(s$misclassification_probabilities$probability,
        paste0("pi_",gsub(" .*","",s$misclassification_probabilities$type),"_",
          s$misclassification_probabilities$pattern)))
  }
  estimate <- quantities(z)
  J <- matrix(NA_real_,length(estimate),length(z),dimnames=list(names(estimate),names(z)))
  for(j in seq_along(z)) {
    h <- rel_step*max(1,abs(z[j])); zp <- zm <- z; zp[j] <- zp[j]+h; zm[j] <- zm[j]-h
    J[,j] <- (quantities(zp)-quantities(zm))/(2*h)
  }
  vcov <- J%*%inference$vcov%*%t(J)
  list(summary=data.frame(quantity=names(estimate),estimate=unname(estimate),
    se=sqrt(pmax(diag(vcov),0))),vcov=vcov,jacobian=J,method="delta_method")
}
