# Tail-safe piecewise-constant duration hazards for the epsilon model.

.piecewise_eps_unpack <- function(z, pi_cap=.49, eps_floor=1e-4,
                                  eps_cap=.95, timegap_contamination=FALSE,
                                  eps_d_floor=1e-6, eps_d_cap=.95) {
  hg <- exp(unname(z[paste0("log_hg",1:5)]))
  hd <- exp(unname(z[paste0("log_hd",1:5)]))
  p_exit0 <- .duration_transition_probability(0,hg,0)
  p_entry0 <- .duration_transition_probability(0,hd,0)
  out <- list(alpha=plogis(unname(z["alpha"])),theta0=p_entry0,
    theta1=1-p_exit0,pi=pi_cap*plogis(unname(z["pi"])),
    eps=eps_floor+(eps_cap-eps_floor)*plogis(unname(z["eps"])),
    lambda_g=hg,beta_g=0,lambda_d=hd,beta_d=0,
    hazard_model="piecewise")
  if (timegap_contamination)
    out$eps_d <- eps_d_floor + (eps_d_cap-eps_d_floor) *
      plogis(unname(z["eps_d"]))
  out
}

.piecewise_eps_pack <- function(params, pi_cap=.49, eps_floor=1e-4,
                                eps_cap=.95, timegap_contamination=FALSE,
                                eps_d_floor=1e-6, eps_d_cap=.95) {
  if (length(params$lambda_g)==1L) params$lambda_g <- rep(params$lambda_g,5L)
  if (length(params$lambda_d)==1L) params$lambda_d <- rep(params$lambda_d,5L)
  out <- c(alpha=qlogis(.bound01(params$alpha)),
    pi=qlogis(.bound01(params$pi/pi_cap)),
    eps=qlogis(.bound01((params$eps-eps_floor)/(eps_cap-eps_floor))),
    setNames(log(params$lambda_g),paste0("log_hg",1:5)),
    setNames(log(params$lambda_d),paste0("log_hd",1:5)))
  if (timegap_contamination) {
    eps_d <- if (is.null(params$eps_d)) .05 else params$eps_d
    out <- append(out,setNames(
      qlogis(.bound01((eps_d-eps_d_floor)/(eps_d_cap-eps_d_floor))),
      "eps_d"),after=3L)
  }
  out
}

piecewise_start_from_power <- function(params,
  representative=c(.125,.625,2,4,7.5)) {
  out <- params
  bg <- if (is.null(params$beta_g)) 0 else params$beta_g
  bd <- if (is.null(params$beta_d)) 0 else params$beta_d
  out$lambda_g <- params$lambda_g[1] * (1+representative)^bg
  out$lambda_d <- params$lambda_d[1] * (1+representative)^bd
  out$beta_g <- 0; out$beta_d <- 0
  out
}

fit_eps_piecewise_hazard <- function(df,start,maxit=500L,reltol=1e-9,
                                     pgtol=1e-7,method="L-BFGS-B",
                                     verbose=0L,
                                     timegap_contamination=FALSE,
                                     timegap_contamination_model="marginal",
                                     timegap_local_decay=1) {
  validate_df_eps(df)
  z0 <- .piecewise_eps_pack(start,
    timegap_contamination=timegap_contamination)
  lower <- rep(-Inf,length(z0)); upper <- rep(Inf,length(z0))
  names(lower) <- names(upper) <- names(z0)
  hz <- c(paste0("log_hg",1:5),paste0("log_hd",1:5))
  lower[hz] <- log(1e-4); upper[hz] <- log(20)
  total_weight <- sum(df$weight)
  objective <- function(z) {
    p <- .piecewise_eps_unpack(z,
      timegap_contamination=timegap_contamination)
    if (timegap_contamination) {
      p$timegap_contamination_model <- timegap_contamination_model
      p$timegap_local_decay <- timegap_local_decay
    }
    value <- tryCatch(e_step_eps(df,p,check_df=FALSE,suff_stats=FALSE)$loglik,
      error=function(e) NA_real_)
    if (!is.finite(value)) return(1e100)
    -value/total_weight
  }
  control <- list(maxit=maxit,reltol=reltol,
    ndeps=rep(1e-3,length(z0)),trace=as.integer(verbose))
  if (identical(method,"L-BFGS-B")) {
    control$factr <- reltol/.Machine$double.eps
    control$pgtol <- pgtol
    opt <- optim(z0,objective,method=method,lower=lower,upper=upper,
      control=control)
  } else {
    opt <- optim(z0,objective,method=method,control=control)
  }
  params <- .piecewise_eps_unpack(opt$par,
    timegap_contamination=timegap_contamination)
  if (timegap_contamination) {
    params$timegap_contamination_model <- timegap_contamination_model
    params$timegap_local_decay <- timegap_local_decay
  }
  estep <- e_step_eps(df,params,check_df=FALSE,suff_stats=FALSE)
  list(params=params,loglik=estep$loglik,gamma=estep$gamma,
    convergence=opt$convergence,message=opt$message,
    iterations=unname(opt$counts["function"]),objective=opt$value,
    gradient_max=NA_real_,par_unconstrained=opt$par,
    objective_function=objective)
}

timegap_contamination_diagnostics <- function(df, fit) {
  eps_d <- fit$params$eps_d
  if (is.null(eps_d)) stop("fit does not contain eps_d")
  hmat <- latent_histories(); gamma <- fit$gamma; wi <- df$weight
  rows <- lapply(1:2, function(t) {
    s_prev <- df[[paste0("y",t)]]
    s_curr <- df[[paste0("y",t+1L)]]
    c_prev <- df[[paste0("timegap_cat",t)]]
    c_curr <- df[[paste0("timegap_cat",t+1L)]]
    observed <- s_prev==0L & s_curr==0L
    numer <- denom <- impossible <- 0
    for (j in seq_len(nrow(hmat))) {
      if (hmat[j,t]!=0L || hmat[j,t+1L]!=0L) next
      use <- observed
      if (!any(use)) next
      lc <- log_emission_transition_d(c_curr[use],c_prev[use],
        fit$params$lambda_d,fit$params$beta_d)
      lp <- log_emission_interval_d(c_curr[use],fit$params$lambda_d,
        fit$params$beta_d)
      omega <- .omega_rho(lc,lp,eps_d)
      w <- wi[use]*gamma[use,j]
      numer <- numer+sum(w*omega)
      denom <- denom+sum(w)
      impossible <- impossible+sum(w*is.infinite(lc))
    }
    data.frame(transition=paste0(t,"->",t+1L),
      posterior_eligible_weight=denom,
      posterior_contaminated_share=numer/denom,
      clock_impossible_share=impossible/denom)
  })
  do.call(rbind,rows)
}

timegap_contamination_decomposition <- function(df, fit) {
  eps_d <- fit$params$eps_d
  if (is.null(eps_d)) stop("fit does not contain eps_d")
  if (nrow(fit$gamma)!=nrow(df)) stop("fit and data have different row counts")
  hmat <- latent_histories(); wi <- df$weight
  detail <- lapply(1:2,function(t) {
    observed <- df[[paste0("y",t)]]==0L & df[[paste0("y",t+1L)]]==0L
    eligible_histories <- which(hmat[,t]==0L & hmat[,t+1L]==0L)
    latent_weight <- wi*rowSums(fit$gamma[,eligible_histories,drop=FALSE])
    prev <- df[[paste0("timegap_cat",t)]]
    curr <- df[[paste0("timegap_cat",t+1L)]]
    lc <- log_emission_transition_d(curr,prev,fit$params$lambda_d,
      fit$params$beta_d)
    lp <- log_emission_interval_d(curr,fit$params$lambda_d,
      fit$params$beta_d)
    omega <- .omega_rho(lc,lp,eps_d)
    clock_feasible <- is.finite(lc)
    mechanism <- ifelse(clock_feasible,"Clock-feasible",
      ifelse(curr==1L,"Reset-compatible","Not reset-compatible"))
    data.frame(transition=paste0(t,"->",t+1L),prev_category=prev,
      curr_category=curr,mechanism=mechanism,posterior_contamination=omega,
      posterior_weight=ifelse(observed,latent_weight,0))
  })
  detail <- do.call(rbind,detail)
  use <- detail$posterior_weight>0
  detail <- detail[use,,drop=FALSE]
  groups <- split(seq_len(nrow(detail)),
    interaction(detail$transition,detail$mechanism,drop=TRUE))
  summary <- do.call(rbind,lapply(groups,function(idx) {
    w <- detail$posterior_weight[idx]
    data.frame(transition=detail$transition[idx[1]],
      mechanism=detail$mechanism[idx[1]],posterior_weight=sum(w),
      posterior_contamination_share=sum(w*detail$posterior_contamination[idx])/sum(w))
  }))
  totals <- aggregate(posterior_weight~transition,detail,sum)
  names(totals)[2] <- "total_weight"
  summary <- merge(summary,totals,by="transition",sort=FALSE)
  summary$share_of_latent_continuations <- summary$posterior_weight/summary$total_weight
  summary <- summary[order(summary$transition,
    match(summary$mechanism,c("Clock-feasible","Reset-compatible",
      "Not reset-compatible"))),]

  cells <- aggregate(cbind(posterior_weight,
    posterior_contamination_mass=posterior_weight*posterior_contamination)~
      transition+prev_category+curr_category+mechanism,detail,sum)
  cells$share_within_transition <- cells$posterior_weight/
    ave(cells$posterior_weight,cells$transition,FUN=sum)
  cells$posterior_contamination_share <-
    cells$posterior_contamination_mass/cells$posterior_weight
  cells <- cells[order(cells$transition,-cells$posterior_weight),]
  list(summary=summary,cells=cells,detail=detail)
}

fit_eps_piecewise_multistart <- function(df, starts, ...) {
  fits <- lapply(starts,function(s) tryCatch(
    fit_eps_piecewise_hazard(df,s,...),
    error=function(e) structure(list(error=conditionMessage(e)),
      class="piecewise_eps_error")))
  tab <- do.call(rbind,lapply(seq_along(fits),function(k) {
    x <- fits[[k]]
    if (inherits(x,"piecewise_eps_error")) return(data.frame(start=k,
      convergence=NA_integer_,loglik=NA_real_,gradient_max=NA_real_,
      min_hazard=NA_real_,max_hazard=NA_real_,error=x$error))
    data.frame(start=k,convergence=x$convergence,loglik=x$loglik,
      gradient_max=x$gradient_max,
      min_hazard=min(c(x$params$lambda_g,x$params$lambda_d)),
      max_hazard=max(c(x$params$lambda_g,x$params$lambda_d)),error=NA_character_)
  }))
  eligible <- which(is.finite(tab$loglik))
  if (!length(eligible)) stop("No piecewise-hazard fit returned a finite likelihood")
  best <- eligible[which.max(tab$loglik[eligible])]
  grad <- .duration_eps_numeric_gradient(fits[[best]]$objective_function,
    fits[[best]]$par_unconstrained,step=1e-3)
  fits[[best]]$gradient_max <- max(abs(grad))
  names(grad) <- names(fits[[best]]$par_unconstrained)
  fits[[best]]$gradient <- grad
  fits[[best]]$gradient_parameter <- names(grad)[which.max(abs(grad))]
  tab$gradient_max[best] <- fits[[best]]$gradient_max
  list(best=fits[[best]],best_start=best,fits=fits,summary=tab)
}

piecewise_hazard_table <- function(fit) {
  data.frame(interval=c("0-3 months","3-12 months","1-3 years",
    "3-5 years","5+ years"),lower_years=.DURATION_PIECEWISE_KNOTS[1:5],
    exit_hazard=fit$params$lambda_g,entry_hazard=fit$params$lambda_d)
}
