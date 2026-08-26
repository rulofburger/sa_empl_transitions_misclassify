# ==============================================================================
# Two-type finite-mixture extension of the four-wave AR(1) model
# ==============================================================================

.fmm_ar1_4w_names <- function(model_type, stationary) {
  out <- c("theta0_A","theta1_A")
  if (!stationary) out <- c(out,"alpha_A")
  out <- c(out,"theta0_B","theta1_B")
  if (!stationary) out <- c(out,"alpha_B")
  out <- c(out,"phi")
  if (model_type=="symmetric") out <- c(out,"pi")
  out
}

pack_fmm_ar1_4w <- function(params,model_type="symmetric",stationary=FALSE) {
  clamp <- function(x,eps=1e-8) pmin(pmax(x,eps),1-eps)
  nms <- .fmm_ar1_4w_names(model_type,stationary)
  out <- setNames(numeric(length(nms)),nms)
  for(nm in setdiff(nms,"pi")) out[nm] <- qlogis(clamp(params[[nm]]))
  if(model_type=="symmetric") out["pi"] <- qlogis(clamp(params$pi/.5))
  out
}

unpack_fmm_ar1_4w <- function(eta,model_type="symmetric",stationary=FALSE) {
  eta <- setNames(as.numeric(eta),.fmm_ar1_4w_names(model_type,stationary))
  p <- as.list(plogis(eta))
  if(stationary) {
    p$alpha_A <- stationary_alpha_ar1(p$theta0_A,p$theta1_A)
    p$alpha_B <- stationary_alpha_ar1(p$theta0_B,p$theta1_B)
  }
  if(model_type=="symmetric") p$pi <- .5*plogis(eta["pi"])
  p
}

resolve_fmm_ar1_4w_labels <- function(p) {
  if(p$theta1_A >= p$theta1_B) return(p)
  out <- p
  for(stem in c("theta0","theta1","alpha")) {
    out[[paste0(stem,"_A")]] <- p[[paste0(stem,"_B")]]
    out[[paste0(stem,"_B")]] <- p[[paste0(stem,"_A")]]
  }
  out$phi <- 1-p$phi
  out
}

fmm_ar1_4w_cell_probabilities <- function(params,model_type="symmetric") {
  pi <- if(model_type=="symmetric") params$pi else 0
  pA <- ar1_4w_cell_probabilities(list(theta0=params$theta0_A,
    theta1=params$theta1_A,alpha=params$alpha_A,pi=pi),model_type)
  pB <- ar1_4w_cell_probabilities(list(theta0=params$theta0_B,
    theta1=params$theta1_B,alpha=params$alpha_B,pi=pi),model_type)
  params$phi*pA+(1-params$phi)*pB
}

initial_fmm_ar1_4w <- function(model_type="symmetric",stationary=FALSE,
                               stable_share=.65) {
  p <- list(theta0_A=.025,theta1_A=.975,alpha_A=.55,
            theta0_B=.16,theta1_B=.82,alpha_B=.38,phi=stable_share)
  if(stationary) {
    p$alpha_A <- stationary_alpha_ar1(p$theta0_A,p$theta1_A)
    p$alpha_B <- stationary_alpha_ar1(p$theta0_B,p$theta1_B)
  }
  if(model_type=="symmetric") p$pi <- .015
  p
}

.random_fmm_ar1_4w_start <- function(model_type,stationary) {
  p <- list(theta0_A=runif(1,.003,.10),theta1_A=runif(1,.90,.997),
            alpha_A=runif(1,.30,.75),theta0_B=runif(1,.04,.40),
            theta1_B=runif(1,.50,.92),alpha_B=runif(1,.20,.70),
            phi=runif(1,.10,.90))
  if(stationary) {
    p$alpha_A <- stationary_alpha_ar1(p$theta0_A,p$theta1_A)
    p$alpha_B <- stationary_alpha_ar1(p$theta0_B,p$theta1_B)
  }
  if(model_type=="symmetric") p$pi <- runif(1,.002,.08)
  p
}

.em_fmm_ar1_4w_cells <- function(cells,start,model_type,stationary,
                                  max_iter=1500L,ll_tol=1e-10,param_tol=1e-7) {
  cache <- .ar1_4w_cache(); h <- cache$latent; w <- cells$weight
  W <- sum(w); p <- start
  update_type <- function(wh,current) {
    mass <- sum(wh); C1 <- sum(wh*h[,1]); C0 <- mass-C1
    D0<-D1<-T01<-T11<-0
    for(tt in 2:4) {
      D0 <- D0+sum(wh*(h[,tt-1]==0)); D1 <- D1+sum(wh*(h[,tt-1]==1))
      T01 <- T01+sum(wh*(h[,tt-1]==0 & h[,tt]==1))
      T11 <- T11+sum(wh*(h[,tt-1]==1 & h[,tt]==1))
    }
    if(!stationary) return(list(theta0=T01/D0,theta1=T11/D1,alpha=C1/mass,mass=mass))
    q <- function(z) {
      t0<-plogis(z[1]);t1<-plogis(z[2]);a<-stationary_alpha_ar1(t0,t1)
      -(C1*log(a)+C0*log(1-a)+T01*log(t0)+(D0-T01)*log(1-t0)+
        T11*log(t1)+(D1-T11)*log(1-t1))/mass
    }
    opt<-optim(qlogis(c(current$theta0,current$theta1)),q,method="BFGS",
               control=list(maxit=100L,reltol=1e-12))
    t0<-plogis(opt$par[1]);t1<-plogis(opt$par[2])
    list(theta0=t0,theta1=t1,alpha=stationary_alpha_ar1(t0,t1),mass=mass)
  }
  ll_prev <- -Inf
  for(iter in seq_len(max_iter)) {
    pi <- if(model_type=="symmetric") p$pi else 0
    emission <- (1-pi)^(4L-cache$mismatch)*pi^cache$mismatch
    prior_one <- function(stem) {
      t0<-p[[paste0("theta0_",stem)]];t1<-p[[paste0("theta1_",stem)]]
      a<-p[[paste0("alpha_",stem)]];pr<-ifelse(h[,1]==1,a,1-a)
      for(tt in 2:4){q<-ifelse(h[,tt-1]==1,t1,t0);pr<-pr*ifelse(h[,tt]==1,q,1-q)}
      pr
    }
    jA<-sweep(emission,2,prior_one("A"),"*")*p$phi
    jB<-sweep(emission,2,prior_one("B"),"*")*(1-p$phi)
    den<-rowSums(jA)+rowSums(jB);ll<-sum(w*log(pmax(den,1e-300)))
    gA<-jA/den;gB<-jB/den
    whA<-colSums(gA*w);whB<-colSums(gB*w)
    A<-update_type(whA,list(theta0=p$theta0_A,theta1=p$theta1_A))
    B<-update_type(whB,list(theta0=p$theta0_B,theta1=p$theta1_B))
    newp<-list(theta0_A=A$theta0,theta1_A=A$theta1,alpha_A=A$alpha,
               theta0_B=B$theta0,theta1_B=B$theta1,alpha_B=B$alpha,phi=A$mass/W)
    if(model_type=="symmetric") {
      mismatches <- sum((gA+gB)*w*cache$mismatch)
      newp$pi <- pmin(pmax(mismatches/(4*W),1e-8),.49)
    }
    newp<-resolve_fmm_ar1_4w_labels(newp)
    oldeta<-pack_fmm_ar1_4w(resolve_fmm_ar1_4w_labels(p),model_type,stationary)
    neweta<-pack_fmm_ar1_4w(newp,model_type,stationary)
    change<-max(abs(neweta-oldeta));rel<-abs(ll-ll_prev)/(abs(ll_prev)+1e-16)
    p<-newp
    if(iter>1 && rel<ll_tol && change<param_tol) break
    ll_prev<-ll
  }
  probs<-fmm_ar1_4w_cell_probabilities(p,model_type)
  list(params=p,eta=pack_fmm_ar1_4w(p,model_type,stationary),
       loglik=sum(w*log(pmax(probs,1e-300))),iterations=iter,
       converged=iter<max_iter,fixed_point_error=change)
}

fit_fmm_ar1_4w <- function(df,model_type="symmetric",stationary=FALSE,
                            starts=NULL,random_starts=60L,seed=20260819L,
                            screen_maxit=500L,refine_top=8L,
                            maxit=4000L,reltol=1e-12,verbose=1L) {
  if(!model_type %in% c("none","symmetric")) stop("invalid model_type")
  cells <- collapse_ar1_4w_cells(df)
  if(is.null(starts)) starts <- list(initial_fmm_ar1_4w(model_type,stationary))
  set.seed(seed)
  if(random_starts>0) for(i in seq_len(random_starts))
    starts[[length(starts)+1L]] <- .random_fmm_ar1_4w_start(model_type,stationary)
  fn <- function(z) {
    p <- unpack_fmm_ar1_4w(z,model_type,stationary)
    probs <- pmax(fmm_ar1_4w_cell_probabilities(p,model_type),1e-300)
    -sum(cells$weight*log(probs))/sum(cells$weight)
  }
  # Multi-start EM cheaply moves each candidate toward a mixture maximum.
  screened <- lapply(starts,function(s) tryCatch(
    .em_fmm_ar1_4w_cells(cells,s,model_type,stationary,
      max_iter=max(200L,screen_maxit),ll_tol=1e-9,param_tol=1e-6),error=function(e)NULL))
  screened <- Filter(function(x)!is.null(x)&&is.finite(x$loglik),screened)
  if(!length(screened)) stop("all FMM starts failed during screening")
  screened <- screened[order(vapply(screened,`[[`,numeric(1),"loglik"),decreasing=TRUE)]
  screened <- screened[seq_len(min(length(screened),refine_top))]
  candidates <- lapply(screened,function(first) tryCatch({
    # Tight EM convergence first, then a short exact observed-likelihood polish.
    em <- .em_fmm_ar1_4w_cells(cells,first$params,model_type,stationary,
      max_iter=maxit,ll_tol=1e-12,param_tol=1e-9)
    opt <- optim(em$eta,fn,method="BFGS",control=list(maxit=1000L,reltol=reltol))
    if(!is.finite(opt$value)) return(NULL)
    p <- resolve_fmm_ar1_4w_labels(unpack_fmm_ar1_4w(opt$par,model_type,stationary))
    eta <- pack_fmm_ar1_4w(p,model_type,stationary)
    list(opt=opt,params=p,eta=eta,loglik=-fn(eta)*sum(cells$weight))
  },error=function(e) NULL))
  candidates <- Filter(Negate(is.null),candidates)
  if(!length(candidates)) stop("all FMM starts failed during refinement")
  best <- candidates[[which.max(vapply(candidates,`[[`,numeric(1),"loglik"))]]
  # One stringent final polish after label normalization.
  opt <- optim(best$eta,fn,method="BFGS",control=list(maxit=maxit,reltol=reltol))
  params <- resolve_fmm_ar1_4w_labels(unpack_fmm_ar1_4w(opt$par,model_type,stationary))
  eta <- pack_fmm_ar1_4w(params,model_type,stationary)
  probs <- fmm_ar1_4w_cell_probabilities(params,model_type)
  loglik <- sum(cells$weight*log(pmax(probs,1e-300)))
  score <- .ar1_4w_gradient(fn,eta)
  H <- optimHess(eta,fn); H <- (H+t(H))/2
  eig <- eigen(H,symmetric=TRUE,only.values=TRUE)$values
  prob_jac <- .ar1_4w_jacobian(function(z) setNames(fmm_ar1_4w_cell_probabilities(
    unpack_fmm_ar1_4w(z,model_type,stationary),model_type)[1:15],paste0("cell",1:15)),eta)
  rank <- qr(prob_jac,tol=1e-8)$rank; K <- length(eta)
  candidate_table <- do.call(rbind,lapply(candidates,function(x) data.frame(
    loglik=x$loglik,theta0_A=x$params$theta0_A,exit_A=1-x$params$theta1_A,
    alpha_A=x$params$alpha_A,theta0_B=x$params$theta0_B,
    exit_B=1-x$params$theta1_B,alpha_B=x$params$alpha_B,
    phi=x$params$phi,pi=x$params$pi %||% NA_real_)))
  candidate_table <- candidate_table[order(candidate_table$loglik,decreasing=TRUE),]
  converged <- opt$convergence==0L && max(abs(score))<1e-6 && min(eig)>1e-9
  if(verbose) cat(sprintf("FMM-AR1-4W [%s, stationary=%s]: ll=%.3f score=%.2e rank=%d/%d minEig=%.2e starts=%d\n",
    model_type,stationary,loglik,max(abs(score)),rank,K,min(eig),nrow(candidate_table)))
  list(params=params,eta=eta,loglik=loglik,converged=converged,
       identified=rank==K,model_type=model_type,stationary=stationary,
       cell_probabilities=probs,candidates=candidate_table,
       diagnostics=list(optimizer_code=opt$convergence,max_abs_score=max(abs(score)),
       information_min_eigenvalue=min(eig),information_condition=max(eig)/min(eig),
       probability_jacobian_rank=rank,parameter_count=K,start_count=nrow(candidate_table)),
       sample=list(n=attr(cells,"n_obs"),weight_sum=attr(cells,"weight_sum")),
       estimator="direct_sixteen_cell_two_type_fmm_mle")
}

.fmm_ar1_4w_quantities <- function(eta,model_type,stationary) {
  p <- unpack_fmm_ar1_4w(eta,model_type,stationary)
  p <- resolve_fmm_ar1_4w_labels(p)
  ssA <- stationary_alpha_ar1(p$theta0_A,p$theta1_A)
  ssB <- stationary_alpha_ar1(p$theta0_B,p$theta1_B)
  agg_ss <- p$phi*ssA+(1-p$phi)*ssB
  agg_entry <- (p$phi*(1-ssA)*p$theta0_A+(1-p$phi)*(1-ssB)*p$theta0_B)/(1-agg_ss)
  agg_exit <- (p$phi*ssA*(1-p$theta1_A)+(1-p$phi)*ssB*(1-p$theta1_B))/agg_ss
  out <- c(type_A_share=p$phi,
    A_initial_employment=p$alpha_A,A_entry=p$theta0_A,A_exit=1-p$theta1_A,
    A_steady_employment=ssA,B_initial_employment=p$alpha_B,
    B_entry=p$theta0_B,B_exit=1-p$theta1_B,B_steady_employment=ssB,
    aggregate_initial_employment=p$phi*p$alpha_A+(1-p$phi)*p$alpha_B,
    aggregate_steady_employment=agg_ss,aggregate_entry=agg_entry,aggregate_exit=agg_exit)
  out <- c(out,aggregate_initial_minus_steady=
    p$phi*p$alpha_A+(1-p$phi)*p$alpha_B-agg_ss)
  if(model_type=="symmetric") out <- c(out,pi=unname(p$pi))
  out
}

analytical_se_fmm_ar1_4w <- function(df,fit,finite_sample=TRUE) {
  if(!fit$identified || !fit$converged) stop("FMM fit must be identified and converged")
  cells <- collapse_ar1_4w_cells(df); eta <- fit$eta
  type <- fit$model_type; stationary <- fit$stationary
  logcells <- function(z) setNames(log(pmax(fmm_ar1_4w_cell_probabilities(
    unpack_fmm_ar1_4w(z,type,stationary),type),1e-300)),paste0("cell",1:16))
  scores <- .ar1_4w_jacobian(logcells,eta)
  avg_nll <- function(z) -sum(cells$weight*logcells(z))/sum(cells$weight)
  bread <- optimHess(eta,avg_nll); bread <- (bread+t(bread))/2
  meat <- matrix(0,length(eta),length(eta),dimnames=list(names(eta),names(eta)))
  for(j in 1:16) meat <- meat+cells$weight_sq[j]/sum(cells$weight)^2*tcrossprod(scores[j,])
  inv <- solve(bread); vcov_eta <- inv%*%meat%*%inv
  n <- attr(cells,"n_obs");K<-length(eta);if(finite_sample&&n>K)vcov_eta<-vcov_eta*n/(n-K)
  qfun <- function(z) .fmm_ar1_4w_quantities(z,type,stationary)
  estimates <- qfun(eta); delta <- .ar1_4w_jacobian(qfun,eta)
  vcov_q <- delta%*%vcov_eta%*%t(delta)
  list(summary=data.frame(quantity=names(estimates),estimate=unname(estimates),
    se=sqrt(pmax(diag(vcov_q),0))),vcov_eta=vcov_eta,vcov_quantities=vcov_q,
    bread=bread,meat=meat,delta_jacobian=delta,
    method="individual_survey_weighted_sandwich_delta")
}
