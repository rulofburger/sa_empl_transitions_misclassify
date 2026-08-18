.known_fmm <- function(stationary=FALSE) {
  p<-list(theta0_A=.03,theta1_A=.97,alpha_A=.52,
          theta0_B=.30,theta1_B=.68,alpha_B=.42,phi=.72,pi=.02)
  if(stationary){p$alpha_A<-stationary_alpha_ar1(p$theta0_A,p$theta1_A)
                 p$alpha_B<-stationary_alpha_ar1(p$theta0_B,p$theta1_B)}
  p
}

test_that("FMM probabilities are valid and invariant to label switching", {
  p<-.known_fmm();q<-p
  for(stem in c("theta0","theta1","alpha")){
    q[[paste0(stem,"_A")]]<-p[[paste0(stem,"_B")]]
    q[[paste0(stem,"_B")]]<-p[[paste0(stem,"_A")]]}
  q$phi<-1-p$phi
  a<-fmm_ar1_4w_cell_probabilities(p);b<-fmm_ar1_4w_cell_probabilities(q)
  expect_equal(sum(a),1,tolerance=1e-12);expect_true(all(a>0))
  expect_equal(a,b,tolerance=1e-12)
})

test_that("free-initial symmetric FMM has full local rank", {
  eta<-pack_fmm_ar1_4w(.known_fmm(),"symmetric",FALSE)
  J<-.ar1_4w_jacobian(function(z)setNames(fmm_ar1_4w_cell_probabilities(
    unpack_fmm_ar1_4w(z,"symmetric",FALSE))[1:15],paste0("cell",1:15)),eta)
  expect_equal(qr(J,tol=1e-8)$rank,length(eta))
})

test_that("EM plus likelihood polish recovers an interior stationary mixture", {
  p<-.known_fmm(TRUE);probs<-fmm_ar1_4w_cell_probabilities(p)
  pattern<-as.data.frame(expand.grid(y1=0:1,y2=0:1,y3=0:1,y4=0:1))
  counts<-pmax(1L,round(50000*probs));idx<-rep(seq_len(16),counts)
  d<-pattern[idx,];d$weight<-1
  fit<-fit_fmm_ar1_4w(d,"symmetric",TRUE,starts=list(p),random_starts=0,
    screen_maxit=400,refine_top=1,verbose=0)
  expect_true(fit$converged);expect_true(fit$identified)
  expect_equal(fit$params$phi,p$phi,tolerance=.02)
  expect_equal(unname(fit$params$pi),p$pi,tolerance=.005)
  inf<-analytical_se_fmm_ar1_4w(d,fit)
  expect_true(all(is.finite(inf$summary$se)))
})
