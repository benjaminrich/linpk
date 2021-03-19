#' Compute Empirical Bayes Estimates of Individual Random Effects 
#'
#' Emprical Bayes (a.k.a. posthoc) estimates of individual random effects are
#' frequently used in population PK analysis to obtain the best-fitting individual-level
#' structural parameters based on a model description and some observed data.
#'
#' @param t.obs Observation times.
#' @param conc Observed concentration values.
#' @param dose Observed dosing information.
#' @param par.model A function that computes individual-level parameters from
#' population-level paramters, individual-level random effect, and covariates.
#' It should accept 3 vectors named \code{th}, \code{eta} and \code{covar}, and
#' return a named list (see Examples).
#' @param ruv.model A function describing the residual unexplained variability
#' component of the model. It should accept 4 vectors: \code{ipred}, \code{th},
#' \code{om}, and \code{sg}, and return a vector of standard deviations of the
#' same length as \code{ipred} (see Examples).
#' @param sim.model A function to simulate from the model. It should accept
#' \code{t.obs}, \code{dose}, and a list \code{param}, and return a vector of
#' the same length as \code{t.obs} (see Examples).
#' @param th A vector of population-level parameters (fixed effect).
#' @param om A variance-covariance matrix for the between-individual random effects.
#' @param sg A variance-covariance matrix for the within-individual random effects.
#' @param covar  A vector (or list) or individual-level covariates.
#' @param eta.start  A vector of starting values for the random effects.
#' @return A list with 2 elements: an vector of random effects (named
#' \code{estimate}), and associated Hessian matrix (named \code{hession})
#' @examples
#' th <- list(ka=0.7, cl=2.4, vc=5.4, vp=5.9, q=3.8, clwt=0.75, vcwt=1, vpwt=1, qwt=0.75)
#' 
#' om <- diag(c(0.55, 0.12, 0.75))
#' rownames(om) <- colnames(om) <- c("nka", "ncl", "nvc")
#' 
#' sg <- matrix(0.16)
#' rownames(sg) <- colnames(sg) <- c("errprop")
#' 
#' par.model <- function(th, eta, covar) {
#' 
#'   with(c(as.list(th), as.list(eta), as.list(covar)), {
#' 
#'     # PK Parameters
#'     ka <- ka*exp(nka)
#'     cl <- cl*exp(ncl)
#'     vc <- vc*exp(nvc)
#'     vc <- vp
#'     q  <- q
#' 
#'     # Allometric scaling
#'     cl <- cl*exp(clwt*(log(wt/70)))
#'     vc <- vc*exp(vcwt*(log(wt/70)))
#'     vp <- vp*exp(vpwt*(log(wt/70)))
#'     q  <- q*exp(qwt*(log(wt/70)))
#' 
#'     # Return the computed individual parameters as a list
#'     list(ka=ka, cl=cl, vc=vc, vp=vp, q=q)
#'   })
#' }
#' 
#' ruv.model <- function(ipred, th, om, sg) {
#'     ipred*sqrt(diag(sg)["errprop"])
#' }
#' 
#' sim.model <- function(t.obs, dose, param) {
#'     with(param, pkprofile(t.obs, ka=ka, cl=cl, vc=vc, vp=vp, q=q, dose=dose))
#' }
#' 
#' dose <- list(t.dose=c(0, 24, 48), amt=10)
#' 
#' # Typical subjects
#' plot(pkprofile(seq(0, 72, 0.1), ka=th$ka, cl=th$cl, vc=th$vc, vp=th$vp, q=th$q, dose=dose))
#' 
#' pk <- data.frame(t.obs=c(3, 5, 10, 20), dv=c(0.3, 0.3, 0.2, 0.1))
#' with(pk, points(t.obs, dv, pch=16, col="red", cex=2))
#' covar <- list(wt=61)
#' 
#' eta <- ebe(t.obs=pk$t.obs, conc=pk$dv, dose=dose, par.model=par.model, ruv.model=ruv.model, sim.model=sim.model, th=th, om=om, sg=sg, covar=covar)$estimate
#' 
#' iparam <- par.model(th=th, eta=eta, covar=covar)
#' ipred <- with(iparam, pkprofile(seq(0, 72, 0.1), ka=ka, cl=cl, vc=vc, vp=vp, q=q, dose=dose))
#' lines(ipred, col="red", lty=2)
#' 
#' @export
ebe <- function(t.obs, conc, dose, par.model, ruv.model, sim.model, th, om, sg, covar, eta.start=0*diag(om)) {

  obj <- function(eta) {
    param <- par.model(th=th, eta=setNames(eta, colnames(om)), covar=covar)
    ipred <- sim.model(t.obs=t.obs, param=param, dose=parm)
    ll <- sum(dnorm(conc, mean=ipred, sd=ruv.model(ipred=ipred, th=th, om=om, sg=sg), log=T)) + dmvnorm(eta, sigma=om, log=T)
    -2*ll
  }

  fit <- nlm(f=obj, p=eta.start, hessian=T)

  names(fit$estimate) <- colnames(om)
  dimnames(fit$hessian) <- dimnames(om)

  return(fit[c("estimate", "hessian")])
}

