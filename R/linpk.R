#' Generate a concentration-time profile.
#'
#' This function generates concentration-time profiles from a linear
#' pharmacokinetic (PK) system, possibly with first-order absoption or
#' zero-order infusion, possibly with one or more peripheral compartments, and
#' possibly under steady-state conditions. Single or multiple doses may be
#' specified.

#' @param t.obs A numeric vector of times at which to observe concentrations.
#' @param cl Central clearance parameter.
#' @param vc Central volume parameter.
#' @param q Inter-compartmental clearance. Can be a vector for more than one
#' peripheral compartment, or empty for none. Must match \code{vp} in length.
#' @param vp Peripheral volume. Can be a vector for more than one
#' peripheral compartment, or empty for none. Must match \code{q} in length.
#' @param ka First-order absorption rate parameter. Set to 0 to indicate
#' that there is no first-order absorption (i.e. bolus or infusion).
#' @param dose A \code{list} or \code{data.frame} containing dose information.
#' May contain the following elements:
#' \describe{
#'   \item{\code{t.dose}}{Dose time (default 0).}
#'   \item{\code{amt}}{Dose amount (default 1).}
#'   \item{\code{rate}}{Rate of zero-order infusion, or 0 to ignore (default 0).
#'   Only one of \code{rate} and \code{dur} should be specified.}
#'   \item{\code{dur}}{Duration of zero-order infusion, or 0 to ignore (default 0).
#'   Only one of \code{rate} and \code{dur} should be specified.}
#'   \item{\code{ii}}{Interdose interval (default 24). Only used if addl or ss are used.}
#'   \item{\code{addl}}{Number of \emph{additional} doses (default 0). The
#'   total number of doses given is \code{addl + 1}.}
#'   \item{\code{ss}}{Indicates that a dose is given under steady-state
#'   conditions (default 0 or FALSE; converted to \code{logical} internally).}
#'   \item{\code{cmt}}{The number of the compartment into which the dose is
#'   administered. The default value is 0, which indicates the depot
#'   compartment for first-order absorption (i.e. \code{ka > 0}), and central
#'   compartment otherwise.}
#'   \item{\code{lag}}{Time lag (default 0).}
#'   \item{\code{f}}{Bioavailable fraction (default 0).}
#' }
#' @return An object of class "pkprofile", which simply a numeric vector of
#' concentration values with some attributes attached to it. These include:
#' \describe{
#'   \item{\code{t.obs}}{Numeric vector of concentration times.}
#'   \item{\code{t.dose}}{Numeric vecotr of dose time.}
#'   \item{\code{secondary}}{A list of derived secondary PK parameters. This includes:
#'     \describe{
#'       \item{\code{HLterm}}{Terminal half-life.}
#'       \item{\code{Ctrough}}{Concentration value at the time the dose was
#'       given (assuming this time is one of those in \code{t.obs}, otherwise
#'       at the most recent previous observation time}.
#'       \item{\code{Cmin}}{The minimum concentration observed in the time
#'       interval between two consecutive doses.}
#'       \item{\code{Cmax}}{The minimum concentration observed in the time
#'       interval between two consecutive doses.}
#'       \item{\code{AUC}}{The area under the concentration-time curve in the time
#'       interval between two consecutive doses, calcuated by the trapezoid rule.}
#'     \code{Ctrough}, \code{Cmin}, \code{Cmax} and \code{AUC} are vectors of
#'     length \code{j - 1} where \code{j} is the number of doses given. It is
#'     recommended that \code{t.obs} be a superset of \code{t.dose}.}
#'   }
#' }
#' This object has its own methods for \code{print}, \code{plot}, \code{lines} and \code{points}.
#' @examples
#' # Default values, a bolus injection
#' y <- pkprofile()
#' plot(y)
#' 
#' t.obs <- seq(0, 24, 0.1)
#' dur <- 1
#' amt <- 1
#' ka <- 1
#' cl <- 0.25
#' vc <- 5
#' q <- 2.5
#' vp <- 10
#' 
#' # One-compartment model with first-order absorption, single dose
#' y <- pkprofile(t.obs, cl=cl, vc=vc, ka=ka, dose=list(amt=amt))
#' plot(y)
#' 
#' # Two-compartment model with first-order absorption, single dose
#' y <- pkprofile(t.obs, cl=cl, vc=vc, vp=vp, q=q, ka=ka, dose=list(amt=amt))
#' plot(y)
#' 
#' # One-compartment model with zero-order infusion, single dose
#' y <- pkprofile(t.obs, cl=cl, vc=vc, dose=list(dur=dur, amt=amt))
#' plot(y)
#' 
#' # Two-compartment model with zero-order infusion, single dose
#' y <- pkprofile(t.obs, cl=cl, vc=vc, vp=vp, q=q, dose=list(dur=dur, amt=amt))
#' plot(y)
#' 
#' # Two-compartment model with bolus injection, single dose
#' y <- pkprofile(t.obs, cl=cl, vc=vc, vp=vp, q=q, dose=list(amt=amt))
#' plot(y)
#' 
#' # Two-compartment model with bolus injection into the peripheral compartment, single dose
#' y <- pkprofile(t.obs, cl=cl, vc=vc, vp=vp, q=q, dose=list(amt=amt, cmt=2))
#' plot(y)
#' 
#' # Two-compartment model with zero-order infusion into the peripheral compartment, single dose
#' y <- pkprofile(t.obs, cl=cl, vc=vc, vp=vp, q=q, dose=list(amt=amt, cmt=2, dur=dur))
#' plot(y)
#' 
#' t.obs <- seq(0, 24*6, 1)
#'
#' # One-compartment model with first-order absorption, multiple doses
#' y <- pkprofile(t.obs, cl=cl, vc=vc, ka=ka, dose=list(t.dose=seq(0, 24*5, 12), amt=amt))
#' plot(y)
#' 
#' # One-compartment model with first-order absorption, multiple doses specified by addl and ii
#' y <- pkprofile(t.obs, cl=cl, vc=vc, ka=ka, dose=list(t.dose=0, amt=amt, addl=9, ii=12))
#' plot(y, type="b")
#' points(y, col="blue")
#' 
#' # One-compartment model with first-order absorption, multiple doses under steady-state conditions
#' yss <- pkprofile(t.obs, cl=cl, vc=vc, ka=ka, dose=list(t.dose=0, amt=amt, addl=9, ii=12, ss=1))
#' lines(yss, col="red")
#' points(yss, col="green")
#' 
#' # One-compartment model with zero-order infusion, multiple doses specified by addl and ii
#' y <- pkprofile(t.obs, cl=cl, vc=vc, dose=list(dur=dur, amt=amt, addl=9, ii=12))
#' plot(y, log="y")
#' 
#' # One-compartment model with zero-order infusion, multiple doses  under steady-state conditions
#' yss <- pkprofile(t.obs, cl=cl, vc=vc, dose=list(dur=dur, amt=amt, addl=9, ii=12, ss=1))
#' lines(yss, col="red")
#' 
#' @export
#' @importFrom utils head tail
pkprofile <- function(t.obs=seq(0, 24, 0.1), cl=1, vc=5, q=numeric(0), vp=numeric(0), ka=0,
    dose=list(t.dose=0, amt=1, rate=0, dur=0, ii=24, addl=0, ss=0, cmt=0, lag=0, f=1)) {

    # Check arguments
    if (!(is.numeric(cl) && length(cl) == 1 && !is.na(cl) && cl > 0)) {
        stop("cl must be a single positive numeric value")
    }
    if (!(is.numeric(vc) && length(vc) == 1 && !is.na(vc) && vc > 0)) {
        stop("vc must be a single positive numeric value")
    }
    if (!(is.numeric(q) && is.numeric(vp) && length(q) == length(vp) &&
            all(!is.na(q) & q > 0) && all(!is.na(vp) & vp > 0))) {
        stop("q and vp must be positive numeric vectors of the same length")
    }
    if (!(is.numeric(ka) && length(ka) == 1 && !is.na(ka) && ka >= 0)) {
        stop("ka must be a single non-negative numeric value (set to 0 for no first-order absorption)")
    }
    oral <- (ka > 0)

    if (!is.list(dose)) {
        stop("dose must be given as a list or data.frame")
    }
    dose <- as.data.frame(dose)
    if (nrow(dose) == 0) {
        stop("No dose given")
    }

    # Defaults
    if (is.null(dose$t.dose)) dose$t.dose <- 0
    if (is.null(dose$amt))    dose$amt    <- 1
    if (is.null(dose$rate))   dose$rate   <- 0
    if (is.null(dose$dur))    dose$dur    <- 0
    if (is.null(dose$ii))     dose$ii     <- 24
    if (is.null(dose$addl))   dose$addl   <- 0
    if (is.null(dose$ss))     dose$ss     <- 0
    if (is.null(dose$cmt))    dose$cmt    <- 0
    if (is.null(dose$lag))    dose$lag    <- 0
    if (is.null(dose$f))      dose$f      <- 1

    dose$t.dose [is.na(dose$t.dose)] <- 0
    dose$amt    [is.na(dose$amt   )] <- 1
    dose$rate   [is.na(dose$rate  )] <- 0
    dose$dur    [is.na(dose$dur   )] <- 0
    dose$ii     [is.na(dose$ii    )] <- 24
    dose$addl   [is.na(dose$addl  )] <- 0
    dose$ss     [is.na(dose$ss    )] <- 0
    dose$cmt    [is.na(dose$cmt   )] <- 0
    dose$lag    [is.na(dose$lag   )] <- 0
    dose$f      [is.na(dose$f     )] <- 1

    dose$amt <- dose$amt * dose$f  # Bioavailable fraction

    dose$ss <- as.logical(dose$ss) # Steady state

    # Zero-order infusion
    if (any(dose$rate > 0 & dose$dur > 0 & dose$rate != (dose$amt / dose$dur))) {
        warning("Both rate and dur specified for infusion; only rate will be considered")
    }
    i <- dose$rate == 0 & dose$dur > 0
    dose$rate[i] <- dose$amt[i] / dose$dur[i]
    i <- dose$rate > 0
    dose$dur[i] <- dose$amt[i] / dose$rate[i]

    # Expand addl
    if (any(dose$addl > 0)) {
        expand.addl <- do.call(rbind, lapply(seq.int(nrow(dose)), function(j) {
                    with(dose[j,],
                        data.frame(t.dose=seq(t.dose, by=ii, length.out=addl+1), j=j))
                }))
        dose <- dose[expand.addl$j,]
        dose$t.dose <- expand.addl$t.dose
        dose$addl <- NULL
    }

    dose <- dose[order(dose$t.dose),] # Time order

    ncomp <- 1 + length(q)
    n <- ncomp + oral
    A <- matrix(0, n, n)
    ke <- cl/vc
    A[1,1] <- -ke
    for (j in seq_along(q)) {
        # Peripheral compartment(s)
        k1. <- q[j]/vc
        k.1 <- q[j]/vp[j]
        A[1,1] <- A[1,1] - k1.
        A[1,j+1] <- k.1
        A[j+1,1] <- k1.
        A[j+1,j+1] <- -k.1
    }
    if (oral) {
        A[1,n] <- ka
        A[n,n] <- -ka
        dose$cmt[dose$cmt == 0] <- n
    } else {
        dose$cmt[dose$cmt == 0] <- 1
    }


    eigenA <- eigen(A)
    L <- eigenA$value
    V <- eigenA$vector
    qrV <- qr(V)
    qrA <- qr(A)

    y <- matrix(0, n, length(t.obs))

    for (j in seq.int(nrow(dose))) {
        t.dose <- dose$t.dose[j]
        amt    <- dose$amt   [j]
        rate   <- dose$rate  [j]
        dur    <- dose$dur   [j]
        ii     <- dose$ii    [j]
        ss     <- dose$ss    [j]
        cmt    <- dose$cmt   [j]
        lag    <- dose$lag   [j]

        tad <- t.obs - t.dose - lag

        t1 <- tad
        if (all(t1 < 0)) break

        if (rate > 0) {
            # Zero-order infusion
            b <- rep(0, n)
            b[cmt] <- rate
            ystat <- -solve(qrA, b)
            if (ss) {
                y[,t1 >= 0] <- 0  # Reset
                Q1 <- exp(L * ii) / (1 - exp(L * ii))
                Q2 <- exp(L * dur) / (1 - exp(L * dur))
                y0 <- drop(V %*% (solve(qrV, ystat) * Q1 / Q2))
                if (lag > 0) {
                    i <- tad >= 0 & tad < lag
                    t1[i] <- t1[i] %% ii
                }
            } else {
                y0 <- rep(0, n)
            }
            Cinf <- solve(qrV, y0 - ystat)
            i <- t1 >= 0 & t1 < dur 
            y[,i] <- y[,i] + V %*% (Cinf * exp(L %o% t1[i])) + ystat
            y0 <- drop(V %*% (Cinf * exp(L * dur)) + ystat) # At EOI
            t1 <- tad - dur    # Advance time to EOI
            if (all(t1 < 0)) break
        } else {
            delta <- rep(0, n)
            delta[cmt] <- amt
            if (ss) {
                y[,t1 >= 0] <- 0  # Reset
                y0 <- drop(V %*% (solve(qrV, delta) / (1 - exp(L * ii))))
                if (lag > 0) {
                    i <- tad >= 0 & tad < lag
                    t1[i] <- t1[i] %% ii
                }
            } else {
                y0 <- delta
            }
        }

        C <- solve(qrV, y0)
        i <- t1 >= 0
        y[,i] <- y[,i] + V %*% (C * exp(L %o% t1[i]))
    }

    conc <- y[1,]/vc

    # Derive secondary parameters
    Ctrough <- numeric(nrow(dose))
    Cmin <- numeric(nrow(dose))
    Cmax <- numeric(nrow(dose))
    Tmax <- numeric(nrow(dose))
    AUC <- numeric(nrow(dose))
    for (j in seq.int(nrow(dose))) {
        i <- t.obs <= dose$t.dose[j]
        Ctrough[j] <- tail(conc[i], 1)
        i <- t.obs >= dose$t.dose[j] & t.obs < ifelse(j < nrow(dose), dose$t.dose[j+1], Inf)
        Cmin[j] <- min(conc[i])
        Cmax[j] <- max(conc[i])
        Tmax[j] <- t.obs[i][which.max(conc[i])]
        i <- t.obs >= dose$t.dose[j]
        AUC[j] <- sum(0.5 * (conc[i][-1] + conc[i][-length(conc[i])]) * diff(t.obs[i]))
    }
    AUC <- c(rev(diff(rev(AUC))), tail(AUC, 1))

    structure(conc,
        class = "pkprofile",
        t.obs = t.obs,
        t.dose = dose$t.dose,
        secondary = list(
            HLterm = log(2)/min(-eigen(A[1:ncomp,1:ncomp])$values),
            Ctrough = Ctrough,
            Cmin = Cmin,
            Cmax = Cmax,
            Tmax = Tmax,
            AUC = AUC))
}

#' print method for class pkprofile
#'
#' @keywords internal
#' @importFrom utils head tail
print.pkprofile <- function(x, ...) {
    t.obs <- attr(x, "t.obs")
    tt <- head(t.obs, 5)
    if (length(t.obs) > 6) {
        tt <- c(tt, "...")
    }
    tt <- c(tt, tail(t.obs, 1))
    tt <- paste0(tt, collapse=", ")
    cat("PK concentration-time profile at times: ", tt, "\n")
    print(as.numeric(x))
}

#' plot method for class pkprofile
#'
#' @keywords internal
#' @importFrom graphics plot.default
plot.pkprofile <- function(x, y, ...) {
    if (!missing(y)) {
        NextMethod()
    } else {
        args <- list(...)
        if (is.null(args$xlab)) {
            args$xlab <- "Time"
        }
        if (is.null(args$ylab)) {
            args$ylab <- "Concentration"
        }
        if (is.null(args$type)) {
            args$type <- "l"
        }
        args$x <- attr(x, "t.obs")
        args$y <- as.numeric(x)
        do.call(plot.default, args)
    }
}

#' lines method for class pkprofile
#'
#' @keywords internal
#' @importFrom graphics lines.default
lines.pkprofile <- function(x, y, ...) {
    if (!missing(y)) {
        NextMethod()
    } else {
        args <- list(...)
        args$x <- attr(x, "t.obs")
        args$y <- as.numeric(x)
        do.call(lines.default, args)
    }
}

#' points method for class pkprofile
#'
#' @keywords internal
#' @importFrom graphics points.default
points.pkprofile <- function(x, y, ...) {
    if (!missing(y)) {
        NextMethod()
    } else {
        args <- list(...)
        args$x <- attr(x, "t.obs")
        args$y <- as.numeric(x)
        do.call(points.default, args)
    }
}
