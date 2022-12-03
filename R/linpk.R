#' Generate a concentration-time profile.
#'
#' This function generates concentration-time profiles from a linear
#' pharmacokinetic (PK) system, possibly with first-order absorption or
#' zero-order infusion, possibly with one or more peripheral compartments, and
#' possibly under steady-state conditions. Single or multiple doses may be
#' specified.
#'
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
#'   Only one of \code{rate} and \code{dur} should be specified unless \code{amt} is missing.}
#'   \item{\code{dur}}{Duration of zero-order infusion, or 0 to ignore (default 0).
#'   Only one of \code{rate} and \code{dur} should be specified unless \code{amt} is missing.}
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
#'   \item{\code{f}}{Bioavailable fraction (default 1).}
#' }
#' @param sc A scaling constant for the central compartment. Concentrations are
#' obtained by dividing amounts by this constant.
#' @param A A matrix of first-order rate constants between the compartments.
#' @param defdose The default dose compartment when the compartment is
#' missing or 0.
#' @param initstate A numeric vector containing values to initialize the
#' compartments.
#' @param ... Further arguments passed along to other methods.
#' @return An object of class "pkprofile", which is simply a numeric vector of
#' concentration values with some attributes attached to it.
#' This object has its own methods for \code{print}, \code{plot}, \code{lines} and \code{points}.
#' @section Warning:
#' Pay attention to the default arguments. They are there for convenience, but
#' may lead to undesired results if one is not careful.
#' @seealso
#' \itemize{
#'   \item \code{\link{halflife}}
#'   \item \code{\link{secondary}}
#'   \item \code{\link{print.pkprofile}}
#'   \item \code{\link{plot.pkprofile}}
#'   \item \code{\link{lines.pkprofile}}
#'   \item \code{\link{points.pkprofile}}
#' }
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
pkprofile <- function(...) UseMethod("pkprofile")


#' Continue an existing concentration-time profile.
#'
#' This method can be used to append to an existing PK profile, for instance to
#' simulate a PK profile with parameters that change over time. Each time the
#' parameters change, a new call to this method is used to advance the system
#' with the new parameter values.
#'
#' @param obj An object returned from a previous call to \code{\link{pkprofile}}.
#' @param t.obs A numeric vector of times at which to observe concentrations.
#' @param ... Further arguments passed along.
#' @param append Should the new profile be appended to the current samples?
#' Otherwise, only the new samples are returned.
#' @return An object of class "pkprofile".
#' @section Warning:
#' The new parameters take effect at the time when the previous profile ends.
#' If the previous profile ends before the new sampling starts, the \emph{new}
#' parameters will be used to advance the system to the start of the new
#' sampling.
#'
#' Any ongoing zero-order infusion at the end of the previous profile is
#' dropped. The remaining infusion amount will NOT be carried forward.
#' @seealso
#' \code{\link{pkprofile}}
#' @examples
#' t.obs <- seq(0, 24, 0.1)
#' amt <- 1
#' ka <- 1
#' cl <- 0.25
#' vc <- 5
#' 
#' # One-compartment model with first-order absorption
#' # First dose at time 0
#' y <- pkprofile(t.obs, cl=cl, vc=vc, ka=ka, dose=list(t.dose=0, amt=amt))
#' 
#' # Second dose at 24h with a lower clearance
#' y <- pkprofile(y, t.obs+24, cl=0.5*cl, vc=vc, ka=ka, dose=list(t.dose=24, amt=amt))
#' 
#' # Third dose at 48h with a higher clearance
#' y <- pkprofile(y, t.obs+48, cl=2*cl, vc=vc, ka=ka, dose=list(t.dose=48, amt=amt))
#' plot(y)
#' 
#' @export
pkprofile.pkprofile <- function(obj, t.obs=finaltime(obj) + seq(0, 24, 0.1), ..., append=TRUE) {
    args <- list(...)
    newargs <- args

    t.last <- finaltime(obj)
    if (length(t.last) == 0) {
        newargs$t.obs <- t.obs
        newargs$initstate <- NULL
        t.last <- 0
    } else {
        newargs$t.obs <- (t.obs - t.last)
        newargs$initstate <- finalstate(obj)
    }

    if (!is.null(args$dose)) {
        if (!is.list(args$dose)) {
            stop("dose must be given as a list or data.frame")
        }
        args$dose <- as.data.frame(args$dose)
        if (!is.null(args$dose$t.dose)) {
            newargs$dose$t.dose <- (args$dose$t.dose - t.last)
        } else {
            warning("No t.dose provided. Assuming last time point of obj.")
            newargs$dose$t.dose <- rep(0, length.out=nrow(args$dose))
        }
        if (any(newargs$dose$t.dose < 0)) {
            warning("Doses before the end of the previous profile. Expect strange behavior.")
        }
    } else {
        # Zero dose to avoid warning
        newargs$dose <- data.frame(t.dose=0, amt=0)
    }

    if (any(newargs$t.obs < 0)) {
        warning("Observations before the end of the previous profile requested. Expect strange behavior.")
    }

    obj2 <- do.call(pkprofile, newargs)
    attr(obj2, "t.obs") <- attr(obj2, "t.obs") + t.last
    attr(obj2, "dose")$t.dose <- attr(obj2, "dose")$t.dose + t.last

    if (append) {
        obj3 <- c(obj, obj2)
        attributes(obj3) <- attributes(obj2)
        attr(obj3, "t.obs") <- c(attr(obj, "t.obs"), attr(obj2, "t.obs"))
        attr(obj3, "state") <- cbind(attr(obj, "state"), attr(obj2, "state"))
        dose1 <- dose.frame(obj)
        dose2 <- dose.frame(obj2)
        dose1[,setdiff(names(dose2), names(dose1))] <- NA
        dose2[,setdiff(names(dose1), names(dose2))] <- NA
        attr(obj3, "dose") <- rbind(dose1, dose2[,names(dose1)])
        obj3
    } else {
        obj2
    }
}

#' @describeIn pkprofile Default method.
#' @export
pkprofile.default <- function(t.obs=seq(0, 24, 0.1), cl=1, vc=5, q=numeric(0), vp=numeric(0), ka=0,
    dose=list(t.dose=0, amt=1, rate=0, dur=0, ii=24, addl=0, ss=0, cmt=0, lag=0, f=1),
    sc=vc, initstate=NULL, ...) {

    if (length(list(...)) > 0) {
        warning(sprintf("Unexpected argument(s): %s",
            paste(names(list(...)), collapse=", ")))
    }

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
        defdose <- n
    } else {
        defdose <- 1
    }

    structure(
        pkprofile.matrix(A, t.obs=t.obs, dose=dose, defdose=defdose, sc=sc, initstate=initstate, call=match.call(), ...),
        pkpar=list(cl=cl, vc=vc, q=q, vp=vp, ka=ka))
}

#' @describeIn pkprofile Matrix method.
#' @export
pkprofile.matrix <- function(A, t.obs=seq(0, 24, 0.1),
    dose=list(t.dose=0, amt=1, rate=0, dur=0, ii=24, addl=0, ss=0, cmt=0, lag=0, f=1),
    defdose=1, sc=1, initstate=NULL, ...) {

    call <- if (!is.null(list(...)$call)) list(...)$call else match.call()

    if (nrow(A) < 1 | nrow(A) != ncol(A)) {
        stop("A must be a square matrix")
    }

    if (!is.list(dose)) {
        stop("dose must be given as a list or data.frame")
    }
    dose <- as.data.frame(dose)
    if (nrow(dose) == 0) {
        warning("No dose given")
        dose <- data.frame(t.dose=0, amt=0)
    }

    if (!is.null(dose$amt) && any(dose$addl < 0)) {
        warning("Negative amt. Expect strange behavior.")
    }
    if (!is.null(dose$addl) && any(dose$addl > 0) && is.null(dose$ii)) {
        stop("addl requires that ii be specified")
    }
    if (!is.null(dose$ss) && any(dose$ss > 0) && is.null(dose$ii)) {
        stop("ss requires that ii be specified")
    }
    if (!is.null(dose$rate) & !is.null(dose$dur) & !is.null(dose$amt)) { 
        if (any(!is.na(dose$rate) & !is.na(dose$dur) & !is.na(dose$amt) &
                dose$rate > 0 & dose$dur > 0 & dose$rate != (dose$amt / dose$dur))) {
            stop("amt, rate and dur are inconsistent for infusion")
        }
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

    dosecols <- c("t.dose", "amt", "rate", "dur", "ii", "addl", "ss", "cmt", "lag", "f")
    if (!all(names(dose) %in% dosecols)) {
        warning(sprintf("Unrecognized dose item(s): %s",
            paste(setdiff(names(dose), dosecols), collapse=", ")))
    }
    dose <- dose[, dosecols]

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

    dose$conc <- NA  # Keep track of the concentration at time of dose (Ctrough)

    # Zero-order infusion
    i1 <- dose$rate == 0 & dose$dur > 0
    i2 <- dose$rate > 0  & dose$dur == 0
    i3 <- dose$rate > 0  & dose$dur > 0
    dose$rate[i1] <- dose$amt[i1] / dose$dur[i1]
    dose$dur[i2] <- dose$amt[i2] / dose$rate[i2]
    dose$amt[i3] <- dose$rate[i3] * dose$dur[i3]

    # Expand addl
    if (any(dose$addl > 0)) {
        expand.addl <- do.call(rbind, lapply(seq_len(nrow(dose)), function(j) {
                    with(dose[j,],
                        data.frame(t.dose=seq(t.dose, by=ii, length.out=addl+1), j=j))
                }))
        dose <- dose[expand.addl$j,]
        dose$t.dose <- expand.addl$t.dose
        dose$addl <- NULL
    }

    dose <- dose[order(dose$t.dose),] # Time order

    n <- nrow(A)
    eigenA <- eigen(A)
    L <- eigenA$value
    if (any(L == 0)) {
        stop("The solution to the system of ODE could not be found (singular matrix)")
    }
    if (length(unique(L)) != length(L)) {
        stop("Repeated eigenvalues. Solving these types of systems is not implemented in `linpk`")
    }
    V <- eigenA$vector
    qrV <- qr(V)
    qrA <- qr(A)
    if (qrV$rank < ncol(V)) {
        stop("The solution to the system of ODE could not be found (singular matrix)")
    }

    t.aug <- c(t.obs, dose$t.dose)
    evid <- c(rep(0, length(t.obs)), rep(1, nrow(dose)))
    if (any(dose$rate > 0)) {
        dose$t.eoi <- dose$t.dose + dose$dur
        t.aug <- c(t.aug, dose$t.eoi)
        evid <- c(evid, rep(2, nrow(dose)))
    }
    y <- matrix(0, n, length(t.aug))

    if (!is.null(initstate)) {
        if (!is.numeric(initstate)) {
            stop("initstate must be a numeric vector")
        }
        if (length(initstate) != n) {
            stop("initstate must contain one value for each compartment")
        }
        C <- solve(qrV, initstate)
        y <- V %*% (C * exp(L %o% t.aug))
    }

    for (j in seq_len(nrow(dose))) {
        t.dose <- dose$t.dose[j]
        amt    <- dose$amt   [j]
        rate   <- dose$rate  [j]
        dur    <- dose$dur   [j]
        ii     <- dose$ii    [j]
        ss     <- dose$ss    [j]
        cmt    <- dose$cmt   [j]
        lag    <- dose$lag   [j]

        # Default dose compartment
        if (cmt == 0) cmt <- defdose

        tad <- t.aug - t.dose - lag

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
                    i <- tad < 0 & tad >= -lag
                    t1[i] <- t1[i] %% ii
                }
            } else {
                y0 <- rep(0, n)
            }
            Cinf <- solve(qrV, y0 - ystat)
            if (any(t1==0)) {
                y[,t1==0] <- y[,t1==0] + y0
            }
            i <- t1 > 0 & t1 < dur 
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
                    i <- tad < 0 & tad >= -lag
                    t1[i] <- t1[i] %% ii
                }
            } else {
                y0 <- delta
            }
        }

        C <- solve(qrV, y0)
        if (any(t1==0)) {
            y[,t1==0] <- y[,t1==0] + y0
        }
        i <- t1 > 0
        y[,i] <- y[,i] + V %*% (C * exp(L %o% t1[i]))
    }

    state.all <- Re(y)
    conc.all <- state.all[1,]/sc[1]
    state <- state.all[, evid==0, drop=FALSE]
    conc <- conc.all[evid==0]

    # Keep track of the concentration at time of dose (Ctrough), and EOI
    dose$conc <- conc.all[evid==1]
    if (any(dose$rate > 0)) {
        dose$conc.eoi <- conc.all[evid==2]
        dose$conc.eoi[dose$rate == 0] <- NA
        dose$t.eoi[dose$rate == 0] <- NA
    }

    structure(conc,
        class      = c("pkprofile", class(conc)),
        call       = call,
        t.obs      = t.obs,
        state      = state,
        dose       = dose,
        defdose    = defdose,
        sc         = sc,
        A          = A,
        L          = L,
        V          = V)
}

#' Get the final state or time of a PK profile.
#' @param x A object of class \code{\link{pkprofile}}.
#' @return A \code{numeric} vector containing the state of each compartment at
#' the final observation time (\code{finalstate}), or the final observation
#' time itself (\code{finaltime}).
#' @seealso
#' \itemize{
#'   \item \code{\link{pkprofile}} for generating a PK profile.
#'   \item \code{\link{pkprofile.pkprofile}} for appending to an existing PK profile.
#' }
#' @examples
#' # Administer a dose at time 0 and a second dose using the final state
#' # from the first dose (at 12h) as the initial state for the second dose.
#' t.obs <- seq(0, 12, 0.1)
#' y <- pkprofile(t.obs, cl=0.25, vc=5, ka=1, dose=list(t.dose=0, amt=1))
#' finalstate(y)
#' y2 <- pkprofile(t.obs, cl=0.25, vc=5, ka=1, dose=list(t.dose=0, amt=1), initstate=finalstate(y))
#' plot(y, xlim=c(0, 24), ylim=c(0, max(y2)), col="blue")  # First dose
#' lines(t.obs+12, y2, col="red")                          # Second dose
#'
#' # Add a vertical line to show where the first profile ends.
#' abline(v=finaltime(y), col="gray75", lty=2)
#' @export
finalstate <- function(x) {
    state <- attr(x, "state")
    finalstate <- state[, ncol(state), drop=FALSE]
    as.numeric(finalstate)
}

#' @rdname finalstate
#' @export
finaltime <- function(x) {
    t.obs <- attr(x, "t.obs")
    t.obs[length(t.obs)]
}

#' Get the doses from a PK profile.
#' @param x A object of class \code{\link{pkprofile}}.
#' @return A \code{data.frame} containing the realized doses, one per row. The
#' \code{data.frame} has all the columns described in \code{\link{pkprofile}},
#' except \code{addl}, since all additional doses have been expanded to
#' individual rows. It also has a \code{conc} column with the simulated
#' concentration at the time of the dose.
#' @seealso
#' \code{\link{pkprofile}}
#' @examples
#' t.obs <- seq(0, 6*24, 0.5)
#' y <- pkprofile(t.obs, cl=0.5, vc=11, ka=1.3,
#'     dose=list(t.dose=c(0, 24*2 + 14), amt=c(100, 50), addl=c(4, 0), ii=24))
#' dose.frame(y)
#' @export
dose.frame <- function(x) {
    attr(x, "dose")
}

#' Derive secondary PK parameters.
#' @param x A object of class \code{\link{pkprofile}}.
#' @param From A vector of interval start times. The defaults is the times of the doses.
#' @param To A vector of interval end times. The defaults is the time of the next dose,
#' or last observation time.
#' @param include.dose.times Should dose times (and end of infusion times) be
#' considered in addition to the simulation times?
#' @return A \code{data.frame} with one row for each time interval and with the
#' following columns:
#' \describe{
#'   \item{\code{From}}{The time of the start of the interval. Can differ from
#'   the specified start time because it always corresponds to an actual data
#'   point.}
#'   \item{\code{To}}{The time of the end of the interval. Can differ from the
#'   specified end time because it always corresponds to an actual data point.}
#'   \item{\code{N}}{The number of distinct data points in the interval used to
#'   derive \code{AUC}, \code{Cmax}, etc.}
#'   \item{\code{Ctrough}}{Concentration at the time of dose (i.e. just prior
#'   to the dose). Only present if the start of the interval corresponds to a
#'   dose time.}
#'   \item{\code{Cmin}}{Minimum concentration over the interval.}
#'   \item{\code{Tmin}}{Time of the minimum concentration over the interval.}
#'   \item{\code{Cmax}}{Maximum concentration over the interval.}
#'   \item{\code{Tmax}}{Time of the maximum concentration over the interval.}
#'   \item{\code{Cave}}{Average concentration over the interval (calculated by
#'   the trapezoid rule).}
#'   \item{\code{AUC}}{Area under the concentration-time curve over the
#'   interval (calculated by the trapezoid rule).}}
#' @examples
#' t.obs <- seq(0, 24*4, 0.1)
#' y <- pkprofile(t.obs, cl=0.25, vc=5, ka=1, dose=list(t.dose=0, amt=1, addl=6, ii=12))
#' secondary(y)
#' secondary(y, 0, 48)
#' secondary(y, 0, Inf)
#' sum(secondary(y)$AUC)  # Same as above
#' plot(y)
#' with(secondary(y), points(Tmax, Cmax, pch=19, col="blue"))
#' with(secondary(y), points(Tmin, Cmin, pch=19, col="red"))
#' with(secondary(y), points(From, Ctrough, pch=19, col="green"))
#' with(secondary(y), points(From + 6, Cave, pch=19, col="purple", cex=2))
#' 
#' @export
secondary <- function(x, From=NULL, To=NULL, include.dose.times=T) {
    conc <- as.numeric(x)
    time <- attr(x, "t.obs")
    dose <- attr(x, "dose")

    if (include.dose.times) {
        time <- c(time, dose$t.dose)
        conc <- c(conc, dose$conc)
        if (!is.null(dose$t.eoi)) {
            i <- !is.na(dose$t.eoi)
            time <- c(time, dose$t.eoi[i])
            conc <- c(conc, dose$conc.eoi[i])
        }
    }
    o <- order(time)
    time <- time[o]
    conc <- conc[o]
    duplicate <- c(FALSE, diff(time) == 0 & diff(conc) == 0)
    time <- time[!duplicate]
    conc <- conc[!duplicate]

    remove.equal <- (is.null(From) && is.null(To))
    if (is.null(From)) {
        From <- dose$t.dose
    }
    if (is.null(To)) {
        To <- c(dose$t.dose[-1], max(time))
    }

    d <- data.frame(From=From, To=To)
    if (remove.equal) {
        d <- subset(d, From != To)
    }
    k <- nrow(d)

    d$N       <- numeric(k)
    d$Ctrough <- numeric(k)
    d$Cmin    <- numeric(k)
    d$Cmax    <- numeric(k)
    d$Cave    <- numeric(k)
    d$AUC     <- numeric(k)
    d$Tmin    <- numeric(k)
    d$Tmax    <- numeric(k)


    for (j in seq_len(k)) {
        d$From[j] <- min(time[time >= d$From[j]])
        d$To[j] <- max(time[time <= d$To[j]])
        i <- time >= d$From[j] & time <= d$To[j]
        d$N[j] <- sum(i)
        d$Ctrough[j] <- dose$conc[dose$t.dose == d$From[j]][1]
        d$Cmin[j] <- min(conc[i])
        d$Cmax[j] <- max(conc[i])
        d$Tmax[j] <- time[i][which.max(conc[i])]
        d$Tmin[j] <- time[i][which.min(conc[i])]
        d$AUC[j] <- AUC.by.trapezoid(time[i], conc[i])
        d$Cave[j] <- if (diff(range(time[i])) == 0) {
            d$Cmin[j]
        } else {
            d$AUC[j]/diff(range(time[i]))
        }
    }
    d
}

AUC.by.trapezoid <- function(x, y) {
    i <- order(x)
    x <- x[i]
    y <- y[i]
    sum(0.5 * (y[-1] + y[-length(y)]) * diff(x))
}

Tmax.oral1cpt <- function(cl, vc, ka, ss=FALSE, ii) {
    ke <- cl/vc
    if (ss) {
        if (missing(ii) || is.null(ii)) {
            stop("ss requires that ii be specified")
        }
        log((ka * (1 - exp(-ke*ii)))/(ke*(1 - exp(-ka*ii))))/(ka - ke)
    } else {
        log(ka/ke)/(ka - ke)
    }
}

#' Half-lives of a linear PK system.
#' @param x A object of class \code{\link{pkprofile}}.
#' @return A \code{numeric} vector containing the half-lives for the different
#' phases of the system. The number of phases generally equal the number of
#' compartments, plus one for the absorption phase if the system has first
#' order absorption (i.e. if \code{ka} is specified). The values are returned
#' sorted in ascending order, so the first corresponds to the alpha phase,
#' the second beta, the third gamma, and so on. The absorption half-life, if
#' present, comes last (it can also be identified by comparing it to the value
#' of \code{log(2)/ka}).
#' @examples
#' y <- pkprofile(0, cl=0.25, vc=5, ka=1.1)
#' halflife(y)
#' log(2)/1.1
#' 
#' y <- pkprofile(0, cl=0.25, vc=5, ka=0.01)  # Flip-flop kinetics
#' halflife(y)
#' log(2)/0.01
#' 
#' # Three-compartment model
#' y <- pkprofile(0, cl=2, vc=10, q=c(0.5, 0.3), vp=c(30, 40))
#' halflife(y)
#' 
#' # The terminal half-life can be used to obtain the terminal slope of the
#' # concentration-time curve on the semi-log scale:
#' t.obs <- seq(0, 36, 0.1)
#' y <- pkprofile(t.obs, cl=0.25, vc=5, ka=1, dose=list(t.dose=0, amt=1))
#' plot(log2(y))
#' abline(-2.247927, -1/halflife(y)[1], col=adjustcolor("blue", 0.2), lwd=12)
#' 
#' @export
halflife <- function(x) {
    L <- attr(x, "L")
    HL <- log(2)/(-Re(L))
    HL[HL < 0] <- Inf
    HL <- sort(HL)
    names(HL) <- paste0("HL.", seq_along(HL))
    ka <- as.list(attr(x, "pkpar"))$ka
    if (!is.null(ka)) {
        if (ka > 0) {
            absorp <- which.min(abs((log(2)/ka) - HL))
            HL <- HL[c(seq_along(HL)[-absorp], absorp)] # Put it last
            names(HL) <- paste0("HL.", seq_along(HL))
            names(HL)[length(HL)] <- "HL.a"
        }
    }
    HL
}

#' Coerce a \code{pkprofile} to a \code{data.frame}
#' @param x An object of class \code{pkprofile}.
#' @param ... Further arguments passed along.
#' @param col.names Character vector of length 2 giving the names for the time and concentration columns.
#' @param .state Include the complete state along with \code{time} and \code{conc}?
#' @return A \code{data.frame} with columns \code{time} and \code{conc} (or the
#' names specified in \code{col.names}). If \code{.state == TRUE}, then the
#' complete state is appended (as a matrix column).
#' @export
as.data.frame.pkprofile <- function(x, ..., col.names=c("time", "conc"), .state=FALSE) {
    df <- data.frame(time=attr(x, "t.obs"), conc=as.numeric(x), ...)
    if (!is.character(col.names) | length(col.names) != 2) {
        warning("There should be 2 column names specified, for time and concentration respectively")
    }
    names(df)[1:2] <- col.names
    if (.state) {
        df$state=t(attr(x, "state"))
    }
    df
}

#' Printing and plotting methods for class \code{pkprofile}.
#' @param x An object of class \code{pkprofile}.
#' @param y Any other object. Specifying \code{y} causes the default method to
#' be called instead (effectively overriding the class-specific behaviour).
#' @keywords internal
#' @name pkprofile-methods
NULL

#' @importFrom utils head tail
#' @rdname pkprofile-methods
#' @export
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
    invisible(x)
}

#' @importFrom graphics plot.default
#' @rdname pkprofile-methods
#' @export
plot.pkprofile <- function(x, y, ...) {
    if (!missing(y)) {
        NextMethod()
    } else {
        args <- list()
        args$x <- attr(x, "t.obs")
        args$y <- as.numeric(x)
        args <- c(args, list(...))
        if (is.null(args$xlab)) {
            args$xlab <- "Time"
        }
        if (is.null(args$ylab)) {
            args$ylab <- "Concentration"
        }
        if (is.null(args$type)) {
            args$type <- "l"
        }
        do.call(plot.default, args)
    }
}

#' @importFrom graphics lines.default
#' @rdname pkprofile-methods
#' @export
lines.pkprofile <- function(x, y, ...) {
    if (!missing(y)) {
        NextMethod()
    } else {
        args <- list()
        args$x <- attr(x, "t.obs")
        args$y <- as.numeric(x)
        args <- c(args, list(...))
        do.call(lines.default, args)
    }
}

#' @importFrom graphics points.default
#' @rdname pkprofile-methods
#' @export
points.pkprofile <- function(x, y, ...) {
    if (!missing(y)) {
        NextMethod()
    } else {
        args <- list()
        args$x <- attr(x, "t.obs")
        args$y <- as.numeric(x)
        args <- c(args, list(...))
        do.call(points.default, args)
    }
}

#' Construct a symmetric matrix from its lower triangle.
#' @param LT A numeric vector giving the elements of the lower triangle of the
#' matrix by row (see examples).
#' @param .names,.colnames,.rownames Optionally, specify row and column names of the resulting matrix.
#' @return A symmetric matrix.
#' @examples
#' LTmat(1:6)
#' @export
LTmat <- function(LT, .names=attr(LT, ".names"), .colnames=.names, .rownames=.names) {
    x <- length(LT)
    p <- (sqrt(8*x + 1) - 1)/2
    m <- matrix(0, p, p)
    m[upper.tri(m, diag=T)] <- LT
    m <- m + t(m) - diag(diag(m))
    colnames(m) <- .colnames; rownames(m) <- .rownames
    m
}

#' Convert from standard deviation and correlation matrix to covariance matrix.
#' @param cor A correlation matrix. If \code{sd} is missing, the diagonal
#' entries are taken to be the standard deviations, otherwise they are ignored.
#' @param sd A vector of standard deviations (optional).
#' @return A covariance matrix.
#' @examples
#' cor2cov(matrix(c(1, 0.5, 0.5, 1), 2, 2), 0.1)
#' @export
cor2cov <- function(cor, sd) {
    if (missing(sd)) {
        sd <- diag(cor)
    }
    diag(cor) <- 1
    n <- nrow(cor)
    diag(sd, n) %*% cor %*% diag(sd, n)
}

#' Construct a block-diagonal matrix.
#' @param ... Any number of square matrices making up the diagonal blocks of
#' the matrix.
#' @return A block-diagonal matrix.
#' @examples
#' blockdiag(matrix(1, 2, 2), 2, matrix(3, 4, 4))
#' @export
blockdiag <- function(...) {
    b <- list(...)
    b <- lapply(b, function(x) {
        if (is.numeric(x) && length(x) == 1) {
            x <- as.matrix(x)
        }
        if (!is.matrix(x) || nrow(x) != ncol(x)) {
            stop("All arguments must be square matrices")
        }
        x
    })
    n <- sapply(b, nrow)
    cn <- c(0, cumsum(n))
    ntot <- sum(n)
    m <- matrix(0, ntot, ntot)
    for (j in seq_along(n)) {
        i <- (1:n[j]) + cn[j]
        m[i,i] <- b[[j]]
    }
    rnam <- unlist(lapply(b, function(x) dimnames(x)[[1]]))
    cnam <- unlist(lapply(b, function(x) dimnames(x)[[2]]))
    if (length(rnam) == nrow(m)) { rownames(m) <- rnam }
    if (length(cnam) == ncol(m)) { colnames(m) <- cnam }
    m
}

#' Generate individual random effects from a multivariate normal distribution.
#' @param n The number of individuals.
#' @param omegaLT A numeric vector giving the elements of the lower triangle
#' of the covariance matrix by row.
#' @param omega The covariance matrix.
#' @param eta.names A character vector of names for each random effect
#' (defaults to the column names of \code{omega}, or if \code{NULL} then to
#' ETA1, ETA2, ...).
#' @return An \eqn{n \times p} matrix, where each row contains the vector of random
#' effects for one individual (\eqn{p} is the size of the covariance matrix).
#' @seealso
#' \code{\link{LTmat}}
#' \code{\link{blockdiag}}
#' @examples
#' omegaLT <- c(0.123, 0.045, 0.678)
#' generateETA(10, omegaLT)
#' @export
generateETA <- function(n, omegaLT, omega=LTmat(omegaLT), eta.names=colnames(omega)) {
    if (is.null(eta.names)) {
        eta.names <- sprintf("ETA%d", 1:nrow(omega))
    }
    x <- mvtnorm::rmvnorm(n, rep(0, nrow(omega)), omega)
    colnames(x) <- eta.names
    x
}

#' Runs the interactive shiny app.
#' @param ... Arguments passed to \code{shiny::runApp()}.
#' @return Called for its side effects.
#' @section Note:
#' The app requires the following packages:
#' \itemize{
#'   \item `shiny`
#'   \item `shinyjs`
#'   \item `shinyAce`
#'   \item `dygraphs`
#' }
#' Make they are installed or the app won't work.
#' @examples
#' \dontrun{
#' linpkApp()
#' }
#' @export
linpkApp <- function(...) {
    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Please install `shiny` before running the app.", call.=F)
    }
    appDir <- system.file("demo-app", package="linpk")
    if (appDir == "") {
        stop("It seems that the app is not installed (or could not be found). Try re-installing `linpk`.", call.=F)
    }
    shiny::runApp(appDir, ...)
    invisible(NULL)
}

