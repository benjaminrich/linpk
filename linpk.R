pkprofile <- function(t.obs=seq(0, 24, 0.1), cl=1, vc=5, q=numeric(0), vp=numeric(0), ka=0,
    dose=list(t.dose=0, amt=1, rate=0, dur=0, ii=0, addl=0, ss=0, cmt=0, lag=0, f=1)) {

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
    if (is.null(dose$ii))     dose$ii     <- 0
    if (is.null(dose$addl))   dose$addl   <- 0
    if (is.null(dose$ss))     dose$ss     <- 0
    if (is.null(dose$cmt))    dose$cmt    <- 0
    if (is.null(dose$lag))    dose$lag    <- 0
    if (is.null(dose$f))      dose$f      <- 1

    dose$t.dose [is.na(dose$t.dose)] <- 0
    dose$amt    [is.na(dose$amt   )] <- 1
    dose$rate   [is.na(dose$rate  )] <- 0
    dose$dur    [is.na(dose$dur   )] <- 0
    dose$ii     [is.na(dose$ii    )] <- 0
    dose$addl   [is.na(dose$addl  )] <- 0
    dose$ss     [is.na(dose$ss    )] <- 0
    dose$cmt    [is.na(dose$cmt   )] <- 0
    dose$lag    [is.na(dose$lag   )] <- 0
    dose$f      [is.na(dose$f     )] <- 1

    if (!oral) {
        dose$cmt[dose$cmt == 0] <- 1
    }

    dose$amt <- dose$amt * dose$f  # Bioavailable fraction

    if (any(dose$rate > 0 & dose$dur > 0 & dose$rate != (dose$amt / dose$dur))) {
        warning("Both rate and dur specified for infusion; only rate will be considered")
    }
    i <- dose$rate == 0 & dose$dur > 0
    dose$rate[i] <- dose$amt[i] / dose$dur[i]
    i <- dose$rate > 0
    dose$dur[i] <- dose$amt[i] / dose$rate[i]

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
    }

    eigenA <- eigen(A)
    L <- eigenA$value
    V <- eigenA$vector
    qrV <- qr(V)
    qrA <- qr(A)

    y <- matrix(0, n, length(t.obs))

    for (j in seq.int(nrow(dose))) {

        tad <- t.obs - dose$t.dose[j] - dose$lag[j]

        t1 <- tad
        if (all(t1 < 0)) break

        y0 <- rep(0, n)
        if (dose$rate[j] > 0) {
            # Zero-order infusion
            b <- rep(0, n)
            b[dose$cmt[j]] <- dose$rate[j]
            ystat <- -solve(qrA, b)
            Cinf <- solve(qrV, y0 - ystat)
            i <- t1 >= 0 & t1 < dur 
            y[,i] <- y[,i] + V %*% (Cinf * exp(L %o% t1[i])) + ystat
            y0 <- drop(V %*% (Cinf * exp(L * dur)) + ystat) # At EOI
            t1 <- tad - dose$dur[j]    # Advance time to EOI
            if (all(t1 < 0)) break
        } else if (oral && dose$cmt[j] == 0) {
            y0[n] <- dose$amt[j]
        } else {
            y0[dose$cmt[j]] <- dose$amt[j] # Bolus
        }

        C <- solve(qrV, y0)
        i <- t1 >= 0
        y[,i] <- y[,i] + V %*% (C * exp(L %o% t1[i]))
    }
    y[1,]/vc
}

