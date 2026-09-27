# Interim monitoring of a solved design (bounds.R).
#
#   boundaries at the information actually observed (Lan & DeMets 1983, 1989);
#   repeated confidence intervals (Jennison & Turnbull 1989);
#   conditional power (Lan, Simon & Halperin 1982) and predictive power
#   (Spiegelhalter, Freedman & Blackburn 1986), both exact over the remaining looks;
#   inference after stopping by the stagewise ordering (Fairbanks & Madsen 1982;
#   Tsiatis, Rosner & Mehta 1984): p-value, median-unbiased estimate, interval.
#
# Looks are described on the design's standardised scale: spending time t (fraction
# of the planned maximum information) and Z oriented so that positive favours H1.
# Estimates, standard errors and information are on the natural scale of the effect.

# Boundaries at looks tt (the last is the final analysis), for a design d.
gs_bounds_at <- function(d, tt) {
    spec <- d$spec
    K <- length(tt)
    if (spec$eff$type == "triangular") {
        tri <- d$triangle
        return(tri_bounds(tt, tri[["Vmax"]], tri[["a"]], tri[["c"]], spec$sided, as.logical(tri[["correction"]])))
    }
    fut <- spec$fut0
    if (spec$eff$type == "PT") fut <- list(type = "pt")
    if (spec$sided == 2 && fut$type == "spending")
        fut <- gs_twosided_beta(tt, spec$alpha, spec$beta, spec$eff, fut, spec$binding)
    spec$fut <- fut
    bK <- NA
    for (it in 1:30) {
        s <- gs_seq(tt, d$eta, spec, d$C, bK)
        if (fut$type != "cp" || (!is.na(bK) && abs(s$bK - bK) < 1e-9)) break
        bK <- s$bK
    }
    # classical shapes (WT, PT) do not spend alpha by themselves: re-solve the final
    # boundary so that exactly alpha is spent at the looks that actually happened
    if (spec$eff$type %in% c("WT", "PT")) {
        f <- function(b) {
            u <- s$up; u[K] <- b
            l <- if (spec$binding) s$lo else if (spec$sided == 2) -s$up else rep(-Inf, K)
            w <- if (spec$binding) s$wedge else rep(0, K)
            l[K] <- if (spec$sided == 2) -b else b
            w[K] <- if (spec$sided == 2) b else 0
            p <- gs_probs(tt, u, l, w, theta = 0)
            sum(p$up) + (if (spec$sided == 2) sum(p$lo) else 0) - spec$alpha
        }
        b <- stats::uniroot(f, c(0, 20), tol = 1e-12)$root
        s$up[K] <- b
        if (spec$sided == 1) s$lo[K] <- b else { s$lo[K] <- -b; s$wedge[K] <- b }
    }
    list(up = s$up, lo = s$lo, wedge = s$wedge)
}

# Two-sided level `level` repeated-CI constants at looks tt (same efficacy family).
gs_rci_bounds <- function(d, tt, level) {
    eff <- d$spec$eff
    if (eff$type == "triangular") return(rep(NA_real_, length(tt)))
    a <- 1 - level
    if (eff$type %in% c("WT", "PT")) {
        # constant of the two-sided WT family at the planned timing, shape at tt
        dd <- gs_design(timing = d$t, sided = 2, alpha = a, beta = d$spec$beta,
                        eff = list(type = "WT", delta = eff$delta))
        return(dd$C * tt^(eff$delta - 0.5))
    }
    spec <- list(sided = 2, alpha = a, beta = d$spec$beta, eff = eff, fut = list(type = "none"), binding = FALSE)
    gs_seq(tt, d$eta, spec)$up
}

# Gauss-Hermite nodes and weights (Golub-Welsch) for predictive power.
gauss_hermite <- function(n = 30) {
    i <- seq_len(n - 1)
    J <- matrix(0, n, n)
    J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
    e <- eigen(J, symmetric = TRUE)
    list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}

# Probability of crossing the efficacy boundary at a later look, given Z = z at time t.
# A binding futility rule stops the paths that reach it (lo, wedge); otherwise pass none.
gs_cp <- function(tt, up, t, z, drift, lo = NULL, wedge = NULL) {
    fut <- tt > t + 1e-12
    if (!any(fut)) return(NA_real_)
    n <- sum(fut)
    l <- if (is.null(lo)) rep(-Inf, n) else lo[fut]
    w <- if (is.null(wedge)) rep(0, n) else wedge[fut]
    sum(gs_probs(tt[fut], up[fut], l, w, theta = drift, start = list(z = z, I = t))$up)
}

gs_pp <- function(tt, up, t, z, gh = gauss_hermite(), lo = NULL, wedge = NULL) {
    m <- z / sqrt(t); s <- 1 / sqrt(t)
    sum(gh$w * vapply(m + s * sqrt(2) * gh$x, function(dr) gs_cp(tt, up, t, z, dr, lo, wedge), 0)) / sqrt(pi)
}

# Monitor observed looks.
#   looks: data.frame with t (spending time), z (oriented), est (natural), se (natural)
#   final: the last row is the final analysis (all remaining alpha is spent there)
gs_monitor <- function(d, looks, final = FALSE, dir = 1) {
    K <- nrow(looks)
    sided <- d$spec$sided
    tobs <- looks$t
    final <- final || tobs[K] >= 1 - 1e-9
    later <- d$t[d$t > tobs[K] + 1e-9]
    tt <- if (final) tobs else c(tobs, if (length(later)) later else 1)
    b <- gs_bounds_at(d, tt)
    level <- 1 - 2 * d$spec$alpha / sided
    crit <- gs_rci_bounds(d, tt, level)
    zq <- stats::qnorm(1 - (1 - level) / 2)
    gh <- gauss_hermite()
    nb <- !isTRUE(d$spec$binding)

    decision <- character(K); stopAt <- NA
    for (k in seq_len(K)) {
        z <- looks$z[k]
        last <- final && k == K
        eff <- if (sided == 2) abs(z) >= b$up[k] else z >= b$up[k]
        fut <- if (sided == 2) abs(z) < b$wedge[k] else z <= b$lo[k]
        decision[k] <- if (!is.na(stopAt)) "after" else if (eff) "efficacy" else if (last) "accept" else
            if (fut) (if (nb) "futilityNB" else "futility") else "continue"
        if (is.na(stopAt) && (eff || last || (fut && !nb))) stopAt <- k
    }
    out <- data.frame(t = tobs, z = looks$z, b = looks$z * sqrt(tobs),
                      up = b$up[seq_len(K)],
                      fut = if (sided == 2) b$wedge[seq_len(K)] else b$lo[seq_len(K)],
                      decision = decision,
                      rciLo = looks$est - crit[seq_len(K)] * looks$se,
                      rciHi = looks$est + crit[seq_len(K)] * looks$se)
    out$fut[out$fut <= -Inf | (sided == 2 & out$fut <= 1e-5)] <- NA
    out$fut[K][final] <- NA
    blo <- if (nb) NULL else b$lo
    bw <- if (nb) NULL else b$wedge
    out$cpH1 <- vapply(seq_len(K), function(k) gs_cp(tt, b$up, tobs[k], looks$z[k], d$eta, blo, bw), 0)
    out$cpTrend <- vapply(seq_len(K), function(k)
        gs_cp(tt, b$up, tobs[k], looks$z[k], looks$z[k] / sqrt(tobs[k]), blo, bw), 0)
    out$pp <- vapply(seq_len(K), function(k)
        if (tobs[k] < max(tt) - 1e-12) gs_pp(tt, b$up, tobs[k], looks$z[k], gh, blo, bw) else NA_real_, 0)
    # no conditional power once the trial has stopped (or at the final analysis)
    stopped <- !is.na(stopAt) & seq_len(K) >= stopAt
    out[stopped | (final & seq_len(K) == K), c("cpH1", "cpTrend", "pp")] <- NA

    res <- list(looks = out, stopAt = stopAt, final = final, tt = tt, bounds = b, level = level)
    if (!is.na(stopAt)) res$inference <- gs_final_inference(d, looks[seq_len(stopAt), ], b, level, dir)
    res
}

# Stagewise-ordering inference at the stopping look k* = nrow(looks).
gs_final_inference <- function(d, looks, b, level, dir = 1) {
    k <- nrow(looks)
    sided <- d$spec$sided
    # canonical correlation sqrt(t_j / t_k) from the information fractions (as rpact);
    # the scale comes from the standard error at the stopping look
    I <- looks$t / looks$t[k] / looks$se[k]^2
    zs <- looks$z[k]
    s <- if (sided == 2 && zs < 0) -1 else 1       # mirror a negative path (symmetric bounds)
    up <- b$up[seq_len(k)]
    lo <- if (isTRUE(d$spec$binding)) b$lo[seq_len(k)] else if (sided == 2) -up else rep(-Inf, k)
    wd <- if (isTRUE(d$spec$binding)) b$wedge[seq_len(k)] else rep(0, k)
    up[k] <- s * zs; lo[k] <- -Inf; wd[k] <- 0
    # P_theta(stop earlier through efficacy, or at k* with Z >= z*), theta oriented
    p <- function(theta) sum(gs_probs(I, up, lo, wd, theta = s * theta)$up)
    th <- dir * looks$est[k]
    se <- looks$se[k]
    root <- function(target) {
        f <- function(x) p(x) - target
        lo_ <- th - 10 * se; hi_ <- th + 10 * se
        if (f(lo_) * f(hi_) > 0) return(NA_real_)
        stats::uniroot(f, c(lo_, hi_), tol = 1e-10)$root
    }
    a2 <- (1 - level) / 2
    pv <- p(0)
    ci <- c(root(a2), root(1 - a2))
    mue <- root(0.5)
    # back to the natural scale: dir = -1 flips the sign and the order of the limits
    nat <- function(x) dir * x
    list(stage = k, p.value = if (sided == 2) min(1, 2 * min(pv, 1 - pv)) else pv,
         mue = nat(mue), lower = min(nat(ci)), upper = max(nat(ci)),
         naive = looks$est[k],
         naiveLower = looks$est[k] - stats::qnorm(1 - a2) * se,
         naiveUpper = looks$est[k] + stats::qnorm(1 - a2) * se)
}

# Interim data rows -> looks for gs_monitor(). n: cumulative size or events; either an
# estimate and its SE (natural scale; log HR for survival) or a Z statistic; optionally
# the information fraction (spending time). dir = sign of the design effect, so that the
# oriented Z is positive when the data favour H1.
monitor_looks <- function(n, est = NULL, se = NULL, z = NULL, frac = NULL, Nmax, Imax, dir = 1) {
    useEst <- !is.null(est) && !is.null(se)
    if (!useEst && is.null(z))
        stop("Interim data: give either an estimate and its standard error, or a Z statistic.")
    ok <- is.finite(n) & (if (useEst) is.finite(est) & is.finite(se) & se > 0 else is.finite(z))
    if (!is.null(frac)) ok <- ok & is.finite(frac)
    if (!any(ok)) stop("Interim data: no row has complete values.")
    o <- order(n[ok])
    n <- n[ok][o]
    if (any(n <= 0) || (length(n) > 1 && any(diff(n) <= 0)))
        stop("Interim data: the cumulative N (or events) must be positive and increase from row to row.")
    t <- if (!is.null(frac)) frac[ok][o] else n / Nmax
    if (any(t <= 0) || (length(t) > 1 && any(diff(t) <= 0)))
        stop("Interim data: the information fractions must be positive and increasing.")
    if (useEst) {
        est <- est[ok][o]; se <- se[ok][o]
        zz <- dir * est / se
    } else {
        z <- z[ok][o]
        se <- 1 / sqrt(t * Imax)
        est <- z * se
        zz <- dir * z
    }
    data.frame(n = n, t = t, z = zz, est = est, se = se)
}
