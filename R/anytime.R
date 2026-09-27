# Anytime-valid monitoring: two-sided Gaussian-mixture e-processes and confidence sequences.
#
# Known SD: Robbins' (1970) normal-mixture likelihood ratio, the mSPRT of Johari et al. (2022):
#   E = (1 + nEff g)^(-1/2) exp(nEff g z^2 / (2 (1 + nEff g))),  z = sqrt(nEff) (mean - mu0) / sigma,
# with the mixing variance g on the standardised scale (effect / sigma).
# Unknown SD: the t-statistic version with a right-Haar prior on sigma (Pérez-Ortiz et al. 2024;
# Grünwald, de Heide & Koolen 2024), in closed form:
#   E = (1 + nEff g)^(-1/2) ((1 + t^2 / nu) / (1 + t^2 / (nu (1 + nEff g))))^((nu + 1) / 2).
# Both equal safestats' eGauss e-values; the default g = delta^2 is safestats' choice. Reject H0
# the first time E >= 1 / alpha (Ville's inequality). Confidence sequences invert the test.

av_logE_z <- function(z, nEff, g) -0.5 * log1p(nEff * g) + nEff * g * z^2 / (2 * (1 + nEff * g))

av_logE_t <- function(t, nEff, nu, g) {
    a <- 1 + nEff * g
    -0.5 * log(a) + ((nu + 1) / 2) * (log1p(t^2 / nu) - log1p(t^2 / (nu * a)))
}

# Squared critical value of the statistic (z or t) at which E reaches 1 / alpha.
av_crit2_z <- function(nEff, g, alpha) {
    a <- 1 + nEff * g
    2 * a / (nEff * g) * (log(1 / alpha) + 0.5 * log(a))
}

av_crit2_t <- function(nEff, nu, g, alpha) {
    a <- 1 + nEff * g
    L <- (2 / (nu + 1)) * (log(1 / alpha) + 0.5 * log(a))
    ifelse(exp(L) < a, nu * (exp(L) - 1) / (1 - exp(L) / a), Inf)
}

# ---- data --------------------------------------------------------------------------------

# Path for data in arrival order: one sample against mu0 (g NULL) or two samples (group 2
# minus group 1). sigma = NA uses the t-statistic version.
av_analyse <- function(x, grp = NULL, mu0 = 0, sigma = NA, g, alpha = 0.05) {
    known <- !is.na(sigma)
    if (is.null(grp)) {
        n <- seq_along(x)
        est <- cumsum(x) / n
        nEff <- n
        v <- (cumsum(x^2) - n * est^2) / pmax(n - 1, 1)
        s <- if (known) rep(sigma, length(x)) else sqrt(pmax(v, 0))
        nu <- n - 1
        ok <- if (known) n >= 1 else n >= 2 & v > 0
        est <- est - mu0
    } else {
        st <- two_sample_stats(x, grp)
        nEff <- ifelse(st$n1 > 0 & st$n2 > 0, st$n1 * st$n2 / (st$n1 + st$n2), 0)
        est <- st$diff
        nu <- st$n1 + st$n2 - 2
        s <- if (known) rep(sigma, length(x)) else st$sp
        ok <- st$n1 >= 1 & st$n2 >= 1 & (known | (nu >= 1 & !is.na(st$sp) & st$sp > 0))
    }
    stat <- ifelse(ok, sqrt(nEff) * est / s, 0)
    logE <- rep(0, length(x))
    logE[ok] <- if (known) av_logE_z(stat[ok], nEff[ok], g) else av_logE_t(stat[ok], nEff[ok], nu[ok], g)
    c2 <- rep(Inf, length(x))
    c2[ok] <- if (known) av_crit2_z(nEff[ok], g, alpha) else av_crit2_t(nEff[ok], nu[ok], g, alpha)
    half <- ifelse(ok, sqrt(c2) * s / sqrt(nEff), Inf)
    hit <- which(logE >= log(1 / alpha))[1]
    list(n = seq_along(x), est = ifelse(ok, est + (if (is.null(grp)) mu0 else 0), NA), stat = stat,
         e = exp(logE), p = pmin(1, exp(-cummax(logE))),
         lower = ifelse(ok, est - half, -Inf) + (if (is.null(grp)) mu0 else 0),
         upper = ifelse(ok, est + half, Inf) + (if (is.null(grp)) mu0 else 0),
         reject = !is.na(hit), at = hit, ok = ok)
}

# ---- design --------------------------------------------------------------------------------

# Known SD: exact probabilities of rejecting H0 by each n up to nMax at the standardised effect
# delta, with a test after every observation (one sample) or every pair (two samples, 1:1).
# The mixture test rejects when |z| >= sqrt(crit2), a two-sided boundary for the gsd engine.
av_exact_z <- function(delta, g, alpha, nMax, two = FALSE) {
    m <- seq_len(if (two) floor(nMax / 2) else nMax)
    nEff <- if (two) m / 2 else m
    c <- sqrt(av_crit2_z(nEff, g, alpha))
    p <- gs_probs(nEff, c, -c, theta = delta)
    step <- if (two) 2 else 1
    list(n = m * step, reject = cumsum(p$up + p$lo))
}

# Unknown SD: seeded simulation, same schedule; returns the stopping n (Inf: not by nMax).
av_sim_t <- function(delta, g, alpha, nMax, two = FALSE, nsim = 1000, seed = 1) {
    set.seed(seed)
    m <- seq_len(if (two) floor(nMax / 2) else nMax)
    out <- numeric(nsim)
    for (i in seq_len(nsim)) {
        if (two) {
            xc <- stats::rnorm(length(m)); xt <- stats::rnorm(length(m), delta)
            m1 <- cumsum(xc) / m; m2 <- cumsum(xt) / m
            ss <- cumsum(xc^2) - m * m1^2 + cumsum(xt^2) - m * m2^2
            nu <- 2 * m - 2; sp <- sqrt(pmax(ss, 0) / pmax(nu, 1))
            tt <- (m2 - m1) / (sp * sqrt(2 / m)); nEff <- m / 2
        } else {
            x <- stats::rnorm(length(m), delta)
            mm <- cumsum(x) / m; v <- (cumsum(x^2) - m * mm^2) / pmax(m - 1, 1)
            tt <- mm / sqrt(v / m); nEff <- m; nu <- m - 1
        }
        ok <- nu >= 1 & is.finite(tt)
        lE <- rep(0, length(m)); lE[ok] <- av_logE_t(tt[ok], nEff[ok], nu[ok], g)
        hit <- which(lE >= log(1 / alpha))[1]
        out[i] <- if (is.na(hit)) Inf else m[hit] * (if (two) 2 else 1)
    }
    out
}

# Design summary at delta: power by nMax, expected n (capped at nMax), median and 90 % point
# of the stopping n, the n needed for the target power, and fixed-sample references.
av_design <- function(sdKnown, two, delta, alpha, power, g = NA, nMax = NA, nsim = 1000, seed = 1) {
    if (delta <= 0) stop("The standardised effect must be positive.")
    if (is.na(g) || g <= 0) g <- delta^2
    za <- stats::qnorm(1 - alpha / 2); zb <- stats::qnorm(power)
    nFix <- (if (two) 4 else 1) * ((za + zb) / delta)^2
    if (!sdKnown) nFix <- nFix + (if (two) 1 else 0.5) * za^2       # Guenther's t-test correction
    nFix <- ceiling(nFix - 1e-9)
    if (two) nFix <- 2 * ceiling(nFix / 2)
    if (is.na(nMax)) nMax <- max(20, ceiling(4 * nFix))
    if (sdKnown) {
        h1 <- av_exact_z(delta, g, alpha, nMax, two)
        h0 <- av_exact_z(0, g, alpha, nMax, two)
        n <- h1$n; F1 <- h1$reject
        # expected n capped at nMax: sum of P(N > n) over the looks, times the step
        step <- if (two) 2 else 1
        asn <- step * (1 + sum(1 - F1[-length(F1)]))
        alpha0 <- utils::tail(h0$reject, 1)
    } else {
        s1 <- av_sim_t(delta, g, alpha, nMax, two, nsim, seed)
        s0 <- av_sim_t(0, g, alpha, nMax, two, nsim, seed + 1)
        n <- seq_len(if (two) floor(nMax / 2) else nMax) * (if (two) 2 else 1)
        F1 <- vapply(n, function(q) mean(s1 <= q), 0)
        asn <- mean(pmin(s1, max(n)))
        alpha0 <- mean(is.finite(s0))
    }
    q <- function(pr) { k <- which(F1 >= pr)[1]; if (is.na(k)) NA else n[k] }
    # a single analysis with the e-value at n (no interim looks): power 1 - beta at the smallest n
    single <- function(nn) {
        nEff <- if (two) nn / 4 else nn
        c <- sqrt(if (sdKnown) av_crit2_z(nEff, g, alpha) else av_crit2_t(nEff, if (two) nn - 2 else nn - 1, g, alpha))
        mu <- delta * sqrt(nEff)
        if (sdKnown) stats::pnorm(mu - c) + stats::pnorm(-mu - c)
        else stats::pt(c, if (two) nn - 2 else nn - 1, ncp = mu, lower.tail = FALSE) + stats::pt(-c, if (two) nn - 2 else nn - 1, ncp = mu)
    }
    ns <- seq(if (two) 4 else 3, 20 * nFix, by = if (two) 2 else 1)
    nSingle <- ns[which(vapply(ns, single, 0) >= power)[1]]
    list(g = g, nFix = nFix, nSingle = nSingle, nMax = max(n), power = utils::tail(F1, 1), alpha = alpha0,
         asn = asn, median = q(0.5), q90 = q(0.9), nPlan = q(power), curve = data.frame(n = n, reject = F1),
         method = if (sdKnown) "exact" else "simulation")
}

# ---- orchestration for the jamovi analysis ---------------------------------------------

av_compute <- function(o) {
    two <- o$samples == "two"
    known <- isTRUE(o$sdKnown)
    if (known && o$sigma <= 0) stop("The known standard deviation must be positive.")
    g <- if (o$mixing == "custom") o$g else NA
    if (!is.na(g) && g <= 0) stop("The mixing variance g must be positive.")
    nMax <- if (o$nMax > 0) o$nMax else NA
    av_design(known, two, abs(o$delta), o$alpha, o$power, g, nMax, o$nsim, o$seed)
}

# Rows to display: about `target` evenly spaced observations, always including the first
# rejection and the last observation (every = 0), or every `every`-th.
av_rows <- function(N, at, every = 0, target = 20) {
    step <- if (every > 0) every else max(1, ceiling(N / target))
    sort(unique(c(seq(step, N, by = step), if (!is.na(at)) at, N)))
}
