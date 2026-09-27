# Maximized SPRT for safety surveillance (Kulldorff et al. 2011; Silva & Kulldorff 2015),
# reimplemented from the published definitions (the Sequential package is GPL-2 only and is
# used only as a test oracle).
#
# Poisson: c observed and mu expected events under H0; LLR = (mu - c) + c log(c / mu) if c > mu.
#   Surveillance ends when mu reaches T; a signal needs at least M events.
# Binomial (self-controlled / matched): of n events, c fall in the risk window;
#   p0 = 1 / (1 + z) with matching ratio z; LLR = c log(c / (n p0)) + (n - c) log((n - c) / (n (1 - p0)))
#   if c / n > p0. Surveillance ends after N events.
# Continuous monitoring tests after every event; group monitoring at given looks.

llr_pois <- function(c, mu) ifelse(c > mu, (mu - c) + c * log(c / mu), 0)

llr_binom <- function(c, n, p0) {
    ph <- c / n
    xlogy <- function(x, y) ifelse(x > 0, x * log(y), 0)
    ifelse(ph > p0, xlogy(c, ph / p0) + xlogy(n - c, (1 - ph) / (1 - p0)), 0)
}

# ---- Poisson, continuous --------------------------------------------------------------

# Expected-count "deadline" u_k for the k-th event: signal iff it arrives before u_k.
maxsprt_u <- function(k, cv) {
    # solved on the log scale: the deadline can be tiny for small k and a large cv
    vapply(k, function(kk) exp(stats::uniroot(function(lm) llr_pois(kk, exp(lm)) - cv,
                                              c(-700, log(kk)), tol = 1e-12)$root), 0)
}

# Exact performance: probability of a signal, expected time to signal given a signal (in
# expected events under H0) and expected length of surveillance, for relative risk rr.
maxsprt_pois_perf <- function(cv, T, M = 1, rr = 1) {
    kmax <- ceiling(T + 20 * sqrt(T) + 50)
    u <- maxsprt_u(M:kmax, cv)
    ks <- M:kmax
    last <- which(u >= T)[1]
    if (is.na(last)) stop("increase kmax")
    tau <- c(u[seq_len(last - 1)], T); ks <- ks[seq_len(last)]
    dist <- 1                          # P(N = j, no signal yet) at the previous checkpoint, j = 0..
    prev <- 0; pSig <- 0; eSig <- 0
    for (i in seq_along(ks)) {
        k <- ks[i]; d <- tau[i] - prev; lam <- rr * d
        J <- seq_along(dist) - 1
        # counts at tau: convolve with Poisson(lam), up to k - 1 (no signal) or >= k (signal at event k)
        need <- k - J                  # events still needed to reach k
        pReach <- stats::ppois(need - 1, lam, lower.tail = FALSE)
        pSig <- pSig + sum(dist * pReach)
        # E[S_k 1{signal}] = sum_j dist_j (prev P(G <= d) + E[G 1{G <= d}]), G ~ Gamma(need, rr)
        eSig <- eSig + sum(dist * (prev * pReach + need / rr * stats::pgamma(d, need + 1, rr)))
        # no signal: fewer than k events by tau (truncated convolution)
        dist <- stats::convolve(dist, rev(stats::dpois(0:(k - 1), lam)), type = "open")[seq_len(k)]
        dist[dist < 0] <- 0
        prev <- tau[i]
    }
    list(power = pSig, timeToSignal = if (pSig > 0) eSig / pSig else NA,
         length = eSig + T * (1 - pSig), deadlines = data.frame(events = ks, mu = tau))
}

maxsprt_pois_cv <- function(T, M = 1, alpha = 0.05) {
    stats::uniroot(function(cv) maxsprt_pois_perf(cv, T, M)$power - alpha, c(0.05, 30), tol = 1e-10)$root
}

# ---- group monitoring (Poisson looks at cumulative expected mu_g; binomial at events n_g) ----

# Smallest count that signals at each look (signal when the LLR exceeds cv).
maxsprt_thresholds <- function(model, cv, looks, M = 1, p0 = NULL) {
    vapply(looks, function(l) {
        c <- max(M, if (model == "poisson") floor(l) else 0):(if (model == "poisson") ceiling(l + 20 * sqrt(l) + 60) else l)
        L <- if (model == "poisson") llr_pois(c, l) else llr_binom(c, l, p0)
        k <- which(L > cv + 1e-10)[1]
        if (is.na(k)) Inf else c[k]
    }, 0)
}

# Exact probability of a signal and expected sizes for signalling thresholds at the looks.
maxsprt_group_perf <- function(model, thr, looks, rr = 1, p0 = NULL) {
    G <- length(looks)
    top <- ceiling(max(looks) + 20 * sqrt(max(looks)) + 60)
    dist <- 1; prev <- 0; pSig <- 0; eSig <- 0
    p <- if (model == "binomial") rr * p0 / (rr * p0 + 1 - p0) else NA
    for (g in seq_len(G)) {
        step <- looks[g] - prev
        inc <- if (model == "poisson") stats::dpois(0:top, rr * step) else stats::dbinom(0:step, step, p)
        dist <- stats::convolve(dist, rev(inc), type = "open"); dist[dist < 0] <- 0
        c <- seq_along(dist) - 1
        s <- c >= thr[g]
        ps <- sum(dist[s])
        pSig <- pSig + ps; eSig <- eSig + looks[g] * ps
        dist[s] <- 0
        if (model == "poisson") dist <- dist[seq_len(min(length(dist), top + 1))]
        prev <- looks[g]
    }
    list(power = pSig, timeToSignal = if (pSig > 0) eSig / pSig else NA, length = eSig + max(looks) * (1 - pSig))
}

# Critical value for group (and binomial continuous) designs. The attained alpha is a step
# function of cv: every cv between two consecutive cell LLRs gives the same design. As the
# Sequential package, report the lower end (the largest LLR that does not signal), with a
# signal when the LLR exceeds it.
maxsprt_group_cv <- function(model, looks, M = 1, alpha = 0.05, p0 = NULL) {
    cand <- sort(unique(unlist(lapply(looks, function(l) {
        c <- max(M, 1):(if (model == "poisson") ceiling(l + 20 * sqrt(l) + 60) else l)
        L <- if (model == "poisson") llr_pois(c, l) else llr_binom(c, l, p0)
        L[L > 0]
    }))))
    att <- function(cv) maxsprt_group_perf(model, maxsprt_thresholds(model, cv, looks, M, p0), looks, 1, p0)$power
    # bisection over cand = c(0, cell LLRs): smallest cv with attained alpha <= nominal
    cand <- c(0, cand)
    lo <- 1; hi <- length(cand)
    while (lo < hi) {
        mid <- (lo + hi) %/% 2
        if (att(cand[mid]) <= alpha) hi <- mid else lo <- mid + 1
    }
    list(cv = cand[lo], alpha = att(cand[lo]))
}

# ---- one-stop design ----------------------------------------------------------------

# model "poisson" or "binomial"; monitoring "continuous" or "group"; size T (Poisson) or N
# (binomial); looks (group) as cumulative expected counts or events; rr = relative risks.
maxsprt_design <- function(model, monitoring, alpha, size, M = 1, looks = NULL, z = 1, rr = c(1.5, 2, 3)) {
    p0 <- 1 / (1 + z)
    if (model == "poisson" && monitoring == "continuous") {
        cv <- maxsprt_pois_cv(size, M, alpha)
        att <- alpha
        thr <- maxsprt_pois_perf(cv, size, M)$deadlines
    } else {
        lk <- if (monitoring == "continuous") seq_len(size) else looks
        g <- maxsprt_group_cv(model, lk, M, alpha, p0)
        cv <- g$cv; att <- g$alpha
        thr <- data.frame(look = lk, threshold = maxsprt_thresholds(model, cv, lk, M, p0))
    }
    d <- list(cv = cv, alpha = att, thresholds = thr, model = model, monitoring = monitoring,
              p0 = p0, M = M, size = size)
    d$performance <- maxsprt_design_perf(d, rr)
    d
}

# Power, expected time to signal and expected length of a solved design at relative risks rr.
maxsprt_design_perf <- function(d, rr) {
    perf <- lapply(rr, function(r)
        if (d$model == "poisson" && d$monitoring == "continuous") maxsprt_pois_perf(d$cv, d$size, d$M, r)
        else maxsprt_group_perf(d$model, d$thresholds$threshold, d$thresholds$look, r, d$p0))
    data.frame(rr = rr, power = vapply(perf, `[[`, 0, "power"),
               timeToSignal = vapply(perf, `[[`, 0, "timeToSignal"),
               length = vapply(perf, `[[`, 0, "length"))
}

# Monitoring data: cumulative observed and expected (Poisson) or cases and events (binomial).
maxsprt_analyse <- function(des, obs, exp) {
    llr <- if (des$model == "poisson") llr_pois(obs, exp) else llr_binom(obs, exp, des$p0)
    ok <- obs >= des$M                 # M: minimum number of cases for a signal (both models)
    sig <- ok & llr > des$cv + 1e-10
    first <- which(sig)[1]
    rr <- if (des$model == "poisson") obs / exp else {
        # risk-window share c/n and matching ratio z give RR = (c / n) z / (1 - c / n)
        ph <- obs / exp; ifelse(ph < 1, ph * (1 - des$p0) / ((1 - ph) * des$p0), Inf)
    }
    data.frame(obs = obs, exp = exp, llr = llr, signal = sig, rr = rr,
               after = !is.na(first) & seq_along(obs) > first)
}

# ---- orchestration for the jamovi analysis ---------------------------------------------

# Looks for group monitoring: custom cumulative values, or nLooks equally spaced up to size.
maxsprt_looks <- function(o, size) {
    v <- parse_numbers(o$looks, "Looks")
    if (!length(v)) v <- size * seq_len(o$nLooks) / o$nLooks
    if (any(v <= 0) || (length(v) > 1 && any(diff(v) <= 0)))
        stop("The looks must be positive and increasing.")
    if (o$model == "binomial") {
        if (any(abs(v - round(v)) > 1e-9)) stop("For the binomial model the looks are numbers of cases (integers).")
        v <- round(v)
    }
    v
}

maxsprt_compute <- function(o) {
    size <- if (o$model == "poisson") o$survLength else o$N
    if (o$model == "binomial" && o$z <= 0) stop("The matching ratio must be positive.")
    rr <- parse_numbers(o$rr, "Relative risks")
    if (!length(rr) || any(rr <= 0)) stop("Give one or more positive relative risks.")
    looks <- if (o$monitoring == "group") maxsprt_looks(o, size) else NULL
    if (!is.null(looks)) size <- max(looks)
    d <- maxsprt_design(o$model, o$monitoring, o$alpha, size, o$M, looks, o$z, rr)
    d$looks <- looks
    d
}

# Boundary curve on the count scale: the number of cases that signals at each expected count
# (Poisson) or number of cases so far (binomial).
maxsprt_boundary <- function(d) {
    th <- d$thresholds
    if (d$model == "poisson" && d$monitoring == "continuous")
        return(data.frame(x = c(0, th$mu), y = c(th$events, utils::tail(th$events, 1))))
    ok <- is.finite(th$threshold)
    data.frame(x = th$look[ok], y = th$threshold[ok])
}

# Data: cumulative observed cases and expected counts (Poisson) or total cases (binomial).
maxsprt_monitor <- function(d, obs, exp) {
    ok <- is.finite(obs) & is.finite(exp)
    if (!any(ok)) stop("Surveillance data: no row has both values.")
    obs <- obs[ok]; exp <- exp[ok]
    if (any(exp <= 0) || any(obs < 0)) stop("Surveillance data: the counts must be non-negative and the expected counts positive.")
    if (length(exp) > 1 && (any(diff(exp) <= 0) || any(diff(obs) < 0)))
        stop("Surveillance data: the cumulative counts must not decrease from row to row.")
    if (d$model == "binomial" && any(obs > exp))
        stop("Surveillance data: the cases in the risk window cannot exceed the total number of cases.")
    m <- maxsprt_analyse(d, obs, exp)
    m$beyond <- exp > d$size * (1 + 1e-9)
    m
}
